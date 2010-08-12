require 'rubygems' 
require 'net/http'
# require 'net/geoip'
require 'open-uri'
require 'erb'        

##
# Add geocoding functionality (via Google) to any object.
#
module Geocoder
  ## using geocoding / reverse geocoding from Google Map API V3 
  GOOGLE_MAP_JSON_BASE_URL = 'http://maps.google.com/maps/api/geocode/json?'
  ##
  # Implementation of 'included' hook method.
  #
  def self.included(base)
    base.extend ClassMethods
    base.class_eval do

      # scope: geocoded objects
      send(Geocoder.scope_method_name, :geocoded,
        :conditions => "#{geocoder_options[:latitude]} IS NOT NULL " +
          "AND #{geocoder_options[:longitude]} IS NOT NULL")

      # scope: not-geocoded objects
      send(Geocoder.scope_method_name, :not_geocoded,
        :conditions => "#{geocoder_options[:latitude]} IS NULL " +
          "OR #{geocoder_options[:longitude]} IS NULL")
      
      ##
      # Find all objects within a radius (in miles) of the given location
      # (address string). Location (the first argument) may be either a string
      # to geocode or an array of coordinates (<tt>[lat,long]</tt>).
      #
      send(Geocoder.scope_method_name, :near, lambda{ |location, *args|
        latitude, longitude = location.is_a?(Array) ?
          location : Geocoder.fetch_coordinates(location)
        if latitude and longitude
          near_scope_options(latitude, longitude, *args)
        else
          {}
        end
      })
    end
  end
    
  ##
  # Methods which will be class methods of the including class.
  #
  module ClassMethods
    ##
    # Get options hash suitable for passing to ActiveRecord.find to get
    # records within a radius (in miles) of the given point.
    # Options hash may include:
    # 
    # +order+     :: column(s) for ORDER BY SQL clause
    # +limit+     :: number of records to return (for LIMIT SQL clause)
    # +offset+    :: number of records to skip (for LIMIT SQL clause)
    #
    def near_scope_options(latitude, longitude, radius = 20, options = {})
      if ActiveRecord::Base.connection.adapter_name == "mysql"
        approx_near_scope_options(latitude, longitude, radius, options)
      else
        full_near_scope_options(latitude, longitude, radius, options)
      end
    end


    private # ----------------------------------------------------------------

    ##
    # Scope options hash for use with a database that supports POWER(),
    # SQRT(), PI(), and trigonometric functions (SIN(), COS(), and ASIN()).
    # 
    # Taken from the excellent tutorial at:
    # http://www.scribd.com/doc/2569355/Geo-Distance-Search-with-MySQL
    #
    def full_near_scope_options(latitude, longitude, radius, options)
      options[:order] ||= 'distance ASC'
      lat_attr = geocoder_options[:latitude]
      lon_attr = geocoder_options[:longitude]
      distance = "3956 * 2 * ASIN(SQRT(" +
        "POWER(SIN((#{latitude} - #{lat_attr}) * " +
        "PI() / 180 / 2), 2) + COS(#{latitude} * PI()/180) * " +
        "COS(#{lat_attr} * PI() / 180) * " +
        "POWER(SIN((#{longitude} - #{lon_attr}) * " +
        "PI() / 180 / 2), 2) ))"
      {
        :select => "*, #{distance} AS distance",
        :conditions => \
          ["#{lat_attr} BETWEEN ? AND ? AND #{lon_attr} BETWEEN ? AND ?"] +
          coordinate_bounds(latitude, longitude, radius),
        :having => "#{distance} <= #{radius}",
        :order  => options[:order],
        :limit  => limit_clause(options)
      }
    end

    ##
    # Scope options hash for use with a database without trigonometric
    # functions, like SQLite. Approach is to find objects within a square
    # rather than a circle, so results are very approximate (will include
    # objects outside the given radius).
    #
    def approx_near_scope_options(latitude, longitude, radius, options)
      lat_attr = geocoder_options[:latitude]
      lon_attr = geocoder_options[:longitude]
      {
        :conditions => \
          ["#{lat_attr} BETWEEN ? AND ? AND #{lon_attr} BETWEEN ? AND ?"] +
          coordinate_bounds(latitude, longitude, radius),
        :order  => options[:order],
        :limit  => limit_clause(options)
      }
    end
    
    ##
    # Get the rough high/low lat/long bounds for a geographic point and
    # radius. Returns an array: <tt>[lat_lo, lat_hi, lon_lo, lon_hi]</tt>.
    # Used to constrain search to a (radius x radius) square.
    #
    def coordinate_bounds(latitude, longitude, radius)
      radius = radius.to_f
      factor = (Math::cos(latitude * Math::PI / 180.0) * 69.0).abs
      [
        latitude  - (radius / 69.0),
        latitude  + (radius / 69.0),
        longitude - (radius / factor),
        longitude + (radius / factor)
      ]
    end
    
    ##
    # Build the limit clause for a query based on the same options hash
    # passed to the x_near_scope_options methods.
    #
    def limit_clause(options)
      if options[:limit] or options[:offset]
        "#{options[:offset].to_i},#{options[:limit].to_i}"
      end
    end
  end

  ##
  # Read the coordinates [lat,lon] of an object. This is not great but it
  # seems cleaner than polluting the instance method namespace.
  #
  def read_coordinates
    [:latitude, :longitude].map{ |i| send self.class.geocoder_options[i] }
  end

  ##
  # Is this object geocoded? (Does it have latitude and longitude?)
  #
  def geocoded?
    read_coordinates.compact.size > 0
  end
  
  ##
  # Calculate the distance from the object to a point (lat,lon).
  # Valid units are defined in <tt>distance_between</tt> class method.
  #
  def distance_to(lat, lon, units = :km)
    return nil unless geocoded?
    mylat,mylon = read_coordinates
    Geocoder.distance_between(mylat, mylon, lat, lon, :units => units)
  end
  
  ##
  # Get other geocoded objects within a given radius (in miles). Takes a
  # radius (in miles) and options for passing to the +near+ scope
  # (<tt>:order</tt>, <tt>:limit</tt>, and <tt>:offset</tt>).
  #
  def nearbys(radius = 20, options = {})
    return [] unless geocoded?
    options = {:conditions => ["id != ?", id]}.merge(options)
    self.class.near(read_coordinates, radius, options) - [self]
  end
  
  ##
  # Fetch coordinates and assign +latitude+ and +longitude+. Also returns
  # coordinates as an array: <tt>[lat, lon]</tt>.
  #
  def fetch_coordinates(save = false)
    coords = Geocoder.fetch_coordinates(
      send(self.class.geocoder_options[:method_name])
    )
    unless coords.blank?
      method = (save ? "update" : "write") + "_attribute"
      send method, self.class.geocoder_options[:latitude],  coords[3][0]
      send method, self.class.geocoder_options[:longitude], coords[3][1]
    end
    coords
  end

  ##
  # Fetch coordinates and update (save) +latitude+ and +longitude+ data.
  #
  def fetch_coordinates!
    fetch_coordinates(true)
  end
  
  ##
  # Fetch all data and assign +formatted_address+ +geometry+ +location_type+ +latitude+ and +longitude+. 
  # Don't use these two methods if you don't have +formatted_address+ +geometry+ +location_type+ set up
  def fetch_all(save = false)
    coords = Geocoder.fetch_coordinates(
      send(self.class.geocoder_options[:method_name])
    )
    # Returned 4 data with in 1 array
    unless coords.blank?
      method = (save ? "update" : "write") + "_attribute"
      send method, self.class.geocoder_options[:formatted_address], coords[0]
      send method, self.class.geocoder_options[:geometry], coords[1]
      # send method, self.class.geocoder_options[:location_type], coords[2]
      send method, self.class.geocoder_options[:latitude], coords[3][0]
      send method, self.class.geocoder_options[:longitude], coords[3][1]
    end
    coords
  end

  ##
  # Fetch coordinates and update (save) +formatted_address+ +geometry+ +location_type+ +latitude+ and +longitude+.
  #
  def fetch_all!
    fetch_all(true)
  end
  ##
  # Query Google for the coordinates of the given phrase.
  # Returns array [lat,lon] if found, nil if not found or if network error.
  #
  def self.fetch_coordinates(query)
    return nil if query.blank?
    return nil unless doc = self.search(query)
    
    # make sure search found a result
    e = doc['status']
    return nil unless (e and e == "OK")

    # isolate the first which suggested by Google of the results
    formatted_address = doc['results'][0]['formatted_address']
    geometry = (doc['results'][0]['geometry']).to_json
    lat = doc['results'][0]['geometry']['location']['lat']
    lng = doc['results'][0]['geometry']['location']['lng']
    # location_type = doc['results'][0]['geometry']['location_type']
    coords = lat.to_f, lng.to_f
    return formatted_address, geometry, coords
  end
  
  ##
  # Calculate the distance between two points on Earth (Haversine formula).
  # Takes two sets of coordinates and an options hash:
  # 
  # <tt>:units</tt> :: <tt>:mi</tt> (default) or <tt>:km</tt>
  #
  def self.distance_between(lat1, lon1, lat2, lon2, options = {})
    
    # set default options to km
    options[:units] ||= :km
    
    # define conversion factors
    units = { :mi => 3956, :km => 6371 }
    
    # convert degrees to radians
    lat1 = to_radians(lat1)
    lon1 = to_radians(lon1)
    lat2 = to_radians(lat2)
    lon2 = to_radians(lon2)
    
    # compute distances
    dlat = (lat1 - lat2).abs
    dlon = (lon1 - lon2).abs
    
    a = (Math.sin(dlat / 2))**2 + Math.cos(lat1) *
        (Math.sin(dlon / 2))**2 * Math.cos(lat2)  
    c = 2 * Math.atan2( Math.sqrt(a), Math.sqrt(1-a))  
    c * units[options[:units]]
  end
  
  ##
  # Compute the geographic center (aka geographic midpoint, center of
  # gravity) for an array of geocoded objects and/or [lat,lon] arrays
  # (can be mixed). Any objects missing coordinates are ignored. Follows
  # the procedure documented at http://www.geomidpoint.com/calculation.html.
  #
  def self.geographic_center(points)
  
    # convert objects to [lat,lon] arrays and remove nils
    points = points.map{ |p|
      p.is_a?(Array) ? p : (p.geocoded?? p.read_coordinates : nil)
    }.compact
    
    # convert degrees to radians
    points.map!{ |p| [to_radians(p[0]), to_radians(p[1])] }
    
    # convert to Cartesian coordinates
    x = []; y = []; z = []
    points.each do |p|
      x << Math.cos(p[0]) * Math.cos(p[1])
      y << Math.cos(p[0]) * Math.sin(p[1])
      z << Math.sin(p[0])
    end

    # compute average coordinate values
    xa, ya, za = [x,y,z].map do |c|
      c.inject(0){ |tot,i| tot += i } / c.size.to_f
    end
    
    # convert back to latitude/longitude
    lon = Math.atan2(ya, xa)
    hyp = Math.sqrt(xa**2 + ya**2)
    lat = Math.atan2(za, hyp)
    
    # return answer in degrees
    [to_degrees(lat), to_degrees(lon)]
  end

  ##
  # Convert degrees to radians.
  #
  def self.to_radians(degrees)
    degrees * (Math::PI / 180)
  end
  
  ##
  # Convert radians to degrees.
  #
  def self.to_degrees(radians)
    (radians * 180.0) / Math::PI
  end
  
  ##
  # Query Google for geographic information about the given phrase.
  #
  def self.search(*array)
    query = array[0] # Fetch the first paramter as address or latlng
    sensor = array[1] || "false" # Set default sensor value to false
    language = array[2] || "zh-CN" # Set default language to Simplified Chinese 
    if res = _fetch_json(query, sensor, language) 
      ActiveSupport::JSON.decode(res) #['results'][0]
    end
  end
  
  ##
  # Request an JSON geo search result from Google.
  # This method is not intended for general use (prefer Geocoder.search).
  #
  def self._fetch_json(query, sensor, language)
    params = {
      :sensor => sensor,
      :language => language
    }
    # determine what service should be used. If query is lat and lng, need use reverse geocoding
    if Geocoder.lat_lng_reg_match(query)
      params[:latlng] = query
    else
      params[:address] = query
    end
    url = Geocoder::GOOGLE_MAP_JSON_BASE_URL + params.to_query
    
    # Query geocoder and make sure it responds quickly.
    begin
      resp = nil
      timeout(3) do
        Net::HTTP.get_response(URI.parse(url)).body
      end
    rescue SocketError, TimeoutError
      return nil
    end
      
  end
  
  ##
  # Name of the ActiveRecord scope method.
  #
  def self.scope_method_name
    Rails.version.starts_with?("3") ? :scope : :named_scope
  end
  
  def self.lat_lng_reg_match(latlng)
    # (Latitude,Longitude) Regex 
    # /^-?([1-8]?[1-9]|[1-9]0)\.{1}\d{1,6},-?([1]?[1-7][1-9]|[1]?[1-8][0]|[1-9]?[0-9])\.{1}\d{1,6}$/
    lat_reg = '-?([1-8]?[1-9]|[1-9]0)\.{1}\d{1,6}'
    lng_reg = '-?([1]?[1-7][1-9]|[1]?[1-8][0]|[1-9]?[0-9])\.{1}\d{1,6}'
    latlng =~ /^#{lat_reg},#{lng_reg}$/ ? true : false
  end
end

##
# Add geocoded_by method to ActiveRecord::Base so Geocoder is accessible.
#
ActiveRecord::Base.class_eval do
  
  ##
  # Set attribute names and include the Geocoder module.
  #
  def self.geocoded_by(method_name = :location, options = {})
    class_inheritable_reader :geocoder_options
    write_inheritable_attribute :geocoder_options, {
      :method_name => method_name,
      :latitude    => options[:latitude]  || :latitude,
      :longitude   => options[:longitude] || :longitude, 
      # Add new option for ActiveRecord
      :formatted_address  => options[:formatted_address] || :formatted_address, 
      :geometry => options[:geometry] || :geometry,
      # :location_type => options[:location_type] || :location_type
    }
    include Geocoder
  end
end
