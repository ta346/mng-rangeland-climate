import ee 
import pandas as pd

def create_reduce_region_function(geometry,
                                  reducer=ee.Reducer.mean(),
                                  scale=1000,
                                  crs='EPSG:4326',
                                  bestEffort=True,
                                  maxPixels=1e13,
                                  tileScale=4):
  """Creates a region reduction function.

  Creates a region reduction function intended to be used as the input function
  to ee.ImageCollection.map() for reducing pixels intersecting a provided region
  to a statistic for each image in a collection. See ee.Image.reduceRegion()
  documentation for more details.

  Args:
    geometry:
      An ee.Geometry that defines the region over which to reduce data.
    reducer:
      Optional; An ee.Reducer that defines the reduction method.
    scale:
      Optional; A number that defines the nominal scale in meters of the
      projection to work in.
    crs:
      Optional; An ee.Projection or EPSG string ('EPSG:5070') that defines
      the projection to work in.
    bestEffort:
      Optional; A Boolean indicator for whether to use a larger scale if the
      geometry contains too many pixels at the given scale for the operation
      to succeed.
    maxPixels:
      Optional; A number specifying the maximum number of pixels to reduce.
    tileScale:
      Optional; A number representing the scaling factor used to reduce
      aggregation tile size; using a larger tileScale (e.g. 2 or 4) may enable
      computations that run out of memory with the default.

  Returns:
    A function that accepts an ee.Image and reduces it by region, according to
    the provided arguments.
  """

  def reduce_region_function(img):
    """Applies the ee.Image.reduceRegion() method.

    Args:
      img:
        An ee.Image to reduce to a statistic by region.

    Returns:
      An ee.Feature that contains properties representing the image region
      reduction results per band and the image timestamp formatted as
      milliseconds from Unix epoch (included to enable time series plotting).
    """

    stat = img.reduceRegion(
        reducer=reducer,
        geometry=geometry,
        scale=scale,
        crs=crs,
        bestEffort=bestEffort,
        maxPixels=maxPixels,
        tileScale=tileScale)

    return ee.Feature(geometry, stat).set({'millis': img.date().millis()})
  return reduce_region_function

def reduce_regions_function(reduction_regions,
                                  reducerAll = True, 
                                  scale=1000,
                                  crs='EPSG:4326',
                                  tileScale=1):
  """Creates multiple regions reduction function.

  Creates multiple regions reduction function intended to be used as the input function
  to ee.ImageCollection.map() for reducing pixels intersecting a provided region
  to a statistic for each image in a collection. See ee.Image.reduceRegions()
  documentation for more details.

  Args:
    collection:
      FeatureCollection - The features to reduce over
    reducer:
      Optional; An ee.Reducer that defines the reduction method.
    scale:
      Optional; A number that defines the nominal scale in meters of the
      projection to work in.
    crs:
      Optional; An ee.Projection or EPSG string ('EPSG:5070') that defines
      the projection to work in.
    tileScale:
      Optional; A number representing the scaling factor used to reduce
      aggregation tile size; using a larger tileScale (e.g. 2 or 4) may enable
      computations that run out of memory with the default.

  Returns:
    A function that accepts an ee.Image and reduces it by region, according to
    the provided arguments.
  """

  

  if reducerAll is True:
    reducerAll = (ee.Reducer.mean().combine(ee.Reducer.median(), '', True) 
              .combine(ee.Reducer.minMax(), '', True) 
              .combine(ee.Reducer.stdDev(), '', True) 
              .combine(ee.Reducer.intervalMean(10, 90).setOutputs(['int_mean_10_90']), '', True) 
              .combine(ee.Reducer.count(),'', True))
  else:
    reducerAll = (ee.Reducer.mean().combine(ee.Reducer.median(), '', True))

  def reduce_regions_function(image):
    """Applies the ee.Image.reduceRegion() method.

    Args:
      img:
        An ee.Image to reduce to a statistic by region.

    Returns:
      An ee.Feature that contains properties representing the image region
      reduction results per band and the image timestamp formatted as
      milliseconds from Unix epoch (included to enable time series plotting).
    """

    fc = image.reduceRegions(
        collection=reduction_regions,
        reducer=reducerAll,
        scale=scale,
        crs = crs,
        tileScale = tileScale
        ).set({'millis': image.date().millis()})
        
    filtered = fc.filter(ee.Filter.notNull(fc.first().propertyNames()))

    return filtered
    
  return reduce_regions_function

  # Define a function to transfer feature properties to a dictionary.
def fc_to_dict(fc):
    """Creates multiple regions reduction function.

        Creates multiple regions reduction function intended to be used as the input function
        to ee.ImageCollection.map() for reducing pixels intersecting a provided region
        to a statistic for each image in a collection. See ee.Image.reduceRegions()
        documentation for more details.

        Args:
            collection:
            FeatureCollection - The features to reduce over
            reducer:
            Optional; An ee.Reducer that defines the reduction method.
            scale:
            Optional; A number that defines the nominal scale in meters of the
            projection to work in.
            crs:
            Optional; An ee.Projection or EPSG string ('EPSG:5070') that defines
            the projection to work in.
            tileScale:
            Optional; A number representing the scaling factor used to reduce
            aggregation tile size; using a larger tileScale (e.g. 2 or 4) may enable
            computations that run out of memory with the default.

        Returns:
            A function that accepts an ee.Image and reduces it by region, according to
            the provided arguments.
    """
    prop_names = fc.first().propertyNames()
    prop_lists = fc.reduceColumns(
        reducer=ee.Reducer.toList().repeat(prop_names.size()),
        selectors=prop_names).get('list')
    return ee.Dictionary.fromLists(prop_names, prop_lists)


# Function to add date variables to DataFrame.
def add_date_info(df):
    df['Timestamp'] = pd.to_datetime(df['millis'], unit='ms')
    df['Year'] = pd.DatetimeIndex(df['Timestamp']).year
    df['Month'] = pd.DatetimeIndex(df['Timestamp']).month
    df['Day'] = pd.DatetimeIndex(df['Timestamp']).day
    df['DOY'] = pd.DatetimeIndex(df['Timestamp']).dayofyear
    return df

def bitwiseExtract(img, fromBit, toBit, new_name):

  """ Function to extract the values from specific bits

  param img: image or number to compute
  type img: ee.Number() or ee.Image()

  param fromBit: start bit
  type fromBit: int

  param toBit: end bit
  type toBit: int

  """
  fromBit = ee.Number(fromBit)
  toBit = ee.Number(toBit)
  new_name = ee.String(new_name)

  maskSize = ee.Number(1).add(toBit).subtract(fromBit)
  mask = ee.Number(1).leftShift(maskSize).subtract(1)
  
  return img.rename([new_name]).rightShift(fromBit).bitwiseAnd(mask)
  
def get_from_dict(a_list, a_dict):
    
    """ Get a list of Dict's values from a list object. Keys must be unique
    :param a_list: list of keys
    :type a_list: ee.List
    :param a_dict: dict to get the values for list's keys
    :type a_dict: ee.Dictionary
    :return: a list of values
    :rtype: ee.List
    """
    a_list = ee.List(a_list) if isinstance(a_list, list) else a_list
    a_dict = ee.Dictionary(a_dict) if isinstance(a_dict, dict) else a_dict

    empty = ee.List([])

    def wrap(el, first):
        f = ee.List(first)
        cond = a_dict.contains(el)
        return ee.Algorithms.If(cond, f.add(a_dict.get(el)), f)

    values = ee.List(a_list.iterate(wrap, empty))
    return values

# Applies scaling factors.
def applyScaleFactors(image):
    opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2).toFloat()
    thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0).toFloat()
    return (image.addBands(opticalBands, overwrite = True)
                .addBands(thermalBands, overwrite = True)
                )
def modis_scale_factor(img):
    bands = img.select('sur_ref1_b.').multiply(0.0001)
    return img.addBands(bands, overwrite = True)

def modis43A_scale_factor(img):
    bands = img.select('Nadir_Reflectance_Band.').multiply(0.0001)
    return img.addBands(bands, overwrite = True)

#6.1 Write harmonization between sensors function 
def harmonizationRoy_fromETM_OLI(img):
    '''
    harmonize TM (landsat 5) to OLI (Landsat 8)
    Slope and intercept from: Roy, D.P., Kovalskyy, V., Zhang, H.K., Vermote, E.F., Yan, L., Kumar, S.S, Egorov, A., 2016, 
    Characterization of Landsat-7 to Landsat-8 reflective wavelength and normalized difference vegetation index continuity, 
    Remote Sensing of Environment, 185, 57-70.(http:#dx.doi.Org/10.1016/j.rse.2015.12.024); Table 2 - reduced major axis (RMA) regression coefficients

    '''
    slopes = ee.Image.constant([0.9785, 0.9542, 0.9825, 1.0073, 1.0171, 0.9949])  # RMA - create an image of slopes per band for L7 TO L8 regression line - David Roy
    constants = ee.Image.constant([-0.0095, -0.0016, -0.0022, -0.0021, -0.0030, 0.0029])  # RMA - create an image of y-intercepts per band for L7 TO L8 regression line - David Roy
    y = img.select(['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7']).multiply(slopes).add(constants)  # select TM bands 1-5,7 and apply slope and intercept to band values
    
    return img.addBands(y, overwrite = True)

def harmonizationRoy_fromETMplus_OLI(img):
    '''
    harmonize ETM+ (Landsat 7) to OLI (Landsat 8)
    Slope and intercept from: Roy, D.P., Kovalskyy, V., Zhang, H.K., Vermote, E.F., Yan, L., Kumar, S.S, Egorov, A., 2016, 
    Characterization of Landsat-7 to Landsat-8 reflective wavelength and normalized difference vegetation index continuity, 
    Remote Sensing of Environment, 185, 57-70.(http:#dx.doi.Org/10.1016/j.rse.2015.12.024); Table 2 - reduced major axis (RMA) regression coefficients

    '''
    slopes = ee.Image.constant([0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071]) 
    constants = ee.Image.constant([0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172]) 
    y = img.select(['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7']).multiply(slopes).add(constants)  
    
    return img.addBands(y, overwrite = True)
  

def image_mask(from_image, mask_parameter = None):

    """
    create mask from ee.Image
    USE CASE
    .. code:: python
        # use in a collection
        col = ee.ImageCollection(ID)
        masked = col.map(mask.mask_from_image("B4", "B3"))
        # use in a single image
        img = ee.Image(ID)
        masked = img.mask.mask_from_image()
    
    
    RETURN: The function to apply a map() over a collection
    :rtype: function
    """

    if not isinstance(from_image, ee.Image): 
        raise ValueError('from_iamge has to be ee.Image')

    if not mask_parameter: 
        raise ValueError('mask_parameter has to be list of pixel values e.g., [10, 20]')

    def apply_mask(image):
        
        mask_value_list = ee.List(mask_parameter)

        def compute(element):
            element = ee.Number(element)
            mask = from_image.neq(element)
            return ee.Image(mask)
    
        landcover_mask_list = ee.List(mask_value_list.map(compute))
        
        def compute2(img, previous):
            previous = ee.Image(previous)
            return previous.And(img)
        
        landcover_mask = ee.Image(landcover_mask_list.iterate(compute2, landcover_mask_list.get(0)))

        final = image.updateMask(landcover_mask).copyProperties(image, ['system:time_start'])
        return final
    
    return apply_mask

def apply_modis_lc_mask(mask_parameter = [11, 12, 13, 14, 15, 17], lc_type = "LC_Type1"):

    """
    Create Function to use it for Map(). Function creates mask from MODIS yearly landcover (2001-2020)
    and applies to each image in image collection with respect to image-year. If image year is before 2001, the landcover
    map for 2001 is be used. Similarly, if image year is after 2020, the 2020 landcover map is used.  

    More: https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MCD12Q1?hl=en#description

    The MCD12Q1 V6 product provides global land cover types at yearly intervals (2001-2016) derived 
    from six different classification schemes. It is derived using supervised classifications of 
    MODIS Terra and Aqua reflectance data. The supervised classifications then undergo additional 
    post-processing that incorporate prior knowledge and ancillary information to further refine specific classes.

    USE CASE
    .. code:: python
        # use in a collection
        col = ee.ImageCollection(ID)
        masked = col.map(mask.apply_modis_lc_mask())
        # use in a single image
        img = ee.Image(ID)
        masked = img.mask.apply_modis_lc_mask()
    
    
    Parameters:
    -----------------------------------------------------------------------

    mask_parameter: 
        (Optional) a list of pixel values to mask out.
    type: list

    lc_type:
        (Optional) a landcover type
    type: str


    RETURN: The function to apply a map() over a collection
    :rtype: function
    """

    # if not isinstance(from_image, ee.Image): 
    #     raise ValueError('from_iamge has to be ee.Image')

    if not mask_parameter: 
        raise ValueError('mask_parameter has to be list of pixel values e.g., [10, 20]')

    def apply_mask(image):
        
        mask_value_list = ee.List(mask_parameter)

        year = ee.Date(image.date().format('YYYY'))

        if year is ee.Date("2021-01-01"):
            year = ee.Date('2020-01-01')
        elif year.millis().lt(ee.Date('2001-01-01').millis()):
            year = ee.Date('2001-01-01')
        
        mask_image = (ee.ImageCollection('MODIS/006/MCD12Q1')
                                    .filter(ee.Filter.date(year))
                                    .select(lc_type)).first()

        def compute(element):
            element = ee.Number(element)
            mask = mask_image.neq(element)
            return ee.Image(mask)
    
        landcover_mask_list = ee.List(mask_value_list.map(compute))
        
        def compute2(img, previous):
            previous = ee.Image(previous)
            return previous.And(img)
        
        landcover_mask = ee.Image(landcover_mask_list.iterate(compute2, landcover_mask_list.get(0)))

        final = image.updateMask(landcover_mask).copyProperties(image, ['system:time_start'])
        
        return final
    
    return apply_mask

def landsat578_cloud(masks = ['dilutedCloud', 'cirrus', 'cloud', 'shadow']):

    """ Function to mask out clouds and shadow in Landsat 5, 7 ,8 TOA C2:
    'LANDSAT/LC08/C02/T1_L2', 'LANDSAT/LE07/C02/T1_L2', 'LANDSAT/LT05/C02/T1_L2'

    param masks: list of mask to compute
    type masks: list
    RETURN: The function to apply a map() over a collection
    rtype: function
    
    """

    def apply_mask(image): 

        options = ee.List(masks)
        
        qa = image.select('QA_PIXEL')

        dilutedCloud = bitwiseExtract(qa, 1, 1, 'dilutedCloud').eq(0) # dilated Cloud
        cirrus = bitwiseExtract(qa, 2, 2, 'cirrus').eq(0) # Cirrus
        cloud = bitwiseExtract(qa, 3, 3, 'cloud').eq(0) # cloud
        cloudShadow = bitwiseExtract(qa, 4, 4, 'cloudShadow').eq(0) # cloud shadow

        # Bits 1 is diluted cloud; 3 cloud, and 5 cloud shadow
        # dilutedCloud = (1 << 1)
        # cirrus = (1 << 2)
        # cloud = (1 << 3)
        # cloudShadow = (1 << 4)

        # # Flags should be set to zero, indicating clear conditions.
        # bad_pixel = (qa.bitwiseAnd(dilutedCloud).eq(0)
        #                     .Or(qa.bitwiseAnd(cirrus).eq(0)) 
        #                     .Or(qa.bitwiseAnd(cloud).eq(0)) 
        #                     .Or(qa.bitwiseAnd(cloudShadow).eq(0)))

        # mask_image = bad_pixel.Not()

        cloud_dict = ee.Dictionary({
            'dilutedCloud' : dilutedCloud,
            'cirrus' : cirrus,
            'cloud' : cloud,
            'cloudShadow' : cloudShadow
        })
        
        masks_list = get_from_dict(options, cloud_dict)

        def compute(img, previous):
            previous = ee.Image(previous)
            return previous.And(img)

        # all pixel with values 0 becomes invalid
        mask_image = ee.Image(masks_list.iterate(compute, masks_list.get(0))) 

        return image.updateMask(mask_image).copyProperties(image, ['system:time_start'])
    
    return apply_mask

def modis43A_cloud(masks = [
        'BRDF_Albedo_Band_Mandatory_Quality_Band1',
        'BRDF_Albedo_Band_Mandatory_Quality_Band2',
        'BRDF_Albedo_Band_Mandatory_Quality_Band3',
        'BRDF_Albedo_Band_Mandatory_Quality_Band4',
        'BRDF_Albedo_Band_Mandatory_Quality_Band5',
        'BRDF_Albedo_Band_Mandatory_Quality_Band6',
        'BRDF_Albedo_Band_Mandatory_Quality_Band7'
    ]):

    """ Function to mask out clouds and shadow in Landsat 5, 7 ,8 TOA C2:
    'LANDSAT/LC08/C02/T1_L2', 'LANDSAT/LE07/C02/T1_L2', 'LANDSAT/LT05/C02/T1_L2'

    param masks: list of mask to compute
    type masks: list
    RETURN: The function to apply a map() over a collection
    rtype: function
    
    """
    masks = [
        'BRDF_Albedo_Band_Mandatory_Quality_Band1',
        'BRDF_Albedo_Band_Mandatory_Quality_Band2',
        'BRDF_Albedo_Band_Mandatory_Quality_Band3',
        'BRDF_Albedo_Band_Mandatory_Quality_Band4',
        'BRDF_Albedo_Band_Mandatory_Quality_Band5',
        'BRDF_Albedo_Band_Mandatory_Quality_Band6',
        'BRDF_Albedo_Band_Mandatory_Quality_Band7'
    ]

    def apply_mask(image): 

        options = ee.List(masks)
        
        band1 = image.select('BRDF_Albedo_Band_Mandatory_Quality_Band1').eq(0)
        band2 = image.select('BRDF_Albedo_Band_Mandatory_Quality_Band2').eq(0)
        band3 = image.select('BRDF_Albedo_Band_Mandatory_Quality_Band3').eq(0)
        band4 = image.select('BRDF_Albedo_Band_Mandatory_Quality_Band4').eq(0)
        band5 = image.select('BRDF_Albedo_Band_Mandatory_Quality_Band5').eq(0)
        band6 = image.select('BRDF_Albedo_Band_Mandatory_Quality_Band6').eq(0)
        band7 = image.select('BRDF_Albedo_Band_Mandatory_Quality_Band7').eq(0)

        cloud_dict = ee.Dictionary({
            'BRDF_Albedo_Band_Mandatory_Quality_Band1': band1,
            'BRDF_Albedo_Band_Mandatory_Quality_Band2': band2,
            'BRDF_Albedo_Band_Mandatory_Quality_Band3': band3,
            'BRDF_Albedo_Band_Mandatory_Quality_Band4': band4,
            'BRDF_Albedo_Band_Mandatory_Quality_Band5': band5,
            'BRDF_Albedo_Band_Mandatory_Quality_Band6': band6,
            'BRDF_Albedo_Band_Mandatory_Quality_Band7': band7
        })
        
        masks_list = get_from_dict(options, cloud_dict)

        def compute(img, previous):
            previous = ee.Image(previous)
            return previous.And(img)

        # all pixel with values 0 becomes invalid
        mask_image = ee.Image(masks_list.iterate(compute, masks_list.get(0))) 

        return image.updateMask(mask_image).copyProperties(image, ['system:time_start'])
    
    return apply_mask