////////code to extract raster data at pts that have been loaded 
//to your account as assets, window is split into 3 panes, 
//with middle pane split into 2 sections, and a lower map pane
///this code gets pasted into the lower middle pane, 
//you can then save it as a script if you want, 
//where it will show up in the left pane under the scripts tab in the "Owner" file tree

//load pts data from assests, on left pane click assests tab, 
//upload new pts shape file first, once it shows up in this list,
// click its entry and select import, 
//it will move to the upper middle pane and be called "table", 
//then you can work with it in the lower middle pane
//next find the raster data you want to work with in the search box 
//at the top of the whole window ("search places and datasets ..."), 
//look for "Image Collection Id" in the resulting window for the data you want and copy it
//to use this, paste this script in the lower middle pane, 
//above the map pane, import your pt assest, copy the right image string to the line 
//below that starts with "var myraster" and click "run" at the top of the middle pane
//in the right pane click the tasks tab, an object called "mylabel" will show up, 
//click run, a window pops up asking you what and where to save the file, 
//I usually select "drive", make a folder in my google drive to store the file, 
//enter that folder name here, then select a name for the file, click run again.
//when it's finished, the result is a file in your google drive with the raster info 
//for all the pts on the dates you selected, which you can download to your computer

var pts = ee.FeatureCollection(table);

/**
 * Function to mask clouds based on the pixel_qa band of Landsat 8 SR data.
 * @param {ee.Image} image input Landsat 8 SR image
 * @return {ee.Image} cloudmasked Landsat 8 image
 */
function cloudMaskL8(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = (1 << 3);
  var cloudsBitMask = (1 << 5);
  // Get the pixel QA band.
  var qa = image.select('pixel_qa');
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                 .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask);
}

/**
 * Function to mask clouds based on the pixel_qa band of Landsat SR data, missions 4,5,7.
 * @param {ee.Image} image Input Landsat SR image
 * @return {ee.Image} Cloudmasked Landsat image
 */
var cloudMaskL457 = function(image) {
  var qa = image.select('pixel_qa');
  // If the cloud bit (5) is set and the cloud confidence (7) is high
  // or the cloud shadow bit is set (3), then it's a bad pixel.
  var cloud = qa.bitwiseAnd(1 << 5)
                  .and(qa.bitwiseAnd(1 << 7))
                  .or(qa.bitwiseAnd(1 << 3));
  // Remove edge pixels that don't occur in all bands
  var mask2 = image.mask().reduce(ee.Reducer.min());
  return image.updateMask(cloud.not()).updateMask(mask2);
};


//IMPORT raster, here we use Landsat 8, you can find the Image Collection Id for other image sources in the searchbox at the top of the screen
//for Landsat 5 use: var myraster = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR') instead
//limit the spatial location with filterBounds(geometry) where geometry is a shape you've map with the drawing tool in the lower map window
//limit the dates with filterDate
//limit cloudy scenes with the 'CLOUD_COVER' scene metadata variable
var myraster = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR') 
.filterBounds(geometry)  
.filterDate('2012-12-01','2018-06-05')
//remove images where whole image has lots of cloud cover
.filterMetadata('CLOUD_COVER' , 'less_than', 30)
//remove cloudy pixels for Landsat 4, 5, or 7
//.map(cloudMaskL457);
//remove cloudy pixels for Landsat 8
.map(cloudMaskL8);



//This code below is a function to extract points called "table" from raster called "myraster"
// Empty Collection to fill
var ft = ee.FeatureCollection(ee.List([]))


//With removal of null values ------------------------------------------
//Function to extract values from image collection based on point file and 
//export as a table 
var fill = function(img, ini) {

  // type cast
  var inift = ee.FeatureCollection(ini);

  // gets the values for the points in the current img
  var ft2 = img.reduceRegions(pts, ee.Reducer.first(),30);

  // gets the date of the img
  var date = img.date().format();

  // writes the date in each feature
  var ft3 = ft2.map(function(f){return f.set("date", date)});

  // merges the FeatureCollections

  var ft3a = ft3.filter(ee.Filter.neq('B1', null));//filter first to remove null values
  return inift.merge(ft3a);
};

// Iterates over the ImageCollection
var newft = ee.FeatureCollection(myraster.iterate(fill, ft))

// Export to a csv file called "myfile.csv" into your Google Drive folder called "Engine" and labeled the task "mylabel"
//this will cause a label called "mylabel" to appear in the task pane, click on this and select "run"
Export.table.toDrive(newft,
"mylabel",
"Engine",
"myfile")


//plot scene, centered on a location and plot sample points ------------------------------------
Map.setCenter(-81.2849, 31.4502,15);
Map.addLayer(myraster, {bands: ['B5', 'B4', 'B3'], min: 0, max: 2000}, 'false-color composite');
Map.addLayer(pts);


