///////////////////////////////////////////////////////////////////////////

var ZOOM_TARGETS = {
  'Indus (Skardu-Gilgit, GB)':    {lon: 75.2,  lat: 35.6,  zoom: 10},
  'Hunza River (GB)':             {lon: 74.5,  lat: 36.1,  zoom: 11},
  'Shigar River (GB)':            {lon: 75.72, lat: 35.38, zoom: 11},
  'Swat River (KPK)':             {lon: 72.36, lat: 34.85, zoom: 11},
  'Kabul River (Nowshera, KPK)':  {lon: 71.95, lat: 34.00, zoom: 11},
  'Chitral River (KPK)':          {lon: 71.78, lat: 35.80, zoom: 11}
};
var riverNames = Object.keys(ZOOM_TARGETS);

var DB_DROP_THRESHOLD = -3;

///////////////////////////////////////////////////////////////////////////
// RIVER AUTO DETECTION
///////////////////////////////////////////////////////////////////////////
function detectRiverMask(domain) {
  var gsw = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('max_extent');
  var permanentWater = gsw.eq(1);

  var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(domain)
    .filterDate('2023-01-01', '2025-12-31')
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30))
    .median();
  var ndwi = s2.normalizedDifference(['B3', 'B8']);
  var recentWater = ndwi.gt(0);

  var riverMask = permanentWater.or(recentWater).unmask(0).clip(domain).rename('river_mask');
  return riverMask;
}

///////////////////////////////////////////////////////////////////////////
// TERRAIN PREDICTORS
///////////////////////////////////////////////////////////////////////////
function getPredictors(domain, riverMask) {
  var dem = ee.ImageCollection('COPERNICUS/DEM/GLO30')
    .filterBounds(domain)
    .select('DEM')
    .mosaic()
    .clip(domain);

  // Fix: mosaicked DEM tiles can carry mixed projections; without an
  // explicit reproject, Terrain.slope() silently returns null everywhere.
  dem = dem.reproject({crs: 'EPSG:4326', scale: 30});

  var slope = ee.Terrain.slope(dem).rename('slope');

  var hand = ee.Image('MERIT/Hydro/v1_0_1').select('hnd').clip(domain).rename('hand');

  var riverBinary = ee.Image(riverMask).unmask(0).gt(0);

  // Fix: use each pixel's real area instead of assuming a fixed 30 m grid.
  var distance = riverBinary.fastDistanceTransform(256).sqrt()
    .multiply(ee.Image.pixelArea().sqrt())
    .clip(domain)
    .rename('distance_from_river');

  return {slope: slope, hand: hand, distance: distance};
}

///////////////////////////////////////////////////////////////////////////
// SENTINEL-1 VH
///////////////////////////////////////////////////////////////////////////
function getS1Collection(domain, start, end) {
  return ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(domain)
    .filterDate(start, end)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'));
}

function getS1Vh(domain, start, end) {
  return getS1Collection(domain, start, end)
    .select('VH')
    .median()
    .focal_median(50, 'circle', 'meters');
}

function getLatestS1Vh(domain, daysBack) {
  var end = ee.Date(Date.now());
  var start = end.advance(-daysBack, 'day');
  return getS1Vh(domain, start, end);
}

// Diagnostic only: how many Sentinel-1 passes actually landed in a window
// over this box. A count of 0 means the "flood" result isn't a real
// reading, it's a data gap.
function getS1ImageCount(domain, daysBack) {
  var end = ee.Date(Date.now());
  var start = end.advance(-daysBack, 'day');
  return getS1Collection(domain, start, end).size();
}

///////////////////////////////////////////////////////////////////////////
// FLOOD MASK
///////////////////////////////////////////////////////////////////////////
function computeFloodMask(domain, preStart, preEnd, duringStart, duringEnd, slopeImg, riverMask) {
  var before = getS1Vh(domain, preStart, preEnd);
  var after = getS1Vh(domain, duringStart, duringEnd);
  var change = after.subtract(before);

  var flood = change.lt(DB_DROP_THRESHOLD)
    .and(riverMask.not())
    .and(slopeImg.lt(15))
    .rename('flood');

  return flood;
}

///////////////////////////////////////////////////////////////////////////
// AREA CALCULATION
///////////////////////////////////////////////////////////////////////////
function areaKm2ForBand(maskImg, bandName, domain) {
  var areaImage = maskImg.selfMask().multiply(ee.Image.pixelArea()).rename(bandName);
  var result = areaImage.reduceRegion({
    reducer: ee.Reducer.sum(), geometry: domain, scale: 10, maxPixels: 1e10, bestEffort: true
  });
  return ee.Number(result.get(bandName)).divide(1000000);
}

///////////////////////////////////////////////////////////////////////////
// FORECAST PRECIPITATION (NOAA GFS 0.25deg, forecast out to 16 days)
///////////////////////////////////////////////////////////////////////////
// This is the only genuine forward-looking data source available inside
// Earth Engine. It gives you predicted rainfall, not river discharge or
// water level. Treat the output as an early-warning signal, not a
// guaranteed flood event.
function getForecastPrecip(domain, daysAhead) {
  var maxHour = daysAhead * 24;

  var gfs = ee.ImageCollection('NOAA/GFS0P25')
    .filterBounds(domain)
    .filter(ee.Filter.lte('forecast_hours', maxHour))
    .sort('creation_time', false);

  var latestRunTime = ee.Image(gfs.first()).get('creation_time');
  var latestRun = gfs.filter(ee.Filter.eq('creation_time', latestRunTime));

  // The GFS0P25 collection no longer exposes 'total_precipitation_surface'.
  // It exposes 'precipitation_rate' instead, an instantaneous rate in
  // kg/m^2/s (roughly mm/s), not an accumulated total. Average the rate
  // across the forecast window, then scale it up by the window length in
  // seconds to get an approximate accumulated depth in mm.
  var avgRate = latestRun.select('precipitation_rate').mean();
  var totalPrecip = avgRate
    .multiply(daysAhead * 24 * 3600)
    .clip(domain)
    .rename('forecast_precip_mm');

  return totalPrecip;
}

///////////////////////////////////////////////////////////////////////////
// FLOOD RISK FORECAST (rule-based susceptibility index)
///////////////////////////////////////////////////////////////////////////
// Risk = weighted mix of forecast rainfall + how low, flat, and close to
// the river the land is. Weights are a starting point. Adjust them once
// you compare a few forecasts against what actually happened.
function computeFloodRisk(domain, riverMask, predictors, daysAhead) {
  var precip = getForecastPrecip(domain, daysAhead);

  // Ceiling set to 60mm, not 150mm. In steep Himalayan catchments, 25-40mm
  // in a short window can trigger flash floods -- runoff funnels straight
  // into the channel instead of spreading out. A 150mm ceiling was tuned
  // for monsoon-plains totals and made the model nearly blind to the rain
  // levels that actually matter here.
  var precipNorm = precip.unitScale(0, 60).clamp(0, 1);
  var handNorm = predictors.hand.unitScale(0, 30).clamp(0, 1).multiply(-1).add(1);
  var slopeNorm = predictors.slope.unitScale(0, 30).clamp(0, 1).multiply(-1).add(1);
  var distNorm = predictors.distance.unitScale(0, 500).clamp(0, 1).multiply(-1).add(1);

  var risk = precipNorm.multiply(0.40)
    .add(handNorm.multiply(0.25))
    .add(slopeNorm.multiply(0.15))
    .add(distNorm.multiply(0.20))
    .rename('flood_risk');

  risk = risk.updateMask(riverMask.not());

  return {risk: risk, precip: precip};
}

///////////////////////////////////////////////////////////////////////////
// DRAWING TOOLS (single definition)
///////////////////////////////////////////////////////////////////////////
var drawingTools = Map.drawingTools();
drawingTools.setShown(false);
while (drawingTools.layers().length() > 0) {
  drawingTools.layers().remove(drawingTools.layers().get(0));
}
drawingTools.layers().add(ui.Map.GeometryLayer({geometries: null, name: 'river_box', color: '23cba7'}));

function clearGeometry() {
  var layer = drawingTools.layers().get(0);
  var geometries = layer.geometries();
  if (geometries.length() > 0) { geometries.remove(geometries.get(0)); }
}
function drawRiverBox() { clearGeometry(); drawingTools.setShape('rectangle'); drawingTools.draw(); }
function hasDrawnBox() { return drawingTools.layers().get(0).geometries().length() > 0; }
function getDrawnBox() { return drawingTools.layers().get(0).getEeObject(); }

///////////////////////////////////////////////////////////////////////////
// DATE CHECK
///////////////////////////////////////////////////////////////////////////
function validDate(date) {
  return typeof date === 'string' && /^\d{4}-\d{2}-\d{2}$/.test(date.trim());
}

///////////////////////////////////////////////////////////////////////////
// RESULTS PANEL (single definition)
///////////////////////////////////////////////////////////////////////////
var resultsPanel = ui.Panel({
  style: {width: '360px', maxHeight: '450px', position: 'bottom-right', shown: false}
});
Map.add(resultsPanel);

///////////////////////////////////////////////////////////////////////////
// MAIN FLOOD ANALYSIS (single definition -- this is the only version now)
///////////////////////////////////////////////////////////////////////////
function runFloodAnalysis() {
  resultsPanel.style().set('shown', true);

  if (!hasDrawnBox()) {
    resultsPanel.widgets().reset([ui.Label('Draw river box first', {color: 'red'})]);
    return;
  }

  var domain = getDrawnBox();
  var preStart = filters.preStart.getValue().trim();
  var preEnd = filters.preEnd.getValue().trim();
  var floodStart = filters.duringStart.getValue().trim();
  var floodEnd = filters.duringEnd.getValue().trim();

  if (!validDate(preStart) || !validDate(preEnd) || !validDate(floodStart) || !validDate(floodEnd)) {
    resultsPanel.widgets().reset([ui.Label('Date format must be YYYY-MM-DD', {color: 'red'})]);
    return;
  }

  resultsPanel.widgets().reset([
    ui.Label('Processing Sentinel-1 + River + HAND + Slope...', {color: 'gray'})
  ]);

  var riverMask = detectRiverMask(domain);
  var predictors = getPredictors(domain, riverMask);
  var floodMask = computeFloodMask(domain, preStart, preEnd, floodStart, floodEnd, predictors.slope, riverMask);
  var riverArea = areaKm2ForBand(riverMask, 'river_mask', domain);
  var floodArea = areaKm2ForBand(floodMask, 'flood', domain);

  Map.layers().reset();
  Map.addLayer(riverMask.selfMask(), {palette: PAL.river}, 'Detected river');
  Map.addLayer(predictors.slope, {min: 0, max: 60, palette: PAL.slope}, 'Slope', false);
  Map.addLayer(predictors.hand, {min: 0, max: 100, palette: PAL.hand}, 'HAND', false);
  Map.addLayer(predictors.distance, {min: 0, max: 1000, palette: PAL.distance}, 'Distance from river', false);
  Map.addLayer(floodMask.selfMask(), {palette: PAL.flood}, 'NEW FLOOD');

  riverArea.evaluate(function(r) {
    floodArea.evaluate(function(f) {
      resultsPanel.widgets().reset([
        ui.Label('LIVE FLOOD RESULTS', {fontWeight: 'bold', fontSize: '18px'}),
        ui.Label('Pre flood: ' + preStart + ' to ' + preEnd),
        ui.Label('Flood period: ' + floodStart + ' to ' + floodEnd),
        ui.Label('Detected permanent river: ' + (r ? r.toFixed(3) : '0') + ' km\u00b2', {color: 'blue'}),
        ui.Label('NEW flooded area: ' + (f ? f.toFixed(2) : '0') + ' km\u00b2', {color: 'cyan', fontWeight: 'bold'}),
        ui.Label('Method: Sentinel-1 VH change + HAND + slope + river exclusion',
          {fontSize: '11px', color: 'gray'})
      ]);
    });
  });
}

drawingTools.onDraw(ui.util.debounce(runFloodAnalysis, 500));
drawingTools.onEdit(ui.util.debounce(runFloodAnalysis, 500));

///////////////////////////////////////////////////////////////////////////
// LIVE SCAN (last-N-days vs last-60-days baseline) -- wired to a button
///////////////////////////////////////////////////////////////////////////
function liveFloodScan() {
  if (!hasDrawnBox()) {
    resultsPanel.style().set('shown', true);
    resultsPanel.widgets().reset([ui.Label('Draw river box first', {color: 'red'})]);
    return;
  }

  resultsPanel.style().set('shown', true);
  resultsPanel.widgets().reset([ui.Label('Scanning last 12 days vs 60-day baseline...', {color: 'gray'})]);

  var domain = getDrawnBox();
  var river = detectRiverMask(domain);
  var predictors = getPredictors(domain, river);

  var old = getLatestS1Vh(domain, 60);
  var now = getLatestS1Vh(domain, 12);
  var change = now.subtract(old);

  var flood = change.lt(-2.5)
    .and(river.not())
    .and(predictors.slope.lt(25))
    .rename('LIVE_FLOOD');

  Map.layers().reset();
  Map.addLayer(river.selfMask(), {palette: PAL.river}, 'Detected river');
  Map.addLayer(flood.selfMask(), {palette: PAL.flood}, 'LIVE FLOOD (last 12 days)');

  var area = areaKm2ForBand(flood, 'LIVE_FLOOD', domain);
  var oldCount = getS1ImageCount(domain, 60);
  var nowCount = getS1ImageCount(domain, 12);

  area.evaluate(function(a) {
    oldCount.evaluate(function(oc) {
      nowCount.evaluate(function(nc) {
        var warning = '';
        if (nc === 0) {
          warning = 'No Sentinel-1 passes landed on this box in the last 12 days. ' +
            'The 0.00 km\u00b2 reading is a data gap, not a confirmed all-clear. Try again in a few days.';
        } else if (oc === 0) {
          warning = 'No Sentinel-1 passes landed in the 60-day baseline window. ' +
            'The comparison has nothing reliable to measure against.';
        } else if (nc === 1) {
          warning = 'Only 1 pass in the last 12 days. A single image can miss partial ' +
            'cloud-free coverage of the box; treat this reading as provisional.';
        }

        var widgets = [
          ui.Label('LIVE CONDITIONS SCAN', {fontWeight: 'bold', fontSize: '18px'}),
          ui.Label('Comparing last 12 days against a 60-day baseline.', {fontSize: '11px', color: '555'}),
          ui.Label('Current flood signal: ' + (a ? a.toFixed(2) : '0') + ' km\u00b2',
            {color: 'cyan', fontWeight: 'bold'}),
          ui.Label('Sentinel-1 passes used: ' + nc + ' (last 12 days), ' + oc + ' (60-day baseline)',
            {fontSize: '11px', color: '333'}),
          ui.Label('Looser thresholds than the historical mode (-2.5 dB, slope < 25\u00b0) -- ' +
            'meant to catch developing conditions, not confirm a specific event.',
            {fontSize: '10px', color: '888'})
        ];

        if (warning) {
          widgets.push(ui.Label(warning, {fontSize: '11px', color: 'red', fontWeight: 'bold'}));
        }

        resultsPanel.widgets().reset(widgets);
      });
    });
  });
}

///////////////////////////////////////////////////////////////////////////
// CONSOLE CHARTS -- diagnostic graphs sent to the Console tab
///////////////////////////////////////////////////////////////////////////
function showConsoleCharts() {
  if (!hasDrawnBox()) {
    resultsPanel.style().set('shown', true);
    resultsPanel.widgets().reset([ui.Label('Draw river box first', {color: 'red'})]);
    return;
  }

  var domain = getDrawnBox();

  // Chart 1: Sentinel-1 VH backscatter over the last 60 days.
  // A real flood shows up as a sustained dip. A single missing dot or a
  // one-day blip is more likely a noisy pass than an actual event.
  var end = ee.Date(Date.now());
  var start = end.advance(-60, 'day');
  var s1Series = getS1Collection(domain, start, end).select('VH');

  var backscatterChart = ui.Chart.image.series({
    imageCollection: s1Series,
    region: domain,
    reducer: ee.Reducer.mean(),
    scale: 30,
    xProperty: 'system:time_start'
  }).setOptions({
    title: 'Sentinel-1 VH backscatter, last 60 days (dB, lower = wetter/flooded)',
    vAxis: {title: 'VH (dB)'},
    hAxis: {title: 'Date'},
    lineWidth: 2,
    pointSize: 4
  });
  print(backscatterChart);

  // Chart 2: GFS forecast rainfall by forecast hour, latest model run.
  // Shows whether the rain arrives as one heavy burst or spreads out,
  // which a single "24.1mm average" number hides completely.
  var gfs = ee.ImageCollection('NOAA/GFS0P25')
    .filterBounds(domain)
    .filter(ee.Filter.lte('forecast_hours', 240))
    .sort('creation_time', false);
  var latestRunTime = ee.Image(gfs.first()).get('creation_time');
  var latestRun = gfs.filter(ee.Filter.eq('creation_time', latestRunTime))
    .select('precipitation_rate');

  var rainChart = ui.Chart.image.series({
    imageCollection: latestRun,
    region: domain,
    reducer: ee.Reducer.mean(),
    scale: 25000,
    xProperty: 'forecast_hours'
  }).setOptions({
    title: 'GFS forecast precipitation rate by hour ahead (kg/m\u00b2/s)',
    vAxis: {title: 'Rate'},
    hAxis: {title: 'Hours from now'},
    lineWidth: 2,
    pointSize: 4
  });
  print(rainChart);

  resultsPanel.style().set('shown', true);
  resultsPanel.widgets().reset([
    ui.Label('Charts sent to the Console tab.', {fontWeight: 'bold'}),
    ui.Label('Open Console (top right of the code editor) to view them.', {fontSize: '11px', color: '555'})
  ]);
}

///////////////////////////////////////////////////////////////////////////
// FORECAST SCAN -- forecast rainfall x terrain susceptibility
///////////////////////////////////////////////////////////////////////////
function runFloodForecast() {
  if (!hasDrawnBox()) {
    resultsPanel.style().set('shown', true);
    resultsPanel.widgets().reset([ui.Label('Draw river box first', {color: 'red'})]);
    return;
  }

  resultsPanel.style().set('shown', true);
  resultsPanel.widgets().reset([ui.Label('Pulling GFS rainfall forecast and scoring risk...', {color: 'gray'})]);

  var domain = getDrawnBox();
  var daysAhead = parseInt(filters.forecastDays.getValue(), 10);
  if (!daysAhead || daysAhead < 1) { daysAhead = 3; }
  if (daysAhead > 10) { daysAhead = 10; }

  var riverMask = detectRiverMask(domain);
  var predictors = getPredictors(domain, riverMask);
  var result = computeFloodRisk(domain, riverMask, predictors, daysAhead);

  Map.layers().reset();
  Map.addLayer(riverMask.selfMask(), {palette: PAL.river}, 'Detected river');
  Map.addLayer(result.precip, {min: 0, max: 150, palette: PAL.precip}, 'Forecast rainfall (mm)', false);
  Map.addLayer(result.risk, {min: 0, max: 1, palette: PAL.risk}, 'Flood risk (next ' + daysAhead + ' days)');

  var highRisk = result.risk.gt(0.45);
  var highRiskArea = areaKm2ForBand(highRisk, 'flood_risk', domain);
  var avgPrecipNum = ee.Number(result.precip.reduceRegion({
    reducer: ee.Reducer.mean(), geometry: domain, scale: 250, maxPixels: 1e9, bestEffort: true
  }).get('forecast_precip_mm'));

  avgPrecipNum.evaluate(function(p) {
    highRiskArea.evaluate(function(a) {
      resultsPanel.widgets().reset([
        ui.Label('FLOOD RISK FORECAST', {fontWeight: 'bold', fontSize: '18px'}),
        ui.Label('Next ' + daysAhead + ' days, from the latest NOAA GFS model run.', {fontSize: '11px', color: '555'}),
        ui.Label('Average forecast rainfall: ' + (p ? p.toFixed(1) : '0') + ' mm', {color: 'blue'}),
        ui.Label('High-risk area: ' + (a ? a.toFixed(2) : '0') + ' km\u00b2', {color: 'orange', fontWeight: 'bold'}),
        ui.Label('This scores where forecast rain overlaps low, flat, river-adjacent land. ' +
          'It does not model river discharge or water level. Treat it as an early-warning ' +
          'layer, not a confirmed prediction.',
          {fontSize: '10px', color: '888'})
      ]);
    });
  });
}

///////////////////////////////////////////////////////////////////////////
// COLORS
///////////////////////////////////////////////////////////////////////////
var PAL = {
  slope: ['#1a9850', '#fee08b', '#d73027'],
  hand: ['#2166ac', '#f7f7f7', '#8c510a'],
  distance: ['#d73027', '#ffffff'],
  river: ['#0000FF'],
  flood: ['#00FFFF'],
  precip: ['#ffffff', '#4292c6', '#08306b'],
  risk: ['#2ca25f', '#fee08b', '#d73027']
};

///////////////////////////////////////////////////////////////////////////
// ZOOM
///////////////////////////////////////////////////////////////////////////
function zoomToRiver() {
  var target = ZOOM_TARGETS[filters.riverSelect.getValue()];
  Map.setCenter(target.lon, target.lat, target.zoom);
}

///////////////////////////////////////////////////////////////////////////
// FILTERS / UI CONTROLS
///////////////////////////////////////////////////////////////////////////
var filters = {
  riverSelect: ui.Select({items: riverNames, value: 'Swat River (KPK)', onChange: zoomToRiver}),
  preStart: ui.Textbox('YYYY-MM-DD', '2022-07-01'),
  preEnd: ui.Textbox('YYYY-MM-DD', '2022-07-15'),
  duringStart: ui.Textbox('YYYY-MM-DD', '2022-08-25'),
  duringEnd: ui.Textbox('YYYY-MM-DD', '2022-09-05'),
  forecastDays: ui.Textbox('Days ahead (1-10)', '3'),
  applyButton: ui.Button({label: 'RUN FLOOD DETECTION', onClick: runFloodAnalysis, style: {stretch: 'horizontal'}}),
  liveButton: ui.Button({label: 'CHECK CURRENT CONDITIONS (live)', onClick: liveFloodScan, style: {stretch: 'horizontal'}}),
  forecastButton: ui.Button({label: 'FORECAST FLOOD RISK', onClick: runFloodForecast, style: {stretch: 'horizontal'}}),
  chartsButton: ui.Button({label: 'SHOW CHARTS (console)', onClick: showConsoleCharts, style: {stretch: 'horizontal'}})
};

var controlPanel = ui.Panel({
  widgets: [
    ui.Label('Himalayan Live Flood Dashboard', {fontWeight: 'bold', fontSize: '20px'}),
    ui.Label('Sentinel-1 + River detection + HAND + Slope', {fontSize: '12px', color: '555'}),

    ui.Label('1. Select river'),
    filters.riverSelect,

    ui.Label('2. Draw box over river and floodplain', {fontWeight: 'bold'}),
    ui.Button({label: 'Draw River Box', onClick: drawRiverBox, style: {stretch: 'horizontal'}}),

    ui.Label('3a. Historical event mode', {fontWeight: 'bold'}),
    ui.Label('Pre flood start'), filters.preStart,
    ui.Label('Pre flood end'), filters.preEnd,
    ui.Label('Flood start'), filters.duringStart,
    ui.Label('Flood end'), filters.duringEnd,
    filters.applyButton,

    ui.Label('3b. Live mode (no dates needed)', {fontWeight: 'bold', margin: '10px 8px 4px 8px'}),
    filters.liveButton,

    ui.Label('3c. Forecast mode (GFS rainfall + terrain)', {fontWeight: 'bold', margin: '10px 8px 4px 8px'}),
    filters.forecastDays,
    filters.forecastButton,
    filters.chartsButton,

    ui.Label('Legend', {fontWeight: 'bold'}),
    ui.Label('Blue = detected permanent river', {color: '0000FF'}),
    ui.Label('Cyan = new flood water', {color: '00AAAA'}),
    ui.Label('Green-Yellow-Red = slope / flood risk')
  ],
  style: {position: 'bottom-left', width: '330px'}
});

ui.root.insert(0, controlPanel);

///////////////////////////////////////////////////////////////////////////
// START
///////////////////////////////////////////////////////////////////////////
Map.setOptions('SATELLITE');
zoomToRiver();
