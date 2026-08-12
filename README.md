# Himalayan Flood Monitoring & Early-Warning Dashboard

A Google Earth Engine-based dashboard for monitoring and assessing flood conditions across selected river systems in northern Pakistan.

## Study Rivers

The dashboard currently includes:

- Indus River — Skardu/Gilgit, Gilgit-Baltistan
- Hunza River — Gilgit-Baltistan
- Shigar River — Gilgit-Baltistan
- Swat River — Khyber Pakhtunkhwa
- Kabul River — Nowshera, Khyber Pakhtunkhwa
- Chitral River — Khyber Pakhtunkhwa

## Main Components

### 1. Historical Flood Detection

Historical flood events are detected using Sentinel-1 SAR VH backscatter change.

The workflow combines:

- Sentinel-1 VH backscatter
- Permanent-water masking
- Slope filtering
- River detection

The objective is to identify newly inundated areas while excluding the existing river channel.

### 2. Live Flood Conditions

The live mode compares recent Sentinel-1 observations with a longer baseline.

It also reports the number of Sentinel-1 observations available during the analysis period. This is important because a lack of satellite observations should not automatically be interpreted as an absence of flooding.

### 3. Flood-Risk Forecast

The forecast mode combines NOAA GFS forecast precipitation with terrain and river-related factors.

The risk index uses:

- Forecast precipitation
- HAND (Height Above Nearest Drainage)
- Slope
- Distance from river

The resulting layer represents a rule-based flood-risk susceptibility/early-warning indicator.

It does not directly model river discharge or water level.

## Data Sources

- Sentinel-1 GRD — SAR backscatter
- Sentinel-2 Surface Reflectance — water/river detection
- JRC Global Surface Water — permanent water
- Copernicus DEM GLO-30 — elevation and terrain
- MERIT Hydro — HAND
- NOAA GFS 0.25° — forecast precipitation

