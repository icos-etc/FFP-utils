# FFP-utils
Functions to manipulate data from the ICOS Carbon Portal, perform a Flux Footprint Prediction (FFP) using the method developed by Kljun et al. (2015), and save the outputs in different formats. 

## CPtoFFPinput
Starting from an input zip file downloaded from the ICOS Carbon Portal, the function returns a dataframe ready to apply the Flux Footprint Prediction (FFP) method developed by Kljun et al. (2015).

| Argument         | Description                                                                             |
| ---------------- | ----------------------------------------------------------------------------------------|
| zip.file	        | Full path of the downloaded ICOS level 2 compressed archive data (either L2 or interim) |
| FFP.input.table	 | As alternative of the zip file, full path of an already existent input CSV table        |
| start.date	      | Start date of measurement needed, in the format YYYY-MM-DD                              |
| end.date	        | End date of measurement needed, in the format YYYY-MM-DD                                |
| save.input.table	| If TRUE, the input table for the FFP model is saved as CSV                              |
| return.input	    | If TRUE, the function output is returned                                                |

This function performs data processing and temporal filtering of an ICOS L2 ARCHIVE or L2 INTERIM zip file. Additionally, it can do a simplified processing and a temporal filtering of an existing FFP input table (e.g. a CSV file previously returned by the same function). 

There is no need to unzip the compressed file of each site. When multiple ZIP files are downloaded together from the Carbon Portal, the only action needed is to unzip the parent ZIP file. 
If the dates are set as NULL, data is not filtered according to a specific time interval. Alongside the FFP input data, the EC tower coordinates are extracted from the metadata contained in the ZIP file, or silently within the function if an existing CSV file is used. When multiple zip files are downloaded, the parent zip folder contains a CSV file named ‘!TOC.csv’ that is used to extract the PID of each level 2 product. 

If ‘return.input’ is set to TRUE, the 4-elements matrix is the only output directly returned by the function within R.
In the R working directory the function creates a folder named as the site ID and a sub-folder named "Input", where are stored the processed FFP inputs in a CSV format. The filename follows the structure: site id_start date_end date_FFP input.csv. 

When specified, a 4-element list is returned, containing: site ID string; vector with EC tower lat/lon coordinates; PID; FFP input data frame. Depending on data and site characteristics, not all the elements can always be returned (e.g. PID, if not available when the ‘!TOC.csv’ file doesn't exist; FFP input table when not all the parameters are available).
When the FFP input table is saved as CSV, ‘site_id’, ‘lat, ‘lon’ and ‘PID’ are added in the first 4 columns of the table.  
 
## doFFP
Starting from the input parameters (e.g. a CSV file or an R object returned by the CPtoFFPinput function), it saves in dedicated folders the daily Flux Footprints Predictions (FFP; Kljun et al., 2015) for a specific Eddy Covariance site. 

| Argument           | Description                                                                        |
| ------------------ | -----------------------------------------------------------------------------------|
| FFP.input.df	      | Full path of an already existent csv, or an existent R dataframe                   |
| EC.tower.coords    | A vector containing the EC tower coordinates in WGS84 (lat, lon)                   |
| site.ID	           | ID of the site                                                                     |
| PID                | PID of the input table                                                             |
| FFP.function.path  | Full path of the FFP function                                                      |
| do.parallel        | if TRUE, the calculation of footprint is parallelized using the furrr package      |
| max.gap	           | Maximum number of consecutive half-hours to fill with LOCF                         |
| FFP.domain	        | Domain size [m]                                                                    |
| dx	                | Cell size of domain [m]                                                            |
| FFP.R              | Density thresholds for isopleths calculation                                       |
| which.FFP.R        | Among the density thresholds, which one export in the netCDF file                  |
| drop.trivial.srcs  | If TRUE, removes point sources with a trivial contribution                         |
| save.FFP.mtx.as    | Save the FFP matrix as csv (type ‘csv’), netCDF (type ‘nc') or no save (type NULL) |
| save.plot.FFP.mtx  | If TRUE, a jpeg plot of FFP (both matrix and isopleths) is saved                   |
| return.isopleth    | If TRUE, the chosen isopleth (which.FFP.R) is saved in the netCDF file             |
| save.log           | If TRUE, log messages are written to a file                                        |
| save.ncplot	       | If TRUE, a sample plot of the nc file produced is saved                            |

For every time interval (e.g. 30 minutes) within each day, the function processes a structured input table and runs the FFP function developed by Kljun et al. (2015), saving the output matrices in csv or netCDF (default) format. Optionally, plots of the footprint matrix and a log file of the entire process are returned. 

Tower coordinates, site ID and PID can be omitted if the dataframe already contains the them (e.g. the CSV file saved by the CPtoFFPinput function), otherwise they must be specified (e.g. list returned by the CPtoFFPinput function). The function relies on the FFP model function (Kljun et al., 2015) which should be placed in the R working directory, preferably inside the subfolder ‘local_functions’ (default). Alternatively, is also possible to specify a different FFP function path. Since the FFP model run is time-consuming, is possible to parallelize the function (do.parallel=TRUE) using furrr with a low number od workers (3), to prevent system crashes. 

Gaps are filled up to a defined threshold (4 by default; max.gap argument). FFP function-related parameters, such as the spatial extent (1000 m by default; FFP.domain argument), resolution (1 m by default; dx argument), the FFP density thresholds for isopleths calculation (c(50, 70, 80, 90) by default; FFP.R argument) and which isopleth save in the netCDF file (which.FFP.R) can be specified. In addition, pixels with a minimal contribution to the total density (threshold of 90%; drop.trivial.srcs argument) are removed by default. 

In the R working directory, a folder named as the site ID with an "Output" sub-folder are created if they don't exist. Depending on the saving options, inside the "Output" folder can be created the following sub-folders: 
*	<ins>nc_files</ins>: if ‘nc’ is chosen as saving option, it contains the daily netCDF files. Each file has a set of dimensions, variables and attributes. 

<ul>

The dimensions are related to the spatial (x, y) and temporal resolution (time), to the coordinates of the EC tower (EC_Coords) and to the isopleth chosen to be saved (instance: total number of polygons; nodes: total number of vertices; nchar: maximum number of characters of each entry in the isopleth attribute table). 

| Variable (Dimensions) | Description|
| - | - |
| x (x), y (y), time (time) | Variables strongly linked to the corresponding dimensions. They place each FFP matrix in the right position in terms of space and time |
| FFP (x, y, time) | Daily array containing the 2D FFP matrices calculated every 30 minutes. Each matrix was compressed and converted to integer values, but providing scale and offset factors to allow software to retrieve original values. When the FFP model fails, a no-data grid (-9999 values) is added to keep the total number of time steps equal (48) for every file generated |
| polygon_id (instance), observation_time (n_char, instance), polygon_area (instance) | Variables strictly related to the 70% isopleth (one for each time interval), containing the ID (polygon_id), the corresponding timestamp (observation_time) and the area in m<sup>2</sup> (polygon_area) |
| quality_flag (instance) | Quality flag of FFP. 0: FFP and isopleths correctly calculated; 1: FFP and isopleths not calculated; 2: FFP correctly calculated, but not cropped because the 90% isopleth could not be retrieved; 3: FFP correctly calculated, but not cropped because multiple isopleths could not be retrieved |
| geometry_container, UTM_Coordinate_System | They define the associations between the variables linked to the geometry (the isopleth chosen to be saved) and the correct Coordinate Reference System (CRS, in WGS84/UTM) for each variable that must be spatially located. The right CRS is automatically detected starting from the WGS84 coordinates of the EC tower |
| nodes_count (instance), nodes_x (nodes), nodes_y (nodes) | Respectively, count of the total number of vertices for the isopleth, X and Y coordinates of the isopleth vertices |
| EC_tower_coordinates (EC_Coords) | Coordinates in WGS84/UTM of the eddy covariance flux tower |

The attributes provide metadata for a better description of the file, according to a set of standards (CF - Climate and Forecast; ATMODAT).

</ul>

*	<ins>Individual_footprint_matrices</ins>: if ‘csv’ is chosen as saving option, it contains the 30-min output FFP matrices saved in CSV format. Each file is named following the structure: site ID_FFP2Dmtx_YYYYMMDDHHmm.csv.
*	<ins>Individual_footprint_plots</ins>: if ‘save.plot.FFP.mtx’ is true, it contains the plots of the 30-min output FFP matrices saved in PNG format. Each file is named following the structure: site ID_FFP2Dmtx_YYYYMMDDHHmm.png.
*	<ins>log</ins>: if ‘save.log’ is true, it contains a log file for each function run. The file is named following the structure: site ID_log_start date_end date.log

## References
Kljun, N., Calanca, P., Rotach, M. W., and Schmid, H. P.: A simple two-dimensional parameterisation for Flux Footprint Prediction (FFP), Geosci. Model Dev., 8, 3695-3713, 2015. doi:10.5194/gmd-8-3695-2015


## Examples
```
# Retrieve the file path
zip_file="ZIP\\My data cart\\ICOSETC_BE-Bra_ARCHIVE_INTERIM_L2.zip"

# Generate the FFP inputs
CPtoFFPinput(zip.file=zip_file, save.input.table = T, return.input = T, start.date = NULL, end.date = '2023-12-31'))
```
```
# Get the path of each ZIP file
zip.list=list.files("ZIP\\My data cart\\", pattern = 'zip', full.names = T)

# Generate the FFP inputs for all the files until 2023
lapply(zip.list, function(x) CPtoFFPinput(zip.file=x, save.input.table = T, return.input = T, 
                                           start.date = NULL, end.date = '2023-12-31'))
```
```
# Filter dates from a specific year for an existing input table and save a new filtered table
CPtoFFPinput(FFP.input.table="IT-SR2\\Input\\IT-SR2_20191028_20231231_FFPinput.csv", save.input.table = T, return.input = F,
             start.date = '2023-01-01', end.date = '2023-12-31')
```

```
# Get the path of each ZIP file
zip.list=list.files("ZIP\\My data cart\\", pattern = 'zip', full.names = T)

# Generate the FFP inputs as R object for only one day
input.list=lapply(zip.list, function(x) CPtoFFPinput(zip.file=x, save.input.table = F, return.input = T, 
                                           start.date = '2023-11-21', end.date = '2023-11-21'))

# Process only the first site
Site_1= input.list[[1]]

# Calculate FFP only for the first site
doFFP(FFP.input.df = Site_1[['FFP_input_parameters']],
      EC.tower.coords = Site_1[['tower_coordinates']],
      site.ID=Site_1[['site_id']],
      PID=Site_1[['PID']],
      FFP.function.path="local_functions\\calc_footprint_FFP_climatology.R",
      max.gap = 4, FFP.domain = 1000, dx = 1, FFP.R = c(50,70,80,90), which.FFP.R=70,
      drop.trivial.srcs = TRUE, save.FFP.mtx.as="nc",          
      save.plot.FFP.mtx = FALSE, return.isopleth=TRUE, save.log=TRUE, save.ncplot=TRUE)
```
```
library("purrr")

# Get the path of each ZIP file
zip.list=list.files("ZIP\\My data cart\\", pattern = 'zip', full.names = T)

# Generate the input table as R object for only one day
input.list=map(zip.list, ~ CPtoFFPinput(zip.file=.x, save.input.table = F, return.input = T, 
                                         start.date = '2023-11-21', end.date = '2023-11-21'))

# Filter out sites with NULL dataframes (where not all the input parameters are available) 
input.list.filtered=discard(input.list, ~ is.null(.x[['FFP_input_parameters']]))

# Calculate FFP and save it as netCDF for all the sites. It could be very time consuming 
walk(input.list.filtered, ~ doFFP(FFP.input.df = .x[['FFP_input_parameters']], 
                                  EC.tower.coords = .x[['tower_coordinates']], 
                                  site.ID=.x[['site_id']], 
                                  PID=.x[['PID']], 
                                  FFP.function.path="local_functions\\calc_footprint_FFP_climatology.R",
                                  do.parallel=TRUE, 
                                  max.gap = 4, FFP.domain = 1000, dx = 1, FFP.R = c(50,70,80,90), which.FFP.R=70,
                                  drop.trivial.srcs = TRUE, save.FFP.mtx.as="nc", 
                                  save.plot.FFP.mtx = FALSE, return.isopleth=TRUE, save.log=TRUE, save.ncplot=TRUE))

```
```
# Processing a single daily file for site using a csv table generated using the CPtoFFPinput function and extract the 90% isopleth
doFFP(FFP.input.df = "BE-Bra\\Input\\BE-Bra_20230701_20230701_FFPinput.csv",
      FFP.function.path="local_functions\\calc_footprint_FFP_climatology.R",
      max.gap = 4, FFP.domain = 1000, dx = 1, FFP.R = c(50,70,80,90), which.FFP.R=90,
      drop.trivial.srcs = TRUE, save.FFP.mtx.as="nc",          
      save.plot.FFP.mtx = FALSE, return.isopleth=TRUE, save.log=TRUE, save.ncplot=TRUE)
```
