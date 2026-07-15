# FFP-utils
Functions to manipulate data from the ICOS Carbon Portal, perform a Flux Footprint Prediction (FFP) using the method developed by Kljun et al. (2015), and save the outputs in several formats temporal resolutions.

## CPtoFFPinput

| Argument         | Description                                                                             |
| ---------------- | ----------------------------------------------------------------------------------------|
| zip.file	       | Full path of the downloaded ICOS level 2 compressed archive data (either L2 or interim) |
| FFP.input.table	 | As alternative of the zip file, full path of an already existent input CSV table        |
| start.date	 | Start date of measurement needed, in the format YYYY-MM-DD                              |
| end.date	       | End date of measurement needed, in the format YYYY-MM-DD                                |
| save.input.table | If TRUE, the input table for the FFP model is saved as CSV                              |
| MDS.file.list    | File list (full path) of MDS gap-filling indices                                        |
| ERA5.file.list   | File list (full path) of yearly ERA5 reanalysis data                                    |
| out.dir	       | Directory where save the output                                                         | 

Starting from the EC data available on the ICOS Carbon Portal, CPtoFFPinput returns a dataframe ready for applying the Flux Footprint Prediction (FFP) method developed by Kljun et al. (2015). The function performs data processing and temporal filtering of an ICOS L2 ARCHIVE or L2 INTERIM zip file. There is no need to unzip the compressed file of each site. When multiple ZIP files are downloaded from the Carbon Portal, the only action needed is to unzip the parent compressed file.

Although CPtoFFPinput takes as its primary input a compressed Level-2 Archive ICOS product, a CSV table (e.g. a CSV file previously returned by the same function) can be used as input (via the FFP.input.table parameter). In this case, the function can only perform a temporal filtering on the CSV file, or add station information (coordinates, UTC offset) if missing.

If the dates are set as NULL, the data is not filtered for a specific time interval. Alongside the FFP input data, the EC tower coordinates and the UTC offset are extracted from the metadata contained in the ZIP file, or silently by the function if an existing CSV file is used. When multiple zip files are downloaded, the parent zip folder contains a CSV file named ‘!TOC.csv’, used to extract the PID of each Level-2 product. 

The function always returns a four-element list containing: the site ID string, a vector with EC tower lat/lon coordinates, the PID, and the FFP input data frame. Depending on the data and site characteristics, not all elements may be returned (e.g., the PID if it is not available because the ‘!TOC.csv’ file doesn't exist, or the FFP input table when not all the parameters needed for calculating the FFP are available). 
When the FFP input table is saved as a CSV file, ‘site_id’, ‘lat’, ‘lon’, 'UTC_offset', and ‘PID’ are added to the first 5 columns of the table. The "TIMESTAMP" returned by the function corresponds to "TIMESTAMP_START".

If save.input.table is set to TRUE, in the R working directory the function creates a folder named as the site ID and a sub-folder named "Input", where the processed the FFP inputs are stored in a CSV format. The filename follows the structure: site id_start date_end date_FFP input.csv. 

MDS.input.list and ERA5.file.list are optional additional parameters used for refining the FFP calculation, by linking the FFP input to additional CSV files. They must contain the absolute paths of the corresponding file lists. The additional CSV files can be placed anywhere.

The MDS.input.list parameter uses one of the MDS (Marginal Distribution Sampling) method outputs to correctly locate in time the observations used for performing the gap-filling. The parameter needs at least one path of the file containing the indices of the half-hours used by MDS in the gap-filling process. The production of this file is currently being implemented in the ICOS ETC internal pipeline. When the file path is specified, two columns are added to the FFP data frame: the index of each observation and the list of indices used by the gap-filling approach, when it is performed. Additionally, when this parameter is used, no temporal filtering must be performed. 

The ERA5.file.list is used to link the ERA5 estimation do the planetary boundary layer height (PBLH) to the FFP input dataframe. It replaces the planetary boundary layer height estimated from data modelled internally. The file list must contain a yearly CSV file with the hourly values of PBLH. The CSV structure should consist of 2 columns named 'time' (timestamp in POSIX format) and 'value' (PBLH value). Half-hourly ERA5-derived observations are then obtained by linearly interpolating the hourly observations.

After the processing, both the MDS and ERA5 files are then moved into the 'Input' sub-folder of the site ID folder. If both MDS.input.list and ERA5.input.list are set to NULL (default), no MDS indices or ERA5 PBLH are collected. For the PBLH, if available, the internally modelled one is then used.  

If the FFP input table needs to be saved in another directory, the new path must be specified with the 'out.dir' parameter. Since no additional parameters are collected, any other data refinement and filtering is left to the end-user. 


 
## doFFP
Starting from the input parameters (e.g. a CSV file or an R object returned by the CPtoFFPinput function), it saves in dedicated folders the daily Flux Footprints Predictions (FFP; Kljun et al., 2015) for a specific Eddy Covariance site. The FFP R function can be downloaded here: https://footprint.kljun.net/download_2.php

| Argument           | Description                                                                        |
| ------------------ | -----------------------------------------------------------------------------------|
| FFP.input.df	   | Full path of an already existent csv, or an existent R dataframe                   |
| EC.tower.coords    | A vector containing the EC tower coordinates in WGS84 (lat, lon)                   |
| site.ID	         | ID of the site                                                                     |
| UTC_offset         | UTC offset of the site                                                             |
| PID                | PID of the input table                                                             |
| MDS.indices        | Performs MDS-style gap-filling for FFP. TRUE only in the ICOS ETC internal pipeline|
| FFP.function.path  | Full path of the FFP function                                                      |
| do.parallel        | If TRUE, the calculation of footprint is parallelized using the furrr package      |
| n.workers          | Number of workers for the parallelized function                                    |
| max.gap	         | Maximum number of consecutive half-hours to fill with LOCF                         |
| FFP.domain	   | Domain size [m]                                                                    |
| dx	               | Cell size of domain [m]                                                            |
| FFP.R              | Density thresholds for isopleths calculation                                       |
| which.FFP.R        | Among the density thresholds, which one export in the netCDF file                  |
| drop.trivial.srcs  | If TRUE, removes point sources with a trivial contribution                         |
| save.FFP.mtx.as    | FFP matrix saved as CSV (‘csv’), NetCDF (‘nc'), GeoTIFF ('gtiff') or no save (NULL)|
| save.plot.FFP.mtx  | If TRUE, a jpeg plot of FFP (both matrix and isopleths) is saved                   |
| return.isopleth    | If TRUE, the chosen isopleth (which.FFP.R) is saved in the netCDF file             |
| save.log           | If TRUE, log messages are written to a file                                        |
| save.ncplot	   | If TRUE, a sample plot of the nc file produced is saved                            |
| do.full.climatology| If TRUE, a single FFP climatology is calculated for the entire dataset provided    | 
| do.daily.average   | If TRUE, do daily FFP average                                                      |
| skip               | A vector for possibly Skip single days from computation                            |    
| feel.lucky         | Calculate footprint for a random day                                               |
| out.dir            | Directory of the out file folder, otherwise the working directory is used          |

For every observation (by default every half-hour) in a single day, the function processes a structured input table and runs the FFP function developed by Kljun et al. (2015), saving the output array in CSV, NetCDF (default) or GeoTIFF format. Plots of the footprint matrix and a log file can be saved. Additionally, FFP can be aggregated by averaging the observations to daily level or over the full dataset provided. 

If not already included in the dataframe (e.g. the CSV file saved by the CPtoFFPinput function), tower coordinates, site ID, UTC offset and PID arguments must be specified (e.g. using the R list returned by the CPtoFFPinput function). The function relies on the FFP model function (Kljun et al., 2015) which should be placed in the R working directory, preferably inside the subfolder ‘local_functions’ (default). Alternatively, is also possible to specify a different FFP function path. Since the FFP model calculation is time-consuming, function parallelization can be enabled (do.parallel=TRUE) using the furrr package. To prevent system crashes, the number of workers should be kept low (default: 3, according to the processing power of the machine). 

Input data gaps are filled up to a defined threshold (4 by default; max.gap argument). FFP function-related parameters, such as the spatial extent (1000 m by default; FFP.domain argument), resolution (1 m by default; dx argument), the FFP density thresholds for isopleth calculation (c(50, 70, 80, 90) by default; FFP.R argument), and which isopleth to save in the NetCDF file (which.FFP.R)—can be specified. In addition, pixels with a minimal contribution to the total density (threshold of 90%; drop.trivial.srcs argument) are removed by default. In the ICOS ETC internal pipeline, resolution and extent are optimized for each site by pre-calculating the FFP 90% isopleth length and using it to place each site in a specific extent-resolution class.

By default, a daily 48-layer NetCDF is produced. Alternatively, it is possible to calculate an aggregated FFP in two ways:
*    <ins>daily average</ins> (do.daily.average=TRUE) calculates a daily FFP average. As result, every daily NetCDF will have only one layer (the daily average one)
*    <ins>full climatology</ins> (do.full.climatology=TRUE) calculates a unique FFP climatology for the entire dataset provided, regardless its length. The date used as reference is the first day in the provided dataset
When MDS.indices is set to TRUE, the FFP average (FFP_climatology function) is automatically calculated for those observations where fluxes were gap-filled (performed using MDS), by taking as inputs the same half-hours used by MDS for the gap-filling.  

Single days can be arbitrarily excluded from the FFP computation (skip: FALSE by default). When needed, a boolean vector (TRUE: skip, FALSE: don't skip) with a length equal to the number of days of the dataset has to be specified. If the skip length is equal to one, the skip value is applied to the entire dataset. If is needed to process only a single random day, 'feel.lucky' must be set to TRUE.


### Outputs

In the R working directory, a folder named as the site ID with an 'Output' sub-folder are created if they don't exist. Depending on the saving options, the following sub-folders can be created inside the 'Output' folder: 
*	<ins>nc_files</ins>: if ‘nc’ is chosen as the saving option, this folder contains the daily netCDF files. Each file includes a set of dimensions, variables, and attributes. 

<ul>

The dimensions are related to the spatial (x, y) and temporal resolution (time), to the coordinates of the EC tower (EC_Coords) or the bounding box (BB_Coords), to the polygon isopleth chosen to be saved (instance: total number of polygons; nodes: total number of vertices) and its eventual multiparts (part: total number of parts).

| Variable (Dimensions) | Description|
| - | - |
| x (x), y (y), time (time) | Variables strongly linked to the corresponding dimensions. They place each FFP matrix in the right position and time |
| FFP (x, y, time) | Daily array containing the 2D FFP matrices calculated every 30 minutes. Each matrix was compressed and converted to integer values, providing scale and offset factors to allow software to retrieve original values. When the FFP model fails, a no-data grid (-9999 values) is added to keep the total number of time steps equal (48) for every file generated. When the FFP is aggregated (daily average or full climatology), the NetCDF array always consists of a single layer |
| polygon_id (instance), observation_time (instance), polygon_area (instance) | Variables strictly related to the 70% isopleth (one for each time interval), containing the ID (polygon_id), the corresponding timestamp (observation_time) and the area in m<sup>2</sup> (polygon_area) |
| quality_flag (instance) | Quality flag of FFP. 0: FFP and isopleths correctly calculated; 1: FFP and isopleths not calculated; 2: FFP correctly calculated, but the 90% isopleth could not be retrieved; 3: FFP correctly calculated, but the 70% isopleth could not be retrieved; 4: FFP calculated using the observations identified by the gap-filling approach|
| geometry_container, UTM_Coordinate_System | They define the associations between the variables linked to the geometry (the isopleth chosen for saving) and the correct Coordinate Reference System (CRS, in WGS84/UTM system) for each variable that must be spatially located. The right CRS is automatically detected starting from the WGS84 coordinates of the EC tower |
| nodes_count (instance), nodes_x (nodes), nodes_y (nodes) | Respectively, the count of the total number of vertices for the isopleth, and the X and Y coordinates of the isopleth vertices |
| part_node_count(part), interior_ring(part) | They are needed when multipart polygons o polygons with holes are produced. Part_node_count identifies the number of nodes for each part, while interior_ring detects which polygon has an internal hole |
| EC_tower_coordinates (EC_Coords) | Coordinates in WGS84/UTM of the eddy covariance flux tower |
| BB_coordinates_x (BB_Coords), BB_coordinates_y (BB_Coords) | Vectors containing the bounding box coordinates in WGS84 of the FFP array |  

The attributes provide metadata for a better description of the file, according to the ATMODAT (Atmospheric Model Data) standards. The netCDF file is structured following the CF (Climate and Forecast) 1.9 standard. 

</ul>

*	<ins>GTiff_files</ins>: if ‘gtiff’ is chosen as saving option, this folder contains the output daily FFP matrices saved in a multi-band GeoTIFF file. The file contains only the FFP matrices (single matrix in case of average or climatology options), each one stored as single layer. Each matrix was compressed and converted to integer values, but scale and offset factors are internally provided to allow software to retrieve original values. The file is named following the structure: site ID_FFP2D_YYYYMMDD.tiff. Additionally, a set of metadata related to file creators, methods used, data standards and file version are provided within the GeoTIFF file. The metadata also includes the QC (Quality Control) values. 
*	<ins>Individual_footprint_matrices</ins>: if ‘csv’ is chosen as the saving option, this folder contains the output FFP matrices saved in CSV format. Each file is named following the structure: site ID_FFP2Dmtx_YYYYMMDDHHmm.csv.
*	<ins>Individual_footprint_plots</ins>: if ‘save.plot.FFP.mtx’ is set to TRUE, this folder contains the plots of the output FFP matrices saved in PNG format. Each file is named following the structure: site ID_FFP2Dmtx_YYYYMMDDHHmm.png.
*	<ins>log</ins>: if ‘save.log’ is set to TRUE, this folder contains a log file for each function run. The file is named following the structure: site ID_log_start date_end date.log

## References
Kljun, N., Calanca, P., Rotach, M. W., and Schmid, H. P.: A simple two-dimensional parameterisation for Flux Footprint Prediction (FFP), Geosci. Model Dev., 8, 3695-3713, 2015. doi:10.5194/gmd-8-3695-2015


## Examples
```
# Retrieve the file path
zip_file="ZIP\\My data cart\\ICOSETC_BE-Bra_ARCHIVE_INTERIM_L2.zip"

# Generate the FFP inputs
CPtoFFPinput(zip.file=zip_file, save.input.table = T, return.input = T, start.date = NULL, end.date = '2023-12-31')
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
CPtoFFPinput(FFP.input.table="IT-SR2\\Input\\IT-SR2_20191028_20231231_FFPinput.csv", save.input.table = T,
             start.date = '2023-01-01', end.date = '2023-12-31')
```

```
# Get the path of each ZIP file
zip.list=list.files("ZIP\\My data cart\\", pattern = 'zip', full.names = T)

# Generate the FFP inputs as R object for only one day
input.list=lapply(zip.list, function(x) CPtoFFPinput(zip.file=x, save.input.table = F, 
                                           start.date = '2023-11-21', end.date = '2023-11-21'))

# Process only the first site
Site_1= input.list[[1]]

# Calculate FFP only for the first site
doFFP(FFP.input.df = Site_1[['FFP_input_parameters']],
      EC.tower.coords = Site_1[['tower_coordinates']],
      site.ID=Site_1[['site_id']],
      UTC.offset=.x[['UTC_offset']],
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
                                  UTC.offset=.x[['UTC_offset']],
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
```
# Get the path of each ZIP file
zip.list=list.files("ZIP\\My data cart\\", pattern = 'zip', full.names = T)

# Generate the FFP inputs as R object for only one day
input.list=lapply(zip.list, function(x) CPtoFFPinput(zip.file=x, save.input.table = F, return.input = T, 
                                           start.date = '2023-11-21', end.date = '2023-11-30'))

# Process only the first site
Site_1= input.list[[1]]

# Calculate FFP full climatology
doFFP(FFP.input.df = Site_1[['FFP_input_parameters']],
      EC.tower.coords = Site_1[['tower_coordinates']],
      site.ID=Site_1[['site_id']],
      UTC.offset=.x[['UTC_offset']],
      PID=Site_1[['PID']],
      FFP.function.path="local_functions\\calc_footprint_FFP_climatology.R",
      max.gap = 4, FFP.domain = 1000, dx = 1, FFP.R = c(50,70,80,90), which.FFP.R=70,
      drop.trivial.srcs = TRUE, save.FFP.mtx.as="nc",          
      save.plot.FFP.mtx = FALSE, return.isopleth=TRUE, save.log=TRUE, save.ncplot=TRUE, 
      do.full.climatology=TRUE)
```
```
# Get the path of each ZIP file
zip.list=list.files("ZIP\\My data cart\\", pattern = 'zip', full.names = T)

# Generate the FFP inputs as R object for only one day
input.list=lapply(zip.list, function(x) CPtoFFPinput(zip.file=x, save.input.table = F, return.input = T, 
                                           start.date = '2023-11-21', end.date = '2023-11-30'))

# Process only the first site
Site_1= input.list[[1]]

# Calculate in parallel FFP daily climatology excluding the days 2, 3, 7 # Set 4 workers #
# Save the output as GeoTiff #
doFFP(FFP.input.df = Site_1[['FFP_input_parameters']],
      EC.tower.coords = Site_1[['tower_coordinates']],
      site.ID=Site_1[['site_id']],
      UTC.offset=.x[['UTC_offset']],
      PID=Site_1[['PID']],
      FFP.function.path="local_functions\\calc_footprint_FFP_climatology.R",
      do.parallel=TRUE, n.workers=4,
      max.gap = 4, FFP.domain = 1000, dx = 1, FFP.R = c(50,70,80,90), which.FFP.R=70,
      drop.trivial.srcs = TRUE, save.FFP.mtx.as="gtiff",          
      save.plot.FFP.mtx = FALSE, return.isopleth=TRUE, save.log=TRUE, save.ncplot=TRUE, 
      do.daily.climatology=TRUE, skip=c(F, T, T, F, F, F, T, F, F, F))
```
```
zip.list=list.files("ZIP\\My data cart\\", pattern = 'zip', full.names = T)

# Select only the first site #
zip_file=zip.list[1]

# MDS fiile list #
MDSlist=list.files("MDS_indices\\", full.names = T, pattern=".csv")

# ERA5 file list # 
ERA5list=list.files("ERA5_PBLH\\", full.names = T, pattern=".csv")

# Generate the input list refining the input DF with ERA5 PBLH and MDS-derived data # 
input = CPtoFFPinput(zip.file=zip_file,
                    save.input.table = FALSE, 
                    start.date = NULL, end.date = NULL, 
                    MDS.file.list = MDSlist, ERA5.file.list = ERA5list, 
                    out.dir = NULL)

# Calculate FFP and save it as netCDF # 
doFFP(FFP.input.df = input$FFP_input_parameters, 
      EC.tower.coords = input$tower_coordinates, 
      UTC.offset = input$UTC_offset,
      site.ID = input$site_id, 
      PID = input$PID, 
      FFP.function.path="local_functions/calc_footprint_FFP_climatology.R",
      do.parallel=F, n.workers = 3, 
      skip=FALSE,       
      do.daily.average = FALSE, 
      do.full.climatology = FALSE,
      max.gap = 4, 
      FFP.domain = 1000, 
      dx = 2, 
      FFP.R = c(50,70,80,90), 
      which.FFP.R = 70, 
      drop.trivial.srcs = TRUE, 
      save.FFP.mtx.as = "nc", 
      save.plot.FFP.mtx = FALSE, 
      return.isopleth = TRUE, 
      save.log = FALSE, 
      save.ncplot=FALSE, 
      MDS.indices = TRUE, 
      feel.lucky = FALSE,
      out.dir = NULL)