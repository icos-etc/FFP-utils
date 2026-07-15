## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
## ---------------------- FFP CALCULATION - SAVE RESULTS AS CSV, GEOTIFF OR NETCDF ------------------------ ##
## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##

doFFP=function(FFP.input.df=NULL,         # input dataframe
               EC.tower.coords=NULL,      # Ec tower coordinates (Lat, Lon)
               UTC.offset=NULL,           # UTC offset of the site 
               site.ID=NULL,              # site.ID
               PID=NULL,                  # site PID
               FFP.function.path="local_functions/calc_footprint_FFP_climatology.R", # If specified, full path of FFP function, otherwise local path
               do.parallel=FALSE,         # parallelize the calculation of footprint
               n.workers=3,               # number of workers for future_map
               max.gap = 4,               # maximum number of consecutive half-hours to fill with LOCF
               FFP.domain = 1000,         # Domain size as an array of (xmin xmax ymin ymax) [m] *1.
               dx = 1,                    # Cell size of domain [m] (default is dx = dy = 1 m)
               FFP.R = c(50,70,80,90),    # Isopleths
               which.FFP.R = 70,          # Which isopleth save
               drop.trivial.srcs = TRUE,  # whether to remove point sources with a trivial contribution
               save.FFP.mtx.as='nc',      # csv, nc (NetCDF), gtiff(GeoTiff), NULL (no save)
               save.plot.FFP.mtx = FALSE, # Whether to save a plot of the FFP2D (matrix and isoplethes, not the spatial polygons)
               return.isopleth=TRUE,      # Return chosen isopleth in the netCDF file  
               save.log=TRUE,             # Save messages files to log
               save.ncplot=TRUE,          # Save a sample plot of nc the file
               do.full.climatology=FALSE, # Do full FFP climatology 
               do.daily.average=FALSE,    # Do daily FFP average
               skip=FALSE,                # Skip single days from computation
               MDS.indices=FALSE,         # Calculates average FFPs where data is gap-filled # TRUE only in ETC internal pipeline    
               feel.lucky=FALSE,          # Calculate footprint for a random day
               out.dir = NULL)            # If not null, the directory of the out file folder, otherwise the working directory is used  
{



  ##  
  ## [0] PRELIMINARY PART -----
  ##

  # Install pacman
  if (!require("pacman", quietly = T)) {install.packages("pacman")} 
  # Install EBImage
  if(!require("EBImage", quietly = T)) {pacman::p_load(BiocManager); BiocManager::install("EBImage")}
  # Load or install required packages via pacman
  pacman::p_load(crayon, data.table, sf, ncdf4, terra, ggplot2, logr, purrr, dplyr, readr, lubridate, EBImage)
  
  space="  "

  # Define message themes
  error <- crayon::red
  warn <- crayon::yellow
  note <- crayon::silver
  prog.mes <- crayon::cyan
  
  # Load local functions if they aren't already loaded 
  if(!exists('calc_footprint_FFP_climatology'))   {source(FFP.function.path)}
  
  
  # Starting message
  # Defining the start message total length
  message_total_length <- 86  
  
  # Message text
  message_text <- paste0("------ ", bold(site.ID), ": FFP calculation and ", save.FFP.mtx.as, " file saving ------")
  message_text_length <- nchar(crayon::strip_style(message_text)) 
  
  # Real margin
  margin <- strrep(" ", max(0, floor((message_total_length - message_text_length) / 2)))
  
  # Final cat
  cat('\n', strrep('*', message_total_length), sep='')
  cat('\n', margin, prog.mes(message_text), sep='')
  cat('\n', strrep('*', message_total_length), '\n', sep='')
  
  
  # If DF is character (and then csv path), read it
  if(is.character(FFP.input.df)) {
    FFP.input.df <- read_csv(FFP.input.df, col_types = cols())
    FFP.input.df$'TIMESTAMP' <- with_tz(FFP.input.df$'TIMESTAMP', tz ='GMT') # Set the correct timezone
  }

  # Warn if the DF TZ is not GMT #
  if(tz(FFP.input.df$TIMESTAMP) != "GMT") {
    cat(warn(paste0("The time zone of the input dataframe is not GMT, it shoud be corrected")))
  }
  
  # Inputs availability check
  if(is.null(FFP.input.df)){
    cat(error(paste0("\nERROR: No FFP input dataframe is specified for the station. FFP calculation is skipped!\n")))
    
  } else if(all(is.null(site.ID), !(any(colnames(FFP.input.df) %in% c("site_id"))))) {
    cat(error(paste0("\nERROR: No site ID is specified for the station. FFP calculation is skipped!\n")))
    
  } else if(all(is.null(EC.tower.coords), !(any(colnames(FFP.input.df) == "lat")))) {
    cat(error(paste0("\nERROR: No latitude is specified for the station. FFP calculation is skipped!\n")))
    
  } else if(all(is.null(EC.tower.coords), !(any(colnames(FFP.input.df) == "lon")))) {
    cat(error(paste0("\nERROR: No longitude is specified for the station. FFP calculation is skipped!\n")))
    
  } else {
    
    
    
    ##
    ## [1] DATAFRAME MANIPULATION -----
    ##
    
    # Starting message 
    cat(prog.mes(paste0('\n', bold('[1]'),' FFP input data manipulation\n')))

    # Move site id, coordinates and PID outside the DF, if necessary....
    # EC tower coordinates
    if(is.null(EC.tower.coords)) {
      
      EC.tower.coords <- c(unique(FFP.input.df$lat), unique(FFP.input.df$lon))
      
      # Remove them from DF
      FFP.input.df$lat <- NULL
      FFP.input.df$lon <- NULL
    }
    
    # Site id
    if(is.null(site.ID)) {
      
      site.ID <- unique(FFP.input.df$site_id)
      
      # Remove it from DF
      FFP.input.df$site_id <- NULL
    }
    
    # PID
    if(all(is.null(PID), colnames(FFP.input.df) %in% 'PID')) {
      
      PID <- unique(FFP.input.df$PID)
      
      # Remove it from DF
      FFP.input.df$PID <- NULL
    }
    
    # Check for NAs in the variables created, and in necessary stop the function # 
    stopifnot('Site ID is NA' = !is.na(site.ID), 'Tower coordinates are NAs' = !is.na(EC.tower.coords))
    
    # Assign coordinates to a variable
    latitude <- EC.tower.coords[1]
    longitude <- EC.tower.coords[2]
    
    # Check for the coordinates format (decimal)
    stopifnot(is.numeric(latitude), is.numeric(longitude), latitude < 90 & latitude > -90)
    
    UTM.zone <- dplyr::case_when( 
                                 latitude>56 & latitude<64 & longitude>3 & longitude<12 ~ 32,
                                 latitude>72 & latitude<84 & longitude>9 & longitude<21 ~ 33,
                                 latitude>72 & latitude<84 & longitude>21 & longitude<33 ~ 35,
                                 T ~ (floor((longitude + 180)/6) %% 60) + 1)
    
    EPSG.code <- dplyr::case_when(
                                  latitude>0 & latitude<84 & UTM.zone<10 ~ as.numeric(paste0(326, 0, UTM.zone)),
                                  latitude>0 & latitude<84 & UTM.zone>=10 ~ as.numeric(paste0(326, UTM.zone)),
                                  latitude<0 & latitude>-80 & UTM.zone<10 ~ as.numeric(paste0(327, 0, UTM.zone)),
                                  latitude<0 & latitude>-80 & UTM.zone>=10 ~ as.numeric(paste0(327, UTM.zone)),
                                  latitude > 84 & latitude < 90 ~ 32661,
                                  latitude < -80 & latitude > -90 ~ 32761)
    
    EC.tower <- sf::st_as_sf(x=data.frame(lat=EC.tower.coords[1], lon=EC.tower.coords[2]), 
                             coords = c('lon', 'lat'), crs=4326)

    # EC coordinates in UTM
    EC.tower.utm <- as.vector(st_coordinates(st_transform(EC.tower, EPSG.code)))
    
    if(is.null(out.dir))  {
      
      # Create main output directory (and eventually site directory)
      Site_dir <- paste0(getwd(), '/', site.ID)
      FFP.output.dir <- paste0(getwd(), '/', site.ID, '/', 'Output')
      
    } else  {
      
      # Create main output directory (and eventually site directory)
      Site_dir <- paste0(out.dir, '/', site.ID)
      FFP.output.dir <- paste0(out.dir, '/', site.ID, '/', 'Output')
      
    }
    
    if(!dir.exists(Site_dir)) {dir.create(Site_dir)}
    if(!dir.exists(FFP.output.dir)) {dir.create(FFP.output.dir)}
    
    # NON sensitive variables # These variables are filled with LOCF independently from the gap length(s)
    is.na_nsv <- lapply(FFP.input.df[c('hm', 'hc', 'd', 'z0', 'zm')], function(x) sum(is.na(x)))

    FFP.input.df[c('hm', 'hc', 'd', 'z0', 'zm')] <- lapply(FFP.input.df[c('hm', 'hc', 'd', 'z0', 'zm')], 
                                                           function(x) data.table::nafill(x, type = 'locf'))
    
    invisible(sapply(names(is.na_nsv), function(x) cat(note(paste0('\n- ', is.na_nsv[[x]], 
                                                                   ' missing values of ', bold(x), 
                                                                   ' have been filled using LOCF.')))))
    
    # LOG file starting messages 
    Log_list <- list()
    
    Log_list[[1]] <- paste0('###### Computing the footprints for the ', site.ID, ' station ######')
    Log_list[[2]] <- '' # Blank space
    Log_list[[3]] <- '------ FFP input data manipulation (gap-filling) ------'
    Log_list[[4]] <- '' # Blank space
    Log_list <- c(Log_list,sapply(names(is.na_nsv), 
                                  function(x) paste0('- ', is.na_nsv[[x]], ' missing values of ', x, 
                                                     ' have been filled using LOCF.'), USE.NAMES = F), '')
    
    # Sensitive variables # These variables are filled with LOCF only if gaps are shorter than 2 hours
    
    for(i in c('umean', 'ol', 'sigmav', 'ustar', 'wind_dir')) {

        isna <- sum(is.na(FFP.input.df[[i]]))                                           # count the NAs
        isna.rle <- rle(is.na(FFP.input.df[[i]]))                                       # Run length encoding on NAs
        isna.rle$'values' <- isna.rle$'values' & isna.rle$'lengths' > max.gap           # Re-write the values when NAs come by more than 4: these will not be filled
        isna.f <- which((is.na(FFP.input.df[[i]]) - inverse.rle(isna.rle)) != 0)        # store gapfilled umean indexes
        # FFP.input.df$'TIMESTAMP'[which((is.na(FFP.input.df$'umean') - inverse.rle(isna.umean.rle)) != 0)]   # See which half-hour will have gapfilled umean
        
        cat(note(paste0('\n- ', as.character(isna), ' missing values of ', bold(i), ' have been found, but only ', 
                        bold(as.character(length(isna.f))) ,' have been filled using LOCF lasting for less than 2 hours.')))
        FFP.input.df[[i]] <- data.table::nafill(FFP.input.df[[i]], type = 'locf')       # Fill every NAs
        FFP.input.df[[i]][inverse.rle(isna.rle)] <- NA                                  #  use inverse.rle to build the vector of indices to re-set to NA (longer gaps
        
        # Add message to log list
        Log_list <- c(Log_list, paste0('- ', as.character(isna), ' missing values of ', i,
                                       ' have been found, but only ', as.character(length(isna.f)), 
                                       ' have been filled using LOCF lasting for less than 2 hours'))
    }
    
    
    
    ##
    ## [1.1] DATAFRAME ADDITIONAL MANIPULATION: CLIMATOLOGY-AVERAGE OPTIONS AND TIME ZONE -----
    ##
    
    # Single FFP #
    # if timestamp is reported as ISO (character), convert it to a time object
    if (!any(grepl('POSIX', class(FFP.input.df$'TIMESTAMP')))) { 
      FFP.input.df$'TIMESTAMP' <- as_datetime(FFP.input.df$'TIMESTAMP', tz ='GMT')
    } 
    
    stopifnot('Choose one between do.climatology and do.daily.average' = 
                !(do.full.climatology == T & do.daily.average == T))

    # Climatology and MDS error #
    stopifnot('Choose one between climatology-average options and MDS' = 
                !(all(any(do.full.climatology, do.daily.average), MDS.indices)))
    
    # Set the correct time format # 
    time_format='%Y-%m-%d %H:%M'

    # If climatology, modify the single timestamp keeping only the start day #
    if(do.full.climatology)  { # Calculate the overall mean FFP # 
      
      if(any(skip))
      {
        # Filter the input dataframe according to the skip parameter #
        if(all(length(skip) != length(unique(as_date(FFP.input.df$TIMESTAMP, tz ='GMT'))), length(skip) != 1)) 
        {stop("Skip vector length is not equal than input DF days length")}
        
        else { # Filter the dataframe according to the skip parameter # Then change the skip parameter #
        FFP.input.df <- FFP.input.df[as_date(FFP.input.df$TIMESTAMP, tz ='GMT') %in% (unique(as_date(FFP.input.df$TIMESTAMP, tz ='GMT'))[!skip]), ]
        
        skip=F
        
        }
      }
      
      # DF filtering # If the MDS obs are available, calculate the climatology FFP using only the valid ones # 
      if(any(colnames(FFP.input.df) %in% "MDS_Idx"))  {
        
        cat(warn("\nFiltering out the Gap-filled observation from the climatology calculation"))
        
        FFP.input.df=FFP.input.df[sapply(FFP.input.df$MDS_Idx, is.null), ]  
        
      }
      
      # Set the first and last day # 
      first_day <- min(as_date(FFP.input.df$'TIMESTAMP', tz ='GMT'))
      last_day <- max(as_date(FFP.input.df$'TIMESTAMP', tz ='GMT'))
      
      # Set the correct time format # 
      time_format <- '%Y-%m-%d'
      
      # Modify the timestamp to keep only the first one # 
      FFP.input.df$'TIMESTAMP' <- as_date(FFP.input.df$'TIMESTAMP', tz ='GMT')[1]
      
      # Drop nas (only for climatology-average options) #
      FFP.input.df <- na.omit(FFP.input.df)
      
    }
    
    # Compute the daily average, if required #
    if(do.daily.average)  { # Calculate the daily mean FFP # 
      
      if(any(skip))
      {
        # Filter the input dataframe according to the skip parameter #
        if(all(length(skip) != length(unique(as_date(FFP.input.df$TIMESTAMP, tz ='GMT'))), length(skip) != 1)) 
        {stop("Skip vector length is not equal than input DF days length")}
        
        else { # Filter the dataframe according to the skip parameter # Then change the skip parameter #
          FFP.input.df <- FFP.input.df[as_date(FFP.input.df$TIMESTAMP, tz ='GMT') %in% (unique(as_date(FFP.input.df$TIMESTAMP, tz ='GMT'))[!skip]), ]
          
          skip=F
          
        }
      }
      
      # DF filtering # If the MDS obs are available, calculate the average using only the valid ones # 
      if(any(colnames(FFP.input.df) %in% "MDS_Idx"))  {
        
        cat("\n")
        cat(warn("\nFiltering out the Gap-filled observation from the average calculation"))
        
        FFP.input.df=FFP.input.df[sapply(FFP.input.df$MDS_Idx, is.null), ]  
        
      }
      
      # Set the correct time format # 
      time_format <- '%Y-%m-%d'
      
      # Modify the timestamp to keep only the daily timestamp # 
      FFP.input.df$'TIMESTAMP' <- as_date(FFP.input.df$'TIMESTAMP', tz ='GMT')
      
      # Drop nas (only for climatology-average options) #
      FFP.input.df <- na.omit(FFP.input.df)
      
    }
    
    # Create a days factor to split the dataset into daily chunks
    days.factor <- as.factor(format(FFP.input.df$'TIMESTAMP', '%Y%m%d', tz = "GMT"))
    
    # skip length check 
    if(all(length(skip) != length(levels(days.factor)), length(skip) != 1)) 
    {error("Skip vector length is not equal than input DF days length")}
    
    
    ## -- TIME ZONE -- ## 
    
    if(is.null(UTC.offset)) {
    
      cat(error("Station UTC offset missing!"))
      
      #cat(warning("Retrieve the UTC offset from the EC tower position. It may differ from the UTC offset of the data"))
      #
      #pacman::p_load(lutz)
      # # Get the correct time zone from tower coordinates #
      # station_tz <- tz_lookup_coords(EC.tower.coords[1], EC.tower.coords[2], method = "accurate")
      # 
      # # Calculate the perpetual solar offset (Take a random winter date)
      # winter_date <- as.Date(paste0(year(FFP.input.df$TIMESTAMP[1]), "-01-01"))            
      # time_offset_hours <- tz_offset(winter_date, station_tz)$utc_offset_h
      # 
      # time_offset_hours <- tz_offset("2024-01-01", station_tz)$utc_offset_h
    
    } else  { # Get the time zone from the input data
      
      time_offset_hours=UTC.offset  
      
    }
    
    # Create the offset string 
    time_offset_string <- sprintf("%+03d:00", time_offset_hours)
    
    # Customize domain # When a global domain is provided 
    if(!is.null(FFP.domain)) {FFP.domain.ext <- FFP.domain * c(-1, 1,-1, 1)}
    
    # Total dys available # 
    lev_days <- levels(days.factor)
    
    # Feel lucky? # Process a single random day # 
    if(feel.lucky)  {
      
      if (feel.lucky && length(skip) > 0 && !all(skip == FALSE)) {
        stop("Choose one between skip and feel.lucky")
      }
      
        lev_days <- sample(lev_days, 1)
        
    }
    
    cat('\n')    ### -- DF MANIPULATION ENDING -- ###
    cat('\n', strrep('*', message_total_length), sep='')
    
    
    
    ##
    ## [2] FFP CALCULATION -----
    ##
    
    for (i.day in lev_days[!skip]) {
      
      # Split the input into days chunks and grab the current one
      FFP.input.df.cur <- data.frame() 
      FFP.input.df.cur <- tibble(split(FFP.input.df, days.factor)[[i.day]])
      
      # Set the current day #
      current_day <- unique(as_date(format(FFP.input.df.cur$TIMESTAMP, '%Y-%m-%d', tz = "GMT")))
      
      # Set the current day string #
      current_day_str <- unique(format(FFP.input.df.cur$TIMESTAMP, '%Y-%m-%d', tz = "GMT"))
      
      ## Warn if the current df legnth is zero ## 
      if(nrow(FFP.input.df.cur) == 0) {
        
        cat(warn(paste0(i.day, ": input data is not available")))
        
        next
      } 
      
      # Retrieving Custom dxtent and resolution # When FFP domanin and res is null, search in the input table # 
      if(all(is.null(FFP.domain), is.null(dx))) {
        
        # If the columns are null, stop # 
        stopifnot('No Extent resolution provided' = !all(is.null(FFP.input.df.cur$Resolution)), 
                  'No Extent extent provided' = !all(is.null(FFP.input.df.cur$Extent)))
        
        # Get the extent and resolution from the columns in the file #
        FFP.domain=unique(FFP.input.df.cur$Extent)
        dx=unique(FFP.input.df.cur$Resolution)
        
        # If multiple resolutions # Error # 
        stopifnot('Multiple extents' = length(FFP.domain)==1, 'Multiple resolutions' = length(dx)==1)
        
      }
      
      # Retrieve the bounding box # 
      FFP.domain.ext <- FFP.domain * c(-1, 1,-1, 1)
      
      # Starting message
      # Defining the start message total length
      message_total_length <- 86  
      
      # Message text
      message_text <- paste0('------ ', bold(current_day_str), ': ', 'Start of calculations ------')
      message_text_length <- nchar(crayon::strip_style(message_text)) 
      
      # Real margin
      margin <- strrep(" ", max(0, floor((message_total_length - message_text_length) / 2)))
      
      # Final cat
      cat('\n', margin, prog.mes(message_text), '\n', sep='')
      
      # Quantify survived NAs from input data (don't drop them!)
      rows.na <- nrow(FFP.input.df.cur[, !(names(FFP.input.df.cur) %in% c("MDS_Idx"))]) - 
        nrow(na.omit(FFP.input.df.cur[, !(names(FFP.input.df.cur) %in% c("MDS_Idx"))]))
      cat(warn(paste0('\n[NOTE] ', #bold(current_day_str), ': ', 
                      bold(as.character(rows.na)), ' rows ','(',as.character(round(rows.na/nrow(FFP.input.df.cur)*100)),'% of data) are data gaps that cannot be filled.\n')))
      
      # Add messages to log
      Log_list <- c(Log_list, '', '', 
                    paste0('------ ', current_day_str, ': ', 'start of calculations ------'), '')
      
      Log_list <- c(Log_list, paste0('[NOTE] ', #current_day_str, ': ', 
                                     as.character(rows.na), ' rows ','(',as.character(round(rows.na/nrow(FFP.input.df.cur)*100)),
                                     '% of data) are data gaps that cannot be filled.'), '')
      
      # Check (here) if there are all NA values and remove them   
      if(all(complete.cases(FFP.input.df.cur[, which(names(FFP.input.df.cur) == "d"):which(names(FFP.input.df.cur) == "PBL")]) == F)) {
          
        cat(warn(paste0('\n [!] ', #bold(current_day_str), ': ',  
                        'no input data is available - model is not computed - no file is produced', '\n')))
        
        # Add messages to log
        Log_list <- c(Log_list, '', paste0('[!] ', #current_day_str, ': ',  
                                           'no input data is available - model is not computed - no file is produced'), '')
        
      } else {
        
        if(all(do.full.climatology == F, do.daily.average==F))  {
          
          cat(prog.mes(paste0('\n', bold('[2]'), #bold(current_day_str), ': ', 
                              ' Computing the half-hourly 2D footprint for all the timestamps at once', '\n')))
          
          Log_list <- c(Log_list, paste0("FFP type: half-hourly"))
          
        } else if(do.daily.average == T) {
          
          cat(prog.mes(paste0('\n', bold('[2]'), #bold(current_day_str),': ', 
                              ' Computing the 2D daily footprint average'), '\n'))
          
          Log_list <- c(Log_list, paste0("FFP type: daily average"))
          
        } else if(do.full.climatology == T) {
          
          cat(prog.mes(paste0('\n', bold('[2]'), #bold(first_day), bold('-'), bold(last_day), ': ',
                              ' Computing the 2D footprint climatology for the full dataset provided', '\n')))
          
          Log_list <- c(Log_list, paste0("FFP type: full climatology"))
          
        }
        
        # MDS conditional #
        if(MDS.indices)
        {
          
          # hybrid list, nested map ## 
          FFP.input.list.cur <- map2(.x = split(FFP.input.df.cur, FFP.input.df.cur$TIMESTAMP), 
                                     .y = map(FFP.input.df.cur$MDS_Idx, function(x) 
                                       FFP.input.df[FFP.input.df$Idx %in% x, ]), 
                                     ~ if(nrow(.y) != 0)  {.y} else {.x})
          
        } else  {
          
          FFP.input.list.cur <- split(FFP.input.df.cur, FFP.input.df.cur$TIMESTAMP)
          
        }
       
        
        ## FFP main function ##
        if(do.parallel) {
          
          # FFP calculation (future parallel)
          pacman::p_load(furrr)
          plan(multisession, workers=n.workers) # By default: low number of workers to avoid PC crashing
          
          # FFP calculation #
          FFP <- future_map(FFP.input.list.cur, 
                            function(x) { 
                              tryCatch(
                                {res <- calc_footprint_FFP_climatology(
                                  zm = x$zm,
                                  z0 = x$z0,
                                  umean = x$umean,
                                  h = x$PBL,
                                  ol = x$ol,
                                  sigmav = x$sigmav,
                                  ustar = x$ustar,
                                  wind_dir = x$wind_dir,
                                  domain = FFP.domain.ext,
                                  dx = dx,
                                  dy = dx,
                                  r = FFP.R,
                                  rslayer = 1,
                                  smooth_data = 1,
                                  pulse = 0)
                                
                                list(result = res, error = NULL)},
                                
                                error = function(e) {
                                  list(result = NULL, 
                                       error = paste0(unique(x$TIMESTAMP), ": ", gsub("\n", "", as.character(e), fixed = TRUE)))}
                              )
                            }, .progress = TRUE)
          
        } else {
          
          # Function not parallelized
          FFP <- lapply(FFP.input.list.cur, 
                        FUN = function(x) {
                          tryCatch(
                            {res <- calc_footprint_FFP_climatology(
                              zm = x$zm,
                              z0 = x$z0,
                              umean = x$umean,
                              h = x$PBL,
                              ol = x$ol,
                              sigmav = x$sigmav,
                              ustar = x$ustar,
                              wind_dir = x$wind_dir,
                              domain = FFP.domain.ext,
                              dx = dx,
                              dy = dx,
                              r = FFP.R,
                              rslayer = 1,
                              smooth_data = 1,
                              pulse = 0)
                            
                            list(result = res, error = NULL)},
                            
                            error = function(e) {
                              list(result = NULL, 
                                   error = paste0(unique(x$TIMESTAMP), ": ", gsub("\n", "", as.character(e), fixed = TRUE)))}
                          )
                        }
          )
        }
        
        
        # Split the results # 
        FFP.ls <- sapply(FFP, "[", "result")
        Error <- unlist(sapply(FFP, "[", "error"), use.names = F)
        rm(FFP)
        
        # Free unused memory
        gc()
        
        # Add error message (if there is any) to log #
        if(!(is.null(Error)))
        {
          cat(error(paste0(Error, '\n')))
          Log_list <- c(Log_list, Error, '')  
        }
        
        # Names are not necessary
        names(FFP.ls) <- NULL
        
        # Format time. It needs to be always a valid value (no gaps allowed in the final product)
        cur.ts <- unique(format(FFP.input.df.cur$'TIMESTAMP', '%Y%m%d%H%M', tz = "GMT"))
        
        ## Extract FFP.R list ## 
        FFP.R.list <- lapply(FFP.ls, function(x) x[['r']])
        
        
        
        ##
        ## [2.1] FFP QUALITY ASSESSMENT -----
        ##
        
        # Computing error index # Index 1: no FFP matrix
        FFP.err.indx_1 <- which(sapply(FFP.R.list, function(x) ifelse(is.null(x), yes = T, no=F)))
        
        # Index 2: FFP matrix yes, but not the higher FFP R (matrix is thus not clipped)
        FFP.err.indx_2 <- which(sapply(FFP.R.list, function(x) !any(x == max(FFP.R)/100, na.rm = TRUE)))

        # Index 3: FFP matrix yes, 70% isopleth no 
        if(any(FFP.R %in% which.FFP.R))  {
          
          FFP.err.indx_3 <- which(sapply(FFP.R.list, function(x) !any(x == which.FFP.R/100, na.rm = TRUE)))
          
        } else  {   # If none is selected (unrealistic), flag if multiple isopleths are missing # 
          
          FFP.err.indx_3 <- which(sapply(FFP.R.list, function(x) sum(is.na(x)) > 1))
          
        }
        
        # Flag 4 # If the FFP derives from a gapfilling approach #
        if(MDS.indices) {
          
          # Check when the climatology-average is calculated # 
          FFP.err.indx_4=which(sapply(FFP.input.df.cur$MDS_Idx, function(x) !is.null(x)))
          
        }
        
        # Re-create the full QC modifying 1, 2, 3 and eventually 4  
        QC <- rep(0, length(FFP.ls))

        # The most permissive first # 
        QC[FFP.err.indx_2] <- 2
        QC[FFP.err.indx_3] <- 3
        
        # Including the MDS flag # 
        if(MDS.indices) {
          
          QC[FFP.err.indx_4] <- 4
          
        }

        # The most sever last # Valid for both climatology-average and non-climatology #          
        QC[FFP.err.indx_1] <- 1
        
        
        # Error messages handling # 
        for(i in seq_along(QC)) {  
          if(QC[i]==1) {
            
            # Print #
            cat(warn(paste0(space, '[!] ', bold(format(FFP.input.df.cur$TIMESTAMP[i], time_format, tz = "GMT")), ': ',
                            'according to the model requirements, FFP was not computed', '\n')))
            
            # Log #
            Log_list <- c(Log_list,  
                          paste0('[!] ', format(FFP.input.df.cur$TIMESTAMP[i], time_format, tz = "GMT"), ': ', 
                                 'according to the model requirements, FFP was not computed'))
            
          } else if(QC[i]==2) {
            
            # Print #
            cat(warn(paste0(space, '[!] ', bold(format(FFP.input.df.cur$TIMESTAMP[i], time_format, tz = "GMT")), ': ', 
                            'according to the model requirements, ', FFP.R[which.max(FFP.R)], 
                            '% isoline was not computed', '\n')))
            
            # Log #
            Log_list <- c(Log_list,  
                          paste0('[!] ', format(FFP.input.df.cur$TIMESTAMP[i], time_format, tz = "GMT"), ': ', 
                                 'according to the model requirements, ', FFP.R[which.max(FFP.R)], 
                                 '% isoline was not computed'))
            
          } else if(QC[i]==3) {
            
            # Print #
            cat(warn(paste0(space, '[!] ', bold(format(FFP.input.df.cur$TIMESTAMP[i], time_format, tz = "GMT")), ': ', 
                            'according to the model requirements, ', 
                            paste0(FFP.R[which(is.na(FFP.R.list[[i]]))], '% '), 
                            'isoline was not computed', '\n')), sep = "")
            
            # Log #
            Log_list <- c(Log_list, 
                          paste0('[!] ', format(FFP.input.df.cur$TIMESTAMP[i], time_format, tz = "GMT"), ': ', 
                                 'according to the model requirements, ', paste0(FFP.R[which(is.na(FFP.R.list[[i]]))], '% '), 
                                 'isoline was not computed')) 
          } else if(QC[i]==4) {
            
            # Print #
            cat(warn(paste0(space, '[!] ', bold(format(FFP.input.df.cur$TIMESTAMP[i], time_format, tz = "GMT")), ': ', 
                            'according to the data filtering, gap-filled FFP is computed', '\n')), sep = "")
            
            # Log #
            Log_list <- c(Log_list, 
                          paste0('[!] ', format(FFP.input.df.cur$TIMESTAMP[i], time_format, tz = "GMT"), ': ', 
                                 'according to the data filtering, gap-filled FFP is computed')) 
          }
        } # FFP calculation ending....
        
        
        # Check if at least one FFP was computed #
        if(length(FFP.err.indx_1) == length(FFP.ls))  {
          
          cat(warn(paste0('\n', space, '[!] ', #bold(current_day_str), ': ', 
                          'According to the model requirements, no FFP was computed. If specified, the nc file is created anyway')))
          
          # Add message to log
          Log_list <- c(Log_list, sapply(FFP.err.indx_1, 
                                         function(x) paste0("[!] ", #current_day_str, ': ', 
                                                            'According to the model requirements, no FFP was computed')), '')  
          
        } else { 
          
          # Plot and matrix saving looping over the list elements 
          for(i in which(QC %in% c(0, 2, 3, 4))) {
            
            # Remove eventual (unrealistic!) negative values # 
            FFP.ls[[i]]$'fclim_2d'[FFP.ls[[i]]$'fclim_2d' < 0] <- 0
            
            if (drop.trivial.srcs) {
              
              # fetch the footprint value at 90% (or the maximum desired) cumulative distribution
              FFP.90trshd <- NULL
              FFP.90trshd <- FFP.ls[[i]]$'fr'[which.max(FFP.R)]
              
              # and remove the FFP values lower that that
              FFP.ls[[i]]$'fclim_2d'[FFP.ls[[i]]$'fclim_2d' < FFP.90trshd] <- NA
              
            }
            
            
            
            ##
            ## [3] OUTPUT: PLOT -----
            ##
            
            if (save.plot.FFP.mtx) {
              # # Two‐dimensional footprint with contour lines of R%
              
              # Do it once
              if(i == min(which(QC %in% c(0, 2, 3)))) {
                cat(prog.mes(paste0('\n', bold('[3]'), #bold(current_day_str), ': ', 
                                    ' Saving FFP matrix plots ...\n')))
                
                # Add message to log 
                Log_list <- c(Log_list, paste0(current_day_str, ': ', 
                                               'FFP matrix plots saved'), '')
              }
              
              # Plots folder #
              FFP.mtx.img.dir <- paste0(FFP.output.dir, '/individual_footprints_plots/')
              
              # Create folder if not exists
              if(!dir.exists(FFP.mtx.img.dir))  {dir.create(FFP.mtx.img.dir)}
              
              # Create a dataframe for the plot 
              DF_grid <- data.frame(X=as.vector(t(FFP.ls[[i]]$'x_2d')) + EC.tower.utm[1], 
                                    Y=as.vector(t(FFP.ls[[i]]$'y_2d')) + EC.tower.utm[2], 
                                    Z=as.vector(FFP.ls[[i]]$'fclim_2d'))
              
              # Create a dataframe for the isopleths 
              names(FFP.ls[[i]]$xr) <- FFP.ls[[i]]$r
              names(FFP.ls[[i]]$yr) <- FFP.ls[[i]]$r
              
              DF_iso <- data.frame(stack(FFP.ls[[i]]$xr), stack(FFP.ls[[i]]$yr)['values'])
              DF_iso$values <- DF_iso$values + EC.tower.utm[1]
              DF_iso$values.1 <- DF_iso$values.1 + EC.tower.utm[2]
              
              # Timestamp label and plot name
              if(all(do.full.climatology == F, do.daily.average==F))  {
                
                timestamp_label <- paste0('Timestamp: ', FFP.input.df.cur[i, 'TIMESTAMP'], ' GMT')
                ffp.plotname <- paste0(FFP.mtx.img.dir, site.ID, '_', 'FFP2Dplot', 
                                       format(FFP.input.df.cur$'TIMESTAMP'[i], '%Y%m%d%H%M'), '.jpeg')
                
              } else if(do.full.climatology == T) {
                
                timestamp_label <- paste0('Timestamp: ', first_day, '-', last_day)
                ffp.plotname <- paste0(FFP.mtx.img.dir, site.ID, "_FFP2Dplot_", "_FFP2D_full_climatology_", 
                                       first_day, '_', last_day, ".jpeg")
                
              } else if(do.daily.average == T) {
                
                timestamp_label <- paste0('Timestamp: ', current_day_str)
                ffp.plotname <- paste0(FFP.mtx.img.dir, site.ID, "_FFP2Dplot_", "_FFP2D_daily_average_", 
                                       current_day_str, ".jpeg")
                
              }
              
              # Plot              
              ggplot(DF_grid)+
                geom_raster(aes(x = X, y = Y, fill = Z), alpha=0.85)+
                geom_path(data=DF_iso, aes(x = values, y = values.1, group = ind), color='grey90')+
                scale_fill_viridis_c(na.value = "transparent", option = 'D')+
                geom_vline(xintercept = EC.tower.utm[1], color='grey70')+
                geom_hline(yintercept = EC.tower.utm[2], color='grey70')+
                annotate(geom = 'point', x=EC.tower.utm[1], y=EC.tower.utm[2], 
                         shape=21, fill='red', size=2)+
                labs(title=paste0('Footprint estimation at ', site.ID), subtitle=timestamp_label, 
                     x= 'X (m)', y= 'Y (m)', fill='Density')+
                # 500 m limit from the tower center #
                xlim(EC.tower.utm[1]-500, EC.tower.utm[1]+500)+
                ylim(EC.tower.utm[2]-500, EC.tower.utm[2]+500)+
                theme(axis.line = element_line(colour = "black", linewidth = 0.1), 
                      panel.grid = element_blank(), 
                      axis.title = element_text(face="bold"),
                      axis.title.x = element_text(margin=margin(0.3,0,0,0, unit="cm")),
                      axis.title.y = element_text(margin=margin(0,0.3,0,0, unit="cm")),
                      panel.background=element_rect(fill="white", linewidth = 1, color="black"), 
                      plot.title = element_text(size=14, hjust = 0.5, face='bold'),  
                      plot.subtitle = element_text(face = "bold"), 
                      legend.title = element_text(face = "bold"), 
                      legend.background = element_rect(colour = 'NA', fill=NA))
              
              # Save the plot
              ggsave(ffp.plotname, width = 23, height = 13, units='cm', dpi = 400)
              
            } # Plot condition ending...
            
          } # Loop over valid matrices ending...
          
          # Add and extra space to the messages
          #cat('\n')
          
          
          
          ##
          ## [4] OUTPUT: CSV MATRIX -----
          ##
          
          if(!is.null(save.FFP.mtx.as)) {
            
            if (save.FFP.mtx.as == 'csv') {
              
              cat(prog.mes(paste0('\n', bold('[4]'), #bold(current_day_str), ': ', 
                                  ' Saving FFP matrices as csv(s)', '\n')))
              
              # Add message to log 
              Log_list <- c(Log_list, paste0('\n', current_day_str, ': ', 
                                             'FFP matrices saved as csv(s)'), '')
              
              # Create folder if not exist
              FFP.mtx.dir <- paste0(FFP.output.dir, '/individual_footprints_matrices/'); dir.create(FFP.mtx.dir)
              
              mapply(function(x, y) fwrite(x[['fclim_2d']], 
                                           paste0(FFP.mtx.dir, site.ID, '_FFP2Dmtx_', y, '.csv'), 
                                           quote = F, row.names = F, col.names = F, verbose = FALSE), 
                     FFP.ls[-FFP.err.indx_1],  format(FFP.input.df.cur[-FFP.err.indx_1 ,'TIMESTAMP'], '%Y%m%d%H%M')) 
              
            }
            
          }
          
        } # Condition if at least one FFP is computed and remaining savings ending...
          
        if(is.null(save.FFP.mtx.as)) {
          
          cat(warn(paste0('\n', space,  '[!] No FFP output matrices are saved')))
          
          # Save output to log
          Log_list <- c(Log_list, paste0('No FFP output matrices are saved'))
          
          
          
          ##
          ## [5] OUTPUT: GEOTIFF -----
          ##
          
        } else if (save.FFP.mtx.as == 'gtiff') {
          
          # Write the raster in a GTiff format # 
          cat(prog.mes(paste0('\n', bold('[5]'), ' Saving the FFP array as multiband GeoTiff')))
          
          # Compute the FFP full Extent # 
          FFP.Xdim <- seq(min(FFP.domain.ext), max(FFP.domain.ext), dx)
          FFP.Ydim <- seq(min(FFP.domain.ext), max(FFP.domain.ext), dx) 
          
          # Adding a dummy NA grid for NULL elements #
          for(i in FFP.err.indx_1) {
            FFP.ls[[i]]$fclim_2d <- array(data = NA, dim = c(length(FFP.Xdim), length(FFP.Ydim)))
          }
          
          FFP.lon <- FFP.Xdim + EC.tower.utm[1]
          FFP.lat <- FFP.Ydim + EC.tower.utm[2]
          
          # Extract the temporal domain
          Time <- as.character(format(unique(FFP.input.df.cur$'TIMESTAMP'), time_format))
          
          # Build the raster grid # 
          if (any(QC==0)){
          
          # Build the lat/lon grid
          FFP.grid <- data.frame(X=as.vector(t(FFP.ls[QC==0][[1]]$'x_2d')) + EC.tower.utm[1], 
                                 Y=as.vector(t(FFP.ls[QC==0][[1]]$'y_2d')) + EC.tower.utm[2])
          
          } else { # If no matrix exist use the extent and resolution parameters
            
            FFP.grid <- data.frame(X = rep(FFP.lon, ((FFP.domain*2)/dx)+1),
                                   Y = as.vector(sapply(FFP.lat, rep, ((FFP.domain*2)/dx)+1)))
          } 
          
          # List of FFP matrices # 
          FFP.mtx <- lapply(FFP.ls, function(x) x[['fclim_2d']])

          # Calculate scale and offset factors only if at least one not-NA matrix is returned
          if(any(QC != 1)) {
            
            # Global minimum
            min_val <- Reduce(min, lapply(FFP.mtx, function(mat) min(mat, na.rm=TRUE)))
            # Global maximum
            max_val <- Reduce(max, lapply(FFP.mtx, function(mat) max(mat, na.rm=TRUE)))
            # Compute the range
            range.FFP <- c(min_val, max_val)
            
            # Define int16 range
            type_min <- 0
            type_max <- 65534
            fillvalue <- 65535 # Automatically set by terra #
            
            # Compute scale factor #
            scale.factor <- (range.FFP[2]-range.FFP[1])/(type_max - type_min)
            
            # Compute offset #
            add.offset <- range.FFP[1]-scale.factor*type_min
            
            # Create a named vector #
            scale.offset <- c(sf = scale.factor, ofst = add.offset)

            # Scaling and clipping #
            FFP.mtx.scaled <- lapply(FFP.mtx, function(x) {
               tmp <- as.integer(round((x - scale.offset['ofst']) / scale.offset['sf']))
               tmp[tmp < type_min] <- type_min
               tmp[tmp > type_max] <- type_max
            
               tmp })

            # r_list <- lapply(FFP.mtx, function(x) {
            #    tmp <- rast(cbind(FFP.grid, as.vector(x)), type = "xyz", crs=crs(paste0('EPSG:', EPSG.code)))
            #    tmp })
                        
            # Free unused memory
            rm(FFP.mtx)
            rm(FFP.ls)
            gc()
            
            # Rasters list # 
            r_list <- lapply(FFP.mtx.scaled, function(x) {
              tmp <- rast(cbind(FFP.grid, as.vector(x)), type = "xyz", crs=crs(paste0('EPSG:', EPSG.code)))
              tmp })
            
                        
            # Clean the scaled matrix #
            rm(FFP.mtx.scaled)
            gc()
            
            # Create the multiband raster # 
            
            r <- rast(r_list)
            
            # Add raster names # 
            names(r) <- paste0("FFP_", 1:length(r_list))
            
            # Add scale and factor # Due to rounding, values can't be exactly the same 
            scoff(r) <- cbind(scale.offset['sf'], scale.offset['ofst'])
            
          } else {
            
            # Bind matrix and coordinates - build the spatraster
            r <- rast(cbind(FFP.grid, FFP.mtx), type = "xyz", crs=crs(paste0('EPSG:', EPSG.code)))
            
          }
          
          # Frequency variable # 
          if(do.full.climatology) {
            
            frequency.char <- "full average"
            
          } else if(do.daily.average) {
            
            frequency.char <- "daily average"
            
          } else {
            
            frequency.char <- "30min"
            
          }
          
          ## Add raster metadata, in any case ##
          #metags(r) <- c(Identifier=paste0("FFP at ", site.ID, " - ", as_date(FFP.input.df.cur$TIMESTAMP)[1]))
          metags(r) <- c(
            TITLE=paste0("Footprint estimation at the ", site.ID, " station"),
            CONTACT=paste0("ICOS Ecosystem Thematic Centre, email: info@icos-etc.eu"),
            CONVENTIONS="CF-1.9",
            CREATION_DATE=as.character(Sys.Date()),
            CONTRIBUTORS="ICOS Carbon Portal",
            CREATOR=paste0("Luca Di Fiore, ICOS ETC, CMCC Viterbo, orcid: 0000-0002-5698-4235", "; ", 
                           "Giacomo Nicolini, ICOS ETC, CMCC Viterbo, orcid: 0000-0002-9800-3351", "; ",
                           "Dario Papale, ICOS ETC, UNITUS, orcid: 0000-0001-5170-8648"),
            FREQUENCY=frequency.char,
            GEOSPATIAL_LAT_RESOLUTION=paste0(dx, " m"),
            GEOSPATIAL_LON_RESOLUTION=paste0(dx, " m"),
            INSTITUTION=paste0("CMCC Foundation - Euro-Mediterranean Center on Climate Change, IAFES Division, Viterbo, Italy", "; ", 
                               "Department for Innovation in Biological, Agro-Food and Forest Systems (DIBAF), University of Tuscia, Viterbo, Italy"),
            KEYWORDS="Flux footprints, Eddy covariance",
            LICENSE="CC-BY 4.0",
            PRODUCT_VERSION="1.0",
            #PROJECT="Integrated Carbon Observation System",
            REFERENCES=paste0("FFP model: Kljun et al. (2015), doi:10.5194/gmd‐8‐3695‐2015" ,"; ",
                              "Source code for generating the NetCDF file: https://github.com/icos-etc/FFP-utils"),
            #HISTORY=paste0("ICOS Ecosystem Thematic Centre", date(), sep=", "),
            SOURCE=paste0("PID code of input ICOS eddy covariance data: ", PID),
            SUBJECTS="Flux footprints, Eddy covariance",
            SUMMARY=paste0("The file contains the footprint estimation and the quality layer of the station ", 
                           site.ID, " at 30-minute temporal resolution."),
            VARIABLES="FFP",
            QC=paste0(QC, collapse =" "))
          
          # Create the path if don't exist # 
          ffp.GTiff.fpath <- paste0(FFP.output.dir, "/GTiff_files/")
          if(!dir.exists(ffp.GTiff.fpath)) {dir.create(ffp.GTiff.fpath)}
          
          # GTiff filename and time 
          if(all(do.full.climatology == F, do.daily.average==F))  {
            
            ffp.gtiffname <- paste0(ffp.GTiff.fpath, site.ID, "_FFP2D_", i.day, ".tiff")
            
            # Specify the time with the correct time zone 
            time_strings <- format(unique(FFP.input.df.cur$TIMESTAMP), format="%Y-%m-%d %H:%M:%S GMT")
            
            # Define the fuse for R (is NOT ISO)
            fuse_R <- paste0("Etc/GMT", ifelse(time_offset_hours >= 0, "-", "+"), abs(time_offset_hours))
            
            # Assign the time and the time zone # 
            time(r) <- as.POSIXct(time_strings, tz = fuse_R)
                        
            
          } else if(do.full.climatology == T) {
            
            ffp.gtiffname <- paste0(ffp.GTiff.fpath, site.ID, "_FFP2D_", "_FFP2D_full_climatology_", 
                                    first_day, '_', last_day, ".tiff")
            time(r) <- first_day
            
          } else if(do.daily.average == T) {
            
            ffp.gtiffname <- paste0(ffp.GTiff.fpath, site.ID, "_FFP2D_", "_FFP2D_daily_average_", i.day, ".tiff")
            time(r) <- as.Date(current_day_str)
            
          }
          
          # Write the GTiff file #
          writeRaster(r, filename = ffp.gtiffname, overwrite=T, progress=F, datatype="FLT4S", 
                      NAflag = as.numeric(fillvalue))
          
          # Save output message to log
          Log_list <- c(Log_list, '', paste0(#unique(format(FFP.input.df.cur$TIMESTAMP, '%Y-%m-%d')), ': ', 
                                             'FFP array saved as multiband GeoTiff'))
          
          cat(prog.mes('\n'))
          
          cat('\n**************************************************************************************')
          
          
          
          ##
          ## [6] OUTPUT: NETCDF -----
          ##
          
        } else if (save.FFP.mtx.as == 'nc') {
          
          
          # Create and write a projected netCDF file message
          cat(prog.mes(paste0('\n', bold('[6]'), #bold(current_day_str), ': ', 
                              ' Creating and writing the netCDF file')))
          
          ##
          ## [6.1] INPUT LIST PARAMETERS -----
          ##
          
          # List to store dimention/attributes inputs
          FFP.input.list <- list()
          
          # Define nc single matrix dimensions 
          if (any(QC==0)) { # Only the last valid time is considered 
            
            # store the X,Y dimensions of the FFP arrays (used for nc file dimension and variables)
            FFP.Xdim <- FFP.ls[which(QC == 0)][[1]]$'x_2d'[1, ]  # take the first row of x_2d, 
            FFP.Ydim <- FFP.ls[which(QC == 0)][[1]]$'y_2d'[, 1]  # take the first column of y_2d
            
          } else { # If no matrix exist use the extent and resolution parameters
            
            # store the X,Y dimensions of the FFP arrays (used for nc file dimension and variables)
            FFP.Xdim <- seq(min(FFP.domain.ext), max(FFP.domain.ext), dx)
            FFP.Ydim <- seq(min(FFP.domain.ext), max(FFP.domain.ext), dx) 
            
          } 
          
          # Store the UTM coordinates of the FFP arrays (used for nc file dimension and variables)
          FFP.lon <- FFP.Xdim + EC.tower.utm[1]
          FFP.lat <- FFP.Ydim + EC.tower.utm[2]
          
          # Create an index of the correct isopleth calculation
          Idx_Iso <- which(sapply(FFP.R.list, function(x) any(x %in% (which.FFP.R/100))))
          Idx_NoIso <- which(sapply(FFP.R.list, function(x) !any(x %in% (which.FFP.R/100))))
          
          # Choose of the right isopleth 
          Pos <- which(sort(FFP.R) %in% which.FFP.R)
          
          # Where the isopleth is not correctly calculated, add 3 zero vertices coords # 
          # Add 3 vertices where the chosen isopleth is NA (for tower center) directly in XR and YR 
          
          if(length(FFP.ls[Idx_NoIso]) > 0)  {
            
            for(i in Idx_NoIso) #FFP.sf[FFP.sf$Area==0, 'Index']
            {
              # Three nodes with zero coords in the second sublist
              FFP.ls[[i]]$xr[[Pos]] <- rep(0, 3)
              FFP.ls[[i]]$yr[[Pos]] <- rep(0, 3)
            
              # If noISO is 48: zero coords to all #   
            }
          } 
          
          
          # When no FFP matrix is calculated at all, zero coords to all # 
          if(length(FFP.ls) == length(FFP.err.indx_1))  {
            
            for(i in 1:length(FFP.ls)) #FFP.sf[FFP.sf$Area==0, 'Index']
            {
              # Three nodes with zero coords in the second sublist
              FFP.ls[[i]]$xr[[Pos]] <- rep(0, 3)
              FFP.ls[[i]]$yr[[Pos]] <- rep(0, 3)
              
            }
          }
          
          
          
          ##
          ## [6.2] FFP ISOPLETH PARAMETERS -----
          ##
          
          # Define parameters linked to isopleths 
          if(return.isopleth) {
            
            cat(prog.mes(paste0('\n', bold('[6.2]'), 
                                ' Processing the FFP output to store the ', which.FFP.R, 
                                '% isopleth as polygon')))
            
            # Area calculation on ST_Polygons # 
            FFP.sf <- data.frame(id=NULL, type=NULL, Area=NULL, Area_TOT=NULL, Index=NULL, nodes=NULL)
            
            if(length(Pos)==0)  {
              stop('FFP isopleth chosen to be returned is not in the isopleths list of the FFP function')
            }
            
            # Create a list of vectors with the coordinates in the right order according to the hole vector #
            nodes_x_tot <- c()
            nodes_y_tot <- c()
            
            
            for(i in 1:length(FFP.ls)) {
                
                DF <- data.frame(Lon=FFP.ls[[i]]$xr[[Pos]] + EC.tower.utm[1], 
                                 Lat=FFP.ls[[i]]$yr[[Pos]] + EC.tower.utm[2])
                
                # Group index generating # 
                pol_idx <- cumsum(is.na(DF$Lon))[!is.na(DF$Lon)] # IMPORTANT # Remove the NA 
                # Calculate the node sum of each type #  
                pol_idx_nodes=as.vector(table(pol_idx))
                
                # Associate coords for regular and holes polygons # 
                DF_idx=cbind(na.omit(DF), pol_idx)
                
                # Split the polygons #
                pol_list <- split(DF_idx, DF_idx$pol_idx)
    
                # Convex hull by polygons # 
                pol <- pol_list %>% map(function(x) 
                                          st_convex_hull(st_union(st_as_sf(na.omit(x), coords=c('Lon', 'Lat'), crs=EPSG.code))))
                pol_union <- do.call(c, pol)
                
                # Relationship matrix among all the polygons
                relation_matrix <- st_contains(pol_union, pol_union, sparse = FALSE)
                
                # Remove the diagonal
                diag(relation_matrix) <- FALSE
                
                # If at least one TRUE, that poligon is contained in something else # 
                is_hole <- colSums(relation_matrix) > 0

                
                nodes_x_list <- c()
                nodes_y_list <- c()
                
                for(j in seq_along(pol_list)) {
                  
                  # Estraiamo i vettori grezzi dal singolo elemento della lista
                  lon_raw <- pol_list[[j]]$Lon
                  lat_raw <- pol_list[[j]]$Lat
                  
                  # If the polygon is standard --> REV # Otherwise no # 
                  if(!is_hole[j]) {
                    
                    nodes_x_list[[j]] <- rev(lon_raw)
                    nodes_y_list[[j]] <- rev(lat_raw)
                    
                  } else {
                    
                    nodes_x_list[[j]] <- lon_raw
                    nodes_y_list[[j]] <- lat_raw
                    
                  }
                }
                
                # AREA calculation #
                
                pol_sf <- st_sf(
                  id = names(pol), 
                  type = ifelse(is_hole, "hole", "standard"),
                  nodes=pol_idx_nodes,
                  geometry = pol_union
                )

                pol_sf$Area=as.vector(st_area(pol_sf$geometry))
                
                # One value for each half-hour. Join full polygons and holes # 
                pol_sf$Area_TOT=sum(ifelse(pol_sf$type == "standard", pol_sf$Area, -pol_sf$Area))
                pol_sf=st_drop_geometry(pol_sf)
                
                
                # Get the polygon index #
                pol_sf$Index <- i
                
                # Join the i half-hour with the others # 
                FFP.sf <- rbind(FFP.sf, pol_sf)
                
                # Bind the coords vectors reversed # 
                nodes_x_tot=c(nodes_x_tot, unlist(nodes_x_list))
                nodes_y_tot=c(nodes_y_tot, unlist(nodes_y_list))
                
            }
            
            # Complex polygons warnings # 
            for(i in unique(FFP.sf$Index))  {
              
              tmp <- FFP.sf[FFP.sf$Index==i, ]
              
              if(nrow(tmp) > 1) {  # two polygons # 
                
                if(any(tmp$type %in% "hole")) { # Any hole? # 
                  
                  cat(warn(paste0('[!] ', bold(unique(FFP.input.df.cur$TIMESTAMP[tmp$Index])), ': ',
                                  which.FFP.R, '% isopleth has a hole', '\n')))

                  Log_list <- c(Log_list, '',
                                paste0('[!] ', unique(FFP.input.df.cur$TIMESTAMP[tmp$Index]), ': ',
                                       which.FFP.R, '% isopleth has a hole'), '')
                  
                } else  {  # Two distinct polygons #
                  
                  cat(warn(paste0('[!] ', bold(unique(FFP.input.df.cur$TIMESTAMP[tmp$Index])), ': ',
                                  which.FFP.R, '% isopleth has multiple sub-polygons', '\n')))
                  
                  Log_list <- c(Log_list, '',
                                paste0('[!] ', unique(FFP.input.df.cur$TIMESTAMP[tmp$Index]), ': ',
                                       which.FFP.R, '% isopleth has multiple sub-polygons'), '')
                }
              }
            }
            
            ## Input list definition (for dimensions and attributes)
            # Adding variables to the list
            FFP.input.list[['polygon_index']] <- unique(FFP.sf$Index)
            # Add the area # 
            FFP.input.list[['polygon_area']] <- unique(FFP.sf[c("Area_TOT", "Index")])[['Area_TOT']]
            
            # Insert time variable and char numbers for the date field in the polygon attribute table
            if(all(do.full.climatology == F, do.daily.average==F))  {

              FFP.input.list[['time']] <- as.numeric(format(unique(FFP.input.df.cur$'TIMESTAMP'), "%Y%m%d%H%M"))
              #FFP.input.list[['char_length']] <- nchar(FFP.input.list$time[1])

            } else if(do.full.climatology == T) {
              
              FFP.input.list[['time']] <- as.numeric(format(first_day, "%Y%m%d"))
              #FFP.input.list[['time']] <- as.character(paste0(first_day , '_', last_day))
              #FFP.input.list[['char_length']] <- nchar(FFP.input.list$time[1])
              
            } else if(do.daily.average == T) {
              
              FFP.input.list[['time']] <- as.numeric(format(current_day, "%Y%m%d"))
              #FFP.input.list[['char_length']] <- nchar(FFP.input.list$time[1])
              
            }
            
            # Selection of the chosen isoplete and remove the  NAs #
            FFP.iso <- lapply(FFP.ls, FUN=function(x) na.omit(x$xr[[Pos]]))
            
            # FFP.input.list build # 
            # Info collection # 
            
            # Real polygons (even zero ones) contained in the input data
            FFP.input.list[['polygon_number']] <- length(FFP.iso) # if half-hours, always 48
            
            # Total nodes for each obs. feature (polygon + eventual holes or parts)
            FFP.input.list[['nodes_number_polygon']] <- lengths(FFP.iso) 
            
            # Nodes for each single element # 
            FFP.input.list[['part_node_count']] <- FFP.sf$nodes
            
            # Add part number # 
            FFP.input.list[['part_number']] <- length(FFP.input.list$part_node_count)
            
            # Holes identifier (0: standard, 1: hole) #
            FFP.input.list[['interior_ring']] <- ifelse(FFP.sf$type=="standard", yes=0, no=1)
            
            # Total number of nodes # 
            FFP.input.list[['nodes_total']] <- sum(lengths(FFP.iso)) 
            
            # Plain coordinate vectors (adjusted for the holes, which are reversed) #
            FFP.input.list[['nodes_x']] <- nodes_x_tot
            FFP.input.list[['nodes_y']] <- nodes_y_tot
            
            # Remove the sf object 
            #rm(FFP.sf)
            
          } # Isopleth input phase ....... ENDING
          
          # Populating the FFP input list with additional information #  
          # Retrieving QC flag and adding it to the list
          FFP.input.list[['QC_Flag']] <- QC
          
          # Add EC tower coords to the input list 
          FFP.input.list[['EC_Tower_x']] <- EC.tower.utm[1]
          FFP.input.list[['EC_Tower_y']] <- EC.tower.utm[2]
          
          
          
          ##
          ## [6.3] FFP BOX COORDINATES -----
          ##
          
          cat(prog.mes(paste0('\n', bold('[6.3]'), 
                              ' Adding bounding box vertices coordinates in Lat/Lon')))
          
          # BB (buonding box) creation #
          x.min <- min(FFP.lon)# - (dx / 2)
          x.max <- max(FFP.lon)# + (dx / 2)
          y.min <- min(FFP.lat)# - (dx / 2)
          y.max <- max(FFP.lat)# + (dx / 2)
          
          # 4 anti-clockwise vertices starting from bottom left # 
          X_vertices <- c(x.min, x.max, x.max, x.min, x.min)
          Y_vertices <- c(y.min, y.min, y.max, y.max, y.min)
          
          # UTM spatial obkect
          BB_matrix <- cbind(X_vertices, Y_vertices)
          BB_polygon  <- st_sfc(st_polygon(list(BB_matrix)), crs=EPSG.code)
          
          # Project the polygon to WGS84 # 
          BB_polygon_WGS84 <- st_transform(BB_polygon, crs = 4326)
          
          # BB coordinates # BL - TL - TR - BR #
          BB_coords_WGS84 <- unique(as.data.frame(st_coordinates(BB_polygon_WGS84)[, c("X", "Y")]))
          
          # Convert to a plain vector # 
          BB_coords_WGS84_X <- unlist(as.vector(BB_coords_WGS84)[1], use.names = F)
          BB_coords_WGS84_Y <- unlist(as.vector(BB_coords_WGS84)[2], use.names = F)

          # Add the coordinates to the input list #
          FFP.input.list[['BB_Coords_X']] <- BB_coords_WGS84_X
          FFP.input.list[['BB_Coords_Y']] <- BB_coords_WGS84_Y
          
          
          ##
          ## [6.4] FFP DATA COMPRESSION -----
          ##
          
          # Adding a dummy NA grid for NULL elements #
          for(i in FFP.err.indx_1)  {
            FFP.ls[[i]]$fclim_2d <- array(data = NA, dim = c(length(FFP.Xdim), length(FFP.Ydim)))
          }
          
          # Extract the FFP array and apply the scale-offset factor # 
          #FFP.mtx <- lapply(FFP.ls, function(x) x[['fclim_2d']]/sum(x[['fclim_2d']], na.rm = T))
          FFP.mtx <- lapply(FFP.ls, function(x) x[['fclim_2d']])
          
          # Calculate scale and offset factors only if at least one not-NA matrix is returned
          if(any(FFP.input.list[['QC_Flag']]!=1)) {
            
            # Add a scaling factor for the FFP values (to reduce object size)
            # Global minimum
            min_val <- Reduce(function(a, b) min(a, b, na.rm = TRUE), 
                              lapply(FFP.mtx, function(mat) {
                                if (all(is.na(mat))) return(NULL) # Skip NA matrices
                                min(mat, na.rm = TRUE)}))
            
            # Global maximum
            max_val <- Reduce(function(a, b) max(a, b, na.rm = TRUE), 
                              lapply(FFP.mtx, function(mat) {
                                if (all(is.na(mat))) return(NULL) # Skip NA matrices
                                max(mat, na.rm = TRUE)}))

            # Compute the range
            range.FFP <- c(min_val, max_val)
            
            # Define int16 range and fillvalue
            type_min <- 0 
            type_max <- 65534
            fillvalue <- 65535
            
            # Compute scale factor #
            scale.factor <- (range.FFP[2]-range.FFP[1])/(type_max - type_min)
            
            # Compute offset # 
            add.offset <- range.FFP[1]-scale.factor*type_min
            
            # Create a named vector # 
            scale.offset <- c(sf = scale.factor, ofst = add.offset)
            
            # Create an array of FFP matrices, using scale and offset factors
            FFP.array <- array(NA_integer_, dim = c(length(FFP.Xdim), length(FFP.Ydim), length(FFP.mtx)))
            for (i in seq_along(FFP.mtx)) {
              FFP.array[,,i] <- round((FFP.mtx[[i]] - add.offset) / scale.factor)
            }

            # Data clipping # Prevents hidden rounding errors # 
            FFP.array[FFP.array < type_min] <- type_min
            FFP.array[FFP.array > type_max] <- type_max
            #FFP.array <- pmin(pmax(FFP.array, type_min), type_max)
            
          } else {
            
            # Extract only the array without calculating scale and offset factors
            FFP.array <- array(unlist(FFP.mtx, use.names = F), 
                               dim = c(length(FFP.Xdim), length(FFP.Ydim), length(FFP.mtx)))
            
          }
          
          # Free unused memory
          rm(FFP.mtx)
          rm(FFP.ls)
          gc()
          
          # Convert all NA's to the fillvalue used in the creation of nc file
          FFP.array[is.na(FFP.array)] <- fillvalue
          
          
          ##
          ## [6.5] NETCDF BUILD -----
          ##
          
          cat(prog.mes(paste0('\n', bold('[6.5]'), 
                              ' Saving the NetCDF file')))
          
          # . * Create the netCDF filename and folder path 
          # path if not exist and file name, set dname
          ffp.ncfpath <- paste0(FFP.output.dir, "/nc_files/")
          
          if(!dir.exists(ffp.ncfpath)) {dir.create(ffp.ncfpath)}
          
          # File ID and name #
          if(all(do.full.climatology == F, do.daily.average==F))  {
            ffp.ncfname <- paste0(ffp.ncfpath, "ICOSETC_", site.ID, "_FFP2D_", i.day, ".nc")
          } else if(do.full.climatology == T) {
            ffp.ncfname <- paste0(ffp.ncfpath, "ICOSETC_", site.ID, "_FFP2D_full_climatology_", 
                                  first_day, '_', last_day, ".nc")
          } else if(do.daily.average == T) {
            ffp.ncfname <- paste0(ffp.ncfpath, "ICOSETC_", site.ID, "_FFP2D_daily_average_", 
                                  i.day, ".nc")
          }

          ffp.dfname <- "FFP"
          
          
          
          ## [4.3] ## Definition of netCDF dimensions and variables
          
          # FFP grid dimension definition (x, y, time)
          FFP_dim_x <- ncdim_def(name='x', 
                                 units='m',
                                 longname="X coordinate of projection",
                                 vals=FFP.lon)
          
          FFP_dim_y <- ncdim_def(name='y', 
                                 units='m',
                                 longname="Y coordinate of projection",
                                 vals=FFP.lat)
          
          # Set the correct utc zone #
          if(any(do.full.climatology, do.daily.average))  {
            
            FFP_dim_time <- ncdim_def(name="time", 
                                      units=paste0("seconds since 1970-01-01 00:00:00 ", time_offset_string),
                                      longname=paste0("Time in seconds since 1970-01-01 00:00:00 ", time_offset_string),
                                      vals=as.integer(as.POSIXct(unique(cur.ts), format='%Y%m%d', tz='GMT')), 
                                      unlim=T, # Time is unlimited (but not for humans!)
                                      calendar="standard")
            
          } else  {
            
            FFP_dim_time <- ncdim_def(name="time", 
                                      units=paste0("seconds since 1970-01-01 00:00:00 ", time_offset_string),
                                      longname=paste0("Time in seconds since 1970-01-01 00:00:00 ", time_offset_string),
                                      vals=as.integer(as.POSIXct(unique(cur.ts), format='%Y%m%d%H%M', tz='GMT')), 
                                      unlim=T, # Time is unlimited (but not for humans!)
                                      calendar="standard")
          }
          
          # Dimension (size=2) for the EC tower coordinates 
          FFP_dim_ECoords <- ncdim_def(name='EC_Coords', 
                                       units='',
                                       longname="Projected (x and y, respectively) coordinates of EC tower",
                                       vals=1:2, 
                                       create_dimvar=F)
          
          # Dimension (size=4) for the BB coordinates 
          FFP_dim_BBcoords <- ncdim_def(name='BB_Coords', 
                                       units='',
                                       longname="WGS84 coordinates (Bottom Left - Bottom Right) of the matrix bounding box",
                                       vals=1:4, 
                                       create_dimvar=F)
          
          # Dimension (size=48) for the QC 
          FFP_dim_QC <- ncdim_def(name='QC_dim', 
                                  units='',
                                  longname="Quality flag of FFP matrix",
                                  vals=1:length(FFP.input.list$QC_Flag), 
                                  create_dimvar=F)
          
          # GRID variables
          # FFP grid #
          FFP_mtx30_def <- ncvar_def(name=ffp.dfname,
                                     units="m-2",
                                     dim=list(FFP_dim_x, FFP_dim_y, FFP_dim_time),
                                     missval=fillvalue,
                                     longname="Matrix of normalised 2D footprint values",
                                     prec="integer",
                                     compression=9, 
                                     chunksizes=c(length(FFP.lon),length(FFP.lat), 1), 
                                     shuffle=TRUE)
          
          # CRS #
          FFP_var_proj <- ncvar_def(name="UTM_coordinate_system",
                                    units='m',
                                    dim=NULL,
                                    missval=NULL,
                                    longname='CRS reference system (Universal Transverse Mercator)',
                                    prec="char")
          
          # Quality flag #
          FFP_var_QC <- ncvar_def(name="quality_flag",
                                  units='1',
                                  dim=FFP_dim_QC, 
                                  longname="Quality control on footprint calculation",
                                  missval=NULL,
                                  prec="integer")
          
          # EC Tower coordinates
          FFP_var_ECoord <- ncvar_def(name="EC_tower_coordinates",
                                      units="m",
                                      dim=FFP_dim_ECoords,
                                      missval=NULL,
                                      longname="Eddy covarinace tower coordinates in UTM",
                                      prec="float")
          
          # BB coordinates
          FFP_var_BBcoord_X <- ncvar_def(name="BB_coordinates_x",
                                      units="degrees_east",
                                      dim=FFP_dim_BBcoords,
                                      missval=NULL,
                                      longname="Bounding box coordinates in WGS84",
                                      prec="float")
          
          FFP_var_BBcoord_Y <- ncvar_def(name="BB_coordinates_y",
                                        units="degrees_north",
                                        dim=FFP_dim_BBcoords,
                                        missval=NULL,
                                        longname="Bounding box coordinates in WGS84",
                                        prec="float")
          
          # If isopleth should be returned  
          
          if(return.isopleth) {
            
            # nodes # Nodes number
            FFP_dim_node <- ncdim_def(name='nodes', 
                                      units='', ## If dimvar F, then units should be empty ## 
                                      longname="Number of nodes",
                                      vals=1:FFP.input.list$nodes_total, 
                                      create_dimvar=F)
            
            # instance # Number of polygons
            FFP_dim_inst <- ncdim_def(name='instance', 
                                      units='',
                                      longname="Number of polygons",
                                      vals=1:FFP.input.list$polygon_number, 
                                      create_dimvar=F)
            
            # Part # Total number of parts across all the obs
            FFP_dim_part <- ncdim_def(name='part', 
                                      units='', 
                                      longname="Total number of geometry parts",
                                      vals=1:FFP.input.list$part_number, 
                                      #vals=1:length(FFP.input.list$part_node_count), 
                                      create_dimvar=F)
            
            # Part node count # Node count for each part
            FFP_var_part_nodes <- ncvar_def(name="part_node_count",
                                            units='1',
                                            dim=FFP_dim_part, 
                                            missval=NULL,
                                            longname='Number of nodes for each polygon part',
                                            prec="integer")
            
            # Internal hole indicator 
            FFP_var_int_ring <- ncvar_def(name="interior_ring",
                                          units='1',
                                          dim=FFP_dim_part, 
                                          missval=NULL,
                                          longname='Interior ring indicator (0=exterior, 1=interior)',
                                          prec="integer")
            
            # n_char # Number of character of each string # DEPRECATED # 
            # Maximum dimension equal to the maximum length of a string
            # FFP_dim_nchar <- ncdim_def(name='n_char', 
            #                            units='',
            #                            longname="Length of each time string",
            #                            vals=1:FFP.input.list$char_length, 
            #                            create_dimvar=F)
            
            # Geometry container #
            FFP_var_geomc <- ncvar_def(name="geometry_container",
                                       units='1',
                                       dim=NULL,
                                       missval=NULL,
                                       longname='Container of the geometry',
                                       prec="char")
            
            # Nodes count #
            FFP_var_nodes <- ncvar_def(name="nodes_count",
                                       units='1',
                                       dim=FFP_dim_inst, 
                                       missval=NULL,
                                       longname='Number of nodes for each feature',
                                       prec="integer")
            
            # Node coordinates (X, Y) # 
            FFP_var_nodeX <- ncvar_def(name="nodes_x",
                                       units='m', 
                                       dim=FFP_dim_node, 
                                       missval=NULL,
                                       longname='All nodes coordinates',
                                       prec="float")
            
            FFP_var_nodeY <- ncvar_def(name="nodes_y",
                                       units='m', 
                                       dim=FFP_dim_node, 
                                       missval=NULL,
                                       longname='All nodes coordinates',
                                       prec="float")
            
            # Time # 
            FFP_var_time <- ncvar_def(name="timestamp_start",
                                      units='1',
                                      dim=FFP_dim_inst,
                                      missval=NULL,
                                      longname='Time of each observation',
                                      prec="double")
            
            # Quality flag # 
            FFP_var_QC <- ncvar_def(name="quality_flag",
                                    units='1',
                                    dim=FFP_dim_inst,
                                    longname="Quality control on footprint calculation",
                                    missval=NULL,
                                    prec="integer")
            
            # Polygon Area #
            FFP_var_area <- ncvar_def(name="polygon_area",
                                      units='m2',
                                      dim=FFP_dim_inst, 
                                      missval=NULL,
                                      longname='Area of the isoplete polygon in squared meters',
                                      prec="float")
            
            # polygon ID # 
            FFP_var_id <- ncvar_def(name="polygon_id",
                                    units='1',
                                    dim=FFP_dim_inst, 
                                    missval=NULL,
                                    longname=paste0('Id of each polygon referred to the ', which.FFP.R, '% isoplete'),
                                    prec="integer")
            
            
            # NC creation if isopleth is desired
            ffp.ncout <- nc_create(filename=ffp.ncfname, 
                                   vars=list(FFP_mtx30_def, 
                                             FFP_var_id,
                                             FFP_var_time,
                                             FFP_var_QC, 
                                             FFP_var_area, 
                                             FFP_var_geomc, 
                                             FFP_var_proj, 
                                             FFP_var_part_nodes, 
                                             FFP_var_int_ring,
                                             FFP_var_nodes, 
                                             FFP_var_nodeX, 
                                             FFP_var_nodeY, 
                                             FFP_var_ECoord, 
                                             FFP_var_BBcoord_X, 
                                             FFP_var_BBcoord_Y), 
                                   force_v4 = TRUE)
          } else {
            
            # NC creation if isopleth is not desired
            ffp.ncout <- nc_create(filename=ffp.ncfname,
                                   vars=list(FFP_mtx30_def,
                                             FFP_var_proj, 
                                             FFP_var_QC,
                                             FFP_var_ECoord, 
                                             FFP_var_BBcoord_X, 
                                             FFP_var_BBcoord_Y), 
                                   force_v4=TRUE)
            
          }
          
          #nc_close(ffp.ncout)
          
          
          ## [4.4] ## Store variables in the NC file
          
          # FFP_mtx30_def
          ncvar_put(nc=ffp.ncout, varid=FFP_mtx30_def, 
                    vals=FFP.array) ## Array
          
          # FFP_var_ECoord
          ncvar_put(nc=ffp.ncout, varid=FFP_var_ECoord, 
                    vals=c(FFP.input.list$EC_Tower_x, FFP.input.list$EC_Tower_y)) ## EC coords
          
          # FFP_var_BBcoord_X
          ncvar_put(nc=ffp.ncout, varid=FFP_var_BBcoord_X, 
                    vals=FFP.input.list$BB_Coords_X) ## BB coords X
          
          # FFP_var_BBcoord_Y
          ncvar_put(nc=ffp.ncout, varid=FFP_var_BBcoord_Y, 
                    vals=FFP.input.list$BB_Coords_Y) ## BB coords Y
          
          # FFP_var_QC 
          ncvar_put(nc=ffp.ncout, varid=FFP_var_QC, 
                    vals=FFP.input.list$QC_Flag) ## Flag
          
          
          # Additional variables related to isopleth
          if(return.isopleth) {  
            
            # FFP_var_id 
            ncvar_put(nc=ffp.ncout, varid=FFP_var_id, vals=FFP.input.list$polygon_index) ## Polygon ID
            
            # FFP_var_nodes 
            ncvar_put(nc=ffp.ncout, varid=FFP_var_nodes, vals=FFP.input.list$nodes_number_polygon) ## Number of nodes
            
            # FFP_var_nodeX
            ncvar_put(nc=ffp.ncout, varid=FFP_var_nodeX, vals=FFP.input.list$nodes_x) ## Nodes coord x
            
            # FFP_var_nodeY 
            ncvar_put(nc=ffp.ncout, varid=FFP_var_nodeY, vals=FFP.input.list$nodes_y) ## Nodes coord y
            
            # FFP_var_area
            ncvar_put(nc=ffp.ncout, varid=FFP_var_area, vals=FFP.input.list$polygon_area) ## Pol area
            
            # FFP_var_time 
            ncvar_put(nc=ffp.ncout, varid=FFP_var_time, vals=FFP.input.list$time) ## Observation time
            
            # FFP_var_part_nodes
            ncvar_put(nc=ffp.ncout, varid=FFP_var_part_nodes, vals=FFP.input.list$part_node_count) ## Nodes for each part
            
            # FFP_var_int_ring
            ncvar_put(nc=ffp.ncout, varid=FFP_var_int_ring, vals=FFP.input.list$interior_ring) ## Holes tag
            
          }
          
          
          # Grid mapping variable
          eastern_boundary <- UTM.zone * 6 - 180 # eastern boundary of the UTM zone, degrees
          western_boundary <- eastern_boundary - 6       # western boundary of the UTM zone, degrees
          central_meridian <- western_boundary + 3
          
          
          # Grid mapping name
          ncatt_put(nc=ffp.ncout, varid="UTM_coordinate_system", attname="grid_mapping_name", 
                    attval='transverse_mercator')
          
          # WKT: SPHEROID[“", , , ...]
          ncatt_put(nc=ffp.ncout, varid="UTM_coordinate_system", attname="semi_major_axis", 
                    attval=6378137)
          
          # WKT: SPHEROID[“", , , ...]
          ncatt_put(nc=ffp.ncout, varid="UTM_coordinate_system", attname="inverse_flattening", 
                    attval=298.257223563)
          
          # WKT: PRIMEM[“", , ...]
          ncatt_put(nc=ffp.ncout, varid="UTM_coordinate_system", attname="longitude_of_prime_meridian", 
                    attval=0)
          
          # WKT: PARAMETER[“Latitude of natural origin”, ]
          ncatt_put(nc=ffp.ncout, varid="UTM_coordinate_system", attname="latitude_of_projection_origin", 
                    attval=0)
          
          # WKT: PARAMETER[“Longitude of natural origin”, ]
          ncatt_put(nc=ffp.ncout, varid="UTM_coordinate_system", attname="longitude_of_central_meridian", 
                    attval=central_meridian)
          
          # WKT: PARAMETER[“Scale factor at natural origin”, ]
          ncatt_put(nc=ffp.ncout, varid="UTM_coordinate_system", attname="scale_factor_at_central_meridian", 
                    attval=0.9996)
          
          # WKT: PARAMETER[“False easting”, ]
          ncatt_put(nc=ffp.ncout, varid="UTM_coordinate_system", attname="false_easting", 
                    attval=500000)
          
          # WKT: PARAMETER[“False northing”, ]
          ncatt_put(nc=ffp.ncout, varid="UTM_coordinate_system", attname="false_northing", 
                    attval=0)
          
          # Add WKT (2) code using sf 
          ncatt_put(nc=ffp.ncout, varid="UTM_coordinate_system", attname="crs_wkt", 
                    attval=sf::st_crs(paste0('epsg:', EPSG.code))[['wkt']])
          
          
          # FFP_mtx30_def
          ncatt_put(nc=ffp.ncout, varid=ffp.dfname, attname="grid_mapping", 
                    attval='UTM_coordinate_system')
          ncatt_put(nc=ffp.ncout, varid=ffp.dfname, attname="coordinates", 
                    attval='x y')
          
          
          # Name of xy 
          ncatt_put(ffp.ncout, varid="x", attname="standard_name", 
                    attval="projection_x_coordinate")
          ncatt_put(ffp.ncout, varid="y", attname="standard_name", 
                    attval="projection_y_coordinate")
          
          if(any(FFP.input.list[['QC_Flag']]!=1)) {
            ncatt_put(nc=ffp.ncout, varid=ffp.dfname, attname="scale_factor", prec = "float",
                      attval=scale.offset['sf'])
            ncatt_put(nc=ffp.ncout, varid = ffp.dfname, attname="add_offset", prec = "float",
                      attval=scale.offset['ofst'])
          }
          
          
          # QC # 
          ncatt_put(ffp.ncout, varid="quality_flag", attname="standard_name",
                    attval='quality_flag')
          
          
          if(return.isopleth) { 
            
            
            ## Geometry container ## 
            ncatt_put(nc=ffp.ncout, varid="geometry_container", attname="geometry_type", 
                      attval='polygon')
            ncatt_put(nc=ffp.ncout, varid="geometry_container", attname="node_count", 
                      attval='nodes_count')
            ncatt_put(nc=ffp.ncout, varid="geometry_container", attname="node_coordinates", 
                      attval='nodes_x nodes_y')
            ncatt_put(nc=ffp.ncout, varid="geometry_container", attname="grid_mapping", 
                      attval='UTM_coordinate_system')
            ncatt_put(nc=ffp.ncout, varid="geometry_container", attname="part_node_count",  
                      attval='part_node_count')
            ncatt_put(nc=ffp.ncout, varid="geometry_container", attname="interior_ring",  
                      attval='interior_ring')
            
            
            # QC # 
            ncatt_put(ffp.ncout, varid="quality_flag", attname="standard_name",
                      attval='quality_flag')
            ncatt_put(ffp.ncout, varid="quality_flag", attname="grid_mapping",
                      attval='UTM_coordinate_system')
            ncatt_put(ffp.ncout, varid="quality_flag", attname="geometry",
                      attval='geometry_container')
            
            
            # polygon_id # Timeseries identifier
            ncatt_put(ffp.ncout, varid="polygon_id", attname="grid_mapping",
                      attval='UTM_coordinate_system')
            ncatt_put(ffp.ncout, varid="polygon_id", attname="geometry",
                      attval='geometry_container')
            ncatt_put(ffp.ncout, varid="polygon_id", attname="cf_role",
                      attval='timeseries_id')
            
            
            # Name of xy 
            ncatt_put(ffp.ncout, varid="x", attname="standard_name", 
                      attval="projection_x_coordinate")
            ncatt_put(ffp.ncout, varid="y", attname="standard_name", 
                      attval="projection_y_coordinate")
            
            
            # FFP_var_time 
            ncatt_put(ffp.ncout, varid = "timestamp_start", attname = "grid_mapping",
                      attval='UTM_coordinate_system')
            ncatt_put(ffp.ncout, varid = "timestamp_start", attname = "geometry",
                      attval='geometry_container')
            
            
            # FFP_var_nodes 
            ncatt_put(ffp.ncout, varid="nodes_count", attname="long_name",
                      attval="number of nodes for each polygon feature")
            
            
            # FFP_var_part_nodes
            ncatt_put(ffp.ncout, varid="part_node_count", attname="long_name",
                      attval="number of nodes for each geometry part")
            
            
            # FFP_var_int_ring
            ncatt_put(ffp.ncout, varid="interior_ring", attname="long_name",
                      attval="interior ring indicator")
            
            
            # FFP_var_nodeX
            ncatt_put(ffp.ncout, varid="nodes_x", attname="axis",
                      attval='X')
            ncatt_put(ffp.ncout, varid="nodes_x", attname="units",
                      attval='m')
            #ncatt_put(ffp.ncout, varid="nodes_x", attname="geometry",  # <-- Come visto prima
            #          attval='geometry_container')
            ncatt_put(ffp.ncout, varid=FFP_var_nodeX, attname="grid_mapping", 
                      attval='UTM_coordinate_system')
            
            
            # FFP_var_nodeY
            ncatt_put(ffp.ncout, varid="nodes_y", attname="axis",
                      attval='Y')
            ncatt_put(ffp.ncout, varid="nodes_y", attname="units",
                      attval='m')
            #ncatt_put(ffp.ncout, varid="nodes_y", attname="geometry",  # <-- Come visto prima
            #          attval='geometry_container')
            ncatt_put(ffp.ncout, varid=FFP_var_nodeY, attname="grid_mapping", 
                      attval='UTM_coordinate_system')
            
            
            
            # FFP_var_area
            ncatt_put(ffp.ncout, varid="polygon_area", attname="grid_mapping",
                      attval='UTM_coordinate_system')
            ncatt_put(ffp.ncout, varid="polygon_area", attname="geometry",
                      attval='geometry_container')
            
            
            # Feature type #
            ncatt_put(ffp.ncout, 0, "featureType", 
                      "timeSeries")
            
            
            # THREDDS type #
            ncatt_put(ffp.ncout, 0, "cdm_data_type", 
                      "Station")
            
          }
          
          
          ## [4.5] ## Store global attributes
          
          # title
          ncatt_put(ffp.ncout, 0, "title", paste0("Footprint estimation at the ", site.ID, " station"))
          
          # Contact
          ncatt_put(ffp.ncout, 0, "contact", 
                    paste0("ICOS Ecosystem Thematic Centre, email: info@icos-etc.eu"))
          
          # Contributors
          ncatt_put(ffp.ncout, 0, "contributors", 
                    "ICOS Carbon Portal")
          
          # Conventions
          ncatt_put(ffp.ncout, 0, "Conventions", 
                    "CF-1.9")
          
          # creation_date
          ncatt_put(ffp.ncout, 0, "creation_date", 
                    as.character(Sys.Date()))
          
          # creator #
          ncatt_put(ffp.ncout, 0, "creator", 
                    paste0("Luca Di Fiore, ICOS ETC, CMCC Viterbo, orcid: 0000-0002-5698-4235", "; ", 
                           "Giacomo Nicolini, ICOS ETC, CMCC Viterbo, orcid: 0000-0002-9800-3351", "; ",
                           "Dario Papale, ICOS ETC, UNITUS, orcid: 0000-0001-5170-8648"))
          
          # crs
          if(FFP.input.list$EC_Tower_y>0) {
            ncatt_put(ffp.ncout, 0, "crs", paste0("WGS84/UTM", UTM.zone, 'N'))
          } else { 
            ncatt_put(ffp.ncout, 0, "crs", paste0("WGS84/UTM", UTM.zone, 'S'))
          }
          
          # Frequency variable # 
          if(do.full.climatology) {
            
            frequency.char <- "full average"
            
          } else if(do.daily.average) {
            
            frequency.char <- "daily average"
            
          } else {
            
            frequency.char <- "30min"
            
          }
          
          # frequency
          ncatt_put(ffp.ncout, 0, "frequency", 
                    frequency.char)
          
          # geospatial_lat_resolution
          ncatt_put(ffp.ncout, 0, "geospatial_lat_resolution", 
                    paste0(dx, " m"))
          
          # geospatial_lon_resolution
          ncatt_put(ffp.ncout, 0, "geospatial_lon_resolution", 
                    paste0(dx, " m"))
          
          # geospatial_vertical_resolution
          #ncatt_put(ffp.ncout, 0, "geospatial_vertical_resolution", 
          #          "To define...")
          
          # institution
          ncatt_put(ffp.ncout, 0, "institution", 
                    paste0("CMCC Foundation - Euro-Mediterranean Center on Climate Change, IAFES Division, Viterbo, Italy", "; ", 
                           "Department for Innovation in Biological, Agro-Food and Forest Systems (DIBAF), University of Tuscia, Viterbo, Italy"))
                    
          # keywords
          ncatt_put(ffp.ncout, 0, "keywords", 
                    "Flux footprints, Eddy covariance")
          
          # license
          ncatt_put(ffp.ncout, 0, "license", 
                    "CC-BY 4.0")
          
          # product_version
          ncatt_put(ffp.ncout, 0, "product_version", 
                    "1.0")
          
          # references
          ncatt_put(ffp.ncout, 0, "references", 
                    paste0("FFP model: Kljun et al. (2015), doi:10.5194/gmd‐8‐3695‐2015" ,"; ",
                           "Source code for generating the NetCDF file: https://github.com/icos-etc/FFP-utils"))
          
          # source
          ncatt_put(ffp.ncout, 0, "source", 
                    paste0("PID code of input ICOS eddy covariance data: ", PID)) 
          
          # Additional metadata ## 
          ncatt_put(ffp.ncout, 0, "Subjects", 
                    "Flux footprints, Eddy covariance")
          
          if(return.isopleth) { 
            # summary
            ncatt_put(ffp.ncout, 0, "summary", 
                      paste0("The file contains the footprint estimation and the quality layer of the station ", 
                             site.ID, " at 30-minute temporal resolution. Additionally, a geometry containing the ", which.FFP.R, "% isopleth boundary is provided."))
          } else {
            
            # summary
            ncatt_put(ffp.ncout, 0, "summary", 
                      paste0("The file contains the footprint estimation and the quality layer of the station ", 
                             site.ID, " at 30-minute temporal resolution."))
            
          }
          
          # ## history
          # ncatt_put(ffp.ncout, 0, "history",
          #           paste("ICOS Ecosystem Thematic Centre", date(), sep=", "))
          
          # FundingReference
          #ncatt_put(ffp.ncout, 0, "FundingReference",
          #          "To define...")
          
          # Documentation of the model
          #ncatt_put(ffp.ncout, 0, "Documentation_of_the_model",
          #          "To define...")
          
          if(return.isopleth) { 
            
            # Variables
            ncatt_put(ffp.ncout, 0, "variables",
                      paste0("FFP, QC, ", which.FFP.R, "% isopleth area, ", 
                             "EC tower coordinates, bounding box coordinates"))
            
          } else {
            
            # Variables
            ncatt_put(ffp.ncout, 0, "variables",
                      paste0("FFP, QC, EC tower coordinates, bounding box coordinates"))
          }
          
          
          ## [4.6] ## Close the file, write it on the disk
          
          nc_close(ffp.ncout)
          
          
          ##
          ## [7] RAST PLOTS -----
          ##
          
          if(save.ncplot) {
            
            cat(prog.mes(paste0('\n', bold('[7]'), ' Saving the NetCDF sample plot')))
            
            ffp.ncplot.fpath <- paste0(FFP.output.dir, "/nc_files/", 'nc_plots/')
            if(!dir.exists(ffp.ncplot.fpath)) {dir.create(ffp.ncplot.fpath)}
          
            # Set the filename #  
            if(do.full.climatology) {
                
              ffp.ncplot.fname <- paste0(ffp.ncplot.fpath, site.ID, "_FFP2D_", 
                                           first_day, "_", last_day, "_full_climatology", '.jpeg')
              
            } else if(do.daily.average) {
              
                ffp.ncplot.fname <- paste0(ffp.ncplot.fpath, site.ID, "_FFP2D_", i.day, "_daily_average", '.jpeg')
              
            } else {
                
                ffp.ncplot.fname <- paste0(ffp.ncplot.fpath, site.ID, "_FFP2D_", i.day, '.jpeg')
                
            }
            
            # Save the image # 
            jpeg(ffp.ncplot.fname, width = 20, height = 13, units = 'cm', res = 300)
            
            if(any(do.full.climatology, do.daily.average))
            {
              # Plot the only layer available #
              terra::plot(terra::rast(ffp.ncfname))
            
              } else {
              
              # Sampling of 16 layers (including missing ones) and plotting #
              terra::plot(terra::rast(ffp.ncfname)[[sort(sample(1:nlyr(terra::rast(ffp.ncfname)), 16))]])  
            }
            
            dev.off()
            
          }
          
          #cat(prog.mes(' done.\n'))
          cat(prog.mes('\n'))
          cat('\n**************************************************************************************')
          
          # Save output to log
          Log_list <- c(Log_list, '', paste0(current_day_str, ': ', 
                                             'FFP matrices saved as netCDF'))
          
        } # NETcdf file creation chunk ending
        
      } # Conditional DF validity check ending
      
    } # Days loop ending
    
    
    
    ##
    ## [8] LOG -----
    ##
    
    # Save messages to a unique log file
    if(save.log)  {
      
      cat(prog.mes(paste0('\n', bold('[8]'), ' Saving the log file\n')))
      
      ffp.log.fpath <- paste0(FFP.output.dir, "/log/")
      if(!dir.exists(ffp.log.fpath)) {dir.create(ffp.log.fpath)}
      
      # Create log file
      if(do.full.climatology) {
        
        log <- file.path(ffp.log.fpath, paste0(site.ID, '_log_', 
                                               first_day, "_", last_day, "_full_climatology",  ".log"))
      } else if(do.daily.average) {
        
        log <- file.path(ffp.log.fpath, paste0(site.ID, '_log_', 
                                               as_date(min(FFP.input.df$'TIMESTAMP')), '_', 
                                               as_date(max(FFP.input.df$'TIMESTAMP')), "_daily_average", ".log"))
      } else {
        
        log <- file.path(ffp.log.fpath, paste0(site.ID, '_log_', 
                                               as_date(min(FFP.input.df$'TIMESTAMP')), '_', 
                                               as_date(max(FFP.input.df$'TIMESTAMP')), ".log"))
      }

      # Open log
      lf <- logr::log_open(log)
      
      # Write messages saved to log file 
      invisible(sapply(Log_list, function(x) logr::log_print(x, console = F, hide_notes = T, 
                                                             blank_after = F, msg = F)))
      
      # Close log
      logr::log_close()
      
    } # Log conditional ending
    
  } # DF, Coords and site.ID validity ending...
  
} # FFP function ending