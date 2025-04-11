## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
## -------------------- FOOTPRINT MODEL RUNNING AND RESULTS SAVING AS CSV OR NC FILE ---------------------- ##
## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##

doFFP=function(FFP.input.df=NULL,        # input dataframe
               EC.tower.coords=NULL,     # Ec tower coordinates (Lat, Lon)
               site.ID=NULL,             # site.ID
               PID=NULL,                 # site PID
               FFP.function.path="local_functions\\calc_footprint_FFP_climatology.R", # If specified, full path of FFP function, otherwise local path
               do.parallel=FALSE,        # parallelize the calculation of footprint
               max.gap = 4,              # maximum number of consecutive half-hours to fill with LOCF
               FFP.domain = 1000,        # Domain size as an array of (xmin xmax ymin ymax) [m] *1.
               dx = 1,                   # Cell size of domain [m] (default is dx = dy = 1 m)
               FFP.R = c(50,70,80,90),   # Isopleths
               which.FFP.R = 70,         # Which isopleth save
               drop.trivial.srcs = TRUE, # whether to remove point sources with a trivial contribution
               save.FFP.mtx.as='nc',     # csv, nc (NetCDF), NULL (no save)
               save.plot.FFP.mtx = FALSE,# Whether to save a plot of the FFP2D (matrix and isoplethes, not the spatial polygons)
               return.isopleth=TRUE,     # Return chosen isopleth in the netCDF file  
               save.log=TRUE,            # Save messages files to log
               save.ncplot=TRUE          # Save a sample plot of nc the file
)     
{

  
  # Install pacman
  if (!require("pacman", quietly = T)) {install.packages("pacman")} 
  # Install EBImage
  if(!require("EBImage", quietly = T)) {pacman::p_load(BiocManager); BiocManager::install("EBImage")}
  # Load or install required packages via pacman
  pacman::p_load(crayon, data.table, sf, ncdf4, terra, ggplot2, logr, purrr, dplyr, readr, lubridate, EBImage)

  # Define message themes
  error <- crayon::red
  warn <- crayon::yellow
  note <- crayon::silver
  prog.mes <- crayon::cyan
  
  # Customize domain 
  FFP.domain.ext <- FFP.domain * c(-1, 1,-1, 1)
  
  # Load local functions if they aren't already loaded 
  if(!exists('calc_footprint_FFP_climatology'))   {source(FFP.function.path)}
  
  # If DF is character (and then csv path), read it
  if(is.character(FFP.input.df)) 
  {
    FFP.input.df=read_csv(FFP.input.df, col_types = cols())
    FFP.input.df$'TIMESTAMP' = with_tz(FFP.input.df$'TIMESTAMP', tz ='GMT') # Set the correct timezone
  }
  
  # INPUTS AVAILABILITY CHECK
  if(is.null(FFP.input.df))
  {
    cat(error(paste0("\nERROR: No FFP input dataframe is specified for the station. FFP calculation is skipped!\n")))
    
  } else if(all(is.null(site.ID), !(any(colnames(FFP.input.df) %in% c("site_id"))))) 
  {
    cat(error(paste0("\nERROR: No site ID is specified for the station. FFP calculation is skipped!\n")))
    
  } else if(all(is.null(EC.tower.coords),
                !(any(colnames(FFP.input.df) == "lat")))) 
  {
    cat(error(paste0("\nERROR: No latitude is specified for the station. FFP calculation is skipped!\n")))
    
  } else if(all(is.null(EC.tower.coords),
                !(any(colnames(FFP.input.df) == "lon")))) 
  {
    cat(error(paste0("\nERROR: No longitude is specified for the station. FFP calculation is skipped!\n")))
    
  } else {
    
    #-- Dataframe manipulation --#
    
    # Move site id, coordinates and PID outside the DF, if necessary....
    # EC tower coordinates
    if(is.null(EC.tower.coords)) 
    {
      EC.tower.coords=c(unique(FFP.input.df$lat), unique(FFP.input.df$lon))
      
      # Remove them from DF
      FFP.input.df$lat=NULL
      FFP.input.df$lon=NULL
    }
    
    # Site id
    if(is.null(site.ID)) 
    {
      site.ID=unique(FFP.input.df$site_id)
      
      # Remove it from DF
      FFP.input.df$site_id=NULL
    }
    
    # PID
    if(all(is.null(PID), colnames(FFP.input.df) %in% 'PID'))
    {
      PID=unique(FFP.input.df$PID)
      
      # Remove it from DF
      FFP.input.df$PID=NULL
    }
    
    # Check for NAs in the variables created, and in case stop the function # 
    stopifnot('Site ID is NA' = !is.na(site.ID), 
              'Tower coordinates are NAs' = !is.na(EC.tower.coords))
    
    # Assign coordinates to a variable
    latitude = EC.tower.coords[1]
    longitude = EC.tower.coords[2]
    
    # Check for the coordinates format (decimal)
    stopifnot(is.numeric(latitude), is.numeric(longitude), latitude < 90 & latitude > -90)
    
    UTM.zone=dplyr::case_when( 
      latitude>56 & latitude<64 & longitude>3 & longitude<12 ~ 32,
      latitude>72 & latitude<84 & longitude>9 & longitude<21 ~ 33,
      latitude>72 & latitude<84 & longitude>21 & longitude<33 ~ 35,
      is.numeric(latitude) & is.numeric(longitude) ~ (floor((longitude + 180)/6) %% 60) + 1)
    
    EPSG.code=dplyr::case_when(
      latitude>0 & latitude<84 & UTM.zone<10 ~ paste0(326, 0, UTM.zone) %>% as.numeric(),
      latitude>0 & latitude<84 & UTM.zone>10 ~ paste0(326, UTM.zone) %>% as.numeric(),
      latitude<0 & latitude>-80 & UTM.zone<10 ~ paste0(327, 0, UTM.zone) %>% as.numeric(),
      latitude<0 & latitude>-80 & UTM.zone>10 ~ paste0(327, UTM.zone) %>% as.numeric(),
      latitude > 84 & latitude < 90 ~ 32661,
      latitude < -80 & latitude > -90 ~ 32761)
    
    EC.tower=sf::st_as_sf(x = data.frame(lat=EC.tower.coords[1], lon=EC.tower.coords[2]), 
                          coords = c('lon', 'lat'), crs=4326)
    
    # EC coordinates in UTM
    EC.tower.utm <- st_transform(EC.tower, EPSG.code) %>% st_coordinates() %>% as.vector()
    
    # Create main output directory (and eventually site directory)
    Site_dir=paste0(getwd(), '\\', site.ID)
    FFP.output.dir=paste0(getwd(), '\\', site.ID, '\\', 'Output')
    
    if(!dir.exists(Site_dir)) {dir.create(Site_dir)}
    if(!dir.exists(FFP.output.dir)) {dir.create(FFP.output.dir)}
    
    # Compute footprints
    cat('\n**************************************************************************************')
    cat(prog.mes(paste0('\nComputing the footprints for the ', bold(site.ID), ' station.')))
    cat('\n**************************************************************************************\n')
    
    # FFP input data manipulation -----------------------------------------------
    cat(prog.mes('\nFFP input data manipulation (gap-filling)'))
    
    # NON sensitive variables ********
    # These variables are filled with LOCF independently from the gap length(s)
    
    is.na_nsv=lapply(FFP.input.df[c('hm', 'hc', 'd', 'z0', 'zm')], function(x) sum(is.na(x)))
    
    FFP.input.df[c('hm', 'hc', 'd', 'z0', 'zm')]=
      as.data.frame(apply(FFP.input.df[c('hm', 'hc', 'd', 'z0', 'zm')], 2,
                          function(x) data.table::nafill(x, type = 'locf')))
    
    invisible(sapply(names(is.na_nsv), function(x) cat(note(paste0('\n- ', is.na_nsv[[x]], 
                                                                   ' missing values of ', bold(x), 
                                                                   ' have been filled using LOCF.')))))
    
    # Preliminary part ---- save log file #
    # Save data to a log file
    Log_list=list()
    
    Log_list[[1]]=paste0('###### Computing the footprints for the ', site.ID, ' station ######')
    Log_list[[2]]='' # Blank space
    Log_list[[3]]='------ FFP input data manipulation (gap-filling) ------'
    Log_list[[4]]='' # Blank space
    Log_list=c(Log_list,sapply(names(is.na_nsv), function(x) 
      paste0('- ', is.na_nsv[[x]], ' missing values of ', x, ' have been filled using LOCF.'), USE.NAMES = F), 
      '')
    
    # Sensitive variables ********
    # These variables are filled with LOCF only if gaps are shorter than 2 hours
    
    for(i in c('umean', 'ol', 'sigmav', 'ustar', 'wind_dir'))
      
    {
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
      Log_list=c(Log_list, paste0('- ', as.character(isna), ' missing values of ', i,
                                  ' have been found, but only ', as.character(length(isna.f)), 
                                  ' have been filled using LOCF lasting for less than 2 hours'))
    }
    
    cat('\n')
    cat('\n**************************************************************************************')
    
    # Single FFP ################################################################
    
    # Create a days factor to split the dataset into daily chunks
    if (!any(grepl('POSIX', class(FFP.input.df$'TIMESTAMP')))) 
    { # if timestamp is reported as ISO (character), convert it to a time object
      FFP.input.df$'TIMESTAMP' <- as_datetime(FFP.input.df$'TIMESTAMP', tz ='GMT')
    } 
    
    days.factor <- as.factor(format(FFP.input.df$'TIMESTAMP', '%Y%m%d'))
    
    i.day <- 1
    
    for (i.day in 1:length(levels(days.factor))) { 
      
      # Split the input into days chunks and grab the current one
      FFP.input.df.cur <- data.frame() 
      FFP.input.df.cur <- as.data.frame(split(FFP.input.df, days.factor)[[levels(days.factor)[i.day]]])
      
      # Calculation start message
      cat(prog.mes(paste0('\n', bold(unique(format(FFP.input.df.cur$TIMESTAMP, '%Y-%m-%d'))), ': ', 
                          'Start of calculations\n')))
      
      # Quantify survived NAs from input data (don't drop them!)
      rows.na <- nrow(FFP.input.df.cur) - nrow(na.omit(FFP.input.df.cur))
      cat(warn(paste0('\n[NOTE] ', bold(unique(format(FFP.input.df.cur$TIMESTAMP, '%Y-%m-%d'))), ': ', 
                      bold(as.character(rows.na)), ' rows ','(',as.character(round(rows.na/nrow(FFP.input.df.cur)*100)),'% of data) are data gaps that cannot be filled.\n')))
      
      # Add messages to log
      Log_list=c(Log_list, '', '', 
                 paste0('------ ', unique(format(FFP.input.df.cur$TIMESTAMP, '%Y-%m-%d')), ': ', 
                        'start of calculations ------'), '')
      
      Log_list=c(Log_list, paste0('[NOTE] ', unique(format(FFP.input.df.cur$TIMESTAMP, '%Y-%m-%d')), ': ', 
                                  as.character(rows.na), ' rows ','(',as.character(round(rows.na/nrow(FFP.input.df.cur)*100)),
                                  '% of data) are data gaps that cannot be filled.'), '')
      
      # Check (here) if there are all NA values and remove them   
      if(all(complete.cases(FFP.input.df.cur[4:ncol(FFP.input.df.cur)])==F)) 
        
      {
        
        cat(warn(paste0('\n [!] ', bold(unique(format(FFP.input.df.cur$TIMESTAMP, '%Y-%m-%d'))), ': ',  
                        'no input data is available - model is not computed - no file is produced', '\n')))
        
        # Add messages to log
        Log_list=c(Log_list, '', paste0('[!] ', unique(format(FFP.input.df.cur$TIMESTAMP, '%Y-%m-%d')), ': ',  
                                       'no input data is available - model is not computed - no file is produced'), '')
        
      } else 
        
      {
        
        ## -1- ## Footprint loop over the half-hours ---------------
        
        cat(prog.mes(paste0('\n', bold(unique(format(FFP.input.df.cur$TIMESTAMP, '%Y-%m-%d'))), ': ', 
                            'Computing the 2D footprint for all the timestamps at once')))
        
        if(do.parallel)
          
          {
           # FFP calculation (future parallel)
           pacman::p_load(furrr)
           plan(multisession, workers = 3) # Set a relatively low number of workers to avoid PC crashing
          
           FFP.ls=future_map(.x = split(FFP.input.df.cur, FFP.input.df.cur$TIMESTAMP), # Split the input to create a list
                             .f = function(x) 
                              tryCatch(expr = calc_footprint_FFP_climatology(
                                zm = x$zm,                  
                                z0 = x$z0,                 
                                umean = x$umean,            
                                h = x$PBL,                  
                                ol = x$ol,                  
                                sigmav = x$sigmav,         
                                ustar = x$ustar,            
                                wind_dir = x$wind_dir,     
                                
                                # optional input
                                domain = FFP.domain.ext, dx = dx, dy = dx,                 
                                r = FFP.R, rslayer =  1, smooth_data = 1,                            
                                pulse = 0),
                                error = function(e) NULL), .progress = T)

          # Free unused memory
          gc()
          
        } else 
        
        {
          # Function not parallelized
          FFP.ls=lapply(split(FFP.input.df.cur, FFP.input.df.cur$TIMESTAMP), # Split the input to create a list
                        FUN = function(x) 
                          tryCatch(expr = calc_footprint_FFP_climatology( 
                            zm = x$zm,                  
                            z0 = x$z0,                 
                            umean = x$umean,            
                            h = x$PBL,                  
                            ol = x$ol,                  
                            sigmav = x$sigmav,         
                            ustar = x$ustar,            
                            wind_dir = x$wind_dir,     
                            
                            # optional input
                            domain = FFP.domain.ext, dx = dx, dy = dx,                 
                            r = FFP.R, rslayer =  1, smooth_data = 1,                            
                            pulse = 0),
                            error = function(e) NULL))
  
          }
        
        # Names are not necessary
        names(FFP.ls)=NULL
        
        # Format time. It needs to be always a valid value (no gaps allowed in the final product)
        cur.ts <- format(FFP.input.df.cur$'TIMESTAMP', '%Y%m%d%H%M')
        
        ## Extract FFP.R list ## 
        FFP.R.list=lapply(FFP.ls, function(x) x[['r']])
        
        # Computing error index # Index 1: no FFP matrix
        FFP.err.indx_1 <- which(sapply(FFP.R.list, function(x) ifelse(is.null(x), yes = T, no=F)))

        # Index 2: FFP matrix yes, but not the higher FFP R (matrix is thus not clipped)
        FFP.err.indx_2=which(sapply(FFP.R.list, function(x) 
          ifelse(all(sum(is.na(x))==1, is.na(x[which.max(FFP.R)])), yes = T, no=F)))

        # Index 3: FFP matrix yes, more than one isopleth no
        FFP.err.indx_3=which(lapply(FFP.R.list, function(x) ifelse(sum(is.na(x))>1, yes = T, no=F))==T)

        FFP.err.indx=sort(c(FFP.err.indx_1, FFP.err.indx_2, FFP.err.indx_3))

        # Re-create the full QC modifying 1, 2 and 3 
        QC=rep(0, nrow(FFP.input.df.cur))
        
        QC[FFP.err.indx_1]=1
        QC[FFP.err.indx_2]=2
        QC[FFP.err.indx_3]=3
        
        # Plot warning where the model was not computed 
        invisible(sapply(FFP.err.indx_1, 
                         function(x) cat(warn(paste0('\n [!] according to the model requirements, FFP was not computed for the timestamp ',
                                                     as.character(format(FFP.input.df.cur[x ,'TIMESTAMP'], '%Y-%m-%d %H:%M'))))))) 
        
        invisible(sapply(FFP.err.indx_2, 
                         function(x) cat(warn(paste0('\n [!] according to the model requirements, ', FFP.R[which.max(FFP.R)], 
                                                     '% isoline was not computed for the timestamp ',
                                                     as.character(format(FFP.input.df.cur[x ,'TIMESTAMP'], '%Y-%m-%d %H:%M'))))))) 
        
        invisible(sapply(FFP.err.indx_3, 
                         function(x) cat(warn(paste0('\n [!] according to the model requirements, ', 
                                                     paste0(FFP.R[which(is.na(FFP.R.list[[x]]))], '% '), 
                                                     'isoline was not computed for the timestamp ',
                                                     as.character(format(FFP.input.df.cur[x ,'TIMESTAMP'], '%Y-%m-%d %H:%M'))))))) 
        
        # Add messages to log
        Log_list=c(Log_list, sapply(FFP.err.indx_1, function(x) paste0('[!] according to the model requirements, FFP was not computed for the timestamp ',
                                                                       as.character(format(FFP.input.df.cur[x ,'TIMESTAMP'], '%Y-%m-%d %H:%M')))))
        
        Log_list=c(Log_list, sapply(FFP.err.indx_2, function(x) paste0('\n [!] according to the model requirements, ', FFP.R[which.max(FFP.R)], 
                                                                       '% isoline was not computed for the timestamp ',
                                                                       as.character(format(FFP.input.df.cur[x ,'TIMESTAMP'], '%Y-%m-%d %H:%M')))))
        
        Log_list=c(Log_list, sapply(FFP.err.indx_3, function(x) paste0('\n [!] according to the model requirements, ', 
                                                                       paste0(FFP.R[which(is.na(FFP.R.list[[x]]))], '% '), 
                                                                       'isoline was not computed for the timestamp ',
                                                                       as.character(format(FFP.input.df.cur[x ,'TIMESTAMP'], '%Y-%m-%d %H:%M')))))
        
        # Check if at least one FFP was computed #
        if(length(FFP.err.indx_1) == length(FFP.ls))
          
        {
          cat(warn(paste0('\n', bold(unique(format(FFP.input.df.cur$TIMESTAMP, '%Y-%m-%d'))), ': ', 
                          '[!] according to the model requirements, no FFP was computed. If specified, the nc file is created anyway')))
          
          # Add message to log
          Log_list=c(Log_list, sapply(FFP.err.indx_1, 
                                      function(x) paste0(unique(format(FFP.input.df.cur$TIMESTAMP, '%Y-%m-%d')), ': ', 
                                                         ' [!] according to the model requirements, no FFP was computed')), '') 
          
        } else ## 
          
        {
          
          ## -2- ## Plot and matrix saving looping over the list elements # Except for matrices not produced #
          
          for(i in which(QC %in% c(0, 2, 3)))
            
          {
            
            if (drop.trivial.srcs) {
              
              # fetch the footprint value at 90% (or the maximum desired) cumulative distribution
              FFP.90trshd <- NULL
              FFP.90trshd <- FFP.ls[[i]]$'fr'[which.max(FFP.R)]
              
              # and remove the FFP values lower that that
              FFP.ls[[i]]$'fclim_2d'[FFP.ls[[i]]$'fclim_2d' < FFP.90trshd] <- NA
              
            }
            
            
            # . Save single plots ---------------
            if (save.plot.FFP.mtx) {
              # # Two‐dimensional footprint with contour lines of R%
              
              # Do it once
              if(i == min(which(QC %in% c(0, 2, 3)))){
                cat(prog.mes(paste0('\n', bold(unique(format(FFP.input.df.cur$TIMESTAMP, '%Y-%m-%d'))), ': ', 
                                    'saving FFP matrix plots ...\n')))
                
                # Add message to log 
                Log_list=c(Log_list, paste0(unique(format(FFP.input.df.cur$TIMESTAMP, '%Y-%m-%d')), ': ', 
                                            'FFP matrix plots saved'), '')
              }
              
              # Create folder if not exists
              if(!dir.exists(paste0(FFP.output.dir, '/individual_footprints_plots/'))){
                FFP.mtx.img.dir <- paste0(FFP.output.dir, '/individual_footprints_plots/'); dir.create(FFP.mtx.img.dir)}
              
              # Create a dataframe for the plot 
              DF_grid=data.frame(X=as.vector(t(FFP.ls[[i]]$'x_2d')) + EC.tower.utm[1], 
                                 Y=as.vector(t(FFP.ls[[i]]$'y_2d')) + EC.tower.utm[2], 
                                 Z=as.vector(FFP.ls[[i]]$'fclim_2d'))
              
              # Create a dataframe for the isopleths 
              names(FFP.ls[[i]]$xr)=FFP.ls[[i]]$r
              names(FFP.ls[[i]]$yr)=FFP.ls[[i]]$r
              
              DF_iso=data.frame(stack(FFP.ls[[i]]$xr), stack(FFP.ls[[i]]$yr)['values'])
              DF_iso$values=DF_iso$values + EC.tower.utm[1]
              DF_iso$values.1=DF_iso$values.1 + EC.tower.utm[2]
              
              # Plot              
              ggplot(DF_grid)+
                geom_raster(aes(x = X, y = Y, fill = Z), alpha=0.85)+
                geom_path(data=DF_iso, aes(x = values, y = values.1, group = ind), color='grey90')+
                scale_fill_viridis_c(na.value = "transparent", option = 'D')+
                geom_vline(xintercept = EC.tower.utm[1], color='grey70')+
                geom_hline(yintercept = EC.tower.utm[2], color='grey70')+
                annotate(geom = 'point', x=EC.tower.utm[1], y=EC.tower.utm[2], 
                         shape=21, fill='red', size=2)+
                labs(title=paste0('Footprint estimation at ', site.ID), 
                     subtitle=paste0('Timestamp: ', FFP.input.df.cur[i, 'TIMESTAMP'], ' GMT'), 
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
              ggsave(paste0(FFP.mtx.img.dir, site.ID, '_', 'FFP2Dmtx_', 
                            format(FFP.input.df.cur$'TIMESTAMP'[i], '%Y%m%d%H%M'),'.png'), 
                     width = 23, height = 13, units='cm', dpi = 400)
              
            } # Plot condition ending...
            
          } # Loop over valid matrices ending...
          
          # Add and extra space to the messages
          cat('\n')
          
          # . Save CSV single matrix ---------------
          
          if(!is.null(save.FFP.mtx.as)){
            
            if (save.FFP.mtx.as == 'csv') {
              
              cat(prog.mes(paste0('\n', bold(unique(format(FFP.input.df.cur$TIMESTAMP, '%Y-%m-%d'))), ': ', 
                                  'saving FFP matrices as csv(s) ...', '\n')))
              
              # Add message to log 
              Log_list=c(Log_list, paste0('\n', unique(format(FFP.input.df.cur$TIMESTAMP, '%Y-%m-%d')), ': ', 
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
        
        
        ## -3- ## Save NC file ---------------
        
        if(is.null(save.FFP.mtx.as)) {
          
          cat(warn(('\n [!] No FFP output matrices are saved')))
          
          # Save output to log
          Log_list=c(Log_list, paste0('No FFP output matrices are saved'))
          
        } else if (save.FFP.mtx.as == 'nc') {
          
          # 3.1 # FFP matrix preparation ---------------
          
          # List to store dimention/attributes inputs
          FFP.input.list=list()
          
          # . Define nc single matrix dimensions ---------------
          
          if (any(QC==0)) { # Only the last valid time is considered 
            
            # store the X,Y dimensions of the FFP arrays (used for nc file dimension and variables)
            FFP.Xdim <- FFP.ls[which(QC == 0)][[1]]$'x_2d'[1, ]  # take the first row of x_2d, 
            FFP.Ydim <- FFP.ls[which(QC == 0)][[1]]$'y_2d'[, 1]  # take the first column of y_2d
            
          } else { # If no matrix exist use the extent and resolution parameters
            
            # store the X,Y dimensions of the FFP arrays (used for nc file dimension and variables)
            FFP.Xdim <- seq(min(FFP.domain.ext), max(FFP.domain.ext), dx)
            FFP.Ydim <- seq(min(FFP.domain.ext), max(FFP.domain.ext), dx) 
            
          } 
          
          # store the UTM coordinates of the FFP arrays (used for nc file dimension and variables)
          FFP.lon <- FFP.Xdim + EC.tower.utm[1]
          FFP.lat <- FFP.Ydim + EC.tower.utm[2]
          
          # Create an index of the correct isopleth calculation
          Idx_Iso=which(sapply(FFP.R.list, function(x) any(x %in% (which.FFP.R/100))))
          
          # . Define parameters linked to isopleths ---------------
          if(return.isopleth)
            
          {
            
            # Area calculation on ST_Polygons # 
            FFP.sf=data.frame(Index=NULL, Area=NULL)
            
            # Choose of the right isopleth 
            Pos=which(sort(FFP.R) %in% which.FFP.R)
            
            if(length(Pos)==0)
              
            {
              stop('FFP isopleth chosen to be returned is not in the isopleths list of the FFP function')
            }
            
            
            for(i in 1:length(FFP.ls))
              
            {
              
              if(i %in% Idx_Iso)
                
              {
                DF=data.frame(Lon=FFP.ls[[i]]$xr[[Pos]]+ EC.tower.utm[2], 
                              Lat=FFP.ls[[i]]$yr[[Pos]]+ EC.tower.utm[1])
                tmp=st_convex_hull(st_union(st_as_sf(na.omit(DF), coords=c('Lon', 'Lat'), crs=EPSG.code)))
                
                # Add a warning in case of NAs 
                if(nrow(DF) != nrow(na.omit(DF)))
                  
                {
                  cat(warn(paste0('\n [!] ', bold(unique(format(FFP.input.df.cur[i ,'TIMESTAMP'], '%Y-%m-%d %H:%M'))), ': ',  
                                  nrow(DF) - nrow(na.omit(DF)), ' NA(s) were found in the ', which.FFP.R, '% isopleth coordinates list. The geometry is still calculated', '\n')))
                  
                  Log_list=c(Log_list, '', 
                             paste0('[!] ', bold(unique(format(FFP.input.df.cur[i ,'TIMESTAMP'], '%Y-%m-%d %H:%M'))), ': ', 
                                    nrow(DF) - nrow(na.omit(DF)), ' NA(s) were found in the ', which.FFP.R, '% isopleth coordinates list. The geometry is still calculated'), '')
                }
                
                tmp=st_sf(geometry=tmp)
                tmp$Area=as.vector(st_area(tmp$geometry))
                tmp$Index=i
                tmp=st_drop_geometry(tmp)
                
                FFP.sf=rbind(FFP.sf, tmp)
                
              } else {
                
                tmp=data.frame(Index=i, Area=0)
                FFP.sf=rbind(FFP.sf, tmp)
                
              }
              
            }
            
            # Input list definition (for dimensions and attributes)
            # Adding variables to the list
            FFP.input.list[['polygon_index']]=FFP.sf$Index
            FFP.input.list[['polygon_area']]=FFP.sf$Area
            
            # Insert time variable 
            FFP.input.list[['time']]=as.character(format(FFP.input.df.cur$'TIMESTAMP', "%Y-%m-%d %H:%M"))
            
            # Define the maximum number of char for the date field in the polygon attribute table
            FFP.input.list[['char_length']]=nchar(FFP.input.list$time[1])
            
            ## Add 3 vertices where the chosen isopleth is NA (for tower center) directly in XR and YR 
            for(i in FFP.sf[FFP.sf$Area==0, 'Index']) #FFP.sf[FFP.sf$Area==0, 'Index']
            {
              # Three nodes with zero coords in the second sublist
              FFP.ls[[i]]$xr[[Pos]]=rep(0, 3)
              FFP.ls[[i]]$yr[[Pos]]=rep(0, 3)
              
            }
            
            # Selection of the chosen isoplete 
            FFP.iso=lapply(FFP.ls, FUN=function(x) na.omit(x$xr[[Pos]]))
            
            # Extraction of total number of nodes, number of nodes for each polygon
            FFP.input.list[['nodes_number_polygon']]=lengths(FFP.iso) 
            FFP.input.list[['nodes_total']]=sum(lengths(FFP.iso))  
            FFP.input.list[['polygon_number']]=length(FFP.iso) 
            
            # Point coordinates in a counterclock-wise order 
            FFP.input.list[['nodes_x']]=as.vector(na.omit(unlist(lapply(FFP.ls, 
                                                                        FUN=function(x) rev(x$xr[[Pos]]) + EC.tower.utm[1]))))
            
            FFP.input.list[['nodes_y']]=as.vector(na.omit(unlist(lapply(FFP.ls, 
                                                                        FUN=function(x) rev(x$yr[[Pos]]) + EC.tower.utm[2]))))
            
            # Remove the sf object 
            rm(FFP.sf)
            
          } # Isopleth input phase ....... ENDING
          
         
          # . Define input parameters (QC - EC Tower Coordinates) ---------------
          
          # Retrieving QC flag and adding it to the list
          FFP.input.list[['QC_Flag']]=QC
          
          # Add EC tower coords to the input list 
          FFP.input.list[['EC_Tower_x']]=EC.tower.utm[1]
          FFP.input.list[['EC_Tower_y']]=EC.tower.utm[2]
          
          
          # . FFP matrix handling ---------------
          
          # Adding a dummy NA grid for NULL elements #
          for(i in FFP.err.indx_1)
          {
            
            FFP.ls[[i]]$fclim_2d=array(data = NA, 
                                       dim = c(length(FFP.Xdim), length(FFP.Ydim)))
            
          }
          
          # . Extract FFP matrix and remove FFP overall list, from the matrix compute the array and normalize 
          #   each FFP matrix to sum 1 (now is around 0.9 as the threshold used)
          FFP.mtx=lapply(FFP.ls, function(x) x[['fclim_2d']]/sum(x[['fclim_2d']], na.rm = T))
          
          # Calculate scale and offset factors only if at least one not-NA matrix is returned
          if(any(FFP.input.list[['QC_Flag']]!=1))
            
          {
            # Add a scaling factor for the FFP values (to reduce object size)
            max.FFP <- max(unlist(FFP.mtx, use.names = F), na.rm=T)
            min.FFP <- min(unlist(FFP.mtx, use.names = F), na.rm=T)
            
            # Compute scale and offset factor 
            # [http://james.hiebert.name/blog/work/2015/04/18/NetCDF-Scale-Factors.html]
            nbits = 16
            
            # stretch/compress data to the available packed range
            scale.factor <- (max.FFP - min.FFP) / (2 ** nbits - 1)
            
            # translate the range to be symmetric about zero
            add.offset <- min.FFP + 2 ** (nbits - 1) * scale.factor
            
            scale.offset = c(sf = scale.factor, ofst = add.offset)
            # FFP.array <- FFP.array * 10^-7
            
            # Create an array of FFP matrices, adding scaling and offset factors
            FFP.array <- array()
            FFP.array <- array(round((unlist(FFP.mtx, use.names = F) - scale.offset['ofst'])/ scale.offset['sf']), 
                               dim = c(length(FFP.Xdim), length(FFP.Ydim), length(FFP.mtx)))
            
          } else {
            
            # Extract only the array without calculating scale and offset factors
            FFP.array <- array(unlist(FFP.mtx, use.names = F), 
                               dim = c(length(FFP.Xdim), length(FFP.Ydim), length(FFP.mtx)))
            
          }
                
          # Free unused memory
          rm(FFP.mtx)
          rm(FFP.ls)
          gc()
          
          # and convert to integers (to reduce object size)
          # format(object.size(FFP.array), 'Gb')
          mode(FFP.array) <- 'integer'
          # format(object.size(FFP.array), 'Gb')
          
          # windows()
          # image(FFP.array[,,1], col=RColorBrewer::brewer.pal(9,"Reds"))
          # sum((FFP.array[,,1]*scale.offset['sf']+scale.offset['ofst'])[!is.na(FFP.array[,,1])])
          
          # Add not scaled values as fill values # They are evaluated by scaling or offset #
          # Convert all NA's to the fillvalue used in the creation of nc file
          fillvalue <- -9999
          FFP.array[is.na(FFP.array)]=fillvalue
          
          # . * Create the netCDF filename and folder path 
          # path if not exist and file name, set dname
          ffp.ncfpath <- paste0(FFP.output.dir, "/nc_files/")
          
          if(!dir.exists(ffp.ncfpath)) {dir.create(ffp.ncfpath)}
          
          # File ID and name #
          ffp.ncfname <- paste0(ffp.ncfpath, site.ID, "_FFP2D_", levels(days.factor)[i.day], ".nc")
          ffp.dfname <- "FFP"
          
          # Create and write a projected netCDF file ++++++++++++++++++++++++++++++
          # Creating and writing (new) netCDF files involves first defining or “laying out” the dimensions and coordinate variables and 
          # the individual variables, and the attributes of each, and then creating the file and “putting” the data into the file, along 
          # with additional attributes or metadata.
          cat(prog.mes(paste0('\n', bold(unique(format(FFP.input.df.cur$TIMESTAMP, '%Y-%m-%d'))), ': ', 
                              'Creating and writing the netCDF file ...')))
          
          # 3.2 # Definition of netCDF dimensions and variables --------------- 
          
          # FFP grid dimension definition (x, y, time)
          FFP_dim_x <- ncdim_def(name = 'x', 
                                 units= 'm',
                                 longname = "x coordinate of projection",
                                 vals = FFP.lon)
          
          FFP_dim_y <- ncdim_def(name = 'y', 
                                 units= 'm',
                                 longname = "y coordinate of projection",
                                 vals = FFP.lat)
          
          FFP_dim_time <- ncdim_def(name = "time", 
                                    units="seconds since 1970-01-01 00:00:00 GMT",
                                    longname = "time in seconds since 1970-01-01 00:00:00 GMT",
                                    vals = as.integer(as.POSIXct(cur.ts, format='%Y%m%d%H%M', tz='GMT')), 
                                    unlim = T, # Time is unlimited (but not for humans!)
                                    calendar="gregorian")
          
          # Dimension (size=2) for the EC tower coordinates 
          FFP_dim_ECoords <- ncdim_def(name = 'EC_Coords', 
                                       units= '',
                                       longname = "projected (x and y, respectively) coordinates of EC tower",
                                       vals = 1:2, 
                                       create_dimvar = F)
          
          # Dimension (size=48) for the QC 
          FFP_dim_QC <- ncdim_def(name = 'QC_dim', 
                                  units= '',
                                  longname = "Quality flag of FFP matrix",
                                  vals = 1:length(FFP.input.list$QC), 
                                  create_dimvar = F)
          
          
          # GRID variables

          # FFP grid #
          FFP_mtx30_def <- ncvar_def(name = ffp.dfname,
                                     units = "m-2",
                                     dim = list(FFP_dim_x, FFP_dim_y, FFP_dim_time),
                                     missval = fillvalue,
                                     longname = "Matrix of normalised 2D footprint values",
                                     prec = "integer",
                                     compression = 9, 
                                     chunksizes = c(length(FFP.lat),length(FFP.lon), 1), 
                                     shuffle = TRUE)
          
          # CRS #
          FFP_var_proj <- ncvar_def(name = "UTM_Coordinate_System",
                                    units='m',
                                    dim = NULL,
                                    missval = NULL,
                                    longname = 'CRS reference system (Universal Transverse Mercator)',
                                    prec = "char")
          
          # Quality flag #
          FFP_var_QC <- ncvar_def(name = "quality_flag",
                                  units='1',
                                  dim = FFP_dim_QC, 
                                  longname = "Quality control on footprint calculation",
                                  missval = NULL,
                                  prec = "integer")
          
          # EC Tower coordinates
          FFP_var_ECoord <- ncvar_def(name = "EC_tower_coordinates",
                                      units = "m",
                                      dim = FFP_dim_ECoords,
                                      missval = NULL,
                                      longname = "Eddy covarinace tower coordinates in UTM",
                                      prec = "float")
          
          # If isopleth should be returned  
          
          if(return.isopleth) {
            
            # nodes # Nodes number
            FFP_dim_node <- ncdim_def(name = 'nodes', 
                                      units= '', ## If dimvar F, then units should be empty ## 
                                      longname = "number of nodes",
                                      vals=1:FFP.input.list$nodes_total, 
                                      create_dimvar = F)
            
            # instance # Number of polygons
            FFP_dim_inst <- ncdim_def(name = 'instance', 
                                      units = '',
                                      longname = "number of polygons",
                                      vals=1:FFP.input.list$polygon_number, 
                                      create_dimvar = F)
            
            # n_char # Number of character of each string #
            # Maximum dimension equal to the maximum length of a string
            FFP_dim_nchar <- ncdim_def(name = 'n_char', 
                                       units= '',
                                       longname = "length of each time string",
                                       vals=1:FFP.input.list$char_length, 
                                       create_dimvar = F)
            
            # Geometry container #
            FFP_var_geomc <- ncvar_def(name = "geometry_container",
                                       units='1',
                                       dim = NULL,
                                       missval = NULL,
                                       longname = 'container of the geometry',
                                       prec = "char")
            
            # Nodes count #
            FFP_var_nodes <- ncvar_def(name = "nodes_count",
                                       units='1',
                                       dim = FFP_dim_inst, 
                                       missval = NULL,
                                       longname = 'number of nodes for each feature',
                                       prec = "integer")
            
            # Node coordinates (X, Y) # 
            FFP_var_nodeX <- ncvar_def(name = "nodes_x",
                                       units='m', 
                                       dim = FFP_dim_node, 
                                       missval = NULL,
                                       longname = 'all nodes coordinates',
                                       prec = "float")
            
            FFP_var_nodeY <- ncvar_def(name = "nodes_y",
                                       units='m', 
                                       dim = FFP_dim_node, 
                                       missval = NULL,
                                       longname = 'all nodes coordinates',
                                       prec = "float")
            
            # Time (char) # 
            FFP_var_time <- ncvar_def(name = "observation_time",
                                      units='',
                                      dim = list(FFP_dim_nchar, FFP_dim_inst), 
                                      missval = NULL,
                                      longname = 'time of each observation',
                                      prec = "char")
            
            # Quality flag # 
            FFP_var_QC <- ncvar_def(name = "quality_flag",
                                    units='1',
                                    dim = FFP_dim_inst,
                                    longname = "Quality control on footprint calculation",
                                    missval = NULL,
                                    prec = "integer")
            
            # Polygon Area #
            FFP_var_area <- ncvar_def(name = "polygon_area",
                                      units='m',
                                      dim = FFP_dim_inst, 
                                      missval = NULL,
                                      longname = 'area of the isoplete polygon in squared meters',
                                      prec = "float")
            
            # polygon ID # 
            FFP_var_id <- ncvar_def(name = "polygon_id",
                                    units='1',
                                    dim = FFP_dim_inst, 
                                    missval = NULL,
                                    longname = paste0('id of each polygon referred to the ', which.FFP.R, '% isoplete'),
                                    prec = "integer")
            
            
            # NC creation if isopleth is desired
            ffp.ncout <- nc_create(filename = ffp.ncfname, 
                                   vars = list(FFP_mtx30_def, 
                                               FFP_var_id,
                                               FFP_var_time,
                                               FFP_var_QC, 
                                               FFP_var_area, 
                                               FFP_var_geomc, 
                                               FFP_var_proj, 
                                               FFP_var_nodes, 
                                               FFP_var_nodeX, 
                                               FFP_var_nodeY, 
                                               FFP_var_ECoord), 
                                   force_v4 = TRUE)
          } else {
            
            # NC creation if isopleth is not desired
            ffp.ncout <- nc_create(filename = ffp.ncfname,
                                   vars = list(FFP_mtx30_def,
                                               FFP_var_proj, 
                                               FFP_var_QC,
                                               FFP_var_ECoord), 
                                   force_v4 = TRUE)
            
          }
          
          #nc_close(ffp.ncout)
          
          # 3.3 # Put variables into the file --------------- 
          
          # FFP_mtx30_def
          ncvar_put(nc = ffp.ncout, varid = FFP_mtx30_def, vals = FFP.array) ## Array
          
          # FFP_var_ECoord
          ncvar_put(nc = ffp.ncout, varid = FFP_var_ECoord, 
                    vals = c(FFP.input.list$EC_Tower_x, FFP.input.list$EC_Tower_y)) ## Observation time
          
          # FFP_var_QC 
          ncvar_put(nc = ffp.ncout, varid = FFP_var_QC, 
                    vals = FFP.input.list$QC_Flag) ## Flag
          
          
          # Additional variables related to isopleth
          if(return.isopleth) {  
            
            # FFP_var_id 
            ncvar_put(nc = ffp.ncout, varid = FFP_var_id, vals = FFP.input.list$polygon_index) ## Polygon ID
            
            # FFP_var_nodes 
            ncvar_put(nc = ffp.ncout, varid = FFP_var_nodes, vals = FFP.input.list$nodes_number_polygon) ## Number of nodes
            
            # FFP_var_nodeX
            ncvar_put(nc = ffp.ncout, varid = FFP_var_nodeX, vals = FFP.input.list$nodes_x) ## Nodes coord x
            
            # FFP_var_nodeY 
            ncvar_put(nc = ffp.ncout, varid = FFP_var_nodeY, vals = FFP.input.list$nodes_y) ## Nodes coord y
            
            # FFP_var_area
            ncvar_put(nc = ffp.ncout, varid = FFP_var_area, vals = FFP.input.list$polygon_area) ## Pol area
            
            # FFP_var_time
            ncvar_put(nc = ffp.ncout, varid = FFP_var_time, vals = FFP.input.list$time) ## Observation time
            
          }
          
          # 3.4 # Put additional attributes into dimension and data variables ---------------
          
          
          # Grid mapping variable
          eastern_boundary <- UTM.zone * 6 - 180 # eastern boundary of the UTM zone, degrees
          western_boundary <- eastern_boundary - 6       # western boundary of the UTM zone, degrees
          central_meridian <- western_boundary + 3
          
          
          # Grid mapping name
          ncatt_put(nc = ffp.ncout, varid = "UTM_Coordinate_System", attname = "grid_mapping_name", 
                    attval = 'transverse_mercator')
          
          # WKT: SPHEROID[“", , , ...]
          ncatt_put(nc = ffp.ncout, varid = "UTM_Coordinate_System", attname = "semi_major_axis", 
                    attval = 6378137)
          
          # WKT: SPHEROID[“", , , ...]
          ncatt_put(nc = ffp.ncout, varid = "UTM_Coordinate_System", attname = "inverse_flattening", 
                    attval = 298.257223563)
          
          # WKT: PRIMEM[“", , ...]
          ncatt_put(nc = ffp.ncout, varid = "UTM_Coordinate_System", attname = "longitude_of_prime_meridian", 
                    attval = 0)
          
          # WKT: PARAMETER[“Latitude of natural origin”, ]
          ncatt_put(nc = ffp.ncout, varid = "UTM_Coordinate_System", attname = "latitude_of_projection_origin", 
                    attval = 0)
          
          # WKT: PARAMETER[“Longitude of natural origin”, ]
          ncatt_put(nc = ffp.ncout, varid = "UTM_Coordinate_System", attname = "longitude_of_central_meridian", 
                    attval = central_meridian)
          
          # WKT: PARAMETER[“Scale factor at natural origin”, ]
          ncatt_put(nc = ffp.ncout, varid = "UTM_Coordinate_System", attname = "scale_factor_at_central_meridian", 
                    attval = 0.9996)
          
          # WKT: PARAMETER[“False easting”, ]
          ncatt_put(nc = ffp.ncout, varid = "UTM_Coordinate_System", attname = "false_easting", 
                    attval = 500000)
          
          # WKT: PARAMETER[“False northing”, ]
          ncatt_put(nc = ffp.ncout, varid = "UTM_Coordinate_System", attname = "false_northing", 
                    attval = 0)
          
          # Add WKT (2) code using sf 
          ncatt_put(nc = ffp.ncout, varid = "UTM_Coordinate_System", attname = "crs_wkt", 
                    attval = sf::st_crs(paste0('epsg:', EPSG.code))[['wkt']])
          
          
          # FFP_mtx30_def
          ncatt_put(nc = ffp.ncout, varid = ffp.dfname, attname = "grid_mapping", 
                           attval = 'UTM_Coordinate_System')
          ncatt_put(nc = ffp.ncout, varid = ffp.dfname, attname = "coordinates", 
                           attval = 'x y')
          
          
          # Name of xy 
          ncatt_put(ffp.ncout, varid = "x", attname = "standard_name", 
                           attval = "projection_x_coordinate")
          ncatt_put(ffp.ncout, varid = "y", attname = "standard_name", 
                           attval = "projection_y_coordinate")

          
          if(any(FFP.input.list[['QC_Flag']]!=1))
            
          {
            ncatt_put(nc = ffp.ncout, varid = ffp.dfname, attname = "scale_factor", 
                             attval = scale.offset['sf'])
            ncatt_put(nc = ffp.ncout, varid = ffp.dfname, attname = "add_offset", 
                             attval = scale.offset['ofst'])
          }

                         
          # QC # 
          ncatt_put(ffp.ncout, varid = "quality_flag", attname = "standard_name",
                    attval='quality_flag')
          
          if(return.isopleth) { 
            
            ## Geometry container ## 
            ncatt_put(nc = ffp.ncout, varid = "geometry_container", attname = "geometry_type", 
                      attval = 'polygon')
            ncatt_put(nc = ffp.ncout, varid = "geometry_container", attname = "node_count", 
                      attval = 'nodes_count')
            ncatt_put(nc = ffp.ncout, varid = "geometry_container", attname = "node_coordinates", 
                      attval = 'nodes_x nodes_y')
            ncatt_put(nc = ffp.ncout, varid = "geometry_container", attname = "grid_mapping", 
                      attval = 'UTM_Coordinate_System')
            
            
            # QC # 
            ncatt_put(ffp.ncout, varid = "quality_flag", attname = "standard_name",
                      attval='quality_flag')
            ncatt_put(ffp.ncout, varid = "quality_flag", attname = "grid_mapping",
                      attval='UTM_Coordinate_System')
            ncatt_put(ffp.ncout, varid = "quality_flag", attname = "geometry",
                      attval='geometry_container')
            
            
            # polygon_id # Timeseries identifier
            ncatt_put(ffp.ncout, varid = "polygon_id", attname = "grid_mapping",
                      attval='UTM_Coordinate_System')
            ncatt_put(ffp.ncout, varid = "polygon_id", attname = "geometry",
                      attval='geometry_container')
            ncatt_put(ffp.ncout, varid = "polygon_id", attname = "cf_role",
                      attval='timeseries_id')
            
            
            # Name of xy 
            ncatt_put(ffp.ncout, varid = "x", attname = "standard_name", 
                      attval = "projection_x_coordinate")
            ncatt_put(ffp.ncout, varid = "y", attname = "standard_name", 
                      attval = "projection_y_coordinate")
            
            
            # FFP_var_time
            ncatt_put(ffp.ncout, varid = "observation_time", attname = "grid_mapping",
                      attval='UTM_Coordinate_System')
            ncatt_put(ffp.ncout, varid = "observation_time", attname = "geometry",
                      attval='geometry_container')
            
            
            # FFP_var_nodes # Additional specifications not required #
            
            
            # FFP_var_nodeX
            ncatt_put(ffp.ncout, varid = "nodes_x", attname = "axis",
                      attval='X')
            ncatt_put(ffp.ncout, varid = "nodes_x", attname = "units",
                      attval='m')
            
            
            # FFP_var_nodeY
            ncatt_put(ffp.ncout, varid = "nodes_y", attname = "axis",
                      attval='Y')
            ncatt_put(ffp.ncout, varid = "nodes_y", attname = "units",
                      attval='m')
            
            
            # FFP_var_area
            ncatt_put(ffp.ncout, varid = "polygon_area", attname = "grid_mapping",
                      attval='UTM_Coordinate_System')
            ncatt_put(ffp.ncout, varid = "polygon_area", attname = "geometry",
                      attval='geometry_container')
            
            
            # Feature type #
            ncatt_put(ffp.ncout, 0, "featureType", 
                      "timeSeries")
            
            
            # THREDDS type #
            ncatt_put(ffp.ncout, 0, "cdm_data_type", 
                      "Station")
            
          }
          
          # 3.5 # Global attributes ---------------
          
          # title
          ncatt_put(ffp.ncout, 0, "title", paste0("Footprint estimation at the ", site.ID, " station"))
          
          # Contact
          ncatt_put(ffp.ncout, 0, "Contact", 
                    paste0("Giacomo Nicolini and Luca Di Fiore, ICOS Ecosystem Thematic Center, ", 
                           "Euro-Mediterranean Center for Climate Change. ", 
                           "Email: giacomo.nicolini@cmcc.it; luca.difiore@cmcc.it"))
          
          # Conventions
          ncatt_put(ffp.ncout, 0, "Conventions", 
                    "CF-1.8")
          
          # creation_date
          ncatt_put(ffp.ncout, 0, "creation_date", 
                    date())
          
          # creator
          ncatt_put(ffp.ncout, 0, "creator", 
                    paste0("Giacomo Nicolini, ICOS Ecosystem Thematic Center, ", 
                           "Euro-Mediterranean Center for Climate Change, Viterbo, Italy; ",
                           "Luca Di Fiore, ICOS Ecosystem Thematic Center, ", 
                           "Euro-Mediterranean Center for Climate Change, Viterbo, Italy"))
          
          # institution
          ncatt_put(ffp.ncout, 0, "institution", 
                    "ICOS Ecosystem Thematic Center, Euro-Mediterranean Center for Climate Change, Viterbo, Italy")
          
          # keywords
          ncatt_put(ffp.ncout, 0, "keywords", 
                    "Flux footprints")
          
          # license
          ncatt_put(ffp.ncout, 0, "license", 
                    "CC-BY 4.0")
          
          # product_version
          ncatt_put(ffp.ncout, 0, "product_version", 
                    "1.0")
          
          # project
          ncatt_put(ffp.ncout, 0, "project", 
                    "Integrated Carbon Observation System")
          
          # references
          ncatt_put(ffp.ncout, 0, "references", 
                    "Kljun et al. (2015), doi:10.5194/gmd‐8‐3695‐2015")
          
          # source
          #ncatt_put(ffp.ncout, 0, "source", 
          #          "To define...")
          
          if(return.isopleth) { 
            # summary
            ncatt_put(ffp.ncout, 0, "summary", 
                      paste0("The file contains the footprint estimation and the quality layer of the station", 
                             site.ID, "at 30-min of temporal resolution. Additionally, a geometry containing the ", which.FFP.R, "% isopleth boundary is provided. ", 
                             "PID code of input data: ", PID))
          } else {
            
            # summary
            ncatt_put(ffp.ncout, 0, "summary", 
                      paste0("The file contains the footprint estimation and the quality layer of the station", 
                             site.ID, "at 30-min of temporal resolution. ", "PID code of input data: ", PID))
            
          }
          
          # frequency
          ncatt_put(ffp.ncout, 0, "frequency", 
                    "30min")
          
          # crs
          ncatt_put(ffp.ncout, 0, "crs", 
                    if(FFP.input.list$EC_Tower_y>0) 
                    {
                      paste0("WGS84/UTM", UTM.zone, 'N')
                    } else 
                    { paste0("WGS84/UTM", UTM.zone, 'S')})
          
          # geospatial_lat_resolution
          ncatt_put(ffp.ncout, 0, "geospatial_lat_resolution", 
                    paste0(dx, " m"))
          
          # geospatial_lon_resolution
          ncatt_put(ffp.ncout, 0, "geospatial_lon_resolution", 
                    paste0(dx, " m"))
          
          # geospatial_vertical_resolution
          #ncatt_put(ffp.ncout, 0, "geospatial_vertical_resolution", 
          #          "To define...")
          
          ## history
          ncatt_put(ffp.ncout, 0, "history",
                    paste("G. Nicolini and L. Di Fiore", date(), sep=", "))
          
          # Additional metadata ## 
          
          ncatt_put(ffp.ncout, 0, "Subjects", 
                    "Flux footprints")
          
          # Contributors
          ncatt_put(ffp.ncout, 0, "Contributors", 
                    paste0("Giacomo Nicolini, ICOS Ecosystem Thematic Center, ", 
                          "Euro-Mediterranean Center for Climate Change, Viterbo, Italy; ",
                          "Luca Di Fiore, ICOS Ecosystem Thematic Center, ", 
                          "Euro-Mediterranean Center for Climate Change, Viterbo, Italy"))
          
          # FundingReference
          #ncatt_put(ffp.ncout, 0, "FundingReference",
          #          "To define...")
          
          # Documentation of the model
          #ncatt_put(ffp.ncout, 0, "Documentation_of_the_model",
          #          "To define...")
          
          
          if(return.isopleth) { 
            
            # Variables
            ncatt_put(ffp.ncout, 0, "Variables",
                      paste0("FFP, QC, ", which.FFP.R, "% isopleth area"))
            
          } else {
            
            # Variables
            ncatt_put(ffp.ncout, 0, "Variables",
                      "FFP, QC")
          }
          
          # 3.6 # Close the file, writing data to disk ---------------
          nc_close(ffp.ncout)
          
          # Save jpeg file #
          
          if(save.ncplot)
            
          {
            
            ffp.ncplot.fpath <- paste0(FFP.output.dir, "/nc_files/", 'nc_plots/')
            if(!dir.exists(ffp.ncplot.fpath)) {dir.create(ffp.ncplot.fpath)}
            
            ffp.ncplot.fname <- paste0(ffp.ncplot.fpath, site.ID, "_FFP2D_", 
                                       levels(days.factor)[i.day], '.jpeg')
            
            jpeg(ffp.ncplot.fname, width = 20, height = 13, units = 'cm', res = 300)
            
            # Sampling of 16 layers (including missing ones) and plotting #
            terra::plot(terra::rast(ffp.ncfname)[[sort(sample(1:nlyr(terra::rast(ffp.ncfname)), 16))]])  
            
            dev.off()
            
          }
          
          cat(prog.mes(' done.\n'))
          
          cat('\n**************************************************************************************')
          
          # Save output to log
          Log_list=c(Log_list, '', paste0(unique(format(FFP.input.df.cur$TIMESTAMP, '%Y-%m-%d')), ': ', 
                                         'FFP matrices saved as netCDF'))
          
          # NETcdf file creation chunk ending
          
        } 
        
      } # Conditional DF validity check ending
      
    } # Days loop ending
    
    # Save messages to a unique log file
    if(save.log)
      
    {
      
      cat(prog.mes('\nSaving the log file.\n'))
      
      ffp.log.fpath <- paste0(FFP.output.dir, "/log/")
      if(!dir.exists(ffp.log.fpath)) {dir.create(ffp.log.fpath)}
      
      # Create log file
      log <- file.path(ffp.log.fpath, paste0(site.ID, '_log_', 
                                             as_date(min(FFP.input.df$'TIMESTAMP')), '_', 
                                             as_date(max(FFP.input.df$'TIMESTAMP')), ".log"))
      
      # Open log
      lf <- logr::log_open(log)
      
      # Write messages saved to log file 
      invisible(sapply(Log_list, function(x) logr::log_print(x, console = F, hide_notes = T, 
                                                             blank_after = F, msg = F)))
      
      # Close log
      logr::log_close()
      
    }
    
    } # DF, Coords and site.ID validity ending...

  } # FFP function ending
