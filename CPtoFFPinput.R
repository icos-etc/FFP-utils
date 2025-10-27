## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
## - MANIPULATE DATA FROM CARBON PORTAL AND BUILD AN OUTPUT READY FOR BEING PROCESSED BY THE FFP FUNCTION - ##
## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##

CPtoFFPinput = function(zip.file = NULL,             # ZIP file path
                        FFP.input.table = NULL,      # Existing input table csv 
                        start.date = NULL,           # Date start in the format YYYY-MM-DD
                        end.date = NULL,             # Date end in the format YYYY-MM-DD
                        save.input.table = TRUE,     # Save the output data as csv
                        return.input = TRUE,         # Return the output data as R list
                        MDS.clim = FALSE)            # Inlcude MDS-derived indices for gap-filled data # TRUE only in internal ETC pipeline

{

  # Stop the function if no input table is provided
  if (is.null(FFP.input.table) & is.null(zip.file))
  {
    stop("No ZIP file or input table were provided, please specify one of them!")
  }
  
  
  ## PREPARATORY PHASE ------
  
  # Install pacman and required packages if they are not already installed
  if (!require("pacman", quietly = T)){install.packages("pacman")}
  pacman::p_load(crayon, dplyr, tidyr, purrr, lubridate, stringr, sf, tibble, readr)
  
  # Error messages
  error <- crayon::red
  warn <- crayon::yellow
  note <- crayon::silver
  prog.mes <- crayon::cyan
  
  
  ## ZIP FILE PROCESSING ------
  
  # -- Build input DF from zip file -- #
  if (!is.null(zip.file))
    
  {
    # List file name within the zip folder #
    file_list <- utils::unzip(zip.file, list = T) %>% pull(Name)
    
    # -- SITE ID -- #
    site.ID <- read_csv(unz(zip.file, # Zip connection
                            file_list[str_detect(file_list, pattern = 'SITEINFO')]), col_types = cols()) %>%
      pull(SITE_ID) %>% unique()
    
    # -- STARTING MESSAGE -- #
    cat('\n********************************************************************\n')
    cat(prog.mes(paste0('Building the FFP input table for the ', site.ID, ' station.')))
    cat('\n********************************************************************\n')
    
    # -- COORDINATES -- #
    Coords <- read_csv(unz(zip.file, # Zip connection
                           file_list[grepl(file_list, pattern = 'SITEINFO')]), col_types = 'ccccc')
    
    Coords.df <- Coords %>% filter(VARIABLE %in% c('LOCATION_LONG', 'LOCATION_LAT')) %>% 
      select(VARIABLE, DATAVALUE) %>%
      mutate(VARIABLE = ifelse(str_detect(string = VARIABLE, pattern = 'LAT'), yes = 'lat', no = 'lon'), 
             DATAVALUE = as.numeric(DATAVALUE))
    
    # Store coordinates as a named vector
    Coords.vec <- Coords.df %>% pull(DATAVALUE, VARIABLE)
    
    # -- PID -- #
    PID <- NULL
    TOC <- list.files(path = dirname(zip.file), pattern = 'TOC', full.names = T)
    
    if (length(TOC) == 1)
      
    {
      # Read the TOC file and extract the PID for the specific site
      PID <- read_csv(TOC, name_repair = function(x) x %>% str_replace(pattern = ' ', replacement = '_'), col_types = cols()) %>%
        rowwise() %>% 
        mutate(SITE_CODE = (str_split(File_name, pattern = '_'))[[1]][2]) %>% 
        filter(SITE_CODE %in% site.ID) %>% 
        pull(PID)
      
      # If more than one TOC is available, PID is set to NA
    } else if (length(TOC) > 1)
    {
      cat(warn(paste0('\n [!] A unique PID cannot be retrieved for the ', site.ID, ' station.\n')))
      
    } else if (length(TOC) == 0)
    {
      cat(warn(paste0('\n [!] PID cannot be retrieved for the ', site.ID, ' station.\n')))
    }
    
    
    
    # -- DATA AVAILABILITY CHECK -- # -- FLUXES VARIABLES -- #
    
    # Create an empty dataframe (to prevent failures of the main dataframe processing part) 
    DF_Input=data.frame()
    
    Fluxes <- read_csv(unz(zip.file, # Zip connection
                           file_list[str_detect(file_list, pattern = 'FLUXES') & 
                                       str_detect(file_list, pattern = 'VARINFO', negate = T)]), col_types = 'ccccc')
    
    Fluxvar <- c('WS', 'MO_LENGTH', 'V_SIGMA', 'USTAR', 'WD', 'PBLH')
    
    # Wind variables check 
    if(!(any(colnames(Fluxes) %in% c('WD', 'WS')))) # If WS/WD are not in the fluxes, check in the meteo file # 
      
    {
      # Meteo # 
      Meteo <- read_csv(unz(zip.file, # Zip connection
                            file_list[str_detect(file_list, pattern = '_METEO_') & 
                                        str_detect(file_list, pattern = 'VARINFO', negate = T)]), col_types = 'ccccc')
    }  
    
    # Warning message #     
    if(exists('Meteo'))
    {
      if (any(!(Fluxvar %in% c(colnames(Fluxes), colnames(Meteo)))))
      {
        if (sum(!(Fluxvar %in% c(colnames(Fluxes), colnames(Meteo)))) > 1)
        {
          cat(warn(paste0('\n [!] ', 
                          bold(str_flatten_comma(Fluxvar[!(Fluxvar %in% c(colnames(Fluxes), colnames(Meteo)))])),
                          ' are not available for the ', site.ID, ' station. FFP input table is not created.\n')))
          
        } else if (sum(!(Fluxvar %in% c(colnames(Fluxes), colnames(Meteo)))) == 1)
        {
          cat(warn(paste0('\n [!] ', 
                          bold(str_flatten_comma(Fluxvar[!(Fluxvar %in% c(colnames(Fluxes), colnames(Meteo)))])),
                          ' is not available for the ', site.ID, ' station. FFP input table is not created.\n')))
        }
        
        DF_Input <- NULL
      }
    } else if (any(!(Fluxvar %in% colnames(Fluxes))))
      {
      
      if (sum(!(Fluxvar %in% colnames(Fluxes))) > 1)
      {
        cat(warn(paste0('\n [!] ', 
                        bold(str_flatten_comma(Fluxvar[!(Fluxvar %in% colnames(Fluxes))])),
                        ' are not available for the ', site.ID, ' station. FFP input table is not created.\n')))
        
      } else if (sum(!(Fluxvar %in% colnames(Fluxes))) == 1)
      {
        cat(warn(paste0('\n [!] ', 
                        bold(str_flatten_comma(Fluxvar[!(Fluxvar %in% colnames(Fluxes))])),
                        ' is not available for the ', site.ID, ' station. FFP input table is not created.\n')))
      }
        DF_Input <- NULL
      }
    
    
    
    # -- DATA AVAILABILITY CHECK -- # -- HC -- #
    
    # First check (and warn) if the ancillary file exists #
    if (is_empty(file_list[str_detect(file_list, pattern = 'ANCILLARY')]))
      
    {
      cat(warn(paste0("\n [!] Canopy height can't be retrieved (ancillary file doesn't exist) for the ",
                      site.ID, " station. FFP input table is not created.\n")))
      
      DF_Input <- NULL
      
    } else {
      
      # Canopy height # Spread ancillary data by group
      Ancillary <- read_csv(unz(zip.file, # Zip connection
                                file_list[str_detect(file_list, pattern = 'ANCILLARY')]), col_types = 'ccccc')
      
      # Check for the height group
      if (all(!(Ancillary$VARIABLE_GROUP %in% "GRP_HEIGHTC"))) # If group doesn't exist give a warning #
        
      {
        cat(warn(paste0("\n [!] Canopy height can't be retrieved from the ancillary file for the ",
                        site.ID, " station. FFP input table is not created.\n")))
        
        DF_Input <- NULL
        
      } # Height availability within ancillary ending...
      
    } # Ancillary availability ending....
    
    
    # -- DATA PROCESSING -- # -- FLUXES VARIABLES -- #
    
    if(!is.null(DF_Input))
      
    {
      # Format the time span to keep only the start
      Fluxes$TIMESTAMP <- strptime(Fluxes$TIMESTAMP_START, format = '%Y%m%d%H%M', tz = 'GMT')
      
      # Select the required columns # c('umean', 'ol', 'sigmav', 'ustar', 'wind_dir', 'h')
      Fluxes.df <- Fluxes %>% select(any_of(c('TIMESTAMP', Fluxvar)))
      
      # Add meteo data if exist # 
      if(exists('Meteo'))
      {
        Meteo$TIMESTAMP <- strptime(Meteo$TIMESTAMP_START, format = '%Y%m%d%H%M', tz = 'GMT')
        
        # Select the required columns # c('umean', 'ol', 'sigmav', 'ustar', 'wind_dir', 'h')
        Meteo.df <- Meteo %>% select(any_of(c('TIMESTAMP', Fluxvar)))
        
        # Join the variables # 
        Fluxes.df <- Fluxes.df %>% left_join(Meteo.df, by='TIMESTAMP')
      }
      
      
      
      # -- DATA PROCESSING -- # -- HC -- #
      
      Ancillary.df <- Ancillary %>%
        split(.$VARIABLE_GROUP) %>%
        map( ~ .x %>% spread(key = VARIABLE, value = DATAVALUE)) #%>%
        #map(~ .x %>% write_csv(paste0(.x$VARIABLE_GROUP %>% unique(), '.csv'))) # CSV write
      
      # Force DF to have date and spp columns, filled with NAs
      cols <- c(HEIGHTC_DATE = NA, HEIGHTC_SPP = NA, HEIGHTC_DATE_START = NA, HEIGHTC_VEG_STATUS = "Alive")
      
      # DF filter to keep only HC
      GRP_HEIGHTC_DF <- Ancillary.df[['GRP_HEIGHTC']] %>% 
        tibble::add_column(., !!!cols[!names(cols) %in% names(.)]) # Add columns if not available # 
      
      # If HC 90th percentile is available
      if(any(GRP_HEIGHTC_DF$HEIGHTC_STATISTIC %in% "90th Percentile")) 
        
      {
        hc.df.tmp <- GRP_HEIGHTC_DF %>%
          filter(HEIGHTC_STATISTIC == '90th Percentile' & is.na(HEIGHTC_SPP) & HEIGHTC_VEG_STATUS == 'Alive') %>% 
          mutate(HEIGHTC_DATE = case_when(is.na(HEIGHTC_DATE) ~ HEIGHTC_DATE_START, # sometimes HEIGHTC_DATE is NA
                                          TRUE ~ HEIGHTC_DATE)) %>%
          select(HEIGHTC_DATE, HEIGHTC) %>%
          rename('TIMESTAMP' = 'HEIGHTC_DATE', 'hc' = 'HEIGHTC') %>%
          mutate(TIMESTAMP = strptime(TIMESTAMP, format = '%Y%m%d', tz = 'GMT'), hc = as.numeric(hc))
        
      } else # HC 90th percentile unavailable (mean is always available) # 
        
      {
        # Build the DF
        hc.df.tmp <- GRP_HEIGHTC_DF %>%
          filter(HEIGHTC_STATISTIC == 'Mean' & is.na(HEIGHTC_SPP) & HEIGHTC_VEG_STATUS == 'Alive') %>% 
          mutate(HEIGHTC_DATE = case_when(is.na(HEIGHTC_DATE) ~ HEIGHTC_DATE_START, # sometimes HEIGHTC_DATE is NA
                                          TRUE ~ HEIGHTC_DATE)) %>%
          select(HEIGHTC_DATE, HEIGHTC) %>%
          rename('TIMESTAMP' = 'HEIGHTC_DATE', 'hc' = 'HEIGHTC') %>% 
          mutate(TIMESTAMP=case_when(nchar(TIMESTAMP)==6 ~ paste0(TIMESTAMP, '15'), # Missing day 
                                     nchar(TIMESTAMP)==4 ~ paste0(TIMESTAMP, '0615'), # Missing month 
                                     T ~ TIMESTAMP))
        
        # Format timestamp and height 
        hc.df.tmp <- hc.df.tmp %>% 
          mutate(TIMESTAMP = strptime(TIMESTAMP, format = '%Y%m%d', tz = 'GMT'), hc = as.numeric(hc))
        
        # HC warning (average is the only measurement available)
        cat(warn(paste0("\n [!] Canopy height value extracted from the site average for the ",
                        site.ID, " station. FFP extent is underestimated.\n")))
        
      }
      
      # If multiple HC obs. exists, take the average # 
      if(any(duplicated(hc.df.tmp$TIMESTAMP))) 
        
      {
        hc.df.tmp <- hc.df.tmp %>% 
          group_by(TIMESTAMP) %>% 
          summarise(hc=mean(hc, na.rm=T))
      }
      
      # Create a dummy DF with the full dates and join it with the first date of measurement
      hc.df <- data.frame(TIMESTAMP = seq.POSIXt(min(hc.df.tmp$TIMESTAMP), max(Fluxes.df$TIMESTAMP), '30 min')) %>%
        left_join(hc.df.tmp, by = 'TIMESTAMP')
      
      # Fill the dataframe
      hc.df <- fill(hc.df, hc, .direction = 'down')
      
      # Correction: if canopy height is zero, correct to 0.2 # Otherwise the FFP fails # 
      hc.df <- hc.df %>% mutate(hc=ifelse(hc == 0, yes=0.2, no=hc))
    
      
      
      # -- DATA PROCESSING -- # -- HM -- #
      
      # Using USTAR as reference variable for HM # 
      VarInfo <- read_csv(unz(zip.file, # Zip connection
                              file_list[str_detect(file_list, pattern = 'FLUXES') &
                                         str_detect(file_list, pattern = 'VARINFO')]), col_types = 'ccccc') %>%
        spread(key = VARIABLE, value = DATAVALUE)
      
      # Extract unique pair(s) of HM and corresponding date (excluding NAs)
      hm.df.tmp <- VarInfo %>%
        filter(VAR_INFO_VARNAME == 'USTAR') %>% 
        select(VAR_INFO_DATE, VAR_INFO_HEIGHT) %>% unique() %>% na.omit() %>%
        rename('TIMESTAMP' = 'VAR_INFO_DATE', 'hm' = 'VAR_INFO_HEIGHT') %>%
        mutate(TIMESTAMP = strptime(TIMESTAMP, format = '%Y%m%d', tz = 'GMT'), hm = as.numeric(hm))
      
      # Create a dummy DF with the full dates and join it with the first date of
      hm.df <- data.frame(TIMESTAMP = seq.POSIXt(min(Fluxes.df$TIMESTAMP), max(Fluxes.df$TIMESTAMP), '30 min')) %>%
        left_join(hm.df.tmp, by = 'TIMESTAMP')
      
      # Fill the dataframe
      hm.df <- fill(hm.df, hm, .direction = 'down')
      
      
      # --- IDX BUILD --- #
      # Get dataset years
      Years=Fluxes.df$TIMESTAMP %>% year() %>% unique()
      DF_Idx=Years %>% map(function(x) 
                      seq(as.POSIXct(paste0(x, "-01-01 00:00:00"), format = "%Y-%m-%d %H:%M", tz="GMT"), 
                          as.POSIXct(paste0(x, "-12-31 23:30:00"), format = "%Y-%m-%d %H:%M", tz="GMT"), 
                          by = "30 mins") %>% tibble(TIMESTAMP=.)) %>% 
        map(function(x) x %>% mutate(Idx=row_number(TIMESTAMP)-1)) %>% 
        bind_rows()
      
      
      # --- FULL JOIN --- #
      DF_Input <- Fluxes.df %>% left_join(hm.df, by = 'TIMESTAMP') %>% 
        left_join(hc.df, by = 'TIMESTAMP') %>% 
        na.omit()
      
      
      # -- D, Z0, Zm -- #
      DF_Input <- DF_Input %>% mutate(d = 2 / 3 * hc, z0 = 0.15 * hc, zm = hm - d)
      
      # Order and rename variables #
      DF_Input <- DF_Input %>% select(TIMESTAMP, hm, hc, d, z0, zm, WS, MO_LENGTH, V_SIGMA, USTAR, WD, PBLH) %>%
        rename('umean' = 'WS', 'ol' = 'MO_LENGTH', 'sigmav' = 'V_SIGMA', 'ustar' = 'USTAR', 'wind_dir' = 'WD',
               'PBL' = 'PBLH')
      
      
      # -- DF REFINING -- #
      
      # Change nodata values (from -9999 to NA)
      DF_Input[DF_Input == -9999] <- NA
      
      
      # -- MDS CHUNK -- #
      
      if(MDS.clim)
        
      {
       
        # Add rownumber (Index) #
        DF_Input=DF_Input %>% left_join(DF_Idx, by = "TIMESTAMP") 
         
        # MDS indice path #
        Idx_dir <- paste0(getwd(), '/', site.ID, '/', 'Input', "\\", "MDS")
        
        # MDS output file reading #
        Idx=(list.files(path = Idx_dir, full.names = T, pattern = "indices"))
        
        # MDS gap-filling indices extraction # 
        Idx_full=Idx %>% map(function(x) x %>% read_csv(col_types = "iic") %>% 
                               rowwise() %>% mutate(MDS_Idx=str_split(ZERO_BASED_INDICES, pattern=",")) %>% 
                               mutate(MDS_Idx=list(as.numeric(MDS_Idx))) %>% 
                               select(-ZERO_BASED_INDICES, -SAMPLES_COUNT) %>% 
                               rename("Idx"="ZERO_BASED_INDEX")) 
        
        # Incorporate MDS gap-filling indices in the dataframe for each half-hour #
        DF_Input=DF_Input %>% mutate(Anno=year(TIMESTAMP)) %>% split(.$Anno) %>% 
          map(function(x) x %>% select(-Anno) %>% 
                mutate(ZERO_BASED_INDEX=row_number(TIMESTAMP)-1)) %>% 
          map2(.x=., .y=Idx_full, ~ .x %>% left_join(.y, by="Idx")) %>% 
          map(function(x) x %>% select(-ZERO_BASED_INDEX)) %>% 
          bind_rows()
        
      }
      
      # -- DF FILTERING -- #
      
      # Filter between start and end dates
      if (!is.null(start.date) & is.null(end.date))
      {
        DF_Input <- DF_Input %>% filter(as_date(TIMESTAMP) >= start.date)
        
      } else if (is.null(start.date) & !is.null(end.date))
      {
        DF_Input <- DF_Input %>% filter(as_date(TIMESTAMP) <= end.date)
        
      } else if (!is.null(start.date) & !is.null(end.date))
      {
        DF_Input <- DF_Input %>% filter(as_date(TIMESTAMP) >= start.date & as_date(TIMESTAMP) <= end.date)
      }
      
      # Check data overlap
      if (nrow(DF_Input) == 0)
        
      {
        cat(warn(paste0("\n [!] No overlap between date window and dataframe timestamp for the ",
                        site.ID, " station. FFP input table is not created.\n")))
        
        DF_Input <- NULL
      }
      
    } # data processing ending....
    
    
    # -- OUTPUT -- #
    
    Model_Input <- list()
    
    Model_Input[['site_id']] <- site.ID
    Model_Input[['tower_coordinates']] <- Coords.vec
    Model_Input[['PID']] <- PID
    Model_Input[['FFP_input_parameters']] <- DF_Input
    
    # Add the data as data frame
    Model_DF <- data.frame(
      site_id = site.ID,
      lat = Coords.vec['lat'] %>% as.numeric(),
      lon = Coords.vec['lon'] %>% as.numeric(),
      PID = PID) %>%
      bind_cols(DF_Input)
    
    # ZIP file chunk ending...
    
  } else if (!is.null(FFP.input.table))
    
    ## EXISTING INPUT TABLE PROCESSING ------
  {
    # Get site code
    site.ID <- str_split(basename(FFP.input.table), pattern = '_')[[1]][1]
    
    cat(warn(paste0('\nNo zip file provided for the ', site.ID, ' station.')))
    cat(prog.mes(paste0('\nUsing FFP input table as input and built-in table for coordinates....')))
    cat(prog.mes(paste0('\nCreating internal list with filtered dataframe, site code, coordinates and, if available, PID....\n')))
    
    # Get built-in coordinates
    Coords <- data.frame(
      site_id = c('DE-SfN', 'FI-Hyy', 'FR-Bil', 'SE-Svb', 'LMP', 'DE-HGT', 'IE-GtF', 'IT-Lsn', 'BE-Dor', 'IT-OXm', 
                  'SE-Nor', 'FI-Sii', 'IT-TrF', 'FI-Ouk', '58US', 'IZO', 'CZ-BK1', 'FKL', 'UK-AMo', 'NO-Hur', 
                  'CUXHAVEN', 'NO-SOOP-Bergen Kirkenes', 'FI-Tvm', 'MLH', 'FI-Sod', 'SSL', 'GF-Guy', 'IT-Ren', 'UTO', 
                  'BE-Wm2', 'BE-Wm1', 'FR-Fon', 'BE-Maa', 'DE-Brs', 'IPR', 'IT-FOS-W1M3A', 'OXK', 'WES', 'IE-JtC', 
                  'ES-LMa', 'GOSars', 'SE-OES', 'FR-MsS', 'SE-TAV', 'FR-Lam', 'SE-Myc', 'IT-Tor', 'BE-Bra', 'SNO', 
                  'ARN', 'FI-Ant', 'FR-Lus', 'FI-TVA', 'FI-Kvr', '1199', 'HEL', 'IT-SR2', 'FI-SER', 'DE-Gri', 'GR-Prt', 
                  '74FS', 'DE-RuW', 'FR-Tou', 'UK-WCO', 'DE-Msr', 'DE-RuS', 'CBW', 'CMN', 'IT-Noe', 'GL-NuF', 'PS-VOS', 
                  'PS-VOS', 'IE-DyL', 'SMR', 'BALTIC-VOS', 'GAT', 'FR-Mej', 'POT', 'TCarrier', 'SVB', 'IT-Cp2', 'BIR', 
                  'DK-RCW', 'DE-BeR', '11SS', 'STE', 'CRP', 'TOH', 'JUE', 'DE-HoH', 'KIT', 'IT-PCm', 'IT-Sas', 'LUT', 
                  'DE-Kli', 'IT-FOS-Lampedusa', 'ZEP', 'ES-SOOP-ESTOC', 'GR-HeK', '11BE', 'GR-HeM', 'RUN', 'IE-Doa', 
                  'DE-FOS-CVOO', 'PUY', 'GL-ZaF', 'CZ-wet', 'WAO', 'KRE', 'CIBA', 'IT-Lpd', 'DE-Geb', 'SE-Sto', 'DE-Har', 
                  'FR-Lqu', 'FR-Aur', 'MHD', 'DE-RuR', 'FR-EM2', 'FR-CLt', 'DE-Amv', 'IE-Cra', 'IT-BCi', 'FI-Var', 
                  'CH-BaK', 'DE-GsB', 'PRS', 'IT-BFt', 'JFJ', 'CH-Dav', 'OPE', 'HUN', 'NL-Loo', 'FI-Kmp', 'BE-Vie', 
                  '687B', 'FR-Gri', 'ZSF', 'CD-Ygb', 'DE-Tha', 'NA-VOS', 'GL-ZaH', 'SE-Htm', 'IT-E2M', 'DE-Okd', 'SAC', 
                  'PUI', 'SE-Deg', 'GL-Dsk', 'CZ-Lnz', 'FR-Pue', 'DE-Hai', 'LIN', 'FR-LGt', 'NOR', 'IT-MBo', 'HPB', 
                  'DE-Gwg', 'FI-Lom', 'PAL', 'VTO', 'FR-Hes', 'FI-Ken', 'DK-Sor', 'DE-Hzd', 'DE-Fen', 'IT-Col', 'PALOMA', 
                  'IT-Niv', 'DE-Kie', 'FR-FBn', 'TRN', 'VOS%20France-Brazil', 'BE-Lon', '26RA', 'FI-Let', 'HTM', '48MB', 
                  'RGL'),
      lat = c(47.80639, 61.84741, 44.493652, 64.25611, 35.5181, 79.006833, 53.040253, 45.740482, 50.311874, 43.774464, 
              60.0865, 61.83265, 45.82376, 66.37752, NA, 28.309, 49.502075, 35.3376, 55.792545, 60.37163, 53.87706, NA, 
              59.8418, 55.3718, 67.36239, 47.9167, 5.2787, 46.58686, 59.7839, 51.29747, 51.296432, 48.476357, 50.97987, 
              52.2966, 45.8147, 43.834516, 50.03, 54.9231, 52.298187, 39.94033, NA, 57.430061, 48.53816, NA, 43.496437, 
              58.36503, 45.84444, 51.30761, 81.3605, 37.104, 63.16358, 46.414246, 59.8418, 61.84662, 51.57989, 54.1804, 
              43.73202, NA, 50.95004, 39.5434, 49, 50.50493, 43.572857, 50.25, 47.80918, 50.865906, 51.9703, 44.1936, 
              40.606174, 64.130936, NA, NA, 53.45647, 61.8474, NA, 53.0657, 48.1184, 40.601, NA, 64.256, 41.704266, 
              58.3886, 55.680683, 52.457233, NA, 53.0431, 52.1856, 51.8088, 50.9102, 52.08656, 49.0915, 40.8741, 
              40.716957, 53.4036, 50.89306, 35.817, 78.9072, 29.166667, 35.33613, NA, 35.32343, -21.0796, 52.94873, 
              17.6, 45.7719, 74.48152, 49.02465, 52.95, 49.572, 41.8154861, 35.5252, 51.09973, 68.35594, 47.933, 
              45.6444, 43.54965, 53.326111, 50.621914, 49.87211, 45.0414, 52.1757, 53.32309, 40.52375, 67.7549, 
              47.56173, 52.029648, 45.93, 45.197754, 46.5475, 46.81533, 48.5619, 46.9558, 52.166447, 60.20289, 
              50.304962, NA, 48.84422, 47.4165, 0.814444, 50.96256, NA, 74.4733, 56.09763, 41.5705, 53.4126, 48.7227, 
              62.9096, 64.18203, 69.25349, 48.68155, 43.7413, 51.079407, 52.1663, 47.322918, 60.0864, 46.01468, 
              47.8011, 47.57083, 67.99724, 67.9733, 51.9397, 48.6741, 67.98721, 55.48587, 50.96381, 47.83292, 41.8494, 
              45.6204, 45.49091, 52.97343, 43.24079, 47.9647, NA, 50.55162, NA, 60.64183, 56.0976, 45.698, 51.9975),
      lon = c(11.3275, 24.29477, -0.956092, 19.7745, 12.6322, 4.334167, -8.001742, 12.750297, 4.968113, 11.255111, 
              17.479504, 24.19285, 7.56089, 29.30768, NA, -16.499, 18.536882, 25.6696, -3.2436917, 11.07949, 8.704878, 
              NA, 23.2503, -7.3395, 26.63859, 7.9166, -52.9248, 11.43369, 21.3672, 4.6624384, 4.6572485, 2.780096, 
              5.631851, 10.4487, 8.636, 9.118163, 11.8083, 8.308, -6.499806, -5.77465, NA, 18.984339, 5.31187, NA, 
              1.237878, 12.1694, 7.578055, 4.51984, -16.3943, -6.734, 27.2347, 0.1206518, 23.249, 24.2804, 2.993217, 
              7.8833, 10.29091, NA, 13.51259, 21.4989, -16.5, 6.3309627, 1.37474, -4.217, 11.456168, 6.4471445, 4.9264,
              10.6999, 8.151694, -51.386066, NA, NA, -9.89284, 24.2947, NA, 11.4429, -1.79635, 15.7237, NA, 19.775, 
              12.357293, 8.2519, 12.101398, 13.315827, NA, 8.4588, -6.33686, 10.535, 6.4096, 11.22235, 8.4249, 14.2504, 
              8.575956, 6.3528, 13.52238, 12.783, 11.8867, -15.5, 25.1328, NA, 25.13017, 55.3841, -7.253945, -24.5, 
              2.9658, -20.555773, 14.77035, 1.121, 15.08, -4.9321167, 12.5365, 10.91463, 19.045208, 7.5981, 2.7349, 
              1.106103, -9.903889, 6.304126, 3.02065, 6.41053, 6.95537, -7.641774, 14.957444, 29.61, 7.58049, 11.10478, 
              7.7, 10.741966, 7.9851, 9.85591, 5.5036, 16.6522, 5.74355, 24.9611, 5.998099, NA, 1.95191, 10.9796, 
              24.502472, 13.56515, NA, -20.550869, 13.41897, 18.7746, 9.09247, 2.142, 27.6549, 19.55654, -53.51413, 
              16.946331, 3.5957, 10.452089, 14.1226, 2.284102, 17.4794, 11.045831, 11.0246, 11.03261, 24.209179, 24.1157,
              -10.2445, 7.06465, 24.24301, 11.644645, 13.48978, 11.06056, 13.5881, 13.5658, 7.13943, 13.64393,
              5.67865, 2.1125, NA, 4.746234, NA, 23.95952, 13.4189, 13.708, -2.5399))
    
    # Extract coordinates of the site code
    Coords <- Coords %>% filter(site_id %in% site.ID)
    Coords.vec <- c(Coords %>% pull(lat), Coords %>% pull(lon))
    
    
    # -- DATE FILTERING -- #
    DF_Input <- read_csv(FFP.input.table, col_types = cols())
    
    if (!is.null(start.date) & is.null(end.date))
    {
      DF_Input <- DF_Input %>% filter(as_date(TIMESTAMP) >= start.date)
      
    } else if (is.null(start.date) & !is.null(end.date))
    {
      DF_Input <- DF_Input %>% filter(as_date(TIMESTAMP) <= end.date)
      
    } else if (!is.null(start.date) & !is.null(end.date))
    {
      DF_Input <- DF_Input %>% filter(as_date(TIMESTAMP) >= start.date & as_date(TIMESTAMP) <= end.date)
    }
    
    # If no overlap is returned after the filtering, set DF equal to NA, otherwise an empty DF is saved
    if (nrow(DF_Input) == 0)
      
    {
      cat(warn(paste0("\n [!] No overlap between date window and dataframe timestamp for the ",
                      site.ID, " station. FFP input table is not created.\n")))
      
      DF_Input <- NULL
    }
    
    
    
    # -- PID -- #
    
    # Extract PID from TOC file, located everywhere inside the project folder
    PID <- NULL
    TOC <- paste0(dirname(list.files(path = getwd(), pattern = paste0('ICOSETC_', site.ID, '_ARCHIVE'), 
                                     full.names = T, recursive = T)), '/!TOC.csv')
    
    if (length(TOC) == 1)
      
    {
      # Read the TOC file and extract the PID for the specific site
      PID <- read_csv(TOC, name_repair = function(x) x %>% str_replace(pattern = ' ', replacement = '_'), 
                      col_types = cols()) %>%
        rowwise() %>%
        mutate(SITE_CODE = (str_split(File_name, pattern = '_'))[[1]][2]) %>%
        filter(SITE_CODE %in% site.ID) %>%
        pull(PID)
      
      # If more than one TOC is available, PID is set to NA
    } else if (length(TOC) > 1)
    {
      cat(warn(paste0('\n [!] A unique PID cannot be retrieved for the ', site.ID,' station.\n')))
      
    } else if (length(TOC) == 0)
    {
      cat(warn(paste0('\n [!] PID cannot be retrieved for the ', site.ID, ' station.\n')))
    }
    
    # Return the output list
    Model_Input <- list()
    
    Model_Input[['site_id']] <- site.ID
    Model_Input[['tower_coordinates']] <- Coords.vec
    Model_Input[['PID']] <- PID
    Model_Input[['FFP_input_parameters']] <- DF_Input %>% select(- any_of(c('site_id', 'lat', 'lon', 'PID')))
    
    # Add the data as data frame
    Model_DF <- data.frame(
      site_id = site.ID,
      lat = Coords.vec[1] %>% as.numeric(),
      lon = Coords.vec[2] %>% as.numeric(),
      PID = PID) %>%
      bind_cols(DF_Input %>% select(- any_of(c('site_id', 'lat', 'lon', 'PID'))))
    
  } else if (is.null(FFP.input.table) & is.null(zip.file))
    
  {
    cat(warn(paste0('\nNo zip file or dataframe provided for the station. Please specify it.')))
  }
  
  
  ## SAVE AND RETURN OPTIONS ------
  
  # Save input tables in the FFP input folder
  if (all(save.input.table & is.data.frame(Model_Input[['FFP_input_parameters']])))
    
  {
    
    # In the project dir, create a folder with site code and input dir inside #
    Site_dir <- paste0(getwd(), '/', site.ID)
    Input_dir <- paste0(getwd(), '/', site.ID, '/', 'Input')
    
    if (!dir.exists(Site_dir)) {dir.create(Site_dir)}
    if (!dir.exists(Input_dir)) {dir.create(Input_dir)}
    
    # Write the CSV file
    write_csv(Model_DF, file = paste0(Input_dir, '/', 
                                      site.ID, '_', min(as_date(Model_DF$TIMESTAMP)) %>% str_replace_all('-', ''), '_',
                                      max(as_date(Model_DF$TIMESTAMP)) %>% str_replace_all('-', ''), '_FFPinput.csv'))
    
    message(prog.mes('\nFFP input table saved'))
    
  }
  
  
  # -- RETURN INPUT TABLE -- #
  if (return.input == TRUE) {return(Model_Input)}
  
} # Function ending