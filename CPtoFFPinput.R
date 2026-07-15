## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
## - MANIPULATE DATA FROM CARBON PORTAL AND BUILD AN OUTPUT READY FOR BEING PROCESSED BY THE FFP FUNCTION - ##
## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##

CPtoFFPinput = function(zip.file = NULL,             # ZIP file path
                        FFP.input.table = NULL,      # Existing input table csv 
                        start.date = NULL,           # Date start in the format YYYY-MM-DD
                        end.date = NULL,             # Date end in the format YYYY-MM-DD
                        save.input.table = TRUE,     # Save the output data as csv
                        MDS.file.list = NULL,        # Inlcude MDS-derived indices for gap-filled data
                        ERA5.file.list = NULL,       # Optional, full path of yearly ERA5 reanalysis data, used internally to retrieve PBLH
                        out.dir = NULL)              # If not null, the directory of the out file folder, otherwise the working directory is used 
{

  
  # Stop the function if no input table is provided
  if (is.null(FFP.input.table) & is.null(zip.file))
  {
    stop("No ZIP file or input table were provided, please specify one of them!")
  }
  
  
  ##
  ## [1] PREPARATORY PHASE -----
  ##
  
  # Install pacman and required packages if they are not already installed
  if (!require("pacman", quietly = T)){install.packages("pacman")}
  pacman::p_load(crayon, dplyr, tidyr, purrr, lubridate, stringr, sf, tibble, readr, zoo)
  
  # Error messages
  error <- crayon::red
  warn <- crayon::yellow
  note <- crayon::silver
  prog.mes <- crayon::cyan
  
  
  # -- Build input DF from zip file -- #
  if (!is.null(zip.file))
    
  {
    
    # List file name within the zip folder #
    file_list <- utils::unzip(zip.file, list = T) %>% pull(Name)
    
    # -- SITE ID -- #
    site.ID <- read_csv(unz(zip.file, # Zip connection
                            file_list[str_detect(file_list, pattern = 'SITEINFO')]), col_types = cols()) %>%
      pull(SITE_ID) %>% unique()

    # Starting message
    # Defining the start message total length
    message_total_length <- 86  
    
    # Message text
    message_text <- paste0('Building the FFP input table for the ', bold(site.ID), ' station')
    message_text_length <- nchar(crayon::strip_style(message_text)) 
    
    # Real margin
    margin <- strrep(" ", max(0, floor((message_total_length - message_text_length) / 2)))
    
    # Final cat
    cat('\n', strrep('*', message_total_length), sep='')
    cat('\n', margin, prog.mes(message_text), sep='')
    cat('\n', strrep('*', message_total_length), '\n', sep='')

    
    
    # -- COORDINATES -- #
    Coords <- read_csv(unz(zip.file, # Zip connection
                           file_list[grepl(file_list, pattern = 'SITEINFO')]), col_types = 'ccccc')
    
    Coords.df <- Coords %>% filter(VARIABLE %in% c('LOCATION_LONG', 'LOCATION_LAT')) %>% 
      select(VARIABLE, DATAVALUE) %>%
      mutate(VARIABLE = ifelse(str_detect(string = VARIABLE, pattern = 'LAT'), yes = 'lat', no = 'lon'), 
             DATAVALUE = as.numeric(DATAVALUE))
    
    # Store coordinates as a named vector
    Coords.vec <- Coords.df %>% pull(DATAVALUE, VARIABLE)
    

    
    # -- UTC OFFSET -- #
    UTC.vec <- Coords %>% filter(VARIABLE %in% c('UTC_OFFSET')) %>% 
      mutate(DATAVALUE=as.numeric(DATAVALUE)) %>% 
      pull(DATAVALUE, VARIABLE) %>% unname()
    
    
    
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
    
    
    
    ##
    ## [2] DATA AVAILABILITY CHECK -----
    ##
    
    cat(prog.mes(paste0('\n', bold('[2]'),' Data availability check\n')))
    
    
    
    # -- FLUXES VARIABLES -- #
    
    # Create an empty dataframe (to prevent failures of the main dataframe processing part) 
    DF_Input=data.frame()
    
    Fluxes <- read_csv(unz(zip.file, # Zip connection
                           file_list[str_detect(file_list, pattern = 'FLUXES') & 
                                       str_detect(file_list, pattern = 'VARINFO', negate = T)]), col_types = 'ccccc')
    
    Fluxvar <- c('WS', 'MO_LENGTH', 'V_SIGMA', 'USTAR', 'WD', 'PBLH')
    
    
    
    ## -- ## ERA5 CHUNK ## -- ## 
    
    if(!(is.null(ERA5.file.list))) { # File list #
      
      # Starting message # 
      cat(prog.mes(paste0('\nIntegrating ERA5 PBLH for the ', site.ID, ' station.\n')))
      
      # Select only the years of interest # 
      Data_years <- str_sub(string=Fluxes$TIMESTAMP_START, start=1, end=4) %>% unique() %>% as.numeric()
      
      # Filter the ERA list according to the data years # 
      ERA5list_filtered <- ERA5.file.list[ERA5.file.list %>% str_detect(paste(Data_years, collapse = "|"))]
            
      # Filter the ERA list according to the SITE_CODE # 
      ERA5list_filtered <- ERA5list_filtered[ERA5list_filtered %>% str_detect(site.ID)]
      
      # Data reading # Select only time and value # 
      PBLH_ERA5 <- read_csv(ERA5list_filtered, col_types = cols()) %>% 
        select(any_of(c("time", "value")))
      
      # Check if the columns are present # 
      if(all(colnames(PBLH_ERA5) %in% c("time", "value"))) { 
        
        # Check the correct format #
        if(all(is.POSIXct(PBLH_ERA5$time), is.numeric(PBLH_ERA5$value)))  {
          
          # - 1 - # TABLE CREATION #
          
          # Create a dummy DF to fill the half hours # add the last obs # The last one is repeated # Not optimal #   
          PBLH_ERA5_30min <- tibble(time=seq.POSIXt(from = min(PBLH_ERA5$time), to = max(PBLH_ERA5$time), by = "30 min"))
          
          # Join data # 
          PBLH_ERA5_full <- PBLH_ERA5_30min %>% left_join(PBLH_ERA5 %>% select(value, time), by="time")
          
          # Approx # 
          PBLH_ERA5_full <- PBLH_ERA5_full %>% mutate(value=na.approx(value))
          
          # If the last obs is 23:00, add one repeated row # 
          if(hour(max(PBLH_ERA5_full$time))==23)  {
            
            PBLH_ERA5_full <- PBLH_ERA5_full %>% 
              bind_rows(tibble(time=max(PBLH_ERA5_full$time)+minutes(30), 
                               value=PBLH_ERA5_full %>% filter(PBLH_ERA5_full$time==max(PBLH_ERA5_full$time)) %>% 
                                 pull(value)))
          }
          
          # Convert to timestamp format # 
          PBLH_ERA5_full <- PBLH_ERA5_full %>% mutate(TIMESTAMP_START=format(time, format = "%Y%m%d%H%M"))
          
          # If the pblh column is null, create one with all NA (-9999) values # 
          
          if(!("PBLH" %in% colnames(Fluxes)))  {
          
          ## Add fluxes 
          Fluxes=Fluxes %>% mutate(PBLH = -9999)
          
          }
          
          # Create a DF PBLH # 
          PBLH_DF <- PBLH_ERA5_full %>% 
            rename("PBLH_ERA"="value") %>% 
            left_join(Fluxes %>% select(TIMESTAMP_START, PBLH), by="TIMESTAMP_START") %>% 
            rename("PBLH_DATA"="PBLH") %>% 
            select(-time)
          
          # Check if ERA is always ok and make a flag # 
          PBLH_DF <- PBLH_DF %>% 
            mutate(PBLH_QC=case_when(PBLH_ERA == -9999 & PBLH_DATA == -9999 ~ 4,   # 4: ERA and model absent
                                     PBLH_ERA != -9999 ~ 2,                        # 2: ERA ok
                                     PBLH_DATA == -9999 & PBLH_DATA != -9999 ~ 3)) # 3: DATA ok 
          
          # PBLH Choose # PBLH measured # 2 # PBLH from ERA # 3 # PBLH modelled #
          PBLH_DF <- PBLH_DF %>% 
            mutate(PBLH=case_when(PBLH_QC == 2 ~ PBLH_ERA,
                                  PBLH_QC == 3 ~ PBLH_DATA,
                                  PBLH_QC == 4 ~ NA))
          
          # 1 # PBLH measured chunk # 
          
          
          ## TODO .... ##
          #if() # Independent mesure of PBLH # 
          
          
          # 2 # PBLH from ERA/modelled # Join the PBLH and associated data flag # 
          # When PBLH measured is not available, use the best option Measured >> ERA >> Model >> NA
          
          if(any(PBLH_DF$PBLH_QC == 2)) { # When there is at least one valid PBLH measure from ERA # 
            
            # Count the number of PBLH estimated with ERA # 
            count_QC <- PBLH_DF %>% count(PBLH_QC) 
            
            # Calculated percentage # 
            count_QC <- count_QC %>% 
              mutate(TOT=nrow(PBLH_DF)) %>% 
              mutate(Perc=n/TOT*100)
            
            # Message # 
            cat(warn(paste0('\n [!] PBLH is estimated using ERA5 data (',
                            count_QC %>% filter(PBLH_QC==2) %>% pull(Perc),
                            '% of records) for the ', site.ID, ' station.\n')))
            
          } 
            
          # Remove PBLH (modelled) from Fluxes table # 
          Fluxes <- Fluxes %>% select(-PBLH)
            
          # Join the final ERA5 PBLH # 
          Fluxes <- Fluxes %>% 
            left_join(PBLH_DF %>% select(TIMESTAMP_START, PBLH, PBLH_QC), by="TIMESTAMP_START")
          
          
          ## Folder path creation ##  
          
          if(is.null(out.dir))  {
            
            # MDS indice path #
            ERA5_dir <- paste0(getwd(), '/', site.ID, '/', 'Input', "/", "ERA5")
            
          } else  {
            
            # MDS indice path #
            ERA5_dir <- paste0(out.dir, '/', site.ID, '/', 'Input', "/", "ERA5")
            
          }
          
        } # Correct format chunk ending #
        
      } else { # No time and value columns #  
          
        # Message # 
        cat(warn(paste0("\n [!] PBLH can't be estimated by ERA5 data (not available or data names don't match the requirements) for the ", site.ID, " station.\n")))
        
      }
           
    } # PBLH chunk END # QC will not exported # 
      
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
    
    
    
    # -- HC -- #
    
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
    
    
    
    ##
    ## [3] DATA PROCESSING -----
    ##
    
    cat(prog.mes(paste0('\n', bold('[3]'),' Data processing\n')))
    
    
    
    # -- FLUXES VARIABLES -- #
    
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
      
      # Database optional chunk # 
      #if(HC_get_DB)  {} # Get the canopy height 
      # else  {}
      
    
        
      # -- HC -- #
      
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
      
      
      
      # -- HM -- #
      
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
      
      
      
      ##
      ## [4] FINAL JOIN AND DF REFINING -----
      ##

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
      
      # Initialize the MDS object #
      MDS.file.site=NULL
      
      if(!is.null(MDS.file.list))
        
      {

        if(is.null(out.dir))  {
          
          # MDS indice path #
          Idx_dir <- paste0(getwd(), '/', site.ID, '/', 'Input', "\\", "MDS")
          
        } else  {
          
          # MDS indice path #
          Idx_dir <- paste0(out.dir, '/', site.ID, '/', 'Input', "\\", "MDS")
            
        }
        
        # Identify the MDS files and their type # 
        MDS.file.site=MDS.file.list[MDS.file.list %>% grepl(pattern=site.ID)]
        MDS.file.site.type=MDS.file.list[MDS.file.list %>% grepl(pattern=site.ID)] %>% str_split(pattern="_") %>% 
          map(function(x) tail(x, n = 1)) %>% str_replace(pattern = ".csv", replacement = "")
        
        # Il the file list is not empty, process the data # Combine CUT and VUT # Prefer VUT to CUT #
        if(length(MDS.file.site)>0)  {
          
          # Starting message # 
          cat(prog.mes(paste0('\nIntegrating MDS gap-filling indices for the ', site.ID, ' station.\n')))
           
          
          # Reading the Fluxnet (get the correct time-span) # 
          Fluxnet <- read_csv(unz(zip.file, # Zip connection
                                 file_list[str_detect(file_list, pattern = 'FLUXNET_HH') & 
                                             str_detect(file_list, pattern = 'VARINFO', negate = T) & 
                                             str_detect(file_list, pattern = 'product', negate = T)]), col_types = 'ccccc',
                              col_select = c("TIMESTAMP_START")) %>% pull(TIMESTAMP_START) %>% 
            as.POSIXct(format = "%Y%m%d%H%M", tz="GMT")
          
          
          # --- IDX BUILD --- #
          # Fluxnet-based # 
          Years=Fluxnet %>% year() %>% unique()
          DF_Idx=seq(as.POSIXct(paste0(min(Years), "-01-01 00:00:00"), format = "%Y-%m-%d %H:%M", tz="GMT"), 
                       as.POSIXct(paste0(max(Years), "-12-31 23:30:00"), format = "%Y-%m-%d %H:%M", tz="GMT"), 
                       by = "30 mins") %>% tibble(TIMESTAMP=.) %>% 
            mutate(Idx=row_number(TIMESTAMP)-1)  # Take zero into account # 
          
          # Add rownumber to fluxes (Index) # Remove where indices are NA # 
          DF_Input=DF_Input %>% left_join(DF_Idx, by = "TIMESTAMP") %>% 
            filter(!is.na(Idx))
          
          # MDS output file reading #
          if(length(MDS.file.site)==2)  {
            
            # File reading and CUT-VUT merging # 
            MDS.file.site.df=map2(.x = MDS.file.site, .y = MDS.file.site.type, .f = ~ .x %>% read_csv(col_types = "iic") %>% 
                                                          rowwise() %>% mutate(MDS_Idx=str_split(ZERO_BASED_INDICES, pattern=",")) %>% 
                                                          mutate(MDS_Idx=list(as.numeric(MDS_Idx))) %>% 
                                                          select(-ZERO_BASED_INDICES, -SAMPLES_COUNT) %>% 
                                                          rename("Idx" = "ZERO_BASED_INDEX") %>% 
                   mutate(type=.y)) %>% 
              ## Prefer VUT to CUT ##
              bind_rows() %>% 
              group_by(Idx) %>%
              slice_min(match(type, c("y", "c")), n = 1) %>%
              ungroup()
            
            # Incorporate MDS gap-filling indices in the dataframe for each half-hour #
            DF_Input=DF_Input %>% 
              left_join(MDS.file.site.df, by = c("Idx"))
            
          } else if(length(MDS.file.site)==1) {
            
            # File reading and CUT-VUT merging # 
            MDS.file.site.df=map2(.x = MDS.file.site, .y = MDS.file.site.type, .f = ~ .x %>% read_csv(col_types = "iic") %>% 
                                    rowwise() %>% mutate(MDS_Idx=str_split(ZERO_BASED_INDICES, pattern=",")) %>% 
                                    mutate(MDS_Idx=list(as.numeric(MDS_Idx))) %>% 
                                    select(-ZERO_BASED_INDICES, -SAMPLES_COUNT) %>% 
                                    rename("Idx" = "ZERO_BASED_INDEX") %>% 
                                    mutate(type=.y)) %>% 
              bind_rows()
            
            # Incorporate MDS gap-filling indices in the dataframe for each half-hour #
            DF_Input=DF_Input %>% 
              left_join(MDS.file.site.df, by = c("Idx"))
            
          } else {
            
            cat(warn(paste0("\n [!] MDS indices files are not available for the ", site.ID, " station.\n")))
            
          }
          
          # Remove the type column # Not needed anymore
          DF_Input=DF_Input %>% select(-type)
          
        }
      
      } ## MDS chunk end 
      
      
      
      ##
      ## [5] DATE FILTERING AND OUTPUT -----
      ##
      
      
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
    Model_Input[['UTC_offset']] <- UTC.vec
    Model_Input[['PID']] <- PID
    Model_Input[['FFP_input_parameters']] <- DF_Input
    
    # Add the data as data frame
    Model_DF <- tibble(
      site_id = site.ID,
      lat = Coords.vec['lat'] %>% as.numeric(),
      lon = Coords.vec['lon'] %>% as.numeric(),
      UTC_offset=UTC.vec,
      PID = PID) %>%
      bind_cols(DF_Input)
    
    # ZIP file chunk ending...

    
    ##
    ## [6] EXISTING TABLE PROCESSING -----
    ##
    
        
  } else if (!is.null(FFP.input.table))
    
  {
    
    cat(prog.mes(paste0('\n', bold('[6]'),' Existing table processing\n')))
    
    # Get site code
    site.ID <- str_split(basename(FFP.input.table), pattern = '_')[[1]][1]
    
    cat(warn(paste0('\nNo zip file provided for the ', site.ID, ' station.')))
    cat(prog.mes(paste0('\nUsing FFP input table as input and built-in table for coordinates....')))
    cat(prog.mes(paste0('\nCreating internal list with filtered dataframe, site code, coordinates and, if available, PID....\n')))
    
    
    # -- DATA READING -- #
    DF_Input <- read_csv(FFP.input.table, col_types = cols())
    
    # Check if the MDS indices is in the column list # 
    if(any(colnames(DF_Input) %in% "MDS_Idx"))  {
      
      # Read the file again # 
      DF_Input <- read_csv(FFP.input.table, col_types = cols(.default = col_guess(), MDS_Idx = col_character()))
      
      # Change the string column in a vector list # 
      DF_Input <- DF_Input %>%
        mutate(MDS_Idx = MDS_Idx %>% str_replace_all("NULL|NA", "") %>% str_split(",") %>% map(as.numeric))
    }
    
    
    # Verify that the site_id is available #
    stopifnot("site_id is missing from the input DF" = "site_id" %in% colnames(DF_Input))
    
    # If one of the input site parameters is missing, retrieve it from a lookup table # 
    if(!(all(c("site_id", "lat", "lon", "UTC_offset") %in% colnames(DF_Input))))  {
      
      # Get built-in coordinates
      Coords <- data.frame(
        site_id = c('BE-Bra', 'BE-Dor', 'BE-Lcr', 'BE-Lon', 'BE-Maa', 'BE-Vie', 'CD-Ygb', 'CH-BaK', 'CH-Dav', 'CZ-BK1', 
                    'CZ-Lnz', 'CZ-wet', 'DE-Amv', 'DE-BeR', 'DE-Brs', 'DE-Fen', 'DE-Geb', 'DE-Gri', 'DE-GsB', 'DE-Gwg', 
                    'DE-Hai', 'DE-Har', 'DE-HoH', 'DE-Hzd', 'DE-Kli', 'DE-Msr', 'DE-RuR', 'DE-RuS', 'DE-RuW', 'DE-Tha', 
                    'DK-Gds', 'DK-RCW', 'DK-Skj', 'DK-Sor', 'DK-Vng', 'ES-LMa', 'FI-Hyy', 'FI-Ken', 'FI-Kmp', 'FI-Kvr', 
                    'FI-Let', 'FI-Lom', 'FI-Sii', 'FI-Sod', 'FI-Tvm', 'FI-Var', 'FR-Aur', 'FR-Bil', 'FR-CLt', 'FR-EM2', 
                    'FR-FBn', 'FR-Fon', 'FR-Gri', 'FR-Hes', 'FR-Lam', 'FR-LGt', 'FR-Lqu', 'FR-Lus', 'FR-Mej', 'FR-Pue', 
                    'FR-Tou', 'GF-Guy', 'GL-Dsk', 'GL-NuF', 'GL-ZaF', 'GL-ZaH', 'GR-HeK', 'GR-HeM', 'IE-Cra', 'IT-BCi', 
                    'IT-BFt', 'IT-Cp2', 'IT-Lsn', 'IT-MBo', 'IT-Niv', 'IT-Noe', 'IT-OXm', 'IT-PCm', 'IT-Ren', 'IT-Sas', 
                    'IT-SR2', 'IT-Tor', 'IT-TrF', 'NL-Loo', 'NO-Hur', 'SE-Deg', 'SE-Htm', 'SE-Myc', 'SE-Nor', 'SE-Oes', 
                    'SE-Sto', 'SE-Svb', 'UK-AMo'),
        lat = c(51.30761, 50.31188, 51.11218, 50.55162, 50.97987, 50.30496, 0.81444, 47.56173, 46.81533, 49.50208, 48.68155, 
                49.02465, 52.1757, 52.45723, 52.29663, 47.8329, 51.09973, 50.95001, 52.02965, 47.57083, 51.07921, 47.93391, 
                52.08656, 50.96381, 50.89304, 47.80918, 50.62191, 50.86591, 50.50493, 50.9626, 56.0737, 55.68068, 55.91273, 
                55.48587, 56.03748, 39.94033, 61.84741, 67.98721, 60.20289, 61.84662, 60.64183, 67.99724, 61.83265, 67.36239, 
                59.8418, 67.7549, 43.54965, 44.49365, 45.0414, 49.87211, 43.24079, 48.47636, 48.84422, 48.6741, 43.49644, 
                47.32292, 45.6444, 46.41425, 48.1184, 43.7413, 43.57285, 5.2787, 69.25349, 64.13093, 74.48152, 74.4733, 
                35.33613, 35.32343, 53.32308, 40.52375, 45.19776, 41.70427, 45.74048, 46.01468, 45.49091, 40.60618, 43.77446, 
                40.87728, 46.58686, 40.71696, 43.73202, 45.84444, 45.82376, 52.16645, 60.37163, 64.18203, 56.09763, 58.36503, 
                60.0865, 57.4301, 68.35594, 64.25611, 55.79255),
        lon = c(4.51984, 4.96811, 3.85043, 4.74623, 5.63185, 5.9981, 24.50247, 7.58049, 9.85591, 18.53688, 16.94633, 14.77035, 
                6.95537, 13.31583, 10.44871, 11.0607, 10.91463, 13.51253, 11.10478, 11.0326, 10.45217, 7.59814, 11.22235, 
                13.48978, 13.5223, 11.45617, 6.30413, 6.44714, 6.33096, 13.56533, 9.3341, 12.1014, 8.40481, 11.64464, 9.16071, 
                -5.77465, 24.29477, 24.24301, 24.9611, 24.2804, 23.95952, 24.20918, 24.19285, 26.63859, 23.2503, 29.61, 1.1061, 
                -0.95609, 6.41053, 3.02065, 5.67865, 2.7801, 1.95191, 7.06465, 1.23788, 2.2841, 2.7349, 0.12065, -1.79635, 
                3.5957, 1.37474, -52.9248, -53.51413, -51.38607, -20.55577, -20.55087, 25.1328, 25.13017, -7.64177, 14.95744, 
                10.74197, 12.35729, 12.7503, 11.04583, 7.13943, 8.15169, 11.25511, 14.2559, 11.43369, 8.57596, 10.29091, 7.57806, 
                7.56089, 5.74355, 11.07949, 19.55654, 13.41897, 12.1694, 17.4795, 18.98415, 19.04521, 19.7745, -3.24369), 
        UTC_offset=c(1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                     1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -3, -3, -3, 0, 0, 2, 
                     2, 0, 1, 1, 1, 1, 1, 1, +1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0))
      
      # Identify the missing cols #
      target_cols <- c("lat", "lon", "UTC_offset")
      missing_cols <- target_cols[!(target_cols %in% colnames(DF_Input))]
      
      # Select only the missing cols # 
      Coords_to_join <- Coords %>% 
        select(all_of(c("site_id", missing_cols)))
      
      # Selective left join using only the missing columns # 
      DF_Input <- DF_Input %>% 
        left_join(Coords_to_join, by = "site_id")

      cat(warn(paste0("\n[!] ", paste(colnames(Coords_to_join %>% select(-site_id)), collapse = ", "), 
                          "added using a lookup table")))
      
    }
    
    # -- DATA FILTERING -- #
    
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
    
    if(is.null(out.dir))  {
      
      # Extract PID from TOC file, located everywhere inside the project folder
      PID <- NULL
      TOC <- paste0(dirname(list.files(path = getwd(), pattern = paste0('ICOSETC_', site.ID, '_ARCHIVE'), 
                                       full.names = T, recursive = T)), '/!TOC.csv')
      
    } else  {
      
      # Extract PID from TOC file, located everywhere inside the defined folder
      PID <- NULL
      TOC <- paste0(dirname(list.files(path = out.dir, pattern = paste0('ICOSETC_', site.ID, '_ARCHIVE'), 
                                       full.names = T, recursive = T)), '/!TOC.csv')
    }
    
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
    Model_Input[['tower_coordinates']] <- c("Lat"=DF_Input$lat %>% unique() %>% as.numeric(), 
                                            "Lon"=DF_Input$lon %>% unique() %>% as.numeric())
    Model_Input[['UTC_offset']] <- DF_Input$UTC_offset %>% unique() %>% as.numeric()
    Model_Input[['PID']] <- PID
    Model_Input[['FFP_input_parameters']] <- DF_Input %>% select(- any_of(c('site_id', 'lat', 'lon', 'PID', 'UTC_offset')))
    
    
    Model_DF <- DF_Input %>% 
      mutate(lat=as.numeric(lat), 
             lon=as.numeric(lon), 
             UTC_offset=as.numeric(UTC_offset), 
             PID=PID) %>% 
      select(site_id, lat, lon, UTC_offset, PID, colnames(DF_Input %>% 
                                                            select(- any_of(c('site_id', 'lat', 'lon', 'UTC_offset', 'PID')))))

  } else if (is.null(FFP.input.table) & is.null(zip.file))
    
  {
    cat(warn(paste0('\nNo zip file or dataframe provided for the station. Please specify it.')))
  }
  
  
  
  ## -- SAVING OPTIONS -- ##
  
  # Save input tables in the FFP input folder
  if (all(save.input.table & is.data.frame(Model_Input[['FFP_input_parameters']])))
    
  {
    
    if(is.null(out.dir))  {
      
      # In the project dir, create a folder with site code and input dir inside #
      Site_dir <- paste0(getwd(), '/', site.ID)
      Input_dir <- paste0(getwd(), '/', site.ID, '/', 'Input')
      
    } else  {
      
      # In the specified dir, create a folder with site code and input dir inside #
      Site_dir <- paste0(out.dir, '/', site.ID)
      Input_dir <- paste0(out.dir, '/', site.ID, '/', 'Input')
      
    }
    
    if (!dir.exists(Site_dir)) {dir.create(Site_dir)}
    if (!dir.exists(Input_dir)) {dir.create(Input_dir)}
    
    ## Write the MDS input files ##
    if(exists("MDS.file.site")) {

      if (!dir.exists(Idx_dir)) {dir.create(Idx_dir)}

      # Copy the input tables in the MDS folder #
      MDS.file.site %>% map(function(x) file.copy(from = x,
                                                  to = paste0(Idx_dir, '/', basename(x))))

      message(prog.mes('\nMDS input files saved'))

      # For saving the table later # Convert the MDS idx column to text before writing #
      Model_DF <- Model_DF %>%
        mutate(MDS_Idx = map_chr(MDS_Idx, function(x) {
          if (is.null(x) || all(is.na(x))) {return(NA_character_) # If NULL or NA, return empty cell
          } else {return(paste(x, collapse = ",")) # collapse the names
          }}))
    }

    # ## Write the ERA5 input files ## 
    if(exists("ERA5list_filtered")) {

      if (!dir.exists(ERA5_dir)) {dir.create(ERA5_dir)}

      # Copy the input tables in the ERA5 folder #
      ERA5list_filtered %>% map(function(x) file.copy(from = x,
                                                  to = paste0(ERA5_dir, '/', basename(x))))

      message(prog.mes('\nERA5 input files saved'))

    }
    
    
    # Final(ly) # Write the CSV file
    write_csv(Model_DF, file = paste0(Input_dir, '/', 
                                      site.ID, '_', min(as_date(Model_DF$TIMESTAMP)) %>% str_replace_all('-', ''), '_',
                                      max(as_date(Model_DF$TIMESTAMP)) %>% str_replace_all('-', ''), '_FFPinput.csv'))
    
    message(prog.mes('\nFFP input table saved'))
    
    
  }
  
  
  
  # -- RETURN INPUT TABLE -- #
  return(Model_Input)
  
} # Function ending