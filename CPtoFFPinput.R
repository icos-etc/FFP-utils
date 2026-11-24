## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##
## - MANIPULATE DATA FROM CARBON PORTAL AND BUILD AN OUTPUT READY FOR BEING PROCESSED BY THE FFP FUNCTION - ##
## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ## -- ##

CPtoFFPinput = function(zip.file = NULL,             # ZIP file path
                        FFP.input.table = NULL,      # Existing input table csv
                        start.date = NULL,           # Date start in the format YYYY-MM-DD
                        end.date = NULL,             # Date end in the format YYYY-MM-DD
                        save.input.table = TRUE,     # Save the output data as csv
                        MDS.file = NULL,             # Include MDS-derived indices for gap-filled data
                        ERA5.file.list = NULL,       # Optional, full path of yearly ERA5 PBLH files (time in local standard time of the site, as ICOS data), used internally to retrieve PBLH
                        out.dir = NULL,              # If not null, the directory of the out file folder, otherwise the working directory is used
                        merge.dataset=NULL)

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
  pacman::p_load(crayon, dplyr, tidyr, purrr, lubridate, stringr, sf, tibble, readr, zoo, httr, jsonlite)

  # Objects created only under some conditions # Initialized here, so that objects with the same name in the session are never used #
  ERA5list_filtered <- NULL

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

    # Meteo table # Read only if WS/WD are not in the fluxes #
    Meteo <- NULL

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


    # -- PID -- # From the TOC file next to the zip, matched on the zip file name (unique product) #
    PID <- NULL
    TOC <- file.path(dirname(zip.file), '!TOC.csv')

    if (file.exists(TOC))

    {
      # Read the TOC file and extract the PID of this zip file
      PID <- read_csv(TOC, name_repair = function(x) x %>% str_replace(pattern = ' ', replacement = '_'), col_types = cols()) %>%
        filter(File_name == basename(zip.file)) %>%
        pull(PID) %>% unique()

      # The zip file must appear exactly once in the TOC file
      if (length(PID) != 1)
      {
        cat(warn(paste0('\n [!] A unique PID cannot be retrieved from the TOC file for the ', basename(zip.file), ' product.\n')))
        PID <- NULL
      }

    } else
    {
      cat(warn(paste0('\n [!] PID cannot be retrieved for the ', site.ID, ' station (no !TOC.csv next to the zip file).\n')))
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

      # Filter the ERA list according to the data years and the SITE_CODE # Year taken from the file name (field "__YYYY__") #
      ERA5_years <- as.numeric(str_extract(basename(ERA5.file.list), "(?<=__)\\d{4}(?=__)"))
      ERA5list_filtered <- ERA5.file.list[ERA5_years %in% Data_years & str_detect(basename(ERA5.file.list), site.ID)]

      # Same year available as ERA5 (final) and ERA5T (preliminary): ERA5 is kept #
      ERA5_is_T <- str_detect(basename(ERA5list_filtered), "__ERA5T__")
      ERA5_yrs <- as.numeric(str_extract(basename(ERA5list_filtered), "(?<=__)\\d{4}(?=__)"))
      ERA5list_filtered <- ERA5list_filtered[!(ERA5_is_T & ERA5_yrs %in% ERA5_yrs[!ERA5_is_T])]

      if(length(ERA5list_filtered) == 0) {
        cat(warn(paste0("\n [!] No ERA5 file matches the ", site.ID, " station and the data years (",
                        paste(Data_years, collapse = ", "), "): the modelled PBLH is used.\n")))
      }

      # Data reading # Select only time and value # Duplicated hours removed #
      PBLH_ERA5 <- read_csv(ERA5list_filtered, col_types = cols()) %>%
        select(any_of(c("time", "value")))
      if(all(c("time", "value") %in% colnames(PBLH_ERA5))) {
        PBLH_ERA5 <- PBLH_ERA5 %>% distinct(time, .keep_all = TRUE) %>% arrange(time)
      }

      # Check if the columns are present #
      if(all(colnames(PBLH_ERA5) %in% c("time", "value"))) {

        # Check the correct format #
        if(all(is.POSIXct(PBLH_ERA5$time), is.numeric(PBLH_ERA5$value)))  {

          # - 1 - # TABLE CREATION #

          # Create a dummy DF to fill the half hours # add the last obs # The last one is repeated # Not optimal #
          PBLH_ERA5_30min <- tibble(time=seq.POSIXt(from = min(PBLH_ERA5$time), to = max(PBLH_ERA5$time), by = "30 min"))

          # Join data #
          PBLH_ERA5_full <- PBLH_ERA5_30min %>% left_join(PBLH_ERA5 %>% select(value, time), by="time")

          # Approx # Only the half-hour between two consecutive hours # Longer gaps stay NA (modelled PBLH used) #
          PBLH_ERA5_full <- PBLH_ERA5_full %>% mutate(value=na.approx(value, maxgap = 1, na.rm = FALSE))

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

          # 1 # PBLH measured chunk #


          ## TODO .... ##
          #if() # Independent measure of PBLH #


          # 2 # PBLH from ERA/modelled # For every half-hour of the fluxes (no half-hour is lost) #
          # When PBLH measured is not available, use the best option Measured >> ERA >> Model >> NA
          # The modelled PBLH (ETC) is kept and the ERA5 PBLH is added next to it #
          Fluxes <- Fluxes %>%
            left_join(PBLH_ERA5_full %>% select(TIMESTAMP_START, PBLH_ERA = value), by = "TIMESTAMP_START") %>%
            mutate(PBLH_QC = case_when(!is.na(PBLH_ERA) & PBLH_ERA != -9999 ~ 2,     # 2: ERA ok
                                       !is.na(PBLH) & PBLH != -9999 ~ 3,             # 3: DATA (modelled) ok
                                       TRUE ~ 4),                                    # 4: ERA and model absent
                   PBLH = case_when(PBLH_QC == 2 ~ PBLH_ERA,
                                    PBLH_QC == 3 ~ PBLH,
                                    TRUE ~ -9999)) %>%                               # -9999: converted to NA later, the row is kept
            select(-PBLH_ERA)

          # Message # Share of half-hours for each PBLH source #
          count_QC <- Fluxes %>% count(PBLH_QC) %>% mutate(Perc = round(n / sum(n) * 100, 1))
          Perc_QC <- function(q) {p <- count_QC %>% filter(PBLH_QC == q) %>% pull(Perc); if(length(p) == 0) 0 else p}

          cat(warn(paste0('\n [!] PBLH for the ', site.ID, ' station: ERA5 ', Perc_QC(2), '% | modelled (ETC) ', Perc_QC(3),
                          '% | not available ', Perc_QC(4), '% of records.\n')))


          ## Folder path creation ##

          if(is.null(out.dir))  {

            # MDS indices path #
            ERA5_dir <- paste0(getwd(), '/', site.ID, '/', 'Input', "/", "ERA5")

          } else  {

            # MDS indices path #
            ERA5_dir <- paste0(out.dir, '/', site.ID, '/', 'Input', "/", "ERA5")

          }

        } else if(nrow(PBLH_ERA5) > 0) { # Wrong format #

          cat(warn(paste0("\n [!] ERA5 time column not recognised as date-time (or value not numeric) for the ", site.ID,
                          " station: the modelled PBLH is used.\n")))

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
    if(!is.null(Meteo))
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

      } else {

        # Canopy height group in wide format (moved here from the data processing part) #
        Ancillary.df <- Ancillary %>%
          split(.$VARIABLE_GROUP) %>%
          map( ~ .x %>% spread(key = VARIABLE, value = DATAVALUE))

        # Force DF to have date and spp columns, filled with NAs
        cols <- c(HEIGHTC_DATE = NA, HEIGHTC_SPP = NA, HEIGHTC_DATE_START = NA, HEIGHTC_VEG_STATUS = "Alive",
                  HEIGHTC_STATISTIC = NA, HEIGHTC = NA)

        # DF filter to keep only HC
        GRP_HEIGHTC_DF <- Ancillary.df[['GRP_HEIGHTC']] %>%
          tibble::add_column(., !!!cols[!names(cols) %in% names(.)]) # Add columns if not available #

        # Which statistic can be used: site level (no species), alive vegetation, with a value #
        hc_available <- function(stat) any(GRP_HEIGHTC_DF$HEIGHTC_STATISTIC %in% stat & is.na(GRP_HEIGHTC_DF$HEIGHTC_SPP) &
                                           GRP_HEIGHTC_DF$HEIGHTC_VEG_STATUS %in% 'Alive' & !is.na(GRP_HEIGHTC_DF$HEIGHTC))

        hc_statistic <- if (hc_available('90th Percentile')) '90th Percentile' else if (hc_available('Mean')) 'Mean' else NA

        # No site-level canopy height: the table is not created (as for the other missing data) #
        if (is.na(hc_statistic)) {
          cat(warn(paste0("\n [!] Canopy height can't be retrieved (no site-level 90th percentile or mean) for the ",
                          site.ID, " station. FFP input table is not created.\n")))
          DF_Input <- NULL
        }

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
      if(!is.null(Meteo))
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

      # If HC 90th percentile is available (site level) # Decided in the data availability check #
      if(hc_statistic == '90th Percentile')

      {
        hc.df.tmp <- GRP_HEIGHTC_DF %>%
          filter(HEIGHTC_STATISTIC == '90th Percentile' & is.na(HEIGHTC_SPP) & HEIGHTC_VEG_STATUS == 'Alive') %>%
          mutate(HEIGHTC_DATE = case_when(is.na(HEIGHTC_DATE) ~ HEIGHTC_DATE_START, # sometimes HEIGHTC_DATE is NA
                                          TRUE ~ HEIGHTC_DATE)) %>%
          select(HEIGHTC_DATE, HEIGHTC) %>%
          rename('TIMESTAMP' = 'HEIGHTC_DATE', 'hc' = 'HEIGHTC') %>%
          mutate(TIMESTAMP=case_when(nchar(TIMESTAMP)==6 ~ paste0(TIMESTAMP, '15'), # Missing day # As in the mean branch
                                     nchar(TIMESTAMP)==4 ~ paste0(TIMESTAMP, '0615'), # Missing month
                                     T ~ TIMESTAMP)) %>%
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
      # The grid starts from the oldest date between the hm configuration and the fluxes (hm set before the first flux is not lost) #
      hm.df <- data.frame(TIMESTAMP = seq.POSIXt(min(c(as.POSIXct(hm.df.tmp$TIMESTAMP), as.POSIXct(Fluxes.df$TIMESTAMP))), max(Fluxes.df$TIMESTAMP), '30 min')) %>%
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

    # ERA5 PBLH is integrated only when the data are built from the zip file #
    if(!is.null(ERA5.file.list)) {
      cat(warn(paste0('\n [!] ERA5.file.list is used only with zip.file: the PBLH of the input table is kept for the ', site.ID, ' station.')))
    }
    cat(prog.mes(paste0('\nUsing FFP input table as input and built-in table for coordinates....')))
    cat(prog.mes(paste0('\nCreating internal list with filtered dataframe, site code, coordinates and, if available, PID....\n')))


    # -- DATA READING -- #
    DF_Input <- read_csv(FFP.input.table, col_types = cols())

    #DF_Input=DF_Input %>% select(-Idx, -MDS_Idx, -type)


    # Check if the MDS indices is in the column list #
    if(any(colnames(DF_Input) %in% "MDS_Idx"))  {

      # Read the file again #
      DF_Input <- read_csv(FFP.input.table, col_types = cols(.default = col_guess(), MDS_Idx = col_character()))

      # Change the string column in a vector list #
      DF_Input <- DF_Input %>%
        mutate(MDS_Idx = map(MDS_Idx, function(x) if (is.na(x) || x == "NULL") NULL else as.numeric(str_split(x, ",")[[1]])))
    }

    # Change nodata values (from -9999 to NA), as in the ZIP route # Only numeric columns (MDS_Idx is a list) #
    DF_Input <- DF_Input %>% mutate(across(where(is.numeric), ~ na_if(.x, -9999)))

    # IF site.ID is not in the DF, add it after a checking #
    if (nchar(site.ID) != 6) {
      stop(paste("Error: site.ID must have exactly 6 characters, but it has", nchar(site.ID), "('", site.ID, "')"))
    }

    if (!"site_id" %in% colnames(DF_Input)) {
      DF_Input$site_id <- site.ID
    }

    # Verify that the site_id is available #
    stopifnot("site_id is missing from the input DF" = "site_id" %in% colnames(DF_Input))

    # If one of the input site parameters is missing, retrieve it from a sparql call #
    if(!(all(c("site_id", "lat", "lon", "UTC_offset") %in% colnames(DF_Input))))  {

      cat(prog.mes(paste0('\nRetrieving site info data from sparql call...')))

          query_sparql <- "
      PREFIX cpmeta: <http://meta.icos-cp.eu/ontologies/cpmeta/>
      PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
      PREFIX cpst: <http://meta.icos-cp.eu/ontologies/stationentry/>
    
      SELECT
        (coalesce(?tcId, ?hoShortName) as ?StationCode)
        (coalesce(?tcName, ?hoName) as ?Name)
        ?Theme
        (coalesce(?tcClass, ?hoClass) as ?Class)
        (coalesce(?tcCountry, ?hoCountry) as ?Country)
        (coalesce(?tcLat, ?hoLat) as ?Latitude)
        (coalesce(?tcLon, ?hoLon) as ?Longitude)
        (coalesce(?tcUtc, ?hoUtc, ?tcTz, ?hoTz) as ?UtcOffset)
      FROM <http://meta.icos-cp.eu/resources/stationentry/>
      FROM <http://meta.icos-cp.eu/ontologies/stationentry/>
      FROM <http://meta.icos-cp.eu/resources/icos/>
      WHERE {
        ?s cpst:hasCountry ?hoCountry ;
           cpst:hasLongName ?hoName ;
           cpst:hasStationClass ?hoClass0 ;
           a/rdfs:label ?Theme .
        
        OPTIONAL { ?s cpst:hasShortName ?hoShortName }
        OPTIONAL { ?s cpst:hasLat ?hoLat }
        OPTIONAL { ?s cpst:hasLon ?hoLon }
        OPTIONAL { ?s cpst:hasTimeZoneOffset ?hoUtc }
        OPTIONAL { ?s cpst:hasTimeZone ?hoTz }
    
        OPTIONAL {
          ?s cpst:hasProductionCounterpart ?prods .
          BIND(iri(?prods) AS ?prod) .
          OPTIONAL { ?prod cpmeta:hasStationId ?tcId }
          OPTIONAL { ?prod cpmeta:hasName ?tcName }
          OPTIONAL { ?prod cpmeta:hasStationClass ?tcClass }
          OPTIONAL { ?prod cpmeta:countryCode ?tcCountry }
          OPTIONAL { ?prod cpmeta:hasLatitude ?tcLat }
          OPTIONAL { ?prod cpmeta:hasLongitude ?tcLon }
          OPTIONAL { ?prod cpmeta:hasTimeZoneOffset ?tcUtc }
          OPTIONAL { ?prod cpmeta:hasTimeZone ?tcTz }
        }
    
        BIND(IF(contains(?hoClass0, 'Ass'), 'Associated', ?hoClass0) AS ?hoClass)
        FILTER EXISTS { ?s cpst:hasShortName [] }
      }
      ORDER BY ?Theme ?StationCode
      "

          # ICOS endpoint call # If the Carbon Portal is not reachable, the built-in table is used #
          Coords_CP <- tryCatch({

            response <- POST(url = "https://meta.icos-cp.eu/sparql",
                             body = list(query = query_sparql),
                             encode = "form",
                             add_headers(Accept = "application/json"))

            stop_for_status(response)

            # Parsing JSON
            json_data <- fromJSON(content(response, "text", encoding = "UTF-8"), simplifyVector = FALSE)
            rows <- json_data$results$bindings

            # Build the DF
            map_dfr(rows, function(row) {
              tibble(site_id = row$StationCode$value %||% NA_character_,
                     Theme       = row$Theme$value %||% NA_character_,
                     lat         = as.numeric(row$Latitude$value %||% NA),
                     lon         = as.numeric(row$Longitude$value %||% NA),
                     UTC_offset  = row$UtcOffset$value %||% NA_character_)
            }) %>%
              filter(Theme=="Ecosystem Station") %>% select(-Theme) %>%
              mutate(across(2:4, function(x) suppressWarnings(as.numeric(x)))) %>%   # A non-numeric UTC offset (time zone name) becomes NA
              distinct(site_id, .keep_all = TRUE)

          }, error = function(e) {

            cat(warn(paste0("\n [!] ICOS Carbon Portal (SPARQL) not reachable (", conditionMessage(e), "): the built-in table is used.\n")))
            NULL
          })

      # Built-in table (to be kept updated by hand) # Used for the values not available from the Carbon Portal #
      Coords_builtin <- data.frame(
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
        lon = c(4.51984, 4.96811, 3.85043, 4.74623, 5.63185, 5.998099, 24.50247, 7.58049, 9.85591, 18.53688, 16.94633, 14.77035,
                6.95537, 13.31583, 10.44871, 11.0607, 10.91463, 13.51253, 11.10478, 11.0326, 10.45217, 7.59814, 11.22235,
                13.48978, 13.5223, 11.45617, 6.30413, 6.44714, 6.33096, 13.56533, 9.3341, 12.1014, 8.40481, 11.64464, 9.16071,
                -5.77465, 24.29477, 24.24301, 24.9611, 24.2804, 23.95952, 24.20918, 24.19285, 26.63859, 23.2503, 29.61, 1.1061,
                -0.95609, 6.41053, 3.02065, 5.67865, 2.7801, 1.95191, 7.06465, 1.23788, 2.2841, 2.7349, 0.12065, -1.79635,
                3.5957, 1.37474, -52.9248, -53.51413, -51.38607, -20.55577, -20.55087, 25.1328, 25.13017, -7.64177, 14.95744,
                10.74197, 12.35729, 12.7503, 11.04583, 7.13943, 8.15169, 11.25511, 14.2559, 11.43369, 8.57596, 10.29091, 7.57806,
                7.56089, 5.74355, 11.07949, 19.55654, 13.41897, 12.1694, 17.4795, 18.98415, 19.04521, 19.7745, -3.24369),
        UTC_offset=c(1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                     1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -3, -3, -3, 0, 0, 2,
                     2, 0, 1, 1, 1, 1, 1, 1, +1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0))

      # Site values: Carbon Portal first, built-in table for the values not available #
      Coords <- tibble(site_id = unique(DF_Input$site_id)) %>%
        left_join(Coords_builtin %>% rename(lat_bi = lat, lon_bi = lon, UTC_offset_bi = UTC_offset), by = "site_id")

      if(!is.null(Coords_CP)) {
        Coords <- Coords %>% left_join(Coords_CP %>% rename(lat_cp = lat, lon_cp = lon, UTC_offset_cp = UTC_offset), by = "site_id")
      } else {
        Coords <- Coords %>% mutate(lat_cp = NA_real_, lon_cp = NA_real_, UTC_offset_cp = NA_real_)
      }

      # Identify the missing cols #
      target_cols <- c("lat", "lon", "UTC_offset")
      missing_cols <- target_cols[!(target_cols %in% colnames(DF_Input))]

      # For each missing value: source used (Carbon Portal >> built-in table) #
      Sources <- list()
      for (v in missing_cols) {
        cp <- Coords[[paste0(v, "_cp")]]
        bi <- Coords[[paste0(v, "_bi")]]
        Coords[[v]] <- ifelse(!is.na(cp), cp, bi)
        Sources[[v]] <- case_when(!is.na(cp) ~ "the ICOS Carbon Portal (SPARQL)", !is.na(bi) ~ "the built-in table", TRUE ~ "none")
      }

      # Message: source of each added value, for each site #
      for (i in seq_len(nrow(Coords))) {
        src_i <- sapply(Sources, function(s) s[i])
        for (s in unique(src_i)) {
          if (s == "none") {
            cat(warn(paste0("\n[!] ", Coords$site_id[i], ": ", paste(names(src_i)[src_i == s], collapse = ", "),
                            " not available (neither from the ICOS Carbon Portal nor from the built-in table).")))
          } else {
            cat(warn(paste0("\n[!] ", Coords$site_id[i], ": ", paste(names(src_i)[src_i == s], collapse = ", "), " added from ", s, ".")))
          }
        }
      }

      # Selective left join using only the missing columns #
      DF_Input <- DF_Input %>%
        left_join(Coords %>% select(all_of(c("site_id", missing_cols))), by = "site_id")
    }

    ## TIMESTAMP CHECK #
    # END timestamps (whole days): first record at 00:30 and last record at 00:00 of the following day #
    # START timestamps (whole days): first record at 00:00 and last record at 23:30 #
    first_hm <- format(min(DF_Input$TIMESTAMP), "%H:%M")
    last_hm  <- format(max(DF_Input$TIMESTAMP), "%H:%M")

    if (first_hm == "00:30" & last_hm == "00:00") {

      DF_Input <- DF_Input %>%
        mutate(TIMESTAMP = TIMESTAMP - (30 * 60)) # Subtract 30 minutes (in seconds) from the whole column

      message("END timestamp detected (first 00:30, last 00:00). Converted to START timestamp (-30 mins).")

    } else if (first_hm == "00:30" | last_hm == "00:00") { # Ambiguous: timestamps are kept as they are #

      cat(warn(paste0("\n [!] Timestamp convention not recognised for the ", site.ID, " station (first record ", first_hm,
                      ", last record ", last_hm, "): timestamps are kept as they are (START assumed).\n")))
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
    # CSV route: the PID can only come from the input table (it is the PID of the data the table was built from) #
    # No TOC search: a TOC found elsewhere is not guaranteed to refer to the same data version #
    PID <- NULL

    if("PID" %in% colnames(DF_Input))  {

      PID_table <- DF_Input$PID %>% na.omit() %>% unique()

      if(length(PID_table) == 1) {
        PID <- PID_table
      } else if(length(PID_table) > 1) {
        cat(warn(paste0('\n [!] More than one PID in the input table for the ', site.ID, ' station. PID is not assigned.\n')))
      }
    }

    if(is.null(PID))  {
      cat(warn(paste0('\n [!] PID not available in the input table for the ', site.ID, ' station.\n')))
    }


  # --- CREATE MODEL_DF AFTER ALL JOINS --- #

    # Return the output list
    Model_Input <- list()

    Model_Input[['site_id']] <- site.ID
    Model_Input[['tower_coordinates']] <- c("lat"=DF_Input$lat %>% unique() %>% as.numeric(),
                                            "lon"=DF_Input$lon %>% unique() %>% as.numeric())
    Model_Input[['UTC_offset']] <- DF_Input$UTC_offset %>% unique() %>% as.numeric()
    Model_Input[['PID']] <- PID
    Model_Input[['FFP_input_parameters']] <- DF_Input %>% select(- any_of(c('site_id', 'lat', 'lon', 'PID', 'UTC_offset')))


    Model_DF <- DF_Input %>%
      mutate(lat=as.numeric(lat),
             lon=as.numeric(lon),
             UTC_offset=as.numeric(UTC_offset),
             PID=.env$PID) %>%        # The PID chosen above, not the table column # NULL -> no PID column, as in the ZIP branch
      select(any_of(c('site_id', 'lat', 'lon', 'UTC_offset', 'PID')), everything())

  } else if (is.null(FFP.input.table) & is.null(zip.file))

  {
    cat(warn(paste0('\nNo zip file or dataframe provided for the station. Please specify it.')))
  }


  ## -- MDS CHUNK -- ## Common to the ZIP and CSV routes # Only on the whole dataset (MDS indices refer to whole years) #

  if(!is.null(MDS.file) & is.data.frame(Model_Input[['FFP_input_parameters']])) {

    if(!is.null(start.date) | !is.null(end.date)) {

      cat(warn(paste0("\n [!] MDS indices are not added for the ", site.ID,
                      " station: they refer to the whole dataset, start.date and end.date must be NULL.\n")))

    } else if(length(MDS.file) != 1) {

      cat(warn(paste0("\n [!] MDS indices are not added for the ", site.ID, " station: MDS.file must be a single file (",
                      length(MDS.file), " provided). Choose one indices file (y = VUT, preferred, or c = CUT).\n")))

    } else if(!file.exists(MDS.file)) {

      cat(warn(paste0("\n [!] MDS indices are not added for the ", site.ID, " station: MDS file not found (", MDS.file, ").\n")))

    } else {

      # Extract dynamic years from the input dataframe
      Data_Years <- Model_Input[['FFP_input_parameters']]$TIMESTAMP %>% year() %>% unique()

      # Read the data # Explicit parsing: each line is "index,count,i1,i2,..." (indices not quoted in the ONEFlux file) #
      MDS_lines <- readLines(MDS.file)
      stopifnot('Unexpected MDS file header: an indices file (ZERO_BASED_INDEX, SAMPLES_COUNT, ZERO_BASED_INDICES) is required' =
                  grepl("^\"?ZERO_BASED_INDEX\"?,\"?SAMPLES_COUNT\"?,\"?ZERO_BASED_INDICES\"?", MDS_lines[1]))
      MDS_lines <- gsub('"', '', MDS_lines[-1])            # Quotes removed, if any
      MDS_lines <- MDS_lines[MDS_lines != ""]              # Empty lines removed, if any

      MDS_indices_DF <- tibble(
        ZERO_BASED_INDEX   = as.numeric(sub(",.*$", "", MDS_lines)),                      # 1st field
        SAMPLES_COUNT      = as.numeric(sub("^[^,]*,([^,]*).*$", "\\1", MDS_lines)),      # 2nd field
        ZERO_BASED_INDICES = map(sub("^[^,]*,[^,]*,?", "", MDS_lines),                   # all the rest
                                 function(x) if (x %in% c("", "NA", "NULL")) NULL else as.numeric(str_split(x, ",")[[1]])))

      # Each row must contain as many indices as declared in SAMPLES_COUNT #
      MDS_bad <- sum(lengths(MDS_indices_DF$ZERO_BASED_INDICES) != MDS_indices_DF$SAMPLES_COUNT, na.rm = TRUE)
      if(MDS_bad > 0) {
        cat(warn(paste0("\n [!] ", MDS_bad, " rows of the MDS file have a number of indices different from SAMPLES_COUNT for the ",
                        site.ID, " station.\n")))
      }

      cat(prog.mes(paste0('\nIntegrating MDS gap-filling indices for the ', site.ID, ' station.\n')))

      # Generate a sequence of timestamps based on the provided years
      timestamps <- lapply(Data_Years, function(yr) {
        start_time <- as.POSIXct(sprintf("%d-01-01 00:00:00", yr), tz = "UTC")
        end_time   <- as.POSIXct(sprintf("%d-12-31 23:30:00", yr), tz = "UTC")
        seq(from = start_time, to = end_time, by = "30 mins")
      }) %>% do.call(c, .)

      n_obs <- length(timestamps)

      # Create dataframe with continuous indices
      DF_Idx <- data.frame(TIMESTAMP_START = format(timestamps, "%Y%m%d%H%M") %>% as.numeric(),
                           ZERO_BASED_INDEX = seq_len(n_obs) - 1)

      # MDS columns, computed once for the rows of the table #
      MDS_cols <- tibble(TIMESTAMP_START = Model_Input[['FFP_input_parameters']]$TIMESTAMP %>% format("%Y%m%d%H%M") %>% as.numeric()) %>%
        left_join(DF_Idx, by = "TIMESTAMP_START") %>%
        left_join(MDS_indices_DF, by = "ZERO_BASED_INDEX") %>%
        transmute(Idx = ZERO_BASED_INDEX, MDS_Idx = ZERO_BASED_INDICES)

      # The list and the saved table must have the same rows (same timestamps, same order) #
      stopifnot('MDS indices: the number of rows changed after the join (duplicated indices in the MDS file?)' =
                  nrow(MDS_cols) == nrow(Model_Input[['FFP_input_parameters']]),
                'MDS indices: the returned list and the saved table have different timestamps' =
                  identical(as.numeric(as.POSIXct(Model_Input[['FFP_input_parameters']]$TIMESTAMP)), as.numeric(as.POSIXct(Model_DF$TIMESTAMP))))

      # MDS columns already in the input table are replaced #
      if(any(c("Idx", "MDS_Idx") %in% colnames(Model_Input[['FFP_input_parameters']]))) {
        cat(warn(paste0("\n [!] MDS indices already in the input table are replaced by the ones of the MDS file for the ", site.ID, " station.\n")))
      }

      # Added to both the returned list and the saved table #
      Model_Input[['FFP_input_parameters']] <- bind_cols(Model_Input[['FFP_input_parameters']] %>% select(-any_of(c("Idx", "MDS_Idx"))), MDS_cols)
      Model_DF <- bind_cols(Model_DF %>% select(-any_of(c("Idx", "MDS_Idx"))), MDS_cols)

      # --- Save the MDS input file in a specific folder ---

      if(is.null(out.dir)) {
        Idx_dir <- file.path(getwd(), site.ID, 'Input', 'MDS')
      } else {
        Idx_dir <- file.path(out.dir, site.ID, 'Input', 'MDS')
      }

      # Create directory if it doesn't exist
      if(!dir.exists(Idx_dir)) { dir.create(Idx_dir, recursive = TRUE) }

      # Copy the files in the specific folder and rename
      copied_idx <- file.copy(from = MDS.file,
                              to = file.path(Idx_dir, basename(MDS.file)),
                              overwrite = TRUE)

      if(copied_idx) {
        cat(prog.mes(paste0("\nMDS indices input file successfully copied into the ", site.ID,  " MDS input folder\n")))
      } else {
        cat(warn(paste0("\n [!] MDS indices input file was not copied into the ", site.ID,  " MDS input folder\n")))
      }

      # Check the MDS file against the index grid built from the data years #
      # Last half-hour referred by the MDS file (gap-filled or used for gap-filling) #
      MDS_last <- max(c(MDS_indices_DF$ZERO_BASED_INDEX, unlist(MDS_indices_DF$ZERO_BASED_INDICES)), na.rm = TRUE)

      # MDS file beyond the grid: different years (e.g. MDS input starting earlier, or a data year missing) #
      if(MDS_last >= n_obs) {
        cat(warn(paste0("\n [!] The MDS file refers to ", MDS_last + 1, " half-hours, but the data years (",
                        paste(Data_Years, collapse = ", "), ") contain ", n_obs,
                        ": MDS indices are probably misaligned for the ", site.ID, " station.\n")))
      }

      # MDS file ending more than one year before the grid: MDS years missing #
      if(MDS_last < n_obs - 48 * 366) {
        cat(warn(paste0("\n [!] The MDS file ends more than one year before the end of the data years (",
                        paste(Data_Years, collapse = ", "), "): check the MDS years for the ", site.ID, " station.\n")))
      }

    }
  } ## MDS chunk end


  ## -- SAVING OPTIONS -- ##

  # merge.dataset only acts on the saved table #
  if(!is.null(merge.dataset) & !save.input.table) {
    cat(warn(paste0("\n[!] merge.dataset is applied only to the saved table: it is ignored because save.input.table = FALSE")))
  }

  # Save input tables in the FFP input folder
  if (all(save.input.table & is.data.frame(Model_Input[['FFP_input_parameters']]))) {

    # Merge an external dataset # Only to update lat, lon, UTC_offset and PID of the saved table # The returned list is not modified #
    if(!is.null(merge.dataset)) {

      stopifnot('merge.dataset must have a single row' = nrow(merge.dataset) == 1)

      # site_id is never replaced # If provided, it must match the site of the table (safety check against mixing sites) #
      if('site_id' %in% colnames(merge.dataset)) {
        if(!identical(as.character(merge.dataset$site_id), as.character(site.ID))) {
          stop(paste0("merge.dataset refers to the site ", merge.dataset$site_id,
                      ", but the table refers to the site ", site.ID, ": merge not applied"))
        }
      }

      Merge_columns <- c('lat', 'lon', 'UTC_offset', 'PID')

      # Any other column is ignored #
      Ignored_columns <- setdiff(colnames(merge.dataset), c('site_id', Merge_columns))
      if(length(Ignored_columns) > 0) {
        cat(warn(paste0("\n[!] ", paste(Ignored_columns, collapse = ", "),
                        " ignored: merge.dataset can only update lat, lon, UTC_offset and PID")))
      }

      merge.dataset <- merge.dataset %>% select(any_of(Merge_columns))

      if(ncol(merge.dataset) > 0) {
        cat(warn(paste0("\n[!] ", paste(colnames(merge.dataset), collapse = ", "),
                        " are replaced in the saved table by using an external dataframe")))

        # Same value on every row # Metadata always in the first columns
        Model_DF <- bind_cols(merge.dataset, Model_DF %>% select(-any_of(colnames(merge.dataset)))) %>%
          select(any_of(c('site_id', 'lat', 'lon', 'UTC_offset', 'PID')), everything())
      }
    }

    if(is.null(out.dir)) {
      Site_dir <- file.path(getwd(), site.ID)
      Input_dir <- file.path(getwd(), site.ID, 'Input')
    } else {
      Site_dir <- file.path(out.dir, site.ID)
      Input_dir <- file.path(out.dir, site.ID, 'Input')
    }

    if (!dir.exists(Site_dir)) {dir.create(Site_dir, recursive = TRUE)}
    if (!dir.exists(Input_dir)) {dir.create(Input_dir, recursive = TRUE)}

    # Convert the MDS list-column to character string before writing the CSV
    if("MDS_Idx" %in% colnames(Model_DF)) {
      Model_DF <- Model_DF %>%
        mutate(MDS_Idx = map_chr(MDS_Idx, function(x) {
          if (is.null(x) || all(is.na(x))) { return(NA_character_) } # If NULL or NA, return empty cell
          else { return(paste(x, collapse = ",")) } # Collapse the vectors into strings
        }))
    }

    # ## Write the ERA5 input files ##
    if(length(ERA5list_filtered) > 0) {

      ERA5_dir <- file.path(Input_dir, "ERA5")
      if (!dir.exists(ERA5_dir)) {dir.create(ERA5_dir, recursive = TRUE)}

      # Copy the input tables in the ERA5 folder
      ERA5list_filtered %>% map(function(x) file.copy(from = x, to = file.path(ERA5_dir, basename(x))))

      message(prog.mes('\nERA5 input files saved'))
    }

    # Final CSV Writing
    write_csv(Model_DF, file = file.path(Input_dir, paste0(
      site.ID, '_', min(as_date(Model_DF$TIMESTAMP)) %>% str_replace_all('-', ''), '_',
      max(as_date(Model_DF$TIMESTAMP)) %>% str_replace_all('-', ''), '_FFPinput.csv')))

    message(prog.mes('\nFFP input table saved\n'))
  }


  # -- RETURN INPUT TABLE -- #
  return(Model_Input)

} # Function ending
