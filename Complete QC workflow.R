library(tidyverse)

## Set Up ####

## Read params file

DataProcess_Params <- read_delim("~/DataProcess_Params.tsv", 
                                 delim = "\t", escape_double = FALSE, 
                                 col_names = FALSE, trim_ws = TRUE) 

## select only the names of the params

parameters <- DataProcess_Params %>% 
  select(X1) %>% 
  as.matrix() %>% 
  as.vector() 

## select the values for each parameter

config <- DataProcess_Params %>% 
  select(X2) %>% 
  as.matrix() %>% 
  as.list()

names(config) <- parameters

## Correctly parse select parameters

config["disk"] <- as.logical(config["disk"])

config["QC.Only"] <- as.logical(config["QC.Only"])

config["valid.value.threshold"] <- as.numeric(config["valid.value.threshold"])

## Individual Functions ####

ID.Free.QC <- function(mzML.file.dir, Annotation.dir, Anno.Out){
  
  ## Locate raw data files (must be .mzML), R cannot parse .raw/.wiff/,d etc
  
  config[["mzML.file.dir"]]
  setwd(config[["mzML.file.dir"]])
  
  Raw.Files.list <- dir(path = config[["mzML.file.dir"]], recursive = T, pattern = "\\.mzML$", full.names = T)
  
  Raw.file.names <- dir(path = config[["mzML.file.dir"]], recursive = T, pattern = "\\.mzML$", full.names = F)
  
  print("Reading mzML files")
  
  ## Load all mzMLs into R
  
  mzML.list <- lapply(1:length(Raw.Files.list), function(x){
    
    mzR::openMSfile(Raw.Files.list[[x]])
    
    
  })
  
  names(mzML.list) <- Raw.file.names
  
  ## Get Acquisition metadata ####
  
  print("Parsing mzML metadata")
  
  ## for each mzML file returns metadata of Instrument and method
  
  runInfo.list <- lapply(1:length(mzML.list), function(x){
    
    mzR::runInfo(object = mzML.list[[x]]) %>% 
      as.data.frame() %>% 
      mutate(file.name = paste(names(mzML.list[x])))
    
  })
  
  names(runInfo.list) <- names(mzML.list)
  
  ## Get Time stamps for each MS-run, to parse runorder
  
  Time.stamps <- lapply(1:length(runInfo.list), function(x){
    
    y <- runInfo.list[[x]][["startTimeStamp"]] %>% 
      unique() %>% 
      as.data.frame(row.names = names(runInfo.list[x])) 
    
  })
  
  Time.stamps.df <- bind_rows(Time.stamps) %>% 
    distinct() %>% 
    rename(time = 1) %>% 
    separate(col = "time", sep = "T|Z", into = c("day", "time")) %>% 
    separate("time", into = c("hours", "minutes", "seconds"), sep = ":") %>% 
    separate("day", into = c("year", "month", "day"), sep = "-")  %>% 
    mutate(stamp = ((as.numeric(year)*365*24*60*60) + 
                      (as.numeric(month)*30*24*60*60) + 
                      (as.numeric(day)*24*60*60) + 
                      (as.numeric(hours)*60*60) + 
                      (as.numeric(minutes)*60) + 
                      (as.numeric(seconds)))) %>% 
    arrange() %>% 
    mutate(RunOrder = 1:length(Raw.file.names)) %>% 
    rownames_to_column(var = "file.name") %>% 
    select(file.name, RunOrder)
  
  ## dataframe for any sample/biological annotations for dataset
  
  Sample.list <- read_delim(Annotation.dir, 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)
  
  Sample.Anno <- inner_join(Sample.list, Time.stamps.df) 
  
  write_delim(x = Sample.Anno, file = paste(Anno.Out, "/Annotation.tsv", sep = ""), delim = "\t", append = F)
  
  runInfo <- inner_join(bind_rows(runInfo.list, .id = "file.name"), Sample.Anno) %>% 
    select(-msLevels) %>% 
    distinct() ## Final Parsed metadata for each MS-run
  
  ## Get mzML file info ####
  
  ## returns information on each MS scan - MS1 and MS/MS scans (MS3 not implemented)
  ## returns scan information: Retention time (s), Injection Time, m/z etc - see ?mzR::header for all info returned
  
  
  print("Parsing Scan Headers")
  
  
  headers.fun <- function(x){
    
    header.list <- lapply(1:unique(runInfo.list[[x]][["scanCount"]]), function(y){
      
      print(paste("Parsing Headers for", names(mzML.list[x])))
      
      header <- mzR::header(object = mzML.list[[x]], scans = y) %>% 
        as.data.frame()
      
    })
    
    header.df <- bind_rows(header.list)
    
  }
  
  ## parallelises function (uses half available cores)
  
  headers <- parallel::mclapply(X = 1:length(mzML.list), 
                                FUN = headers.fun, 
                                mc.cores = (detectCores()/2), 
                                mc.cleanup = T)
  
  
  names(headers) <- names(mzML.list)
  
  ## Separate msLevels, precursor and fragment
  
  ##data for each MS scan
  
  prec.headers <- inner_join(bind_rows(headers, .id = "file.name"), Sample.Anno) %>% 
    filter(msLevel == 1)
  
  ##data for each MS/MS scan
  
  frag.headers <- inner_join(bind_rows(headers, .id = "file.name"), Sample.Anno) %>%
    filter(msLevel == 2)
  
  rm(headers)
  gc()
  
  ## Scan Counts ####
  
  print("Extracting Scan Counts")
  
  ## Returns the number of scans per file (MS and MS/MS)
  
  Scan.Counts <- inner_join(inner_join((prec.headers %>%
                                          group_by(RunOrder, label) %>% 
                                          summarise(`Precursor Scans` = n_distinct(seqNum)/1000)), 
                                       (frag.headers %>% 
                                          group_by(RunOrder, label) %>% 
                                          summarise(`MS/MS Scans` = n_distinct(seqNum)/1000))
  ), (runInfo %>% 
        dplyr::select(RunOrder, scanCount) %>% 
        mutate(scanCount = scanCount/1000) %>% 
        rename(`Total Scans` = scanCount))) %>% 
    pivot_longer(cols = contains("Scans"), names_to = "Scan.Level", values_to = "Scans") %>% 
    mutate(`Scan Type` = factor(Scan.Level, levels = c("Precursor Scans", "MS/MS Scans", "Total Scans"))) %>% 
    arrange("RunOrder") %>% 
    ggpubr::ggbarplot(x = "label",
                      y = "Scans",
                      ylab = "Number of scans (x 10^3)",
                      color = "Scan Type", 
                      label = T, 
                      lab.pos = "out",
                      position = position_dodge(0.9), 
                      title = "Number of Scans per Raw file",
                      legend = "right",
                      lab.nb.digits = 1) +
    ggpubr::rotate_x_text(45) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ## Sensitivity ####
  
  ## InjectionTime
  
  ## returns Injection Times, for all files, as well as individual files
  ## accross retention time
  
  print("Extracting Injection Times")
  
  overlayed.IT <- prec.headers %>%  
    ggpubr::ggline(x = "retentionTime", 
                   xlab = "Rention Time (s)", 
                   y = "injectionTime", 
                   ylab = "Injection Time (ms)",
                   plot_type = "l", 
                   size = .5, 
                   color = "label", 
                   numeric.x.axis = T, 
                   legend = "right", 
                   title = "Injection time of all Precursor Scans, per Retention Time") +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  prec.list <- split.data.frame(prec.headers, f = prec.headers$file.name) 
  
  single.IT <- lapply(1:length(prec.list), function(x){
    
    prec.list[[x]] %>%  
      ggpubr::ggline(x = "retentionTime", 
                     xlab = "Rention Time (s)", 
                     y = "injectionTime", 
                     ylab = "Injection Time (ms)",
                     plot_type = "l", 
                     size = .5, 
                     color = "label", 
                     numeric.x.axis = T, 
                     legend = "right",
                     title = paste("Injection time of Precursor Scans in", names(prec.list[x]), "- per Retention Time", sep = " ")) +
      theme(plot.title = element_text(hjust = 0.5))
    
  })
  
  names(single.IT) <- names(mzML.list)
  
  
  ## Total Ion Current ####
  
  print("Extracting TICs")
  
  ##plots the MS1 TIC for each file (overlayed)
  
  TIC <- prec.headers %>%  
    ggpubr::ggline(x = "retentionTime", 
                   xlab = "Retention Time (s)",
                   y = "totIonCurrent",
                   ylab = "Total Ion Current",
                   numeric.x.axis = T, 
                   color = "label",
                   legend = "right",
                   size = 0.5, 
                   plot_type = "l",
                   title = "Total Ion Current of all Precursor Scans, per Retention Time") +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ## MS and MS/MS summmaries ####
  
  print("Summarising Intensities")
  
  ## Distribution of MS1 IT
  
  ms1.IT.histo <- prec.headers %>% 
    ggpubr::gghistogram(x = "injectionTime",
                        xlab = "Injection time (ms)",
                        y = "..count..",
                        binwidth = 0.5, 
                        color = "label", 
                        add_density = T, 
                        legend = "right", 
                        facet.by = "condition", 
                        title = "Distribution of Injection times") +
    theme(plot.title = element_text(hjust = 0.5))
  
  ms1.IT.Box <- prec.headers %>%
    arrange(RunOrder) %>% 
    ggpubr::ggviolin(x = "label",
                     y = "injectionTime",
                     ylab = "Injection time (ms)",
                     color = "condition",
                     legend = "right", 
                     title = "Distribution of Injection times", 
                     add = "boxplot") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggpubr::rotate_x_text(45)
  
  ## Distribution of MS2 IT
  
  ms2.IT.histo <- frag.headers %>% 
    ggpubr::gghistogram(x = "injectionTime",
                        xlab = "Injection time (ms)",
                        y = "..count..",
                        binwidth = 0.5, 
                        color = "label", 
                        add_density = T, 
                        legend = "right", 
                        facet.by = "condition", 
                        title = "Distribution of Injection times") +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  
  ms2.IT.Box <- frag.headers %>%
    arrange(RunOrder) %>% 
    ggpubr::ggviolin(x = "label",
                     y = "injectionTime",
                     ylab = "Injection time (ms)",
                     color = "condition",
                     legend = "right", 
                     title = "Distribution of Injection times", 
                     add = "boxplot") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggpubr::rotate_x_text(45)
  
  ## Distribution of Frament Intensities
  
  Frag.TIC.histo <- frag.headers %>% 
    mutate(prec.int = log2(precursorIntensity)) %>% 
    mutate(RunOrder = as.character(RunOrder)) %>%
    ggpubr::gghistogram(x = "prec.int",
                        y = "..count..",
                        density = T,
                        legend = "right", 
                        add_density = T,
                        color = "label", 
                        binwidth = 0.1, 
                        add = "median", 
                        facet.by = "condition",
                        title = "Distribution of Precursor Ion Intensities Selected for Fragmentation") +
    theme(plot.title = element_text(hjust = 0.5))
  
  Frag.TIC.Box <- frag.headers %>%
    mutate(prec.int = log2(precursorIntensity)) %>%
    arrange(RunOrder) %>% 
    ggpubr::ggviolin(x = "label",
                     y = "prec.int",
                     ylab = "log2 Intensity of MS/MS precursors",
                     color = "condition",
                     legend = "right", 
                     add = "boxplot",
                     title = "Distribution of log2 Intensity of MS/MS precursors") +
    theme(plot.title = element_text(hjust = 0.5)) +
    rotate_x_text(45)
  
  ## Distribution of Precursor Intensities
  
  TIC.histo <- prec.headers %>% 
    mutate(prec.int = log2(totIonCurrent)) %>% 
    mutate(RunOrder = as.character(RunOrder)) %>%
    ggpubr::gghistogram(x = "prec.int",
                        y = "..count..",
                        add_density = T,
                        legend = "right",
                        color = "label", 
                        binwidth = 0.1,
                        facet.by = "condition",
                        title = "Distribution of all Precursor Ion Intensities") +
    theme(plot.title = element_text(hjust = 0.5))
  
  TIC.Box <- prec.headers %>%
    mutate(prec.int = log2(totIonCurrent)) %>%
    arrange(RunOrder) %>% 
    ggpubr::ggviolin(x = "label",
                     y = "prec.int",
                     ylab = "log2 Intensity of Precursor Ions ",
                     color = "condition",
                     legend = "right", 
                     title = "Distribution of log2 Intensity of all precursors", 
                     add = "boxplot") +
    theme(plot.title = element_text(hjust = 0.5)) +
    rotate_x_text(45)
  
  
  ## MS/MS Performance ####
  
  
  frag.list <- split.data.frame(frag.headers, f = frag.headers$file.name) 
  
  ## TIC of MS and MS/MS scans accross RT
  
  single.Full.TIC <- lapply(1:length(frag.list), function(x){
    
    Scans <- bind_rows(frag.list[[x]], prec.list[[x]])
    
    
    ggplot(data = Scans, mapping = aes(x = retentionTime, color = factor(msLevel))) +
      geom_line(data = (Scans %>% 
                          filter(msLevel == 1)), mapping = aes(y = totIonCurrent)) +
      geom_line(data = (Scans %>% 
                          filter(msLevel == 2)), mapping = aes(y = totIonCurrent)) +
      labs(color = "msLevel") +
      ggtitle(paste("Total Ion Current of MS1 and MS2 of - ", names(prec.list[x]))) +
      xlab("Retention Time (s)") +
      ylab("Total Ion Current") +
      theme(plot.title = element_text(hjust = 0.5))
    
  })
  
  names(single.Full.TIC) <- names(mzML.list)
  
  print("Calculating MS/MS metrics")
  
  ## TopN performance, how often the MS is saturated for MS/MS scans
  
  TopN <- frag.headers %>%
    group_by(RunOrder, label, precursorScanNum) %>% 
    summarise(TopN = n_distinct(precursorMZ),
              retentionTime = mean(retentionTime)) %>%
    mutate(roll.mean = rollmean(x = TopN, k = 100, fill = NA)) %>% 
    distinct() %>% 
    ggpubr::ggline(x = "retentionTime",
                   xlab = "Retention Time (s)",
                   y = "roll.mean", 
                   ylab = "Rolling Average of Number of MS/MS events per Precursor Scan",
                   legend = "right",
                   color = "label", 
                   plot_type = "l",  
                   numeric.x.axis = T, 
                   title = "Rolling average of MS/MS events per Retention Time") +
    theme(plot.title = element_text(hjust = 0.5))
  
  TopN.list <- split.data.frame(TopN[["data"]], f = TopN[["data"]]$label)
  
  single.TopN <- lapply(1:length(TopN.list), function(x){
    
    ggplot2::ggplot(data = TopN.list[[x]], aes(x = retentionTime)) +
      geom_line(aes(y = roll.mean), color = "blue") +
      geom_point(aes(y = TopN), color = "red") +
      ggtitle(paste("Number of Sequential MS/MS events per Retention Time - ", names(prec.list[x]))) +
      xlab("Retention Time (s)") +
      ylab("Sequential MS/MS events") +
      theme(plot.title = element_text(hjust = 0.5))
    
  })
  
  names(single.TopN) <- names(mzML.list)
  
  ## Angle - Score ####
  
  print("Parsing fragmentation scans")
  
  ## Calculate angle score according to: 
  ##"viQC: Visual and Intuitive Quality Control for Mass Spectrometry-Based Proteome Analysis"
  ## "https://link.springer.com/article/10.1134/S1061934819140119"
  ## https://github.com/lisavetasol/viQC
  
  frag.scans <- lapply(1:length(frag.list), function(x){
    
    frag.list[[x]] %>% 
      select(seqNum) %>% 
      distinct() %>% 
      as.matrix() %>% 
      as.vector()
    
  })
  
  angle.score.fun <- function(x){
    
    scans <- frag.scans[[x]]
    
    scan.mat <- lapply(1:length(scans), function(y){
      
      test <- mzR::peaks(object = mzML.list[[x]], scans = scans[[y]]) %>%
        as.data.frame() %>% 
        rename(`m/z` = V1,
               intensity = V2)
      
      peak.count <- length(test$`m/z`)
      
      avg.int <- test %>% 
        summarise(avg.int = log2(sum(intensity)/peak.count), 
                  peaks = peak.count,
                  seqNum = scans[[y]])
    })
    
    print("Plotting fragmentation performance")
    
    angle.score.mat <- bind_rows(scan.mat) %>% 
      ggpubr::ggscatter(x = "peaks",
                        xlab = "Number of Fragmentation Peaks",
                        y = "avg.int", 
                        ylab = "Average Intensity of Peaks",
                        size = 0.1, 
                        title = paste("Number of Peaks vs Average Intensity of Peaks for -", names(mzML.list[[x]]))) +
      theme(plot.title = element_text(hjust = 0.5))
    
  }
  
  
  angle.score.list <- parallel::mclapply(X = 1:length(frag.scans), 
                                         FUN = angle.score.fun, 
                                         mc.cores = (detectCores()/2), 
                                         mc.cleanup = T)
  
  names(angle.score.list) <- names(mzML.list)
  
  
  ## Dynamic Exclusion
  
  print("Plotting Dynamic Exlusion Performance")
  
  MS2.Exclusion <- lapply(1:length(frag.list), function(x){
    
    single <- split.data.frame(frag.list[[x]] , f = frag.list[[x]]$precursorCharge)
    
    image <- lapply(1:length(single), function(y){
      
      single[[y]] %>% 
        dplyr::select(retentionTime, precursorMZ, precursorCharge) %>% 
        distinct() %>% 
        ggpubr::ggscatter(x = "retentionTime",
                          xlab = "Retention time (s)",
                          y = "precursorMZ",
                          ylab = "Precursor m/z",
                          size = 0.5,
                          title = paste(names(frag.list[x]), "Charge State:", unique(single[[y]]$precursorCharge), sep = " ")) +
        theme(plot.title = element_text(hjust = 0.5))
      
    })
    
    
  })
  
  names(MS2.Exclusion) <- names(mzML.list)
  
  ## Outputs ####
  
  print("Collecting Outputs")
  
  Overlayed.Image.list <- list(Scan.Counts = Scan.Counts,
                               overlayed.IT = overlayed.IT,
                               TIC = TIC,
                               ms1.IT.histo = ms1.IT.histo,
                               ms1.IT.Box = ms1.IT.Box,
                               ms2.IT.Box = ms2.IT.Box,
                               TIC.histo = TIC.histo,
                               TIC.Box = TIC.Box,
                               Frag.TIC.histo = Frag.TIC.histo,
                               Frag.TIC.Box, Frag.TIC.Box,
                               TopN = TopN)
  
  print("Printing ID free metrics, per mzML to disk")
  
  ggsave(
    filename = "IdFreeQC-metrics.pdf", 
    plot = marrangeGrob(Overlayed.Image.list, nrow=1, ncol=1), 
    width = 16, height = 9
  )
  
  Single.Image.List <- lapply(1:length(Raw.file.names), function(x){
    
    name <- Raw.file.names[x]
    
    Image.list <- append(list(TIC = single.Full.TIC[[name]],
                              IT = single.IT[[name]],
                              TopN = single.TopN[[name]],
                              Angle.Score = angle.score.list[[name]]),MS2.Exclusion[[name]])   
    
    print("Printing summarised ID free metricsto disk")
    
    ggsave(
      filename = paste(name, "IdFreeQC-metrics.pdf"), 
      plot = marrangeGrob(Image.list, nrow=1, ncol=1), 
      width = 16, height = 9
    )
    
  })
  
  Outputs <- list(annotation = Sample.Anno)
  
  
}


Peptide.level.QCs <- function(Run.Specific.files, disk){
  
  disk = disk ## To save outputs to disk or not
  
  ## Retention Time QCs ####
  
  ## Run.Specific.files - output of msFragger data import
  
  ## Collect RT information and calculate peak-widths
  
  RT <- Run.Specific.files[["quant"]]  %>% 
    select(Run, label, condition, RunOrder, peptide, modified_peptide, charge, contains("retention")) %>%
    select(-retention_time) %>% 
    filter(!is.na(apex_retention_time)) %>% 
    mutate(peak_width.s = (retention_time_end - retention_time_begin),
           retention_time.min = apex_retention_time/60) %>%
    arrange(RunOrder) %>% 
    distinct()
  
  ## Peakwidth distribution
  
  PeakWidth_BoxPlot <- RT %>% 
    filter(peak_width.s <= 80) %>% 
    ggpubr::ggboxplot(x = "label",
                      color = "condition",
                      y = "peak_width.s", 
                      xlab = "Run",
                      ylab = "Peakwidth in Seconds", 
                      title = "Distribution of Peakwidths per Sample", 
                      numeric.x.axis = TRUE,
                      legend = "right") +
    rotate_x_text(45) +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  ## apex retention time variability
  
  RT_Total_CVs <- RT %>% 
    group_by(modified_peptide,  charge) %>% 
    summarise(meanRT = mean(apex_retention_time),
              sdRT = sd(apex_retention_time),
              medianRT = median(apex_retention_time)) %>% 
    dplyr::mutate(cvRT = sdRT/meanRT*100) %>% 
    ggpubr::gghistogram(x = "cvRT", 
                        y = "..density..", 
                        binwidth = 1, 
                        add_density = TRUE, 
                        title = "Distribution of apex retention time CVs",
                        add = "median", 
                        xlab = "Apex Retention Time CV", 
                        ylab = "Density") + 
    xlim(0, 20) +
    geom_vline(xintercept = 5, color = "red") +
    theme(plot.title = element_text(hjust = 0.5))
  
  ## apex retention time variability - per condition (from annotation file)
  
  RT_condition_CVs <- RT %>% 
    group_by(modified_peptide,  charge, condition) %>% 
    summarise(meanRT = mean(apex_retention_time),
              sdRT = sd(apex_retention_time),
              medianRT = median(apex_retention_time)) %>% 
    dplyr::mutate(cvRT = sdRT/meanRT*100) %>% 
    ggpubr::gghistogram(x = "cvRT", 
                        y = "..density..", 
                        color = "condition",
                        binwidth = 1, 
                        add_density = TRUE, 
                        title = "Distribution of apex retention time CVs per condition",
                        add = "median", 
                        xlab = "Apex Retention Time CV", 
                        ylab = "Density", 
                        legend = "right") + 
    xlim(0, 20) +
    geom_vline(xintercept = 5, color = "red") +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  ## Distribution of PSMs accross RT
  
  ID_Density <- RT %>% 
    ggpubr::ggdensity(x = "apex_retention_time",
                      y = "..density..",
                      ylab = "Density of PSMs",
                      color = "Run",
                      add_density = TRUE,
                      title = "PSM density per Retention Time", 
                      xlab = "Apex Retention Time in Seconds", 
                      legend = "right") +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  Total.Samples <- length(RT$Run %>% unique()) 
  
  set.seed(541)
  
  ## RT variability accross selected peptides (Identified in all samples)
  
  Constant_Peptides_RT <- inner_join((RT %>% 
                                        group_by(modified_peptide, label) %>% 
                                        summarise(meanRT = mean(apex_retention_time)) %>% 
                                        summarise(Valid.Values = sum(!is.na(meanRT))) %>% 
                                        filter(Valid.Values == Total.Samples) %>% 
                                        dplyr::select(modified_peptide) %>% 
                                        distinct() %>% 
                                        sample_n(60)), RT) %>% 
    group_by(label, modified_peptide, RunOrder) %>% 
    summarise(RT = median(apex_retention_time)) %>%
    arrange(RunOrder) %>%
    ggpubr::ggline(x = "label", 
                   y = "RT",  
                   color = "modified_peptide", 
                   ylab = "Retention time in seconds", 
                   xlab = "Acquisition Order", 
                   legend = "none",  
                   title = "Retention TIme Stability") +
    theme(plot.title = element_text(hjust = 0.5)) +
    rotate_x_text(45)
  
  
  ## Mass accuracy ####
  
  ## collect accuracy metrics
  
  accuracy <- Run.Specific.files[["quant"]] %>% 
    select(Run, label, condition, RunOrder, modified_peptide, contains("ppm"), hyperscore)
  
  ## Distribution of mass accuracy (in ppm) - per file
  
  mass.accuracy <- accuracy %>%
    pivot_longer(cols = contains("cal"), names_to = "ppm_diff", values_to = "ppm") %>%
    group_by(Run, label, RunOrder) %>%
    arrange(RunOrder) %>% 
    ggpubr::ggboxplot(x = "label", 
                      color = "condition",
                      xlab = "ppm difference",
                      y = "ppm", 
                      facet.by = "ppm_diff", 
                      legend = "right",
                      nrow = 2, 
                      title = "Distribution of the differences between observed m/z and calculated m/z in ppm") +
    ylim(-20, 20) +
    rotate_x_text(45) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ## peptide characteristics ####
  
  ## collect identified peptide characteristics
  
  pep.char <- Run.Specific.files[["psm"]]  
  
  ## Distribution of peptide length per sample
  
  pep.len <- pep.char %>% 
    select(Run, label, condition, RunOrder, Peptide, `Peptide Length`) %>% 
    distinct() %>% 
    ggpubr::ggboxplot(x = "label",
                      y = "Peptide Length",
                      color = "condition", 
                      legend = "right", 
                      title = "Distribution of Peptide Lengths") +
    rotate_x_text(45) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ## Distribution of peptide charge state per sample
  
  pep.charge <- pep.char %>% 
    select(Run, label, condition, RunOrder, Peptide, Charge) %>% 
    distinct() %>% 
    group_by(label) %>% 
    summarise(peptide.count = n_distinct(Peptide),
              condition = condition,
              RunOrder = RunOrder,
              Peptide = Peptide,
              Charge = Charge) %>% 
    ungroup() %>% 
    group_by(label, Charge) %>% 
    summarise(Charge.count = n_distinct(Peptide),
              condition = condition,
              RunOrder = RunOrder,
              peptide.count = peptide.count, 
              Charge.freq = Charge.count/peptide.count*100) %>% 
    distinct() %>% 
    arrange(RunOrder) %>% 
    ggpubr::ggbarplot(x = "label",
                      color = "condition",
                      y = "Charge.freq", 
                      ylab = "Frequency",
                      facet.by = "Charge", 
                      nrow = 5, 
                      legend = "right", 
                      title = "Frequency of peptides by charge") +
    rotate_x_text(45) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ## Number of missed cleavages per sample
  
  pep.term <- pep.char %>% 
    select(Run, label, condition, RunOrder, Peptide, `Number of Missed Cleavages`) %>% 
    distinct() %>% 
    group_by(label) %>% 
    summarise(peptide.count = n_distinct(Peptide),
              Peptide = Peptide,
              condition = condition,
              RunOrder = RunOrder,
              `Number of Missed Cleavages` = `Number of Missed Cleavages`) %>% 
    ungroup() %>% 
    group_by(label, `Number of Missed Cleavages`) %>% 
    summarise(MC.count = n_distinct(Peptide), 
              peptide.count = peptide.count,
              condition = condition,
              RunOrder = RunOrder,
              MC.freq = MC.count/peptide.count*100) %>%
    rename(Num.MC = `Number of Missed Cleavages`) %>% 
    distinct() %>% 
    arrange(RunOrder) %>% 
    ggpubr::ggbarplot(x = "label", 
                      y = "MC.freq",
                      ylab = "Frequency of Missed Cleavages",
                      facet.by = "Num.MC", 
                      color = "condition",
                      legend = "right",
                      nrow = 3, title = "Frequency of Missed Cleavages per Sample") +
    rotate_x_text(45) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ## counts of identified peptides per sample ####
  
  peptides <- Run.Specific.files[["ion"]]
  
  peptide.counts <- peptides  %>% 
    group_by(label) %>% 
    summarise(stripped.seq = n_distinct(`Peptide Sequence`), 
              `Peptide Sequence` = `Peptide Sequence`,
              `Modified Sequence` = `Modified Sequence`,
              `Assigned Modifications` = `Assigned Modifications`,
              Charge = Charge, 
              condition = condition,
              RunOrder  = RunOrder) %>% 
    ungroup() %>% 
    group_by(label) %>%
    summarise(stripped.seq = stripped.seq, 
              condition = condition,
              RunOrder  = RunOrder,
              dist.seq = n_distinct(`Modified Sequence`), 
              Mod.Peps = sum(!is.na(`Assigned Modifications`)),
              Cysteine = str_detect(`Assigned Modifications`, "57.0215"),
              Methionine = str_detect(`Assigned Modifications`, "15.9949"),
              `Modified Sequence` = `Modified Sequence`) %>% 
    pivot_longer(cols = c("Methionine", "Cysteine"), names_to = "Residue", values_to = "Modified") %>% 
    filter(Modified == TRUE) %>% 
    ungroup() %>% 
    group_by(label, Residue) %>% 
    summarise(stripped.seq = stripped.seq,
              dist.seq = dist.seq,
              Mod.Peps = Mod.Peps,
              Mod.Res = sum(Modified),
              condition = condition,
              RunOrder  = RunOrder) %>%
    arrange(RunOrder) %>% 
    distinct() %>% 
    pivot_wider(id_cols = c("label", "stripped.seq", "dist.seq", "Mod.Peps", "condition", "RunOrder"), names_from = "Residue", values_from = "Mod.Res")
  
  psm.count <- pep.char %>% 
    select(label, Spectrum) %>% 
    group_by(label) %>% 
    summarise(PSMs = n_distinct(Spectrum))
  
  Scan.counts <- inner_join(psm.count, peptide.counts) %>% 
    pivot_longer(cols = c("PSMs", "stripped.seq", "dist.seq", "Mod.Peps", "Cysteine", "Methionine"), names_to = "Property", values_to = "count") %>% 
    mutate(Property = factor(x = Property, levels = c("PSMs", "dist.seq", "stripped.seq", "Mod.Peps", "Cysteine", "Methionine"))) %>% 
    ggpubr::ggbarplot(x = "label", 
                      y = "count", 
                      color = "condition",
                      legend = "right",
                      facet.by = "Property", 
                      nrow = 3, 
                      scales = "free", 
                      panel.labs = list(Property = c("PSMs", "Distinct Sequences", "Stripped Sequences", "Modified Peptides", "Cysteine Carbidomidomethylation", "Methionine Oxidation")),
                      title = "Peptide and PSM counts") +
    rotate_x_text(45) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ## Outputs ####
  
  if(disk == F) {
    
    annotation <- Run.Specific.files[["ion"]] %>% 
      select(Run, file.name, condition, label, Tech.Rep, replicate, RunOrder) %>% 
      distinct()
    
    Output.list <- list(Scan.counts = Scan.counts,
                        mass.accuracy = mass.accuracy,
                        pep.len = pep.len,
                        pep.charge = pep.charge,
                        pep.term = pep.term,
                        RT_Total_CVs = RT_Total_CVs,
                        RT_condition_CVs = RT_condition_CVs,
                        PeakWidth_BoxPlot = PeakWidth_BoxPlot,
                        Constant_Peptides_RT = Constant_Peptides_RT,
                        ID_Density = ID_Density, 
                        annotation = annotation
    )
    
  } else {
    
    Image.list <- list(Scan.counts = Scan.counts,
                       mass.accuracy = mass.accuracy,
                       pep.len = pep.len,
                       pep.charge = pep.charge,
                       pep.term = pep.term,
                       RT_Total_CVs = RT_Total_CVs,
                       RT_condition_CVs = RT_condition_CVs,
                       PeakWidth_BoxPlot = PeakWidth_BoxPlot,
                       Constant_Peptides_RT = Constant_Peptides_RT,
                       ID_Density = ID_Density
    )
    
    ggsave(
      filename = "Peptide_Level_QCs.pdf", 
      plot = marrangeGrob(Image.list, nrow=1, ncol=1), 
      width = 15, height = 9
    )
    
    annotation <- Run.Specific.files[["ion"]] %>% 
      select(Run, file.name, condition, label, Tech.Rep, replicate, RunOrder) %>% 
      distinct()
    
    Output <- list(annotation = annotation)
    
  }
  
} ##

Protein.level.QCs <- function(Search.Output.dir, valid.value.threshold, Comparisons, disk){
  
  ## Import Protein level data ####
  
  setwd(Search.Output.dir)
  
  annotation <- Peptide.QCs[["annotation"]] ## output from peptide level QCs
  
  combined.proteins <- read_delim("combined_protein.tsv", 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE) %>% 
    filter(str_detect(Protein, "contam|rev") == F) %>% 
    select(Protein, `Protein ID`, `Entry Name`, Gene, `Protein Length`, Coverage, Description, contains("Total Intensity")) %>% 
    pivot_longer(cols = contains("Total Intensity"), names_to = "File.Name", values_to = "Intensity") %>% 
    separate(File.Name, into = c("Run", "Quant.type"), sep = " ", extra = "merge") %>% 
    mutate(Quant.type = case_when(str_detect(Quant.type, "MaxLFQ") ~ "MaxLFQ.Intensity",
                                  str_detect(Quant.type, "MaxLFQ") == F ~ "Freequant.Intensity"))
  
  Annotated.data <- inner_join(combined.proteins, annotation)
  
  ## Quantified Protein Counts ####
  
  print("Calculating Protein IDs per Run")
  
  Protein.count.Run <- Annotated.data %>%
    filter(Intensity > 0) %>% 
    group_by(label, Quant.type) %>% 
    summarise(PG.Count = n_distinct(Protein)) %>% 
    ggpubr::ggbarplot(x = "label", 
                      y = "PG.Count", 
                      fill = "Quant.type",
                      legend = "right",
                      position = position_dodge(0.9), 
                      label = TRUE, lab.size = 3) +
    rotate_x_text(45)  
  
  print("Calculating Protein IDs per condition")
  
  Protein.count.condition <- Annotated.data %>%
    filter(Intensity > 0) %>% 
    group_by(condition, Quant.type) %>% 
    summarise(PG.Count = n_distinct(Protein)) %>% 
    ggpubr::ggbarplot(x = "condition", 
                      y = "PG.Count", 
                      fill = "Quant.type",
                      legend = "right",
                      position = position_dodge(0.9), 
                      label = TRUE, lab.size = 3) +
    rotate_x_text(45)
  
  ## Missingness metrics ####
  
  print("Calculating Protein ID overlap per Run")
  
  ## Rates of Identification
  
  Missing.Proportions.Total <- Annotated.data %>%
    dplyr::select(Run, Intensity, Protein, Quant.type) %>% 
    distinct() %>% 
    mutate(PG.Normalised = case_when(Intensity == 0 ~ 0,
                                     Intensity != 0 ~ 1)) %>% 
    dplyr::group_by(Protein, Quant.type) %>% 
    summarise(Valid.Values = sum(PG.Normalised)) %>%   
    distinct() %>% 
    group_by(Quant.type, Valid.Values) %>% 
    summarise(Valid.Counts = n_distinct(Protein)) %>%
    mutate(Valid.Proportion = (Valid.Counts/(Annotated.data$Protein %>% n_distinct()))*100,
           Valid.Values.axis = as.numeric(as.character(Valid.Values))) %>% 
    distinct()  %>%
    ggpubr::ggbarplot(x = "Valid.Values", 
                      xlab = "Number of Samples",
                      y = "Valid.Proportion", 
                      ylab = "Proportion of identified proteins",
                      fill = "Valid.Values.axis",
                      facet.by = "Quant.type",
                      label = T, 
                      lab.nb.digits = 1, 
                      legend = "right", numeric.x.axis = TRUE,
                      title = "Identification frequency as a proportion of all identifications") +
    theme(plot.title = element_text(hjust = 0.5))
  
  print("Calculating Protein ID overlap per condition")
  
  ## How many PGs pass valid value criteria (can be used for statistics later)
  
  Valid.Values.condition <- inner_join(Annotated.data %>%
                                         dplyr::select(Run, condition, Intensity, Protein, Quant.type) %>% 
                                         distinct() %>% 
                                         mutate(PG.Normalised = case_when(Intensity == 0 ~ 0,
                                                                          Intensity != 0 ~ 1)) %>% 
                                         dplyr::group_by(Protein, condition, Quant.type) %>% 
                                         summarise(Valid.Values = sum(PG.Normalised)) %>%   
                                         distinct() %>% 
                                         group_by(condition, Quant.type, Valid.Values) %>% 
                                         summarise(Valid.Counts = n_distinct(Protein)), Protein.count.condition[["data"]]) %>%
    mutate(Valid.Proportion = Valid.Counts/PG.Count,
           valid.value.axis = as.numeric(as.character(Valid.Values))) %>% 
    ggpubr::ggbarplot(x = "Valid.Values", 
                      xlab = "Number of Identifications",
                      y = "Valid.Proportion",
                      ylab = "Proportion of total identifications",
                      fill = "valid.value.axis",
                      facet.by = c("Quant.type", "condition" ), 
                      label = T, 
                      lab.nb.digits = 2, 
                      title = "Frequency of Identifications per Condition",
                      numeric.x.axis = TRUE, legend = "right") +
    theme(plot.title = element_text(hjust = 0.5))
  
  ## Filter data to only include potentially comparable proteins
  
  Valid.value.filter <- Annotated.data %>%
    dplyr::select(Run, condition, Intensity, Protein, Quant.type) %>% 
    distinct() %>% 
    mutate(PG.Normalised = case_when(Intensity == 0 ~ 0,
                                     Intensity != 0 ~ 1)) %>% 
    dplyr::group_by(Protein, condition, Quant.type) %>% 
    summarise(Valid.Values = sum(PG.Normalised)) %>%   
    distinct()
  
  print("Calculating overlaps per condition of potentially comparable protein IDs")
  
  valid.value.threshold <- valid.value.threshold
  
  ## Calculate overlaps of valid values accross conditions
  
  Valid.Values.condition.overlap <- Annotated.data %>%
    dplyr::select(Run, condition, Intensity, Protein, Quant.type) %>% 
    distinct() %>% 
    mutate(PG.Normalised = case_when(Intensity == 0 ~ 0,
                                     Intensity != 0 ~ 1)) %>% 
    dplyr::group_by(Protein, condition, Quant.type) %>% 
    summarise(Valid.Values = sum(PG.Normalised)) %>% 
    filter(Valid.Values >= valid.value.threshold) %>% 
    mutate(Valid = TRUE) %>%
    pivot_wider(id_cols = c("Protein",  "Quant.type"), names_from = "condition", values_from = "Valid", values_fill = FALSE)
  
  overlap.list <- split.data.frame(Valid.Values.condition.overlap, Valid.Values.condition.overlap$Quant.type)
  
  Overlaps <- lapply(1:length(overlap.list), function(x){
    
    Total.Overlap <- ComplexUpset::upset(data = (overlap.list[[x]] %>% 
                                                   select(-Quant.type) %>% 
                                                   column_to_rownames("Protein")), intersect = colnames(overlap.list[[x]] %>%                                             
                                                                                                          select(-Quant.type) %>% 
                                                                                                          column_to_rownames("Protein")), n_intersections = 20)
    
  })
  
  names(Overlaps) <- names(overlap.list)
  
  Comparisons <- read_csv("/media/naadir/Outputs/Workflow.Test/SearchOutputs/Comparisons.txt", 
                          col_names = FALSE) %>% 
    as.matrix() %>% 
    as.vector() ## read comparisons of interest
  
  ## Determine which proteins are comparible
  
  Valid.Comparisons <- lapply(1:length(Comparisons), function(x){
    
    print(paste("Parsing Valid Comparisons for", Comparisons[x]))
    
    contrast <- Comparisons[[x]] %>% 
      as.data.frame() %>% 
      separate_rows(".", sep = "_vs_") %>% 
      rename(condition = 1)
    
    Valid.values.long <- Valid.Values.condition.overlap %>% 
      pivot_longer(cols = 3:(length(colnames(Valid.Values.condition.overlap))), names_to = "condition", values_to = "Valid")
    
    
    Valid.Comparisons <- inner_join(contrast, Valid.values.long) %>% 
      mutate(Con = paste(condition, Quant.type)) %>% 
      pivot_wider(id_cols = c("Protein"), names_from = "Con", values_from = "Valid") %>%
      select(Protein, contains("MaxLFQ"), contains("Freequant")) %>% 
      rename(LFQ_1 = 2,
             LFQ_2 = 3,
             Int_1 = 4,
             Int_2 = 5) %>% 
      mutate(Valid.MaxLFQ.comp = case_when(LFQ_1 == T & LFQ_2 == T ~ TRUE),
             Valid.Freequant.comp = case_when(Int_1 == T & Int_2 == T ~ TRUE),
             Comparison = paste(Comparisons[[x]])) %>% 
      select(Protein, contains("comp")) %>% 
      replace(is.na(.), FALSE)
    
    
  })
  
  print("Calculating overlaps per comparison of comparable protein IDs")
  
  Valid.Comparisons.df <- bind_rows(Valid.Comparisons)
  
  Int.Upset.matrix <- Valid.Comparisons.df %>% 
    select(Protein, Valid.Freequant.comp, Comparison) %>% 
    pivot_wider(id_cols = "Protein", names_from = "Comparison", values_from = "Valid.Freequant.comp") %>% 
    column_to_rownames("Protein")
  
  Int.Upset <- upset(Int.Upset.matrix, intersect = colnames(Int.Upset.matrix), n_intersections = 10) +
    ggtitle("Top 10 Intersections between Freequant Intensity comparisons")
  
  LFQ.Upset.matrix <- Valid.Comparisons.df %>% 
    select(Protein, Valid.MaxLFQ.comp, Comparison) %>% 
    pivot_wider(id_cols = "Protein", names_from = "Comparison", values_from = "Valid.MaxLFQ.comp") %>% 
    column_to_rownames("Protein")
  
  
  LFQ.Upset <- upset(LFQ.Upset.matrix, intersect = colnames(LFQ.Upset.matrix), n_intersections = 10) +
    ggtitle("Top 10 Intersections between MaxLFQ comparisons")
  
  ## Missingness ####
  
  print("Calculating missig data metrics")
  
  Missing.value.data <- Annotated.data %>% 
    as.data.frame() %>% 
    mutate(missing = Intensity == 0,
           Intensity = case_when((missing == FALSE) ~ Intensity))
  
  Missing.Values.stats <- function(data, title, histo){
    
    Missing.data <- data %>% 
      mutate(Intensity = case_when(Intensity != 0 ~ Intensity)) %>% 
      group_by(Protein, condition, Quant.type) %>%  
      summarise(Averagelog2Intensity = mean(log2(Intensity), na.rm = TRUE),
                MissingValue = any(is.na(Intensity)))
    
    Missing.data.list <- split.data.frame(x = Missing.data, f = Missing.data$condition)
    
    
    fig.list <- lapply(1:length(Missing.data.list), function(z){
      
      a <- Missing.data.list[[z]] %>% 
        ggpubr::gghistogram(x = "Averagelog2Intensity",
                            xlab = "Average log2 Intensity",
                            y = histo, 
                            color = "MissingValue", 
                            facet.by = c("Quant.type"), 
                            nrow = 1, 
                            add_density = TRUE, 
                            title = paste(names(Missing.data.list[z]), "Intensity Distribution of Missing data", sep = ": "), 
                            add = "median", 
                            binwidth = 0.5, legend = "right") +
        theme(plot.title = element_text(hjust = 0.5))
      
      y <-  Missing.data.list[[z]] %>% 
        group_by(`MissingValue`) %>% 
        arrange(`Averagelog2Intensity`) %>% 
        mutate(num = 1, cs = cumsum(num), cs_frac = cs/n()) 
      
      b <- y %>% 
        ggpubr::ggline(x = "Averagelog2Intensity", 
                       xlab = "Average log2 Intensity",
                       y = "cs_frac", 
                       color = "MissingValue",
                       facet.by = c("Quant.type"), 
                       nrow = 1, 
                       numeric.x.axis = TRUE, 
                       plot_type = "l", 
                       ylab = " Cumulative Fraction",
                       title = paste(names(Missing.data.list[z]), "Cumulative Fraction of Missing data per Intensity", sep = ": "),
                       legend = "right") +
        theme(plot.title = element_text(hjust = 0.5))
      
      figure <- ggarrange(a, b , nrow = 2)
      print(figure)
    })
    
    names(fig.list) <- names(Missing.data.list)
    
    out <- fig.list
  }  
  
  Missing.data.fig <- Missing.Values.stats(data = Annotated.data, title = "Missing data distributions", histo = "..count..")  
  
  ##Quant QCs ####
  
  ## Distribution of intensities of identified proteins
  
  Box.Plots <- Annotated.data %>% 
    filter(Intensity > 0) %>% 
    mutate(Log2.Intensity = log2(Intensity)) %>%
    ggpubr::ggboxplot(x = "label", 
                      y = "Log2.Intensity", 
                      fill = "condition", 
                      facet.by = "Quant.type", 
                      nrow = 2, 
                      legend = "right") +
    rotate_x_text(45)
  
  Histograms <- Annotated.data %>% 
    filter(Intensity > 0) %>%  
    mutate(Log2.Intensity = log2(Intensity)) %>%
    ggpubr::ggdensity(x = "Log2.Intensity",
                      xlab = "Log2 Intensity",
                      y = "..density..",
                      add_density = TRUE,
                      color = "condition", 
                      facet.by = "Quant.type", 
                      nrow = 2, 
                      legend = "right", add = "median")
  
  ## Variability of Identified proteins
  
  CVs <- Annotated.data %>%  
    filter(Intensity > 0) %>% 
    mutate(Log2.Intensity = log2(Intensity)) %>% 
    group_by(Protein, condition, Quant.type) %>% 
    summarise(mean = mean(Log2.Intensity, na.rm = TRUE),
              sd = sd(Log2.Intensity, na.rm = TRUE),
              CV  = sd/mean*100) %>% 
    distinct() %>% 
    ggpubr::ggdensity(x = "CV",
                      xlab = "Coeffcient of Variation",
                      y = "..count..",
                      add_density = TRUE,
                      color = "condition", 
                      facet.by = "Quant.type", 
                      nrow = 2, 
                      legend = "right", 
                      add = "median")
  
  
  if(disk == F) {
    
    Image.list <- append(list(Protein.count.Run = Protein.count.Run,
                              Protein.count.condition = Protein.count.condition,
                              Missing.Proportions.Total = Missing.Proportions.Total,
                              Valid.Values.condition = Valid.Values.condition,
                              Int.Upset = Int.Upset,
                              LFQ.Upset = LFQ.Upset,
                              Box.Plots = Box.Plots,
                              Histograms = Histograms,
                              CVs = CVs,
                              Valid.Comparisons.df = Valid.Comparisons.df,
                              Valid.Values.filer = Valid.value.filter,
                              Annotated.data = Annotated.data,
                              
                              
    ), Missing.data.fig)
    
  } else {
    
    Image.list <- append(list(Protein.count.Run = Protein.count.Run,
                              Protein.count.condition = Protein.count.condition,
                              Missing.Proportions.Total = Missing.Proportions.Total,
                              Valid.Values.condition = Valid.Values.condition,
                              Int.Upset = Int.Upset,
                              LFQ.Upset = LFQ.Upset,
                              Box.Plots = Box.Plots,
                              Histograms = Histograms,
                              CVs = CVs
    ), Missing.data.fig)
    
    dir.create("RawData")
    dir.create("Filters")
    
    write_delim(x = Valid.Comparisons.df, file = "Filters/ValidComparisons.tsv", delim = "\t", append = F)
    
    write_delim(x = Valid.value.filter, file = "Filters/ValidValues.tsv", delim = "\t", append = F)
    
    write_delim(x = Annotated.data, file = "RawData/Annotated.data.tsv", delim = "\t", append = F)
    
    ggsave(
      filename = "Protein_Level_QCs.pdf", 
      plot = marrangeGrob(Image.list, nrow=1, ncol=1), 
      width = 15, height = 9
    )
    
    Output.dfs <- list(Valid.Comparisons.df = Valid.Comparisons.df,
                       Valid.Values.filer = Valid.value.filter,
                       Annotated.data = Annotated.data)
    
  }
  
  
} ## End of Protein QCs

## Process data for differential expression analysis
## filter to only proteins comparable in comparisons of interest
## perform imputation - Mark proteins that have been imputed
## analyse the effects of imputation accross the data set

Process.data.DE <- function(Valid.value.threshold, Missing.value.threshold, disk){
  
  ## Make Summarized experiment object ####
  
  Annotated.data <- Protein.QCs[["Annotated.data"]]
  
  Valid.value.filter <- Protein.QCs[["Valid.Values.filer"]] %>%  
    filter(Valid.Values >= Valid.value.threshold)
  
  #Protein.QCs[["Valid.Values.filer"]] 
  
  Valid.Comparisons.filter <- inner_join((Protein.QCs[["Valid.Comparisons.df"]]) %>% 
                                           separate_rows(Comparison, sep = "_vs_") %>% 
                                           distinct() %>% 
                                           rename(condition = Comparison,
                                                  Freequant.Intensity = Valid.Freequant.comp,
                                                  MaxLFQ.Intensity = Valid.MaxLFQ.comp) %>% 
                                           pivot_longer(cols = c("Freequant.Intensity", "MaxLFQ.Intensity"), names_to = "Quant.type", values_to = "Valid.Comparison") %>% 
                                           filter(Valid.Comparison == TRUE), (Annotated.data %>% 
                                                                                select(Protein, `Protein ID`) %>% 
                                                                                distinct())) %>% 
    rename(Protein.ID = `Protein ID`)
  
  #Protein.QCs[["Valid.Comparisons.df"]]
  
  True.Comparison.filter <- inner_join((Protein.QCs[["Valid.Comparisons.df"]]), (Annotated.data %>% 
                                                                                   select(Protein, `Protein ID`) %>% 
                                                                                   rename(Protein.ID = `Protein ID`) %>% 
                                                                                   distinct())) %>% 
    rename(Freequant.Intensity = Valid.Freequant.comp,
           MaxLFQ.Intensity = Valid.MaxLFQ.comp) %>% 
    pivot_longer(cols = c("Freequant.Intensity", "MaxLFQ.Intensity"), names_to = "Quant.type", values_to = "Valid.Comparison") %>% 
    filter(Valid.Comparison == TRUE)
  
  Filtered.Annotated.data <- inner_join(Annotated.data, Valid.value.filter)
  
  Annotation <- Annotated.data %>% 
    select(label, replicate, condition) %>% 
    distinct()
  
  
  Intensity.DEP.df <- Filtered.Annotated.data %>% 
    filter(str_detect(Quant.type, "Free") & Intensity > 0) %>% 
    select(`Protein ID`, Gene, label, Intensity) %>% 
    pivot_wider(id_cols = c("Protein ID", "Gene"), names_from = "label", values_from = "Intensity")
  
  MaxLFQ.DEP.df <- Filtered.Annotated.data %>% 
    filter(str_detect(Quant.type, "Max") & Intensity > 0) %>% 
    select(`Protein ID`, Gene, label, Intensity) %>% 
    pivot_wider(id_cols = c("Protein ID", "Gene"), names_from = "label", values_from = "Intensity")
  
  
  DEP.Input.Int <- DEP::make_unique(Intensity.DEP.df, names = "Protein ID", ids = "Gene")
  
  DEP.Input.LFQ <- DEP::make_unique(MaxLFQ.DEP.df, names = "Protein ID", ids = "Gene")
  
  conditions <- Annotated.data$condition %>% unique()
  
  Sample.names <- lapply(1:length(conditions), function(x){
    
    grep(conditions[[x]], colnames(DEP.Input.Int))
    
  })
  
  Sample.names.vec <- unlist(Sample.names)
  
  ExpDesign <- inner_join(DEP.Input.Int %>% 
                            pivot_longer(cols = all_of(Sample.names.vec), 
                                         names_to = "label",
                                         values_to = "Intensity"), Annotation) %>% 
    dplyr::select(label, condition, replicate) %>% 
    distinct()
  
  Before.Imputation.Int <- make_se(proteins_unique = DEP.Input.Int, columns = Sample.names.vec, expdesign = ExpDesign)
  Before.Imputation.LFQ <- make_se(proteins_unique = DEP.Input.LFQ, columns = Sample.names.vec, expdesign = ExpDesign)
  
  ## Data Imputation ####
  
  ## List of Summarised Experiments
  
  Before.Imputation <- list(Freequant.Intensity = Before.Imputation.Int,
                            MaxLFQ.Intensity = Before.Imputation.LFQ)
  
  Imputation.Options <- c("bpca", "knn", "QRILC", "MinDet", "MinProb",
                          "min", "nbavg") ## TODO, make a vector user defined
  
  Imputed.se <- lapply(1:length(Before.Imputation), function(x){
    
    non.imputed.df <- Before.Imputation[[x]]
    
    lapply(1:length(Imputation.Options), function(y){
      
      impute(se = non.imputed.df, fun = Imputation.Options[[y]])
      
    })
  })## End of Imputation function
  
  names(Imputed.se) <- names(Before.Imputation)
  
  names(Imputed.se[[1]]) <- Imputation.Options
  
  names(Imputed.se[[2]]) <- Imputation.Options
  
  
  ## Get DFs of Imputed data
  
  Imputed.df.list <- lapply(1:length(Imputed.se), function(x){
    
    Quant.type.se.list <- Imputed.se[[x]]
    
    lapply(1:length(Imputation.Options), function(y){
      
      a <- get_df_long(se = Quant.type.se.list[[y]])
      
      filtered <- inner_join(Valid.Comparisons.filter %>% 
                               filter(str_detect(Quant.type, names(Imputed.se[x]))), a)  
      
    })
    
  })## End of get df function
  
  names(Imputed.df.list) <- names(Before.Imputation)
  
  names(Imputed.df.list[[1]]) <- Imputation.Options
  
  names(Imputed.df.list[[2]]) <- Imputation.Options
  
  ## Intensity distribution before and after imputation
  
  
  Imputation.Effects.Histogram <- lapply(1:length(Imputed.df.list), function(x){
    
    Quant.type.df.list <- Imputed.df.list[[x]]
    
    lapply(1:length(Imputation.Options), function(y){
      
      Imputation.Distro <- left_join((Quant.type.df.list[[y]] %>% 
                                        select(Protein.ID, intensity, condition, replicate) %>% 
                                        rename(After = intensity)), (Annotated.data %>%
                                                                       filter(str_detect(Quant.type, names(Imputed.df.list[x]))) %>% 
                                                                       select(`Protein ID`, Intensity, condition, replicate) %>% 
                                                                       rename(Before = Intensity,
                                                                              Protein.ID = `Protein ID`) %>% 
                                                                       mutate(Before = log2(Before)))) %>% 
        pivot_longer(cols = c("After", "Before"), names_to = "Imputation Status", values_to = "Intensity") %>% 
        group_by(Protein.ID, condition, `Imputation Status`) %>% 
        summarise(Intensity = mean(Intensity, na.rm = TRUE)) %>%  
        mutate(`Imputation Status` = factor(`Imputation Status`, levels = c("Before", "After"))) %>% 
        ggpubr::ggdensity(x = "Intensity",
                          xlab = "Log2 Intensity", 
                          add = "median", 
                          add_density = TRUE, 
                          title = paste("Log2", names(Imputed.df.list[x]), "before and after", names(Quant.type.df.list[y]), "imputation, per Condition", sep = " "), 
                          y ="..density..",
                          ylab = "Density",
                          color = "Imputation Status",
                          facet.by = c("condition"), 
                          legend = "right", nrow = (length(conditions))%/%3) +
        theme(plot.title = element_text(hjust = 0.5)) 
      
    })
    
  }) ## End of Intensity effects function
  
  
  Histogram.Image.list <- append(Imputation.Effects.Histogram[[1]], Imputation.Effects.Histogram[[2]])
  
  ## Imputation Effects, global mean and CV delta mean and CV
  
  Imputation.effects.global.mean <- lapply(1:length(Imputed.df.list), function(x){
    
    Quant.type.df.list <- Imputed.df.list[[x]]
    
    lapply(1:length(Imputation.Options), function(y){
      
      Imputation.Effect <-  inner_join(inner_join((Quant.type.df.list[[y]] %>% 
                                                     group_by(Protein.ID, condition) %>% 
                                                     summarise(Imputed.mean = mean(intensity, na.rm = TRUE),
                                                               Imputed.sd = sd(intensity),
                                                               Imputed.CV = Imputed.sd/Imputed.mean*100)), (Annotated.data %>% 
                                                                                                              filter(Intensity>0 & str_detect(Quant.type, names(Imputed.df.list[x]))) %>% 
                                                                                                              group_by(`Protein ID`, condition) %>%
                                                                                                              rename(Protein.ID = `Protein ID`) %>% 
                                                                                                              summarise(mean = mean(log2(Intensity)),
                                                                                                                        sd = sd(log2(Intensity)),
                                                                                                                        CV = sd/mean*100))), (Annotated.data %>% 
                                                                                                                                                group_by(`Protein ID`, condition, Quant.type) %>%
                                                                                                                                                filter(Intensity > 0 & str_detect(Quant.type, names(Imputed.df.list[x]))) %>% 
                                                                                                                                                rename(Protein.ID = `Protein ID`) %>% 
                                                                                                                                                summarise(val.Count = n_distinct(Intensity, na.rm = TRUE)))) %>% 
        mutate(missval = case_when(val.Count < 4 ~ TRUE,
                                   val.Count == 4 ~ FALSE)) %>% 
        filter(missval == TRUE)##TODO, vector for valid value cut-off
      
      Imputation.Effect.mean <- Imputation.Effect %>% 
        select(Protein.ID, condition, Imputed.mean, mean) %>% 
        pivot_longer(cols = c("Imputed.mean", "mean"), names_to = "Type", values_to = "mean") %>% 
        ggpubr::ggdensity(x = "mean", 
                          xlab = "Mean Log2 Intensity",
                          facet.by = c("condition"),
                          y = "..density..", 
                          ylab = "Density",
                          color = "Type", 
                          add_density = TRUE, 
                          add = "median",
                          legend = "right",
                          title =  paste("Mean of Log2", names(Imputed.df.list[x]), " of all Protein Groups with missing Values Before and After Imputation with the", names(Quant.type.df.list[y]), "Method", sep = " ")) +
        theme(plot.title = element_text(hjust = 0.5))
      
      
      
    })
    
  })## End of global imputation effects - mean function
  
  Global.mean.Image.list <- append(Imputation.effects.global.mean[[1]], Imputation.effects.global.mean[[2]])
  
  Imputation.effects.global.CV <- lapply(1:length(Imputed.df.list), function(x){
    
    Quant.type.df.list <- Imputed.df.list[[x]]
    
    lapply(1:length(Imputation.Options), function(y){
      
      Imputation.Effect <-  inner_join(inner_join((Quant.type.df.list[[y]] %>% 
                                                     group_by(Protein.ID, condition) %>% 
                                                     summarise(Imputed.mean = mean(intensity, na.rm = TRUE),
                                                               Imputed.sd = sd(intensity),
                                                               Imputed.CV = Imputed.sd/Imputed.mean*100)), (Annotated.data %>% 
                                                                                                              filter(Intensity>0 & str_detect(Quant.type, names(Imputed.df.list[x]))) %>% 
                                                                                                              group_by(`Protein ID`, condition) %>%
                                                                                                              rename(Protein.ID = `Protein ID`) %>% 
                                                                                                              summarise(mean = mean(log2(Intensity)),
                                                                                                                        sd = sd(log2(Intensity)),
                                                                                                                        CV = sd/mean*100))), (Annotated.data %>% 
                                                                                                                                                group_by(`Protein ID`, condition, Quant.type) %>%
                                                                                                                                                filter(Intensity > 0 & str_detect(Quant.type, names(Imputed.df.list[x]))) %>% 
                                                                                                                                                rename(Protein.ID = `Protein ID`) %>% 
                                                                                                                                                summarise(val.Count = n_distinct(Intensity, na.rm = TRUE)))) %>% 
        mutate(missval = case_when(val.Count < 4 ~ TRUE,
                                   val.Count == 4 ~ FALSE)) %>% 
        filter(missval == TRUE)##TODO, vector for valid value cut-off
      
      
      Imputation.Effect.CV <- Imputation.Effect %>% 
        select(Protein.ID, condition, Imputed.CV, CV) %>% 
        pivot_longer(cols = c("Imputed.CV", "CV"), names_to = "Type", values_to = "CV") %>% 
        ggpubr::gghistogram(x = "CV", 
                            xlab = "Coefficient of Variation",
                            facet.by = "condition",
                            y = "..density..", 
                            ylab = "Density",
                            color = "Type", 
                            add_density = TRUE, 
                            add = "median", 
                            legend = "right",
                            title = paste("Coefficient of Variation of Log2", names(Imputed.df.list[x]), "of all Protein Groups with missing Values Before and After Imputation with the", names(Quant.type.df.list[y]), "Method, per Condition", sep = " ")) +
        theme(plot.title = element_text(hjust = 0.5))         
      
      
    })
    
  })## End of global imputation effects - CV function
  
  Global.CV.Image.list <- append(Imputation.effects.global.CV[[1]], Imputation.effects.global.CV[[2]])
  
  Global.mean.CV.list <- append(Global.CV.Image.list, Global.mean.Image.list)
  
  Global.Imputation.Effects.Image.list <- append(Histogram.Image.list, Global.mean.CV.list)
  
  ## Imputation Effects, delta mean and CV
  
  Imputation.Effect.delta <- lapply(1:length(Imputed.df.list), function(x){
    
    bound.res <- bind_rows(Imputed.df.list[[x]], .id = "Imputation.Method")
    
    Imputation.Effect <-  inner_join(inner_join((bound.res %>% 
                                                   group_by(Protein.ID, condition, Imputation.Method) %>% 
                                                   summarise(Imputed.mean = mean(intensity, na.rm = TRUE),
                                                             Imputed.sd = sd(intensity),
                                                             Imputed.CV = Imputed.sd/Imputed.mean*100)), (Annotated.data %>% 
                                                                                                            filter(Intensity>0 & str_detect(Quant.type, names(Imputed.df.list[x]))) %>% 
                                                                                                            group_by(`Protein ID`, condition) %>%
                                                                                                            rename(Protein.ID = `Protein ID`) %>% 
                                                                                                            summarise(mean = mean(log2(Intensity)),
                                                                                                                      sd = sd(log2(Intensity)),
                                                                                                                      CV = sd/mean*100))), (Annotated.data %>% 
                                                                                                                                              group_by(`Protein ID`, condition, Quant.type) %>%
                                                                                                                                              filter(Intensity > 0 & str_detect(Quant.type, names(Imputed.df.list[x]))) %>% 
                                                                                                                                              rename(Protein.ID = `Protein ID`) %>% 
                                                                                                                                              summarise(val.Count = n_distinct(Intensity, na.rm = TRUE)))) %>% 
      mutate(missval = case_when(val.Count < 4 ~ TRUE,
                                 val.Count == 4 ~ FALSE)) %>% 
      filter(missval == TRUE)
    
    Imputation.Effect.DeltaMean <- Imputation.Effect %>% 
      mutate(Delta.mean = Imputed.mean - mean) %>% 
      select(Protein.ID, condition, Delta.mean, Imputation.Method) %>% 
      distinct() %>% 
      ggpubr::ggviolin(x = "condition", 
                       y = "Delta.mean", 
                       color = "Imputation.Method",
                       add = "boxplot", 
                       title = paste("Difference between", names(Imputed.df.list[x]), "means Before and After Imputation, Only Imputed Proteins"), 
                       legend = "right") +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_hline(yintercept = 0, color = 'red')
    
    Imputation.Effect.DeltaCV <- Imputation.Effect %>% 
      mutate(Delta.CV = Imputed.CV - CV) %>% 
      select(Protein.ID, condition, Delta.CV, Imputation.Method) %>% 
      distinct() %>% 
      ggpubr::ggviolin(x = "condition", 
                       y = "Delta.CV", 
                       color = "Imputation.Method",
                       add = "boxplot",
                       legend = "right",
                       title = paste("Difference between", names(Imputed.df.list[x]), "Coefficient of Variation Before and After Imputation, Only Imputed Proteins")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_hline(yintercept = 0, color = 'red')
    
    Image.list <- list(Imputation.Effect.DeltaMean, Imputation.Effect.DeltaCV)
    
  })## End of function for delta stats
  
  Imputation.Delta.Effects.Image.list <- append(Imputation.Effect.delta[[1]], Imputation.Effect.delta[[2]])
  
  ## Save selected data #####
  
  
  if(disk == F) {
    
    print.list <- append(Global.Imputation.Effects.Image.list, Imputation.Delta.Effects.Image.list)
    
  } else {
    
    dir.create("Imputed.data")
    
    lapply(1:length(Imputed.df.list), function(x){
      
      df.list <- Imputed.df.list[[x]]
      
      dir.create(paste("Imputed.data/", names(Imputed.df.list[x]), sep = ""))
      
      lapply(1:length(Imputation.Options), function(y){
        
        dir.create(Imputation.Options[y])
        
        imputed.df <- df.list[[y]]
        
        write_delim(x = imputed.df, file = paste("Imputed.data/", names(Imputed.df.list[x]), "/", Imputation.Options[y], ".tsv", sep = ""), delim = "\t", append = F)
        
      })
      
    })
    
    dir.create("RawData")
    
    write_delim(x = Filtered.Annotated.data, file = "RawData/Filtered.Annotated.Data.tsv", delim = "\t", append = F)
    
    print.list <- append(Global.Imputation.Effects.Image.list, Imputation.Delta.Effects.Image.list)
    
    ggsave(
      filename = "Imputations.pdf", 
      plot = marrangeGrob(print.list, nrow=1, ncol=1), 
      width = 15, height = 9
    )
    
    Output.list <- append(list(Filtered.Annotated.data = Filtered.Annotated.data,
                               Intensity.DEP.df = Intensity.DEP.df,
                               MaxLFQ.DEP.df = MaxLFQ.DEP.df), Imputed.df.list)
    
  }  
  
} ## End of Function



## Wrapped Function ####

## Perform entire workflow in one step (mzML import and ID-free QCs are relatively slow:( )

Workflow.function <- function(mzML.file.dir, Annotation.dir, Anno.Out, Search.Output.dir, 
                              valid.value.threshold, Comparisons, disk, QC.Only){
  ## Load Required Libraries ####
  
  
  library(mzR)
  library(tidyverse)
  library(zoo)
  library(parallel)
  library(ggpubr)
  library(gridExtra)
  library(ComplexUpset)
  library(DEP)
  library(ComplexHeatmap)
  library(circlize)
  library(gridExtra)
  
  ## Load All Functions ####
  
  ## ID-Free ####
  
  ## Function for ID Free QC metrics 
  ## Mostly run-run consistency and instrument performance
  
  ID.Free.QC <- function(mzML.file.dir, Annotation.dir, Anno.Out){
    
    ## Locate raw data files (must be .mzML), R cannot parse .raw/.wiff/,d etc
    
    config[["mzML.file.dir"]]
    setwd(config[["mzML.file.dir"]])
    
    Raw.Files.list <- dir(path = config[["mzML.file.dir"]], recursive = T, pattern = "\\.mzML$", full.names = T)
    
    Raw.file.names <- dir(path = config[["mzML.file.dir"]], recursive = T, pattern = "\\.mzML$", full.names = F)
    
    print("Reading mzML files")
    
    ## Load all mzMLs into R
    
    mzML.list <- lapply(1:length(Raw.Files.list), function(x){
      
      mzR::openMSfile(Raw.Files.list[[x]])
      
      
    })
    
    names(mzML.list) <- Raw.file.names
    
    ## Get Acquisition metadata ####
    
    print("Parsing mzML metadata")
    
    ## for each mzML file returns metadata of Instrument and method
    
    runInfo.list <- lapply(1:length(mzML.list), function(x){
      
      mzR::runInfo(object = mzML.list[[x]]) %>% 
        as.data.frame() %>% 
        mutate(file.name = paste(names(mzML.list[x])))
      
    })
    
    names(runInfo.list) <- names(mzML.list)
    
    ## Get Time stamps for each MS-run, to parse runorder
    
    Time.stamps <- lapply(1:length(runInfo.list), function(x){
      
      y <- runInfo.list[[x]][["startTimeStamp"]] %>% 
        unique() %>% 
        as.data.frame(row.names = names(runInfo.list[x])) 
      
    })
    
    Time.stamps.df <- bind_rows(Time.stamps) %>% 
      distinct() %>% 
      rename(time = 1) %>% 
      separate(col = "time", sep = "T|Z", into = c("day", "time")) %>% 
      separate("time", into = c("hours", "minutes", "seconds"), sep = ":") %>% 
      separate("day", into = c("year", "month", "day"), sep = "-")  %>% 
      mutate(stamp = ((as.numeric(year)*365*24*60*60) + 
                        (as.numeric(month)*30*24*60*60) + 
                        (as.numeric(day)*24*60*60) + 
                        (as.numeric(hours)*60*60) + 
                        (as.numeric(minutes)*60) + 
                        (as.numeric(seconds)))) %>% 
      arrange() %>% 
      mutate(RunOrder = 1:length(Raw.file.names)) %>% 
      rownames_to_column(var = "file.name") %>% 
      select(file.name, RunOrder)
    
    ## dataframe for any sample/biological annotations for dataset
    
    Sample.list <- read_delim(Annotation.dir, 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
    
    Sample.Anno <- inner_join(Sample.list, Time.stamps.df) 
    
    write_delim(x = Sample.Anno, file = paste(Anno.Out, "/Annotation.tsv", sep = ""), delim = "\t", append = F)
    
    runInfo <- inner_join(bind_rows(runInfo.list, .id = "file.name"), Sample.Anno) %>% 
      select(-msLevels) %>% 
      distinct() ## Final Parsed metadata for each MS-run
    
    ## Get mzML file info ####
    
    ## returns information on each MS scan - MS1 and MS/MS scans (MS3 not implemented)
    ## returns scan information: Retention time (s), Injection Time, m/z etc - see ?mzR::header for all info returned
    
    
    print("Parsing Scan Headers")
    
    
    headers.fun <- function(x){
      
      header.list <- lapply(1:unique(runInfo.list[[x]][["scanCount"]]), function(y){
        
        print(paste("Parsing Headers for", names(mzML.list[x])))
        
        header <- mzR::header(object = mzML.list[[x]], scans = y) %>% 
          as.data.frame()
        
      })
      
      header.df <- bind_rows(header.list)
      
    }
    
    ## parallelises function (uses half available cores)
    
    headers <- parallel::mclapply(X = 1:length(mzML.list), 
                                  FUN = headers.fun, 
                                  mc.cores = (detectCores()/2), 
                                  mc.cleanup = T)
    
    
    names(headers) <- names(mzML.list)
    
    ## Separate msLevels, precursor and fragment
    
    ##data for each MS scan
    
    prec.headers <- inner_join(bind_rows(headers, .id = "file.name"), Sample.Anno) %>% 
      filter(msLevel == 1)
    
    ##data for each MS/MS scan
    
    frag.headers <- inner_join(bind_rows(headers, .id = "file.name"), Sample.Anno) %>%
      filter(msLevel == 2)
    
    rm(headers)
    gc()
    
    ## Scan Counts ####
    
    print("Extracting Scan Counts")
    
    ## Returns the number of scans per file (MS and MS/MS)
    
    Scan.Counts <- inner_join(inner_join((prec.headers %>%
                                            group_by(RunOrder, label) %>% 
                                            summarise(`Precursor Scans` = n_distinct(seqNum)/1000)), 
                                         (frag.headers %>% 
                                            group_by(RunOrder, label) %>% 
                                            summarise(`MS/MS Scans` = n_distinct(seqNum)/1000))
    ), (runInfo %>% 
          dplyr::select(RunOrder, scanCount) %>% 
          mutate(scanCount = scanCount/1000) %>% 
          rename(`Total Scans` = scanCount))) %>% 
      pivot_longer(cols = contains("Scans"), names_to = "Scan.Level", values_to = "Scans") %>% 
      mutate(`Scan Type` = factor(Scan.Level, levels = c("Precursor Scans", "MS/MS Scans", "Total Scans"))) %>% 
      arrange("RunOrder") %>% 
      ggpubr::ggbarplot(x = "label",
                        y = "Scans",
                        ylab = "Number of scans (x 10^3)",
                        color = "Scan Type", 
                        label = T, 
                        lab.pos = "out",
                        position = position_dodge(0.9), 
                        title = "Number of Scans per Raw file",
                        legend = "right",
                        lab.nb.digits = 1) +
      ggpubr::rotate_x_text(45) +
      theme(plot.title = element_text(hjust = 0.5))
    
    ## Sensitivity ####
    
    ## InjectionTime
    
    ## returns Injection Times, for all files, as well as individual files
    ## accross retention time
    
    print("Extracting Injection Times")
    
    overlayed.IT <- prec.headers %>%  
      ggpubr::ggline(x = "retentionTime", 
                     xlab = "Rention Time (s)", 
                     y = "injectionTime", 
                     ylab = "Injection Time (ms)",
                     plot_type = "l", 
                     size = .5, 
                     color = "label", 
                     numeric.x.axis = T, 
                     legend = "right", 
                     title = "Injection time of all Precursor Scans, per Retention Time") +
      theme(plot.title = element_text(hjust = 0.5))
    
    
    prec.list <- split.data.frame(prec.headers, f = prec.headers$file.name) 
    
    single.IT <- lapply(1:length(prec.list), function(x){
      
      prec.list[[x]] %>%  
        ggpubr::ggline(x = "retentionTime", 
                       xlab = "Rention Time (s)", 
                       y = "injectionTime", 
                       ylab = "Injection Time (ms)",
                       plot_type = "l", 
                       size = .5, 
                       color = "label", 
                       numeric.x.axis = T, 
                       legend = "right",
                       title = paste("Injection time of Precursor Scans in", names(prec.list[x]), "- per Retention Time", sep = " ")) +
        theme(plot.title = element_text(hjust = 0.5))
      
    })
    
    names(single.IT) <- names(mzML.list)
    
    
    ## Total Ion Current ####
    
    print("Extracting TICs")
    
    ##plots the MS1 TIC for each file (overlayed)
    
    TIC <- prec.headers %>%  
      ggpubr::ggline(x = "retentionTime", 
                     xlab = "Retention Time (s)",
                     y = "totIonCurrent",
                     ylab = "Total Ion Current",
                     numeric.x.axis = T, 
                     color = "label",
                     legend = "right",
                     size = 0.5, 
                     plot_type = "l",
                     title = "Total Ion Current of all Precursor Scans, per Retention Time") +
      theme(plot.title = element_text(hjust = 0.5))
    
    
    ## MS and MS/MS summmaries ####
    
    print("Summarising Intensities")
    
    ## Distribution of MS1 IT
    
    ms1.IT.histo <- prec.headers %>% 
      ggpubr::gghistogram(x = "injectionTime",
                          xlab = "Injection time (ms)",
                          y = "..count..",
                          binwidth = 0.5, 
                          color = "label", 
                          add_density = T, 
                          legend = "right", 
                          facet.by = "condition", 
                          title = "Distribution of Injection times") +
      theme(plot.title = element_text(hjust = 0.5))
    
    ms1.IT.Box <- prec.headers %>%
      arrange(RunOrder) %>% 
      ggpubr::ggviolin(x = "label",
                       y = "injectionTime",
                       ylab = "Injection time (ms)",
                       color = "condition",
                       legend = "right", 
                       title = "Distribution of Injection times", 
                       add = "boxplot") +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggpubr::rotate_x_text(45)
    
    ## Distribution of MS2 IT
    
    ms2.IT.histo <- frag.headers %>% 
      ggpubr::gghistogram(x = "injectionTime",
                          xlab = "Injection time (ms)",
                          y = "..count..",
                          binwidth = 0.5, 
                          color = "label", 
                          add_density = T, 
                          legend = "right", 
                          facet.by = "condition", 
                          title = "Distribution of Injection times") +
      theme(plot.title = element_text(hjust = 0.5))
    
    
    
    ms2.IT.Box <- frag.headers %>%
      arrange(RunOrder) %>% 
      ggpubr::ggviolin(x = "label",
                       y = "injectionTime",
                       ylab = "Injection time (ms)",
                       color = "condition",
                       legend = "right", 
                       title = "Distribution of Injection times", 
                       add = "boxplot") +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggpubr::rotate_x_text(45)
    
    ## Distribution of Frament Intensities
    
    Frag.TIC.histo <- frag.headers %>% 
      mutate(prec.int = log2(precursorIntensity)) %>% 
      mutate(RunOrder = as.character(RunOrder)) %>%
      ggpubr::gghistogram(x = "prec.int",
                          y = "..count..",
                          density = T,
                          legend = "right", 
                          add_density = T,
                          color = "label", 
                          binwidth = 0.1, 
                          add = "median", 
                          facet.by = "condition",
                          title = "Distribution of Precursor Ion Intensities Selected for Fragmentation") +
      theme(plot.title = element_text(hjust = 0.5))
    
    Frag.TIC.Box <- frag.headers %>%
      mutate(prec.int = log2(precursorIntensity)) %>%
      arrange(RunOrder) %>% 
      ggpubr::ggviolin(x = "label",
                       y = "prec.int",
                       ylab = "log2 Intensity of MS/MS precursors",
                       color = "condition",
                       legend = "right", 
                       add = "boxplot",
                       title = "Distribution of log2 Intensity of MS/MS precursors") +
      theme(plot.title = element_text(hjust = 0.5)) +
      rotate_x_text(45)
    
    ## Distribution of Precursor Intensities
    
    TIC.histo <- prec.headers %>% 
      mutate(prec.int = log2(totIonCurrent)) %>% 
      mutate(RunOrder = as.character(RunOrder)) %>%
      ggpubr::gghistogram(x = "prec.int",
                          y = "..count..",
                          add_density = T,
                          legend = "right",
                          color = "label", 
                          binwidth = 0.1,
                          facet.by = "condition",
                          title = "Distribution of all Precursor Ion Intensities") +
      theme(plot.title = element_text(hjust = 0.5))
    
    TIC.Box <- prec.headers %>%
      mutate(prec.int = log2(totIonCurrent)) %>%
      arrange(RunOrder) %>% 
      ggpubr::ggviolin(x = "label",
                       y = "prec.int",
                       ylab = "log2 Intensity of Precursor Ions ",
                       color = "condition",
                       legend = "right", 
                       title = "Distribution of log2 Intensity of all precursors", 
                       add = "boxplot") +
      theme(plot.title = element_text(hjust = 0.5)) +
      rotate_x_text(45)
    
    
    ## MS/MS Performance ####
    
    
    frag.list <- split.data.frame(frag.headers, f = frag.headers$file.name) 
    
    ## TIC of MS and MS/MS scans accross RT
    
    single.Full.TIC <- lapply(1:length(frag.list), function(x){
      
      Scans <- bind_rows(frag.list[[x]], prec.list[[x]])
      
      
      ggplot(data = Scans, mapping = aes(x = retentionTime, color = factor(msLevel))) +
        geom_line(data = (Scans %>% 
                            filter(msLevel == 1)), mapping = aes(y = totIonCurrent)) +
        geom_line(data = (Scans %>% 
                            filter(msLevel == 2)), mapping = aes(y = totIonCurrent)) +
        labs(color = "msLevel") +
        ggtitle(paste("Total Ion Current of MS1 and MS2 of - ", names(prec.list[x]))) +
        xlab("Retention Time (s)") +
        ylab("Total Ion Current") +
        theme(plot.title = element_text(hjust = 0.5))
      
    })
    
    names(single.Full.TIC) <- names(mzML.list)
    
    print("Calculating MS/MS metrics")
    
    ## TopN performance, how often the MS is saturated for MS/MS scans
    
    TopN <- frag.headers %>%
      group_by(RunOrder, label, precursorScanNum) %>% 
      summarise(TopN = n_distinct(precursorMZ),
                retentionTime = mean(retentionTime)) %>%
      mutate(roll.mean = rollmean(x = TopN, k = 100, fill = NA)) %>% 
      distinct() %>% 
      ggpubr::ggline(x = "retentionTime",
                     xlab = "Retention Time (s)",
                     y = "roll.mean", 
                     ylab = "Rolling Average of Number of MS/MS events per Precursor Scan",
                     legend = "right",
                     color = "label", 
                     plot_type = "l",  
                     numeric.x.axis = T, 
                     title = "Rolling average of MS/MS events per Retention Time") +
      theme(plot.title = element_text(hjust = 0.5))
    
    TopN.list <- split.data.frame(TopN[["data"]], f = TopN[["data"]]$label)
    
    single.TopN <- lapply(1:length(TopN.list), function(x){
      
      ggplot2::ggplot(data = TopN.list[[x]], aes(x = retentionTime)) +
        geom_line(aes(y = roll.mean), color = "blue") +
        geom_point(aes(y = TopN), color = "red") +
        ggtitle(paste("Number of Sequential MS/MS events per Retention Time - ", names(prec.list[x]))) +
        xlab("Retention Time (s)") +
        ylab("Sequential MS/MS events") +
        theme(plot.title = element_text(hjust = 0.5))
      
    })
    
    names(single.TopN) <- names(mzML.list)
    
    ## Angle - Score ####
    
    print("Parsing fragmentation scans")
    
    ## Calculate angle score according to: 
    ##"viQC: Visual and Intuitive Quality Control for Mass Spectrometry-Based Proteome Analysis"
    ## "https://link.springer.com/article/10.1134/S1061934819140119"
    ## https://github.com/lisavetasol/viQC
    
    frag.scans <- lapply(1:length(frag.list), function(x){
      
      frag.list[[x]] %>% 
        select(seqNum) %>% 
        distinct() %>% 
        as.matrix() %>% 
        as.vector()
      
    })
    
    angle.score.fun <- function(x){
      
      scans <- frag.scans[[x]]
      
      scan.mat <- lapply(1:length(scans), function(y){
        
        test <- mzR::peaks(object = mzML.list[[x]], scans = scans[[y]]) %>%
          as.data.frame() %>% 
          rename(`m/z` = V1,
                 intensity = V2)
        
        peak.count <- length(test$`m/z`)
        
        avg.int <- test %>% 
          summarise(avg.int = log2(sum(intensity)/peak.count), 
                    peaks = peak.count,
                    seqNum = scans[[y]])
      })
      
      print("Plotting fragmentation performance")
      
      angle.score.mat <- bind_rows(scan.mat) %>% 
        ggpubr::ggscatter(x = "peaks",
                          xlab = "Number of Fragmentation Peaks",
                          y = "avg.int", 
                          ylab = "Average Intensity of Peaks",
                          size = 0.1, 
                          title = paste("Number of Peaks vs Average Intensity of Peaks for -", names(mzML.list[[x]]))) +
        theme(plot.title = element_text(hjust = 0.5))
      
    }
    
    
    angle.score.list <- parallel::mclapply(X = 1:length(frag.scans), 
                                           FUN = angle.score.fun, 
                                           mc.cores = (detectCores()/2), 
                                           mc.cleanup = T)
    
    names(angle.score.list) <- names(mzML.list)
    
    
    ## Dynamic Exclusion
    
    print("Plotting Dynamic Exlusion Performance")
    
    MS2.Exclusion <- lapply(1:length(frag.list), function(x){
      
      single <- split.data.frame(frag.list[[x]] , f = frag.list[[x]]$precursorCharge)
      
      image <- lapply(1:length(single), function(y){
        
        single[[y]] %>% 
          dplyr::select(retentionTime, precursorMZ, precursorCharge) %>% 
          distinct() %>% 
          ggpubr::ggscatter(x = "retentionTime",
                            xlab = "Retention time (s)",
                            y = "precursorMZ",
                            ylab = "Precursor m/z",
                            size = 0.5,
                            title = paste(names(frag.list[x]), "Charge State:", unique(single[[y]]$precursorCharge), sep = " ")) +
          theme(plot.title = element_text(hjust = 0.5))
        
      })
      
      
    })
    
    names(MS2.Exclusion) <- names(mzML.list)
    
    ## Outputs ####
    
    print("Collecting Outputs")
    
    Overlayed.Image.list <- list(Scan.Counts = Scan.Counts,
                                 overlayed.IT = overlayed.IT,
                                 TIC = TIC,
                                 ms1.IT.histo = ms1.IT.histo,
                                 ms1.IT.Box = ms1.IT.Box,
                                 ms2.IT.Box = ms2.IT.Box,
                                 TIC.histo = TIC.histo,
                                 TIC.Box = TIC.Box,
                                 Frag.TIC.histo = Frag.TIC.histo,
                                 Frag.TIC.Box, Frag.TIC.Box,
                                 TopN = TopN)
    
    print("Printing ID free metrics, per mzML to disk")
    
    ggsave(
      filename = "IdFreeQC-metrics.pdf", 
      plot = marrangeGrob(Overlayed.Image.list, nrow=1, ncol=1), 
      width = 16, height = 9
    )
    
    Single.Image.List <- lapply(1:length(Raw.file.names), function(x){
      
      name <- Raw.file.names[x]
      
      Image.list <- append(list(TIC = single.Full.TIC[[name]],
                                IT = single.IT[[name]],
                                TopN = single.TopN[[name]],
                                Angle.Score = angle.score.list[[name]]),MS2.Exclusion[[name]])   
      
      print("Printing summarised ID free metricsto disk")
      
      ggsave(
        filename = paste(name, "IdFreeQC-metrics.pdf"), 
        plot = marrangeGrob(Image.list, nrow=1, ncol=1), 
        width = 16, height = 9
      )
      
    })
    
    Outputs <- list(annotation = Sample.Anno)
    
    
  } ## End of ID free QC
  
  ## msFragger Output - Import ####
  
  msFragger.run_specific_data_import <- function(pattern, Search.Output.dir){
    
    setwd(Search.Output.dir)
    
    annotation <- ID.Free.QCs[["annotation"]]
    
    file.paths <- dir(path = Search.Output.dir, recursive = T, pattern = pattern, full.names = T)
    
    file.names <- dir(path = Search.Output.dir, recursive = T, pattern = pattern, full.names = F)
    
    file.names <- bind_cols(file.names) %>% 
      separate(`...1`, into = "Run", sep = "/") %>% 
      as.matrix() %>% 
      as.vector()
    
    files <- lapply(1:length(file.paths), function(x){
      
      tsv <- grep(x = file.paths[[x]], pattern = "tsv", value = T)
      
      print(paste("Reading file", file.names[x], sep = " "))
      
      if(length(tsv != 0)) {
        
        read_delim(file.paths[[x]], 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE)
        
      } else {read_csv(file = file.paths[[x]])
        
      }
      
    })
    
    names(files) <- file.names
    
    ## Annotate combined files
    
    combined.file <- bind_rows(files, .id = "Run")
    
    inner_join(combined.file, (annotation %>%
                                 mutate(Run = file.name) %>% 
                                 separate(Run, into = "Run", sep = ".mzML")))
    
  } ## End of Function
  
  ## Peptide Level QCs ####
  
  Peptide.level.QCs <- function(Run.Specific.files, disk){
    
    disk = disk ## To save outputs to disk or not
    
    ## Retention Time QCs ####
    
    ## Run.Specific.files - output of msFragger data import
    
    ## Collect RT information and calculate peak-widths
    
    RT <- Run.Specific.files[["quant"]]  %>% 
      select(Run, label, condition, RunOrder, peptide, modified_peptide, charge, contains("retention")) %>%
      select(-retention_time) %>% 
      filter(!is.na(apex_retention_time)) %>% 
      mutate(peak_width.s = (retention_time_end - retention_time_begin),
             retention_time.min = apex_retention_time/60) %>%
      arrange(RunOrder) %>% 
      distinct()
    
    ## Peakwidth distribution
    
    PeakWidth_BoxPlot <- RT %>% 
      filter(peak_width.s <= 80) %>% 
      ggpubr::ggboxplot(x = "label",
                        color = "condition",
                        y = "peak_width.s", 
                        xlab = "Run",
                        ylab = "Peakwidth in Seconds", 
                        title = "Distribution of Peakwidths per Sample", 
                        numeric.x.axis = TRUE,
                        legend = "right") +
      rotate_x_text(45) +
      theme(plot.title = element_text(hjust = 0.5)) 
    
    ## apex retention time variability
    
    RT_Total_CVs <- RT %>% 
      group_by(modified_peptide,  charge) %>% 
      summarise(meanRT = mean(apex_retention_time),
                sdRT = sd(apex_retention_time),
                medianRT = median(apex_retention_time)) %>% 
      dplyr::mutate(cvRT = sdRT/meanRT*100) %>% 
      ggpubr::gghistogram(x = "cvRT", 
                          y = "..density..", 
                          binwidth = 1, 
                          add_density = TRUE, 
                          title = "Distribution of apex retention time CVs",
                          add = "median", 
                          xlab = "Apex Retention Time CV", 
                          ylab = "Density") + 
      xlim(0, 20) +
      geom_vline(xintercept = 5, color = "red") +
      theme(plot.title = element_text(hjust = 0.5))
    
    ## apex retention time variability - per condition (from annotation file)
    
    RT_condition_CVs <- RT %>% 
      group_by(modified_peptide,  charge, condition) %>% 
      summarise(meanRT = mean(apex_retention_time),
                sdRT = sd(apex_retention_time),
                medianRT = median(apex_retention_time)) %>% 
      dplyr::mutate(cvRT = sdRT/meanRT*100) %>% 
      ggpubr::gghistogram(x = "cvRT", 
                          y = "..density..", 
                          color = "condition",
                          binwidth = 1, 
                          add_density = TRUE, 
                          title = "Distribution of apex retention time CVs per condition",
                          add = "median", 
                          xlab = "Apex Retention Time CV", 
                          ylab = "Density", 
                          legend = "right") + 
      xlim(0, 20) +
      geom_vline(xintercept = 5, color = "red") +
      theme(plot.title = element_text(hjust = 0.5)) 
    
    ## Distribution of PSMs accross RT
    
    ID_Density <- RT %>% 
      ggpubr::ggdensity(x = "apex_retention_time",
                        y = "..density..",
                        ylab = "Density of PSMs",
                        color = "Run",
                        add_density = TRUE,
                        title = "PSM density per Retention Time", 
                        xlab = "Apex Retention Time in Seconds", 
                        legend = "right") +
      theme(plot.title = element_text(hjust = 0.5)) 
    
    Total.Samples <- length(RT$Run %>% unique()) 
    
    set.seed(541)
    
    ## RT variability accross selected peptides (Identified in all samples)
    
    Constant_Peptides_RT <- inner_join((RT %>% 
                                          group_by(modified_peptide, label) %>% 
                                          summarise(meanRT = mean(apex_retention_time)) %>% 
                                          summarise(Valid.Values = sum(!is.na(meanRT))) %>% 
                                          filter(Valid.Values == Total.Samples) %>% 
                                          dplyr::select(modified_peptide) %>% 
                                          distinct() %>% 
                                          sample_n(60)), RT) %>% 
      group_by(label, modified_peptide, RunOrder) %>% 
      summarise(RT = median(apex_retention_time)) %>%
      arrange(RunOrder) %>%
      ggpubr::ggline(x = "label", 
                     y = "RT",  
                     color = "modified_peptide", 
                     ylab = "Retention time in seconds", 
                     xlab = "Acquisition Order", 
                     legend = "none",  
                     title = "Retention TIme Stability") +
      theme(plot.title = element_text(hjust = 0.5)) +
      rotate_x_text(45)
    
    
    ## Mass accuracy ####
    
    ## collect accuracy metrics
    
    accuracy <- Run.Specific.files[["quant"]] %>% 
      select(Run, label, condition, RunOrder, modified_peptide, contains("ppm"), hyperscore)
    
    ## Distribution of mass accuracy (in ppm) - per file
    
    mass.accuracy <- accuracy %>%
      pivot_longer(cols = contains("cal"), names_to = "ppm_diff", values_to = "ppm") %>%
      group_by(Run, label, RunOrder) %>%
      arrange(RunOrder) %>% 
      ggpubr::ggboxplot(x = "label", 
                        color = "condition",
                        xlab = "ppm difference",
                        y = "ppm", 
                        facet.by = "ppm_diff", 
                        legend = "right",
                        nrow = 2, 
                        title = "Distribution of the differences between observed m/z and calculated m/z in ppm") +
      ylim(-20, 20) +
      rotate_x_text(45) +
      theme(plot.title = element_text(hjust = 0.5))
    
    
    ## peptide characteristics ####
    
    ## collect identified peptide characteristics
    
    pep.char <- Run.Specific.files[["psm"]]  
    
    ## Distribution of peptide length per sample
    
    pep.len <- pep.char %>% 
      select(Run, label, condition, RunOrder, Peptide, `Peptide Length`) %>% 
      distinct() %>% 
      ggpubr::ggboxplot(x = "label",
                        y = "Peptide Length",
                        color = "condition", 
                        legend = "right", 
                        title = "Distribution of Peptide Lengths") +
      rotate_x_text(45) +
      theme(plot.title = element_text(hjust = 0.5))
    
    ## Distribution of peptide charge state per sample
    
    pep.charge <- pep.char %>% 
      select(Run, label, condition, RunOrder, Peptide, Charge) %>% 
      distinct() %>% 
      group_by(label) %>% 
      summarise(peptide.count = n_distinct(Peptide),
                condition = condition,
                RunOrder = RunOrder,
                Peptide = Peptide,
                Charge = Charge) %>% 
      ungroup() %>% 
      group_by(label, Charge) %>% 
      summarise(Charge.count = n_distinct(Peptide),
                condition = condition,
                RunOrder = RunOrder,
                peptide.count = peptide.count, 
                Charge.freq = Charge.count/peptide.count*100) %>% 
      distinct() %>% 
      arrange(RunOrder) %>% 
      ggpubr::ggbarplot(x = "label",
                        color = "condition",
                        y = "Charge.freq", 
                        ylab = "Frequency",
                        facet.by = "Charge", 
                        nrow = 5, 
                        legend = "right", 
                        title = "Frequency of peptides by charge") +
      rotate_x_text(45) +
      theme(plot.title = element_text(hjust = 0.5))
    
    ## Number of missed cleavages per sample
    
    pep.term <- pep.char %>% 
      select(Run, label, condition, RunOrder, Peptide, `Number of Missed Cleavages`) %>% 
      distinct() %>% 
      group_by(label) %>% 
      summarise(peptide.count = n_distinct(Peptide),
                Peptide = Peptide,
                condition = condition,
                RunOrder = RunOrder,
                `Number of Missed Cleavages` = `Number of Missed Cleavages`) %>% 
      ungroup() %>% 
      group_by(label, `Number of Missed Cleavages`) %>% 
      summarise(MC.count = n_distinct(Peptide), 
                peptide.count = peptide.count,
                condition = condition,
                RunOrder = RunOrder,
                MC.freq = MC.count/peptide.count*100) %>%
      rename(Num.MC = `Number of Missed Cleavages`) %>% 
      distinct() %>% 
      arrange(RunOrder) %>% 
      ggpubr::ggbarplot(x = "label", 
                        y = "MC.freq",
                        ylab = "Frequency of Missed Cleavages",
                        facet.by = "Num.MC", 
                        color = "condition",
                        legend = "right",
                        nrow = 3, title = "Frequency of Missed Cleavages per Sample") +
      rotate_x_text(45) +
      theme(plot.title = element_text(hjust = 0.5))
    
    ## counts of identified peptides per sample ####
    
    peptides <- Run.Specific.files[["ion"]]
    
    peptide.counts <- peptides  %>% 
      group_by(label) %>% 
      summarise(stripped.seq = n_distinct(`Peptide Sequence`), 
                `Peptide Sequence` = `Peptide Sequence`,
                `Modified Sequence` = `Modified Sequence`,
                `Assigned Modifications` = `Assigned Modifications`,
                Charge = Charge, 
                condition = condition,
                RunOrder  = RunOrder) %>% 
      ungroup() %>% 
      group_by(label) %>%
      summarise(stripped.seq = stripped.seq, 
                condition = condition,
                RunOrder  = RunOrder,
                dist.seq = n_distinct(`Modified Sequence`), 
                Mod.Peps = sum(!is.na(`Assigned Modifications`)),
                Cysteine = str_detect(`Assigned Modifications`, "57.0215"),
                Methionine = str_detect(`Assigned Modifications`, "15.9949"),
                `Modified Sequence` = `Modified Sequence`) %>% 
      pivot_longer(cols = c("Methionine", "Cysteine"), names_to = "Residue", values_to = "Modified") %>% 
      filter(Modified == TRUE) %>% 
      ungroup() %>% 
      group_by(label, Residue) %>% 
      summarise(stripped.seq = stripped.seq,
                dist.seq = dist.seq,
                Mod.Peps = Mod.Peps,
                Mod.Res = sum(Modified),
                condition = condition,
                RunOrder  = RunOrder) %>%
      arrange(RunOrder) %>% 
      distinct() %>% 
      pivot_wider(id_cols = c("label", "stripped.seq", "dist.seq", "Mod.Peps", "condition", "RunOrder"), names_from = "Residue", values_from = "Mod.Res")
    
    psm.count <- pep.char %>% 
      select(label, Spectrum) %>% 
      group_by(label) %>% 
      summarise(PSMs = n_distinct(Spectrum))
    
    Scan.counts <- inner_join(psm.count, peptide.counts) %>% 
      pivot_longer(cols = c("PSMs", "stripped.seq", "dist.seq", "Mod.Peps", "Cysteine", "Methionine"), names_to = "Property", values_to = "count") %>% 
      mutate(Property = factor(x = Property, levels = c("PSMs", "dist.seq", "stripped.seq", "Mod.Peps", "Cysteine", "Methionine"))) %>% 
      ggpubr::ggbarplot(x = "label", 
                        y = "count", 
                        color = "condition",
                        legend = "right",
                        facet.by = "Property", 
                        nrow = 3, 
                        scales = "free", 
                        panel.labs = list(Property = c("PSMs", "Distinct Sequences", "Stripped Sequences", "Modified Peptides", "Cysteine Carbidomidomethylation", "Methionine Oxidation")),
                        title = "Peptide and PSM counts") +
      rotate_x_text(45) +
      theme(plot.title = element_text(hjust = 0.5))
    
    
    ## Outputs ####
    
    if(disk == F) {
      
      annotation <- Run.Specific.files[["ion"]] %>% 
        select(Run, file.name, condition, label, Tech.Rep, replicate, RunOrder) %>% 
        distinct()
      
      Output.list <- list(Scan.counts = Scan.counts,
                          mass.accuracy = mass.accuracy,
                          pep.len = pep.len,
                          pep.charge = pep.charge,
                          pep.term = pep.term,
                          RT_Total_CVs = RT_Total_CVs,
                          RT_condition_CVs = RT_condition_CVs,
                          PeakWidth_BoxPlot = PeakWidth_BoxPlot,
                          Constant_Peptides_RT = Constant_Peptides_RT,
                          ID_Density = ID_Density, 
                          annotation = annotation
      )
      
    } else {
      
      Image.list <- list(Scan.counts = Scan.counts,
                         mass.accuracy = mass.accuracy,
                         pep.len = pep.len,
                         pep.charge = pep.charge,
                         pep.term = pep.term,
                         RT_Total_CVs = RT_Total_CVs,
                         RT_condition_CVs = RT_condition_CVs,
                         PeakWidth_BoxPlot = PeakWidth_BoxPlot,
                         Constant_Peptides_RT = Constant_Peptides_RT,
                         ID_Density = ID_Density
      )
      
      ggsave(
        filename = "Peptide_Level_QCs.pdf", 
        plot = marrangeGrob(Image.list, nrow=1, ncol=1), 
        width = 15, height = 9
      )
      
      annotation <- Run.Specific.files[["ion"]] %>% 
        select(Run, file.name, condition, label, Tech.Rep, replicate, RunOrder) %>% 
        distinct()
      
      Output <- list(annotation = annotation)
      
    }
    
  } ## End of Peptide Level QCs
  
  ## Protein level QCs
  
  Protein.level.QCs <- function(Search.Output.dir, valid.value.threshold, Comparisons, disk){
    
    ## Import Protein level data ####
    
    setwd(Search.Output.dir)
    
    annotation <- Peptide.QCs[["annotation"]] ## output from peptide level QCs
    
    combined.proteins <- read_delim("combined_protein.tsv", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE) %>% 
      filter(str_detect(Protein, "contam|rev") == F) %>% 
      select(Protein, `Protein ID`, `Entry Name`, Gene, `Protein Length`, Coverage, Description, contains("Total Intensity")) %>% 
      pivot_longer(cols = contains("Total Intensity"), names_to = "File.Name", values_to = "Intensity") %>% 
      separate(File.Name, into = c("Run", "Quant.type"), sep = " ", extra = "merge") %>% 
      mutate(Quant.type = case_when(str_detect(Quant.type, "MaxLFQ") ~ "MaxLFQ.Intensity",
                                    str_detect(Quant.type, "MaxLFQ") == F ~ "Freequant.Intensity"))
    
    Annotated.data <- inner_join(combined.proteins, annotation)
    
    ## Quantified Protein Counts ####
    
    print("Calculating Protein IDs per Run")
    
    Protein.count.Run <- Annotated.data %>%
      filter(Intensity > 0) %>% 
      group_by(label, Quant.type) %>% 
      summarise(PG.Count = n_distinct(Protein)) %>% 
      ggpubr::ggbarplot(x = "label", 
                        y = "PG.Count", 
                        fill = "Quant.type",
                        legend = "right",
                        position = position_dodge(0.9), 
                        label = TRUE, lab.size = 3) +
      rotate_x_text(45)  
    
    print("Calculating Protein IDs per condition")
    
    Protein.count.condition <- Annotated.data %>%
      filter(Intensity > 0) %>% 
      group_by(condition, Quant.type) %>% 
      summarise(PG.Count = n_distinct(Protein)) %>% 
      ggpubr::ggbarplot(x = "condition", 
                        y = "PG.Count", 
                        fill = "Quant.type",
                        legend = "right",
                        position = position_dodge(0.9), 
                        label = TRUE, lab.size = 3) +
      rotate_x_text(45)
    
    ## Missingness metrics ####
    
    print("Calculating Protein ID overlap per Run")
    
    ## Rates of Identification
    
    Missing.Proportions.Total <- Annotated.data %>%
      dplyr::select(Run, Intensity, Protein, Quant.type) %>% 
      distinct() %>% 
      mutate(PG.Normalised = case_when(Intensity == 0 ~ 0,
                                       Intensity != 0 ~ 1)) %>% 
      dplyr::group_by(Protein, Quant.type) %>% 
      summarise(Valid.Values = sum(PG.Normalised)) %>%   
      distinct() %>% 
      group_by(Quant.type, Valid.Values) %>% 
      summarise(Valid.Counts = n_distinct(Protein)) %>%
      mutate(Valid.Proportion = (Valid.Counts/(Annotated.data$Protein %>% n_distinct()))*100,
             Valid.Values.axis = as.numeric(as.character(Valid.Values))) %>% 
      distinct()  %>%
      ggpubr::ggbarplot(x = "Valid.Values", 
                        xlab = "Number of Samples",
                        y = "Valid.Proportion", 
                        ylab = "Proportion of identified proteins",
                        fill = "Valid.Values.axis",
                        facet.by = "Quant.type",
                        label = T, 
                        lab.nb.digits = 1, 
                        legend = "right", numeric.x.axis = TRUE,
                        title = "Identification frequency as a proportion of all identifications") +
      theme(plot.title = element_text(hjust = 0.5))
    
    print("Calculating Protein ID overlap per condition")
    
    ## How many PGs pass valid value criteria (can be used for statistics later)
    
    Valid.Values.condition <- inner_join(Annotated.data %>%
                                           dplyr::select(Run, condition, Intensity, Protein, Quant.type) %>% 
                                           distinct() %>% 
                                           mutate(PG.Normalised = case_when(Intensity == 0 ~ 0,
                                                                            Intensity != 0 ~ 1)) %>% 
                                           dplyr::group_by(Protein, condition, Quant.type) %>% 
                                           summarise(Valid.Values = sum(PG.Normalised)) %>%   
                                           distinct() %>% 
                                           group_by(condition, Quant.type, Valid.Values) %>% 
                                           summarise(Valid.Counts = n_distinct(Protein)), Protein.count.condition[["data"]]) %>%
      mutate(Valid.Proportion = Valid.Counts/PG.Count,
             valid.value.axis = as.numeric(as.character(Valid.Values))) %>% 
      ggpubr::ggbarplot(x = "Valid.Values", 
                        xlab = "Number of Identifications",
                        y = "Valid.Proportion",
                        ylab = "Proportion of total identifications",
                        fill = "valid.value.axis",
                        facet.by = c("Quant.type", "condition" ), 
                        label = T, 
                        lab.nb.digits = 2, 
                        title = "Frequency of Identifications per Condition",
                        numeric.x.axis = TRUE, legend = "right") +
      theme(plot.title = element_text(hjust = 0.5))
    
    ## Filter data to only include potentially comparable proteins
    
    Valid.value.filter <- Annotated.data %>%
      dplyr::select(Run, condition, Intensity, Protein, Quant.type) %>% 
      distinct() %>% 
      mutate(PG.Normalised = case_when(Intensity == 0 ~ 0,
                                       Intensity != 0 ~ 1)) %>% 
      dplyr::group_by(Protein, condition, Quant.type) %>% 
      summarise(Valid.Values = sum(PG.Normalised)) %>%   
      distinct()
    
    print("Calculating overlaps per condition of potentially comparable protein IDs")
    
    valid.value.threshold <- valid.value.threshold
    
    ## Calculate overlaps of valid values accross conditions
    
    Valid.Values.condition.overlap <- Annotated.data %>%
      dplyr::select(Run, condition, Intensity, Protein, Quant.type) %>% 
      distinct() %>% 
      mutate(PG.Normalised = case_when(Intensity == 0 ~ 0,
                                       Intensity != 0 ~ 1)) %>% 
      dplyr::group_by(Protein, condition, Quant.type) %>% 
      summarise(Valid.Values = sum(PG.Normalised)) %>% 
      filter(Valid.Values >= valid.value.threshold) %>% 
      mutate(Valid = TRUE) %>%
      pivot_wider(id_cols = c("Protein",  "Quant.type"), names_from = "condition", values_from = "Valid", values_fill = FALSE)
    
    overlap.list <- split.data.frame(Valid.Values.condition.overlap, Valid.Values.condition.overlap$Quant.type)
    
    Overlaps <- lapply(1:length(overlap.list), function(x){
      
      Total.Overlap <- ComplexUpset::upset(data = (overlap.list[[x]] %>% 
                                                     select(-Quant.type) %>% 
                                                     column_to_rownames("Protein")), intersect = colnames(overlap.list[[x]] %>%                                             
                                                                                                            select(-Quant.type) %>% 
                                                                                                            column_to_rownames("Protein")), n_intersections = 20)
      
    })
    
    names(Overlaps) <- names(overlap.list)
    
    Comparisons <- read_csv("/media/naadir/Outputs/Workflow.Test/SearchOutputs/Comparisons.txt", 
                            col_names = FALSE) %>% 
      as.matrix() %>% 
      as.vector() ## read comparisons of interest
    
    ## Determine which proteins are comparible
    
    Valid.Comparisons <- lapply(1:length(Comparisons), function(x){
      
      print(paste("Parsing Valid Comparisons for", Comparisons[x]))
      
      contrast <- Comparisons[[x]] %>% 
        as.data.frame() %>% 
        separate_rows(".", sep = "_vs_") %>% 
        rename(condition = 1)
      
      Valid.values.long <- Valid.Values.condition.overlap %>% 
        pivot_longer(cols = 3:(length(colnames(Valid.Values.condition.overlap))), names_to = "condition", values_to = "Valid")
      
      
      Valid.Comparisons <- inner_join(contrast, Valid.values.long) %>% 
        mutate(Con = paste(condition, Quant.type)) %>% 
        pivot_wider(id_cols = c("Protein"), names_from = "Con", values_from = "Valid") %>%
        select(Protein, contains("MaxLFQ"), contains("Freequant")) %>% 
        rename(LFQ_1 = 2,
               LFQ_2 = 3,
               Int_1 = 4,
               Int_2 = 5) %>% 
        mutate(Valid.MaxLFQ.comp = case_when(LFQ_1 == T & LFQ_2 == T ~ TRUE),
               Valid.Freequant.comp = case_when(Int_1 == T & Int_2 == T ~ TRUE),
               Comparison = paste(Comparisons[[x]])) %>% 
        select(Protein, contains("comp")) %>% 
        replace(is.na(.), FALSE)
      
      
    })
    
    print("Calculating overlaps per comparison of comparable protein IDs")
    
    Valid.Comparisons.df <- bind_rows(Valid.Comparisons)
    
    Int.Upset.matrix <- Valid.Comparisons.df %>% 
      select(Protein, Valid.Freequant.comp, Comparison) %>% 
      pivot_wider(id_cols = "Protein", names_from = "Comparison", values_from = "Valid.Freequant.comp") %>% 
      column_to_rownames("Protein")
    
    Int.Upset <- upset(Int.Upset.matrix, intersect = colnames(Int.Upset.matrix), n_intersections = 10) +
      ggtitle("Top 10 Intersections between Freequant Intensity comparisons")
    
    LFQ.Upset.matrix <- Valid.Comparisons.df %>% 
      select(Protein, Valid.MaxLFQ.comp, Comparison) %>% 
      pivot_wider(id_cols = "Protein", names_from = "Comparison", values_from = "Valid.MaxLFQ.comp") %>% 
      column_to_rownames("Protein")
    
    
    LFQ.Upset <- upset(LFQ.Upset.matrix, intersect = colnames(LFQ.Upset.matrix), n_intersections = 10) +
      ggtitle("Top 10 Intersections between MaxLFQ comparisons")
    
    ## Missingness ####
    
    print("Calculating missig data metrics")
    
    Missing.value.data <- Annotated.data %>% 
      as.data.frame() %>% 
      mutate(missing = Intensity == 0,
             Intensity = case_when((missing == FALSE) ~ Intensity))
    
    Missing.Values.stats <- function(data, title, histo){
      
      Missing.data <- data %>% 
        mutate(Intensity = case_when(Intensity != 0 ~ Intensity)) %>% 
        group_by(Protein, condition, Quant.type) %>%  
        summarise(Averagelog2Intensity = mean(log2(Intensity), na.rm = TRUE),
                  MissingValue = any(is.na(Intensity)))
      
      Missing.data.list <- split.data.frame(x = Missing.data, f = Missing.data$condition)
      
      
      fig.list <- lapply(1:length(Missing.data.list), function(z){
        
        a <- Missing.data.list[[z]] %>% 
          ggpubr::gghistogram(x = "Averagelog2Intensity",
                              xlab = "Average log2 Intensity",
                              y = histo, 
                              color = "MissingValue", 
                              facet.by = c("Quant.type"), 
                              nrow = 1, 
                              add_density = TRUE, 
                              title = paste(names(Missing.data.list[z]), "Intensity Distribution of Missing data", sep = ": "), 
                              add = "median", 
                              binwidth = 0.5, legend = "right") +
          theme(plot.title = element_text(hjust = 0.5))
        
        y <-  Missing.data.list[[z]] %>% 
          group_by(`MissingValue`) %>% 
          arrange(`Averagelog2Intensity`) %>% 
          mutate(num = 1, cs = cumsum(num), cs_frac = cs/n()) 
        
        b <- y %>% 
          ggpubr::ggline(x = "Averagelog2Intensity", 
                         xlab = "Average log2 Intensity",
                         y = "cs_frac", 
                         color = "MissingValue",
                         facet.by = c("Quant.type"), 
                         nrow = 1, 
                         numeric.x.axis = TRUE, 
                         plot_type = "l", 
                         ylab = " Cumulative Fraction",
                         title = paste(names(Missing.data.list[z]), "Cumulative Fraction of Missing data per Intensity", sep = ": "),
                         legend = "right") +
          theme(plot.title = element_text(hjust = 0.5))
        
        figure <- ggarrange(a, b , nrow = 2)
        print(figure)
      })
      
      names(fig.list) <- names(Missing.data.list)
      
      out <- fig.list
    }  
    
    Missing.data.fig <- Missing.Values.stats(data = Annotated.data, title = "Missing data distributions", histo = "..count..")  
    
    ##Quant QCs ####
    
    ## Distribution of intensities of identified proteins
    
    Box.Plots <- Annotated.data %>% 
      filter(Intensity > 0) %>% 
      mutate(Log2.Intensity = log2(Intensity)) %>%
      ggpubr::ggboxplot(x = "label", 
                        y = "Log2.Intensity", 
                        fill = "condition", 
                        facet.by = "Quant.type", 
                        nrow = 2, 
                        legend = "right") +
      rotate_x_text(45)
    
    Histograms <- Annotated.data %>% 
      filter(Intensity > 0) %>%  
      mutate(Log2.Intensity = log2(Intensity)) %>%
      ggpubr::ggdensity(x = "Log2.Intensity",
                        xlab = "Log2 Intensity",
                        y = "..density..",
                        add_density = TRUE,
                        color = "condition", 
                        facet.by = "Quant.type", 
                        nrow = 2, 
                        legend = "right", add = "median")
    
    ## Variability of Identified proteins
    
    CVs <- Annotated.data %>%  
      filter(Intensity > 0) %>% 
      mutate(Log2.Intensity = log2(Intensity)) %>% 
      group_by(Protein, condition, Quant.type) %>% 
      summarise(mean = mean(Log2.Intensity, na.rm = TRUE),
                sd = sd(Log2.Intensity, na.rm = TRUE),
                CV  = sd/mean*100) %>% 
      distinct() %>% 
      ggpubr::ggdensity(x = "CV",
                        xlab = "Coeffcient of Variation",
                        y = "..count..",
                        add_density = TRUE,
                        color = "condition", 
                        facet.by = "Quant.type", 
                        nrow = 2, 
                        legend = "right", 
                        add = "median")
    
    
    if(disk == F) {
      
      Image.list <- append(list(Protein.count.Run = Protein.count.Run,
                                Protein.count.condition = Protein.count.condition,
                                Missing.Proportions.Total = Missing.Proportions.Total,
                                Valid.Values.condition = Valid.Values.condition,
                                Int.Upset = Int.Upset,
                                LFQ.Upset = LFQ.Upset,
                                Box.Plots = Box.Plots,
                                Histograms = Histograms,
                                CVs = CVs,
                                Valid.Comparisons.df = Valid.Comparisons.df,
                                Valid.Values.filer = Valid.value.filter,
                                Annotated.data = Annotated.data,
                                
                                
      ), Missing.data.fig)
      
    } else {
      
      Image.list <- append(list(Protein.count.Run = Protein.count.Run,
                                Protein.count.condition = Protein.count.condition,
                                Missing.Proportions.Total = Missing.Proportions.Total,
                                Valid.Values.condition = Valid.Values.condition,
                                Int.Upset = Int.Upset,
                                LFQ.Upset = LFQ.Upset,
                                Box.Plots = Box.Plots,
                                Histograms = Histograms,
                                CVs = CVs
      ), Missing.data.fig)
      
      dir.create("RawData")
      dir.create("Filters")
      
      write_delim(x = Valid.Comparisons.df, file = "Filters/ValidComparisons.tsv", delim = "\t", append = F)
      
      write_delim(x = Valid.value.filter, file = "Filters/ValidValues.tsv", delim = "\t", append = F)
      
      write_delim(x = Annotated.data, file = "RawData/Annotated.data.tsv", delim = "\t", append = F)
      
      ggsave(
        filename = "Protein_Level_QCs.pdf", 
        plot = marrangeGrob(Image.list, nrow=1, ncol=1), 
        width = 15, height = 9
      )
      
      Output.dfs <- list(Valid.Comparisons.df = Valid.Comparisons.df,
                         Valid.Values.filer = Valid.value.filter,
                         Annotated.data = Annotated.data)
      
    }
    
    
  } ## End of Protein QCs
  
  ## Process data for differential expression analysis
  
  Process.data.DE <- function(Valid.value.threshold, Missing.value.threshold, disk){
    
    ## Make Summarized experiment object ####
    
    Annotated.data <- Protein.QCs[["Annotated.data"]]
    
    Valid.value.filter <- Protein.QCs[["Valid.Values.filer"]] %>%  
      filter(Valid.Values >= Valid.value.threshold)
    
    #Protein.QCs[["Valid.Values.filer"]] 
    
    Valid.Comparisons.filter <- inner_join((Protein.QCs[["Valid.Comparisons.df"]]) %>% 
                                             separate_rows(Comparison, sep = "_vs_") %>% 
                                             distinct() %>% 
                                             rename(condition = Comparison,
                                                    Freequant.Intensity = Valid.Freequant.comp,
                                                    MaxLFQ.Intensity = Valid.MaxLFQ.comp) %>% 
                                             pivot_longer(cols = c("Freequant.Intensity", "MaxLFQ.Intensity"), names_to = "Quant.type", values_to = "Valid.Comparison") %>% 
                                             filter(Valid.Comparison == TRUE), (Annotated.data %>% 
                                                                                  select(Protein, `Protein ID`) %>% 
                                                                                  distinct())) %>% 
      rename(Protein.ID = `Protein ID`)
    
    #Protein.QCs[["Valid.Comparisons.df"]]
    
    True.Comparison.filter <- inner_join((Protein.QCs[["Valid.Comparisons.df"]]), (Annotated.data %>% 
                                                                                     select(Protein, `Protein ID`) %>% 
                                                                                     rename(Protein.ID = `Protein ID`) %>% 
                                                                                     distinct())) %>% 
      rename(Freequant.Intensity = Valid.Freequant.comp,
             MaxLFQ.Intensity = Valid.MaxLFQ.comp) %>% 
      pivot_longer(cols = c("Freequant.Intensity", "MaxLFQ.Intensity"), names_to = "Quant.type", values_to = "Valid.Comparison") %>% 
      filter(Valid.Comparison == TRUE)
    
    Filtered.Annotated.data <- inner_join(Annotated.data, Valid.value.filter)
    
    Annotation <- Annotated.data %>% 
      select(label, replicate, condition) %>% 
      distinct()
    
    
    Intensity.DEP.df <- Filtered.Annotated.data %>% 
      filter(str_detect(Quant.type, "Free") & Intensity > 0) %>% 
      select(`Protein ID`, Gene, label, Intensity) %>% 
      pivot_wider(id_cols = c("Protein ID", "Gene"), names_from = "label", values_from = "Intensity")
    
    MaxLFQ.DEP.df <- Filtered.Annotated.data %>% 
      filter(str_detect(Quant.type, "Max") & Intensity > 0) %>% 
      select(`Protein ID`, Gene, label, Intensity) %>% 
      pivot_wider(id_cols = c("Protein ID", "Gene"), names_from = "label", values_from = "Intensity")
    
    
    DEP.Input.Int <- DEP::make_unique(Intensity.DEP.df, names = "Protein ID", ids = "Gene")
    
    DEP.Input.LFQ <- DEP::make_unique(MaxLFQ.DEP.df, names = "Protein ID", ids = "Gene")
    
    conditions <- Annotated.data$condition %>% unique()
    
    Sample.names <- lapply(1:length(conditions), function(x){
      
      grep(conditions[[x]], colnames(DEP.Input.Int))
      
    })
    
    Sample.names.vec <- unlist(Sample.names)
    
    ExpDesign <- inner_join(DEP.Input.Int %>% 
                              pivot_longer(cols = all_of(Sample.names.vec), 
                                           names_to = "label",
                                           values_to = "Intensity"), Annotation) %>% 
      dplyr::select(label, condition, replicate) %>% 
      distinct()
    
    Before.Imputation.Int <- make_se(proteins_unique = DEP.Input.Int, columns = Sample.names.vec, expdesign = ExpDesign)
    Before.Imputation.LFQ <- make_se(proteins_unique = DEP.Input.LFQ, columns = Sample.names.vec, expdesign = ExpDesign)
    
    ## Data Imputation ####
    
    ## List of Summarised Experiments
    
    Before.Imputation <- list(Freequant.Intensity = Before.Imputation.Int,
                              MaxLFQ.Intensity = Before.Imputation.LFQ)
    
    Imputation.Options <- c("bpca", "knn", "QRILC", "MinDet", "MinProb",
                            "min", "nbavg") ## TODO, make a vector user defined
    
    Imputed.se <- lapply(1:length(Before.Imputation), function(x){
      
      non.imputed.df <- Before.Imputation[[x]]
      
      lapply(1:length(Imputation.Options), function(y){
        
        impute(se = non.imputed.df, fun = Imputation.Options[[y]])
        
      })
    })## End of Imputation function
    
    names(Imputed.se) <- names(Before.Imputation)
    
    names(Imputed.se[[1]]) <- Imputation.Options
    
    names(Imputed.se[[2]]) <- Imputation.Options
    
    
    ## Get DFs of Imputed data
    
    Imputed.df.list <- lapply(1:length(Imputed.se), function(x){
      
      Quant.type.se.list <- Imputed.se[[x]]
      
      lapply(1:length(Imputation.Options), function(y){
        
        a <- get_df_long(se = Quant.type.se.list[[y]])
        
        filtered <- inner_join(Valid.Comparisons.filter %>% 
                                 filter(str_detect(Quant.type, names(Imputed.se[x]))), a)  
        
      })
      
    })## End of get df function
    
    names(Imputed.df.list) <- names(Before.Imputation)
    
    names(Imputed.df.list[[1]]) <- Imputation.Options
    
    names(Imputed.df.list[[2]]) <- Imputation.Options
    
    ## Intensity distribution before and after imputation
    
    
    Imputation.Effects.Histogram <- lapply(1:length(Imputed.df.list), function(x){
      
      Quant.type.df.list <- Imputed.df.list[[x]]
      
      lapply(1:length(Imputation.Options), function(y){
        
        Imputation.Distro <- left_join((Quant.type.df.list[[y]] %>% 
                                          select(Protein.ID, intensity, condition, replicate) %>% 
                                          rename(After = intensity)), (Annotated.data %>%
                                                                         filter(str_detect(Quant.type, names(Imputed.df.list[x]))) %>% 
                                                                         select(`Protein ID`, Intensity, condition, replicate) %>% 
                                                                         rename(Before = Intensity,
                                                                                Protein.ID = `Protein ID`) %>% 
                                                                         mutate(Before = log2(Before)))) %>% 
          pivot_longer(cols = c("After", "Before"), names_to = "Imputation Status", values_to = "Intensity") %>% 
          group_by(Protein.ID, condition, `Imputation Status`) %>% 
          summarise(Intensity = mean(Intensity, na.rm = TRUE)) %>%  
          mutate(`Imputation Status` = factor(`Imputation Status`, levels = c("Before", "After"))) %>% 
          ggpubr::ggdensity(x = "Intensity",
                            xlab = "Log2 Intensity", 
                            add = "median", 
                            add_density = TRUE, 
                            title = paste("Log2", names(Imputed.df.list[x]), "before and after", names(Quant.type.df.list[y]), "imputation, per Condition", sep = " "), 
                            y ="..density..",
                            ylab = "Density",
                            color = "Imputation Status",
                            facet.by = c("condition"), 
                            legend = "right", nrow = (length(conditions))%/%3) +
          theme(plot.title = element_text(hjust = 0.5)) 
        
      })
      
    }) ## End of Intensity effects function
    
    
    Histogram.Image.list <- append(Imputation.Effects.Histogram[[1]], Imputation.Effects.Histogram[[2]])
    
    ## Imputation Effects, global mean and CV delta mean and CV
    
    Imputation.effects.global.mean <- lapply(1:length(Imputed.df.list), function(x){
      
      Quant.type.df.list <- Imputed.df.list[[x]]
      
      lapply(1:length(Imputation.Options), function(y){
        
        Imputation.Effect <-  inner_join(inner_join((Quant.type.df.list[[y]] %>% 
                                                       group_by(Protein.ID, condition) %>% 
                                                       summarise(Imputed.mean = mean(intensity, na.rm = TRUE),
                                                                 Imputed.sd = sd(intensity),
                                                                 Imputed.CV = Imputed.sd/Imputed.mean*100)), (Annotated.data %>% 
                                                                                                                filter(Intensity>0 & str_detect(Quant.type, names(Imputed.df.list[x]))) %>% 
                                                                                                                group_by(`Protein ID`, condition) %>%
                                                                                                                rename(Protein.ID = `Protein ID`) %>% 
                                                                                                                summarise(mean = mean(log2(Intensity)),
                                                                                                                          sd = sd(log2(Intensity)),
                                                                                                                          CV = sd/mean*100))), (Annotated.data %>% 
                                                                                                                                                  group_by(`Protein ID`, condition, Quant.type) %>%
                                                                                                                                                  filter(Intensity > 0 & str_detect(Quant.type, names(Imputed.df.list[x]))) %>% 
                                                                                                                                                  rename(Protein.ID = `Protein ID`) %>% 
                                                                                                                                                  summarise(val.Count = n_distinct(Intensity, na.rm = TRUE)))) %>% 
          mutate(missval = case_when(val.Count < 4 ~ TRUE,
                                     val.Count == 4 ~ FALSE)) %>% 
          filter(missval == TRUE)##TODO, vector for valid value cut-off
        
        Imputation.Effect.mean <- Imputation.Effect %>% 
          select(Protein.ID, condition, Imputed.mean, mean) %>% 
          pivot_longer(cols = c("Imputed.mean", "mean"), names_to = "Type", values_to = "mean") %>% 
          ggpubr::ggdensity(x = "mean", 
                            xlab = "Mean Log2 Intensity",
                            facet.by = c("condition"),
                            y = "..density..", 
                            ylab = "Density",
                            color = "Type", 
                            add_density = TRUE, 
                            add = "median",
                            legend = "right",
                            title =  paste("Mean of Log2", names(Imputed.df.list[x]), " of all Protein Groups with missing Values Before and After Imputation with the", names(Quant.type.df.list[y]), "Method", sep = " ")) +
          theme(plot.title = element_text(hjust = 0.5))
        
        
        
      })
      
    })## End of global imputation effects - mean function
    
    Global.mean.Image.list <- append(Imputation.effects.global.mean[[1]], Imputation.effects.global.mean[[2]])
    
    Imputation.effects.global.CV <- lapply(1:length(Imputed.df.list), function(x){
      
      Quant.type.df.list <- Imputed.df.list[[x]]
      
      lapply(1:length(Imputation.Options), function(y){
        
        Imputation.Effect <-  inner_join(inner_join((Quant.type.df.list[[y]] %>% 
                                                       group_by(Protein.ID, condition) %>% 
                                                       summarise(Imputed.mean = mean(intensity, na.rm = TRUE),
                                                                 Imputed.sd = sd(intensity),
                                                                 Imputed.CV = Imputed.sd/Imputed.mean*100)), (Annotated.data %>% 
                                                                                                                filter(Intensity>0 & str_detect(Quant.type, names(Imputed.df.list[x]))) %>% 
                                                                                                                group_by(`Protein ID`, condition) %>%
                                                                                                                rename(Protein.ID = `Protein ID`) %>% 
                                                                                                                summarise(mean = mean(log2(Intensity)),
                                                                                                                          sd = sd(log2(Intensity)),
                                                                                                                          CV = sd/mean*100))), (Annotated.data %>% 
                                                                                                                                                  group_by(`Protein ID`, condition, Quant.type) %>%
                                                                                                                                                  filter(Intensity > 0 & str_detect(Quant.type, names(Imputed.df.list[x]))) %>% 
                                                                                                                                                  rename(Protein.ID = `Protein ID`) %>% 
                                                                                                                                                  summarise(val.Count = n_distinct(Intensity, na.rm = TRUE)))) %>% 
          mutate(missval = case_when(val.Count < 4 ~ TRUE,
                                     val.Count == 4 ~ FALSE)) %>% 
          filter(missval == TRUE)##TODO, vector for valid value cut-off
        
        
        Imputation.Effect.CV <- Imputation.Effect %>% 
          select(Protein.ID, condition, Imputed.CV, CV) %>% 
          pivot_longer(cols = c("Imputed.CV", "CV"), names_to = "Type", values_to = "CV") %>% 
          ggpubr::gghistogram(x = "CV", 
                              xlab = "Coefficient of Variation",
                              facet.by = "condition",
                              y = "..density..", 
                              ylab = "Density",
                              color = "Type", 
                              add_density = TRUE, 
                              add = "median", 
                              legend = "right",
                              title = paste("Coefficient of Variation of Log2", names(Imputed.df.list[x]), "of all Protein Groups with missing Values Before and After Imputation with the", names(Quant.type.df.list[y]), "Method, per Condition", sep = " ")) +
          theme(plot.title = element_text(hjust = 0.5))         
        
        
      })
      
    })## End of global imputation effects - CV function
    
    Global.CV.Image.list <- append(Imputation.effects.global.CV[[1]], Imputation.effects.global.CV[[2]])
    
    Global.mean.CV.list <- append(Global.CV.Image.list, Global.mean.Image.list)
    
    Global.Imputation.Effects.Image.list <- append(Histogram.Image.list, Global.mean.CV.list)
    
    ## Imputation Effects, delta mean and CV
    
    Imputation.Effect.delta <- lapply(1:length(Imputed.df.list), function(x){
      
      bound.res <- bind_rows(Imputed.df.list[[x]], .id = "Imputation.Method")
      
      Imputation.Effect <-  inner_join(inner_join((bound.res %>% 
                                                     group_by(Protein.ID, condition, Imputation.Method) %>% 
                                                     summarise(Imputed.mean = mean(intensity, na.rm = TRUE),
                                                               Imputed.sd = sd(intensity),
                                                               Imputed.CV = Imputed.sd/Imputed.mean*100)), (Annotated.data %>% 
                                                                                                              filter(Intensity>0 & str_detect(Quant.type, names(Imputed.df.list[x]))) %>% 
                                                                                                              group_by(`Protein ID`, condition) %>%
                                                                                                              rename(Protein.ID = `Protein ID`) %>% 
                                                                                                              summarise(mean = mean(log2(Intensity)),
                                                                                                                        sd = sd(log2(Intensity)),
                                                                                                                        CV = sd/mean*100))), (Annotated.data %>% 
                                                                                                                                                group_by(`Protein ID`, condition, Quant.type) %>%
                                                                                                                                                filter(Intensity > 0 & str_detect(Quant.type, names(Imputed.df.list[x]))) %>% 
                                                                                                                                                rename(Protein.ID = `Protein ID`) %>% 
                                                                                                                                                summarise(val.Count = n_distinct(Intensity, na.rm = TRUE)))) %>% 
        mutate(missval = case_when(val.Count < 4 ~ TRUE,
                                   val.Count == 4 ~ FALSE)) %>% 
        filter(missval == TRUE)
      
      Imputation.Effect.DeltaMean <- Imputation.Effect %>% 
        mutate(Delta.mean = Imputed.mean - mean) %>% 
        select(Protein.ID, condition, Delta.mean, Imputation.Method) %>% 
        distinct() %>% 
        ggpubr::ggviolin(x = "condition", 
                         y = "Delta.mean", 
                         color = "Imputation.Method",
                         add = "boxplot", 
                         title = paste("Difference between", names(Imputed.df.list[x]), "means Before and After Imputation, Only Imputed Proteins"), 
                         legend = "right") +
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_hline(yintercept = 0, color = 'red')
      
      Imputation.Effect.DeltaCV <- Imputation.Effect %>% 
        mutate(Delta.CV = Imputed.CV - CV) %>% 
        select(Protein.ID, condition, Delta.CV, Imputation.Method) %>% 
        distinct() %>% 
        ggpubr::ggviolin(x = "condition", 
                         y = "Delta.CV", 
                         color = "Imputation.Method",
                         add = "boxplot",
                         legend = "right",
                         title = paste("Difference between", names(Imputed.df.list[x]), "Coefficient of Variation Before and After Imputation, Only Imputed Proteins")) +
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_hline(yintercept = 0, color = 'red')
      
      Image.list <- list(Imputation.Effect.DeltaMean, Imputation.Effect.DeltaCV)
      
    })## End of function for delta stats
    
    Imputation.Delta.Effects.Image.list <- append(Imputation.Effect.delta[[1]], Imputation.Effect.delta[[2]])
    
    ## Save selected data #####
    
    
    if(disk == F) {
      
      print.list <- append(Global.Imputation.Effects.Image.list, Imputation.Delta.Effects.Image.list)
      
    } else {
      
      dir.create("Imputed.data")
      
      lapply(1:length(Imputed.df.list), function(x){
        
        df.list <- Imputed.df.list[[x]]
        
        dir.create(paste("Imputed.data/", names(Imputed.df.list[x]), sep = ""))
        
        lapply(1:length(Imputation.Options), function(y){
          
          dir.create(Imputation.Options[y])
          
          imputed.df <- df.list[[y]]
          
          write_delim(x = imputed.df, file = paste("Imputed.data/", names(Imputed.df.list[x]), "/", Imputation.Options[y], ".tsv", sep = ""), delim = "\t", append = F)
          
        })
        
      })
      
      dir.create("RawData")
      
      write_delim(x = Filtered.Annotated.data, file = "RawData/Filtered.Annotated.Data.tsv", delim = "\t", append = F)
      
      print.list <- append(Global.Imputation.Effects.Image.list, Imputation.Delta.Effects.Image.list)
      
      ggsave(
        filename = "Imputations.pdf", 
        plot = marrangeGrob(print.list, nrow=1, ncol=1), 
        width = 15, height = 9
      )
      
      Output.list <- append(list(Filtered.Annotated.data = Filtered.Annotated.data,
                                 Intensity.DEP.df = Intensity.DEP.df,
                                 MaxLFQ.DEP.df = MaxLFQ.DEP.df), Imputed.df.list)
      
    }  
    
  } ## End of Function
  
  
  
  ## ID Free ####
  
  if (QC.Only == T) {
    
    ID.Free.QCs <- ID.Free.QC(mzML.file.dir = config[["mzML.file.dir"]],
                              Annotation.dir = config[["Annotation.dir"]], 
                              Anno.Out = config[["Anno.Out"]])
    
    fragger.file.patterns <- c("_quant.csv", "^ion\\.tsv$", "^peptide\\.tsv$", "^protein\\.tsv$", "^psm\\.tsv$")
    fragger.file.types <- c("quant", "ion", "peptide", "protein", "psm")
    
    Run.Specific.files <-  lapply(1:length(fragger.file.patterns), function(x){
      
      msFragger.run_specific_data_import(pattern = fragger.file.patterns[[x]], 
                                         Search.Output.dir = config[["Search.Output.dir"]])
      
    }) 
    
    names(Run.Specific.files) <- fragger.file.types
    
    Peptide.QCs <- Peptide.level.QCs(Run.Specific.files = Run.Specific.files, 
                                     disk = T)
    
    
    Protein.QCs <- Protein.level.QCs(Search.Output.dir = config[["Search.Output.dir"]], 
                                     valid.value.threshold = config[["valid.value.threshold"]], 
                                     Comparisons = config[["Comparisons"]], 
                                     disk = config[["disk"]])
    
    Outputs <- list(Protein.QCs = Protein.QCs, 
                    Peptide.QCs = Peptide.QCs,
                    Comparisons = config[["Comparisons"]])
    
  } else {
    
    
    
    ID.Free.QCs <- ID.Free.QC(mzML.file.dir = config[["mzML.file.dir"]],
                              Annotation.dir = config[["Annotation.dir"]], 
                              Anno.Out = config[["Anno.Out"]])
    
    fragger.file.patterns <- c("_quant.csv", "^ion\\.tsv$", "^peptide\\.tsv$", "^protein\\.tsv$", "^psm\\.tsv$")
    fragger.file.types <- c("quant", "ion", "peptide", "protein", "psm")
    
    Run.Specific.files <-  lapply(1:length(fragger.file.patterns), function(x){
      
      msFragger.run_specific_data_import(pattern = fragger.file.patterns[[x]], 
                                         Search.Output.dir = config[["Search.Output.dir"]])
      
    }) 
    
    names(Run.Specific.files) <- fragger.file.types
    
    Peptide.QCs <- Peptide.level.QCs(Run.Specific.files = Run.Specific.files, 
                                     disk = T)
    
    
    Protein.QCs <- Protein.level.QCs(Search.Output.dir = config[["Search.Output.dir"]], 
                                     valid.value.threshold = config[["valid.value.threshold"]], 
                                     Comparisons = config[["Comparisons"]], 
                                     disk = config[["disk"]])
    
    Processed.Data <- Process.data.DE(Valid.value.threshold =  config[["valid.value.threshold"]], 
                                      Missing.value.threshold = 3, 
                                      disk = config[["disk"]])
    
    
    Outputs <- list(Protein.QCs = Protein.QCs, 
                    Processed.Data = Processed.Data,
                    Peptide.QCs = Peptide.QCs,
                    Comparisons = config[["Comparisons"]])
    
  }
  
  
}

## Total Workflow ####


Processed.Data <- Workflow.function(mzML.file.dir = config[["mzML.file.dir"]], 
                                    Annotation.dir = config[["Annotation.dir"]], 
                                    Anno.Out = config[["Anno.Out"]], 
                                    Search.Output.dir = config[["Search.Output.dir"]], 
                                    valid.value.threshold = config[["valid.value.threshold"]], 
                                    Comparisons = config[["Comparisons"]], 
                                    disk = config[["disk"]], 
                                    QC.Only = config[["QC.Only"]]) ## Change names of imputed df lists to Quant type
## TODO, update angle score plot, with PSMs, by colour

names(Processed.Data[["Processed.Data"]]) <- c("Filtered.Annotated.data", "Intensity.DEP.df", "MaxLFQ.DEP.df", "Freequant", "MaxLFQ")

saveRDS(object = Processed.Data, file = "Processed.Data.RDS")