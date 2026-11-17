#Author: Asher Ferguson
#Last Updated: 11/10/2026
#Email: asher-ferguson@uiowa.edu
#Script Name: Mt. Bigelow Dendrometer Data Analysis
#Script Description: 

#instructions:
#make sure all necessary packages are installed
#set working directory under "file setup"
#run everything until "BELOW THIS LINE ARE TASKS..."

#imports----
library(tidyverse)
library(ggh4x)
library(leaps)
library(caret)
library(springpheno)
library(lme4)
library(nlme)
library(lmerTest)
library(car)

#file setup ----
#replace the working directory with your own
#Replace the working directory with the path on your computer holding the dendrometer, sap flow, and flux data
#setwd("C:\\Users\\asferguson\\OneDrive - University of Iowa\\RA\\Shared Folder w Dannenberg\\data")
#setwd("C:\\Users\\mdannenberg\\OneDrive - University of Iowa\\Ferguson - shared folder\\data")
setwd("C:\\Users\\fergu\\OneDrive - University of Iowa\\RA\\Shared Folder w Dannenberg\\data")

dendro <- read.csv("MtB_DendrometerData.csv")
flux <- read.csv("MtB_FluxData.csv")
sapflow <- read.csv("MtB_SapFlow_SpeciesAverage_2025.csv")

#useful functions----
movingAverage <- function(vector, window, centered = TRUE){
  output <- c()
  for(i in 1:length(vector)){
    if(centered==FALSE){
      start = max(1, i-window+1)
      end = i
    }
    else{
      half = trunc(window/2)
      start = max(1, i-half)
      end = min(length(vector), i+half)
    }
    output[i] <- mean(vector[start:end], na.rm = TRUE)
  }
  return(output)
}

#dendro data processing----
formatDendroData <- function(dendro){
  dendro <- dendro[!is.na(dendro$value), ] #remove NA values
  
  #explicitly add 0:00 to midnight time stamps
  null_time <- nchar(as.character(dendro$ts)) < 11
  dendro$ts[null_time] <- paste(dendro$ts[null_time], "0:00")
  
  #turn ts into date column
  dendro$datetime <- parse_date_time(dendro$ts, orders = "mdy HM")
  
  #get year and doy
  dendro$year <- year(dendro$datetime)
  dendro$doy <- yday(dendro$datetime)
  
  #add decimal hours
  dendro$minute <- dendro$datetime |>
    substr(15,16) |>
    as.numeric()
  dendro$hours <- dendro$datetime |> 
    substr(12,13) |>
    as.numeric()
  dendro$hours <- replace(dendro$hours, is.na(dendro$hours), 0)
  dendro$minute <- replace(dendro$minute, is.na(dendro$minute), 0)
  dendro$hours <- dendro$hours + (dendro$minute/60)
  dendro$minute <- NULL
  
  #create single date variable
  dendro$date <- as.Date(dendro$datetime)
  
  #remove dying tree (18842)
  dendro <- dendro[dendro$series != "x18842", ]
  
  #add species data
  dendro$species <- ifelse(
    dendro$series == "x18835" | dendro$series == "x18836" | dendro$series == "x18837" | dendro$series == "x18838",
    "PSME",
    "PIPO"
  )
  
  #decimalDOY
  dendro$doyDec <- round(dendro$doy + (dendro$hours / 24), 4)
  
  #net growth columns
  dendro$growth <- dendro$value - lag(dendro$value)
  
  #new irreversible growth (above a previous maximum)
  dendro$newGrowth <- dendro$max - lag(dendro$max)
  
  dendro$isNewGrowth <- ifelse(dendro$newGrowth > 0, TRUE, FALSE)
  
  return(dendro)
}

dendro <- formatDendroData(dendro)
dendro <- dendro[c(2,3:7,11:20)]

#summarize half-hourly data by species
dendroSppSums <- dendro |> summarize(value = mean(value),
                                     max = mean(max),
                                     twd = mean(twd),
                                     gro_yr = mean(gro_yr),
                                     doy = first(doy),
                                     hours = first(hours),
                                     date = first(date),
                                     growth = mean(growth),
                                     n_newGrowth = sum(isNewGrowth == TRUE, na.rm = TRUE),
                                     n_trees = sum(!is.na(isNewGrowth)),
                                     prop_newGrowth = mean(isNewGrowth),
                                     .by = c(species, doyDec, year))

#create dlySums dataset
dlySums <- dendro |> summarize(max_value = max(value),
                               min_value = min(value),
                               mean_value = mean(value),
                               max_twd = max(twd),
                               min_twd = min(twd),
                               mean_twd = mean(twd),
                               predawn_twd = ifelse(
                                 any(!is.na(twd[hours>=4 & hours<=8])),
                                 min(twd[hours>=4 & hours<=8], na.rm=TRUE),
                                 min(twd)
                               ),#predawn twd is calculated as the minimum value observed between 4:00 and 8:00. Or the min twd for the day if somehow it would be NA.
                               max_value_time = hours[which.max(value)],
                               min_value_time = hours[which.min(value)],
                               max_twd_time = hours[which.max(twd)],
                               min_twd_time = hours[which.min(twd)],
                               prevMax = max(max),
                               year = first(year),
                               doy = first(doy),
                               species = first(species),
                               totalGrowth = sum(growth, na.rm=TRUE),
                               .by=c("date", "series"))
dlySums$newGrowth <- dlySums$prevMax - lag(dlySums$prevMax)
dlySums$mds <- dlySums$max_twd-dlySums$min_twd
dlySums$mdsNORM <- dlySums$mds/ ave(dlySums$mds, dlySums$series, FUN = max)
dlySums$twdNORM <- dlySums$predawn_twd/ ave(dlySums$mds, dlySums$series, FUN = max)

#daily summarized dendro by species
dlySumsSpp <- dlySums |> summarize(max_value = mean(max_value),
                                   min_value = mean(min_value),
                                   mean_value = mean(mean_value),
                                   max_twd = mean(max_twd),
                                   min_twd = mean(min_twd),
                                   mean_twd = mean(mean_twd),
                                   predawn_twd = mean(predawn_twd),#predawn twd is calculated as the maximum value observed between 4:00 and 8:00
                                   year = first(year),
                                   doy = first(doy),
                                   mds = mean(mds),
                                   mdsNORM = mean(mdsNORM),
                                   twdNORM = mean(twdNORM),
                                   totalGrowth = mean(totalGrowth),
                                   newGrowth = mean(newGrowth),
                                   .by=c("date", "species"))

#flux data processing----
flux$doy <- flux$DoY
flux$DoY <- NULL
flux$date <- as.POSIXct(
  paste0(flux$Year, "-01-01")
) + (flux$doy-1) * 86400
flux$doyDec <- round(flux$doy + flux$Hour/24, 4)

flxDlySums <- summarize(flux,
                        Ustar = mean(Ustar, na.rm = TRUE),
                        NEE = mean(NEE_f, na.rm = TRUE),
                        GPP = mean(GPP_f, na.rm = TRUE),
                        Reco = mean(Reco_f, na.rm = TRUE),
                        LE = mean(LE_f, na.rm = TRUE),
                        H = mean(H_f, na.rm = TRUE),
                        Rg = mean(Rg_f, na.rm = TRUE),
                        Tair = mean(Tair_f, na.rm = TRUE),
                        rH = mean(rH, na.rm = TRUE),
                        VPD = mean(VPD_f, na.rm = TRUE),
                        PA = mean(PA, na.rm = TRUE),
                        .by = c("doy", "Year"))
flxDlySums$date <- as.Date(paste(flxDlySums$Year, 1, 1, sep = "-")) + flxDlySums$doy - 1

#daily average of flux for all years
AverageYearFlux <- summarize(flux, 
                             Ustar = mean(Ustar, na.rm = TRUE),
                             NEE = mean(NEE_f, na.rm = TRUE),
                             GPP = mean(GPP_f, na.rm = TRUE),
                             Reco = mean(Reco_f, na.rm = TRUE),
                             LE = mean(LE_f, na.rm = TRUE),
                             H = mean(H_f, na.rm = TRUE),
                             Rg = mean(Rg_f, na.rm = TRUE),
                             Tair = mean(Tair_f, na.rm = TRUE),
                             rH = mean(rH, na.rm = TRUE),
                             VPD = mean(VPD_f, na.rm = TRUE),
                             .by = "doy"
)

#sap flow data processing----
names(sapflow)[1] <- "doy"
sapflowSppLong <- sapflow |> pivot_longer(cols = c("Average_PIPO", "Average_PISF", "Average_PSME"),
                                          names_to = "spp",
                                          values_to = "value")
sapflowSppLong <- sapflowSppLong[,c(1, 13, 14)]
sapflowSppLong$spp <- sapflowSppLong$spp |> substr(9,12)
sapflowSppLong$spp <- as.factor(sapflowSppLong$spp)
names(sapflowSppLong)[3] <- "sapflow"

sapflowSppLong$day <- trunc(sapflowSppLong$doy, 0)

#create sap flow daily summaries
sapflowDlySums <- sapflowSppLong |> summarize(sapflow = mean(sapflow),
                                              .by = c(day, spp))

#getting day length----
daylength(365, 32.416)  #function is from the springpheno package
daylength <- data.frame(doy = 1:365, daylength = daylength(365, 32.416))

#gathering climate variables----
sapJoinData <- sapflow[,c(1,4,6)]
sapJoinData$date <- as.Date(sapJoinData$doy-1, origin = "2025-01-01")
names(sapJoinData) <- c("doy", "sapPIPO", "sapPSME", "date")
climVars <- left_join(flux, sapJoinData, by = c("date", "doyDec" = "doy"))
climVars$doy.y <- NULL
names(climVars)[names(climVars) == "doy.x"] <- "doy"
climVars$SWC <- rowMeans(
  climVars[,c("SWC_1","SWC_2","SWC_3")],
  na.rm = TRUE
)

#climate variables and dendro data in the same sheet
climGrowth <- left_join(climVars,
                        dendroSppSums,
                        by = c("Year" = "year", "doyDec"))


#BELOW THIS LINE ARE TASKS THAT CAN BE RUN INDIVIDUALLY----

#exploratory graphs----
#dendro: cumulative growth by tree, 2025
dendro[dendro$year == 2025,] |> ggplot(aes(doyDec, gro_yr, group = series, color = series)) + 
  geom_line() +
  theme_bw() +
  labs(
    x = "Date (DOY) of 2025",
    y = "Cumulative growth (μm)",
    color = "Tree ID"
  )
ggsave(
  "outputs/cumGrowthbytree.png",
  width = 12,
  height = 6,
  dpi = 300
)

#dendro: total cumulative growth by DOY, 2025
totalGrowthbyDOY <- dendro |> summarize(growth = sum(gro_yr),
                                        .by = doyDec)
totalGrowthbyDOY |> ggplot(aes(doyDec, growth)) +
  geom_line() +
  geom_hline(yintercept = max(totalGrowthbyDOY$growth) * 0.05, color = "red") +
  geom_hline(yintercept = max(totalGrowthbyDOY$growth) * 0.95, color = "red") +
  theme_bw() +
  labs(
    x = "Date (DOY)",
    y = "Total cumulative growth (μm)",
  ) + 
  annotate("text",
           x = 50,
           y = max(totalGrowthbyDOY$growth) * 0.06,
           label = "5% of growth",
           color = "black") +
  annotate("text",
           x = 50,
           y = max(totalGrowthbyDOY$growth) * 0.94,
           label = "95% of growth",
           color = "black")
ggsave("outputs/totalCumGrowth.png",
       width = 12,
       height = 6,
       dpi = 300)

#dendro: relative growth by tree
growthbytree2025 <- dendro[dendro$year == 2025, ] |> pivot_wider(id_cols = doyDec, names_from = series, values_from = gro_yr)
relGrowthbytree <- growthbytree2025
relGrowthbytree[,-1] <- lapply(relGrowthbytree[,-1], function(x){
  x/max(x,na.rm = TRUE)
})
relGrowthbytreeLong <- relGrowthbytree |> pivot_longer(cols = !"doyDec",
                                                       names_to = "series",
                                                       values_to = "cumNormGrowth")
relGrowthbytreeLong$species <- ifelse(relGrowthbytreeLong$series == "x18839" | relGrowthbytreeLong$series == "x18840" | relGrowthbytreeLong$series == "x18841",
                                      "PIPO",
                                      "PSME")
relGrowthbySpp <- relGrowthbytreeLong |> summarize(cumNormGrowth = mean(cumNormGrowth, na.rm = TRUE),
                                                   .by = c(doyDec, species)
)
relGrowthbySpp |> ggplot(aes(doyDec, cumNormGrowth, color = species)) +
  geom_line() +
  labs(
    x = "Date (DOY)",
    y = "Cumulative Normalized Growth"
  )
ggsave("outputs/cumulativeGrowthbySpp.png",
       width = 12,
       height = 6,
       dpi = 300)

#flux: climate vs 2025 weather
plotdata1 <- left_join(AverageYearFlux,
                       flxDlySums,
                       by = "doy"
)
plotdata1$VPD.y <- movingAverage(plotdata1$VPD.y, 10)
plotdata1$Tair.y <- movingAverage(plotdata1$Tair.y, 10)
plotdata1$rH.y <- movingAverage(plotdata1$rH.y, 10)
plotdata1 <- plotdata1 |> pivot_longer(cols = c(VPD.x, Tair.x, rH.x, VPD.y, Tair.y, rH.y),
                                       names_to = c("variable", "category"),
                                       names_sep = "\\.",
                                       values_to = "value"
)
plotdata1$category[plotdata1$category == "x"] <- "average"
plotdata1$category[plotdata1$category == "y"] <- "2025 (10-day smoothed)"

min(plotdata1$value[plotdata1$variable == "VPD"], na.rm = TRUE)

ggplot(plotdata1, aes(doy, value)) + 
  geom_line() + 
  facet_grid2(category ~ variable, 
              scales = "free_y", 
              independent = "y",
              labeller = labeller(
                variable = c(
                  rH = "Relative Humidity (%)",
                  Tair = "Air temperature (°C)",
                  VPD = "Vapor pressure deficit (hPa)"
                )
              )
  ) +
  facetted_pos_scales(
    y = list(
      variable == "rH" ~ scale_y_continuous(limits = c(0,100)),
      variable == "Tair" ~ scale_y_continuous(limits = c(-10, 35)),
      variable == "VPD" ~ scale_y_continuous(limits = c(0,30))
    )
  ) +
  labs(
    x = "Date (DOY)",
    y = NULL
  ) +
  theme_bw()

ggsave("outputs/2025vsclimate.png",
       width = 12,
       height = 6,
       dpi = 300)

#sap flow
sapflowDlySums |> ggplot(aes(day, sapflow, color = spp)) + 
  geom_line() +
  labs(
    y = "sap flow"
  )

#correlation analysis (moving windows)----
#the x axis is the month of growth data considered
#the y axis is the window of climate data considered
#growthVar must be a column name of the dlySumsSpp dataframe
#climVar must be a column name of the climVars dataframe
#right now I'm considering the spp together. I could swap dlySums for dlySumsSpp to get separate results by spp
#startMonth and endMonth MUST be in format YYYYMM.
movingWindows <- function(startMonth, endMonth, spp, climVar, growthVar){
  sppDlySums <- dlySumsSpp[dlySumsSpp$"species" == spp, ]
  corrDat <- data.frame(date = sppDlySums$date,
                        month = month(sppDlySums$date),
                        growthVar = sppDlySums[[growthVar]],
                        OneDay = NA,
                        ThreeDay = NA,
                        OneWeek = NA,
                        TwoWeek = NA,
                        OneMonth = NA,
                        ThreeMonth = NA,
                        SixMonth = NA,
                        OneYear = NA)
  longMonth <- ifelse(corrDat$month < 10,
                      paste0(0, as.character(corrDat$month)),
                      corrDat$month)
  corrDat$monthYear = as.numeric(paste0(as.character(year(corrDat$date)), longMonth))
  
  if(startMonth < min(corrDat$monthYear)){
    stop("Error: start month too early or wrong format")
  }
  if(endMonth > max(corrDat$monthYear)){
    stop("Error: end date too late or wrong format")
  }
  if(endMonth <= startMonth){
    stop("Error: start date after end date")
  }
  
  corrDat <- corrDat[corrDat$monthYear >= startMonth & corrDat$monthYear <= endMonth, ]
  
  for(row in 1:nrow(corrDat)){
    endDate = corrDat$date[row]
    
    #one day window
    startDate = corrDat$date[row] - 1
    corrDat$OneDay[row] = mean(climVars[[climVar]][climVars$date >= startDate & climVars$date < endDate], na.rm = TRUE)
    
    #three day window
    startDate = corrDat$date[row] - 3
    corrDat$ThreeDay[row] = mean(climVars[[climVar]][climVars$date >= startDate & climVars$date < endDate], na.rm = TRUE)
    
    #one week window
    startDate = corrDat$date[row] - 7
    corrDat$OneWeek[row] = mean(climVars[[climVar]][climVars$date >= startDate & climVars$date < endDate], na.rm = TRUE)
    
    #two week window
    startDate = corrDat$date[row] - 14
    corrDat$TwoWeek[row] = mean(climVars[[climVar]][climVars$date >= startDate & climVars$date < endDate], na.rm = TRUE)
    
    #one month window
    startDate = corrDat$date[row] - 30
    corrDat$OneMonth[row] = mean(climVars[[climVar]][climVars$date >= startDate & climVars$date < endDate], na.rm = TRUE)
    
    #three month (defined as ninety days)
    startDate = corrDat$date[row] - 90
    corrDat$ThreeMonth[row] = mean(climVars[[climVar]][climVars$date >= startDate & climVars$date < endDate], na.rm = TRUE)
    
    #six month (defined as 180 days)
    startDate = corrDat$date[row] - 180
    corrDat$SixMonth[row] = mean(climVars[[climVar]][climVars$date >= startDate & climVars$date < endDate], na.rm = TRUE)
    
    #one year (defined as 365 days)
    startDate = corrDat$date[row] - 365
    corrDat$OneYear[row] = mean(climVars[[climVar]][climVars$date >= startDate & climVars$date < endDate], na.rm = TRUE)
  }
  
  movingCorrs <- expand.grid(
    month = 1:12,
    year = substr(startMonth, 1, 4):substr(endMonth, 1, 4)
  )
  
  longMonth <- ifelse(movingCorrs$month < 10,
                      paste0(0, as.character(movingCorrs$month)),
                      movingCorrs$month)
  movingCorrs$monthYear = as.numeric(paste0(as.character(movingCorrs$year), longMonth))
  movingCorrs <- movingCorrs[movingCorrs$monthYear >= startMonth & movingCorrs$monthYear <= endMonth, ]
  
  safe_cor <- function(x,y){
    complete <- complete.cases(x,y)
    if(sum(complete) < 28){
      return(NA)
    }
    return(cor(x[complete], y[complete]))
  }
  
  for(row in 1:nrow(movingCorrs)){
    monthYear = movingCorrs$monthYear[row]
    corrDatNOW <- corrDat[corrDat$monthYear == monthYear, ]
    
    movingCorrs$OneDay[row] <- safe_cor(corrDatNOW$growthVar, corrDatNOW$OneDay)
    movingCorrs$ThreeDay[row] <- safe_cor(corrDatNOW$growthVar, corrDatNOW$ThreeDay)
    movingCorrs$OneWeek[row] <- safe_cor(corrDatNOW$growthVar, corrDatNOW$OneWeek)
    movingCorrs$TwoWeek[row] <- safe_cor(corrDatNOW$growthVar, corrDatNOW$TwoWeek)
    movingCorrs$OneMonth[row] <- safe_cor(corrDatNOW$growthVar, corrDatNOW$OneMonth)
    movingCorrs$ThreeMonth[row] <- safe_cor(corrDatNOW$growthVar, corrDatNOW$ThreeMonth)
    movingCorrs$SixMonth[row] <- safe_cor(corrDatNOW$growthVar, corrDatNOW$SixMonth)
    movingCorrs$OneYear[row] <- safe_cor(corrDatNOW$growthVar, corrDatNOW$OneYear)
  }
  
  movingCorrsLong <- movingCorrs |> pivot_longer(cols = 4:11,
                                                 names_to = "window_length",
                                                 values_to = "value")
  windowLevels <- c("OneDay", "ThreeDay", "OneWeek", "TwoWeek", "OneMonth", "ThreeMonth", "SixMonth", "OneYear")
  windowLabs = c("1 Day", "3 Days", "1 Week", "2 Weeks", "1 Month", "Three Month", "Six Month", "One Year")
  movingCorrsLong$window_length <- factor(movingCorrsLong$window_length, levels = windowLevels)
  movingCorrsLong$windowType <- c("OneDay" = 1,
                                  "ThreeDay" = 2,
                                  "OneWeek" = 3,
                                  "TwoWeek" = 4,
                                  "OneMonth" = 5,
                                  "ThreeMonth" = 6,
                                  "SixMonth" = 7,
                                  "OneYear" = 8)[movingCorrsLong$window_length]
  movingCorrsLong <- movingCorrsLong[movingCorrsLong$month > 3 & movingCorrsLong$month < 11, ]
  movingCorrsLong$monthYear <- as.factor(movingCorrsLong$monthYear)
  
  startDate <- as.Date(paste0(substr(startMonth, 1, 4), "-", substr(startMonth, 5, 6), "-01"))
  endDate <- as.Date(paste0(substr(endMonth, 1, 4), "-", substr(endMonth, 5, 6), "-01"))
  monthDates <- seq(startDate, endDate, by = "month")
  monthLabs <- format(
    monthDates[lubridate::month(monthDates) %in% 4:9],
    "%b %Y"
  )
  #monthLabs <- format(
  # seq(startDate, endDate, by = "month"),
  #"%b %Y"
  #)
  
  movingCorrsLong <- movingCorrsLong[movingCorrsLong$month %in% 4:9, ]
  
  movingCorrsLong |> ggplot(aes(monthYear, window_length, fill = value)) +
    geom_tile() + 
    scale_x_discrete(labels = monthLabs) +
    geom_text(aes(label = round(value, 2)), size = 6) +
    scale_y_discrete(labels = windowLabs) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(angle = 45, hjust = 1)) + 
    scale_fill_gradient2(low = "orange", mid = "white", high = "green", midpoint = 0) + 
    labs(
      x = "Month (only shows growing season)",
      y = "Window Length of Climate Data",
      fill = "r",
      title = paste("Correlations between monthly",
                    spp,
                    growthVar,
                    "and",
                    climVar,
                    "over various intervals")
    ) + 
    theme_bw()
  
  filepath <- paste0("outputs/",
                     spp,
                     growthVar,
                     climVar,
                     "MovingCorrs.png")
  ggsave(filepath,
         height = 6,
         width = 12,
         dpi = 300)
  
  print(paste("File saved to", filepath))
  
  return(movingCorrsLong)
}

#i haven't bothered writing all of these out yet. There are many more potential calculations here.
movingWindows(202501, 202512, "PIPO", "sapPIPO", "totalGrowth")
movingWindows(202501, 202512, "PIPO", "Tair", "totalGrowth")

#predicting probability of new growth----
#IV must be a variable name in climGrowth
NewGrowthProbModel <- function(IV, onlySummer = TRUE, quadratic = FALSE){
  dat <- climGrowth
  SummerData <- climGrowth[climGrowth$doy.x >= 79 & climGrowth$doy.x <= 312, ]     #what are these dates based on?
  if(onlySummer == TRUE){
    dat <- SummerData 
  }
  
  if(quadratic == TRUE){
    modPIPO <- glm(
      as.formula(paste0("cbind(n_newGrowth, n_trees - n_newGrowth) ~ ", IV, " + I(", IV, "^2)")),
      data = dat[dat$species == "PIPO", ],
      family = binomial
    )
    modPSME <- glm(
      as.formula(paste0("cbind(n_newGrowth, n_trees - n_newGrowth) ~ ", IV, " + I(", IV, "^2)")),
      data = dat[dat$species == "PSME", ],
      family = binomial
    )
  }
  else{
    modPIPO <- glm(
      as.formula(paste0("cbind(n_newGrowth, n_trees - n_newGrowth) ~ ", IV)),
      data = dat[dat$species == "PIPO", ],
      family = binomial
    )
    modPSME <- glm(
      as.formula(paste0("cbind(n_newGrowth, n_trees - n_newGrowth) ~ ", IV)),
      data = dat[dat$species == "PSME", ],
      family = binomial
    )
  }
  
  Pred <- data.frame(
    seq(
      from = min(dat[[IV]], na.rm = TRUE),
      to = max(dat[[IV]], na.rm = TRUE),
      length.out = 200
    )
  )
  names(Pred) <- IV
  
  Pred$probPIPO <- predict(
    modPIPO,
    newdata = Pred,
    type = "response"
  )
  Pred$probPSME <- predict(
    modPSME,
    newdata = Pred,
    type = "response"
  )
  
  xlabs <- c(
    "Hour" = "Time of day (hours)",
    "Tair_f" = "Air temperature (°C)",
    "SWC" = "Soil water content (%, averaged 0-30cm depth)",
    "doyDec" = "Date of year",
    "VPD" = "Vapor pressure deficit (hPa)",
    "sapflow" = "Sap flow velocity (cm/hr)",
    "Rg_f" = "Global radiation (W m^-2)"
  )
  
  Pred |> pivot_longer(cols = c("probPIPO", "probPSME"),
                       names_to = "spp",
                       values_to = "probability") |>
    ggplot(aes(.data[[IV]], probability, color = spp)) +
    geom_line() +
    labs(
      x = xlabs[IV],
      y = ifelse(onlySummer,
                 "Probability of new irreversible growth (DOY 79-312)",
                 "Probability of new irreversible growth (whole year)"
      )
    ) +
    scale_color_discrete(
      labels = c("probPIPO" = "PIPO", "probPSME" = "PSME")
    ) +
    theme_bw()
  
  filepath <- paste0(
    "outputs/GrowthProbability/",
    IV,
    "NewGrowthCurve",
    ifelse(onlySummer, "Summer", "AllYear"),
    ifelse(quadratic, "quadreg", ""),
    ".png"
  )
  
  ggsave(filepath,
         height = 4,
         width = 6,
         dpi = 300)
  
  print(paste0("Graph saved at: ", filepath))
  
  output <- list(
    PIPO = modPIPO,
    PSME = modPSME
  )
  
  return(output)
}

NewGrowthProbModel("VPD")
NewGrowthProbModel("doyDec", FALSE, TRUE)
NewGrowthProbModel("doyDec", TRUE, TRUE)
NewGrowthProbModel("Hour", TRUE, TRUE)
NewGrowthProbModel("Tair_f", TRUE, TRUE)
NewGrowthProbModel("Tair_f")
NewGrowthProbModel("Tair_f", FALSE, TRUE)
NewGrowthProbModel("SWC")
NewGrowthProbModel("SWC", TRUE, TRUE)
NewGrowthProbModel("sapflow")
NewGrowthProbModel("Rg_f")
NewGrowthProbModel("Rg_f", TRUE, TRUE)
NewGrowthProbModel("Rg_f", FALSE, FALSE)
NewGrowthProbModel("Rg_f", FALSE, TRUE)


#growing season calculations----
#as in Javadian et al: GPP season is DOY 72 to 312
flxDlySums$GPP_smoothed <- movingAverage(flxDlySums$GPP, 30)

AverageYearGPP <- summarize(flxDlySums, 
                            GPP_smoothed = mean(GPP_smoothed, na.rm = TRUE),
                            .by = "doy")
amplitude <- max(AverageYearGPP$GPP_smoothed)-min(AverageYearGPP$GPP_smoothed)
halfthreshold <- min(AverageYearGPP$GPP_smoothed + amplitude*0.5)
posthreshold <- min(AverageYearGPP$GPP_smoothed + amplitude*0.8)
AverageYearGPP$season <- ifelse(AverageYearGPP$GPP_smoothed >= halfthreshold, TRUE, FALSE)
AverageYearGPP$POS <- ifelse(AverageYearGPP$GPP_smoothed >= posthreshold, TRUE, FALSE)
AverageYearGPP |> ggplot(aes(doy, GPP_smoothed)) +
  geom_rect(
    aes(xmin=79, xmax=312,
        ymin = -Inf, ymax = Inf),
    fill = "lightgreen",
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(xmin=96, xmax=157,
        ymin = -Inf, ymax = Inf),
    fill = "darkgreen",
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_rect(
    aes(xmin=201, xmax=271,
        ymin = -Inf, ymax = Inf),
    fill = "darkgreen",
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_line()

#as in Etzold et al: 5-95% of annual growth. DOY 139-286
climGrowthPIPO <- climGrowth[climGrowth$species == "PIPO",]
climGrowthPSME <- climGrowth[climGrowth$species == "PSME",]

PIPO5p <- max(climGrowthPIPO$gro_yr, na.rm = TRUE) * 0.05
PIPO95p <- max(climGrowthPIPO$gro_yr, na.rm = TRUE) * 0.95
pipoStart <- climGrowthPIPO$doy.x[which.min(abs(climGrowthPIPO$gro_yr - PIPO5p))]
pipoEnd <- climGrowthPIPO$doy.x[which.min(abs(climGrowthPIPO$gro_yr - PIPO95p))]
PSME5p <- max(climGrowthPSME$gro_yr, na.rm = TRUE) * 0.05
PSME95p <- max(climGrowthPSME$gro_yr, na.rm = TRUE) * 0.95
psmeStart <- climGrowthPSME$doy.x[which.min(abs(climGrowthPSME$gro_yr - PSME5p))]
psmeEnd <- climGrowthPSME$doy.x[which.min(abs(climGrowthPSME$gro_yr - PSME95p))]

#modeling as in Etzold et al----
#note: creating the new datasets is not redundant. This preserves "series" so it can be made into a random effect
df <- left_join(dendro[dendro$year==2025, ],
                climVars[climVars$Year==2025, ],
                by="doyDec"
)

df <- df[,c("date.x","year","doy.x","doyDec","species","series","value","max","twd","gro_yr","growth","newGrowth","isNewGrowth",
            "Ustar","NEE_f","GPP_f","Reco_f","LE_f","H_f","Rg_f","Tair_f","rH","VPD_f","P","PA","SWC","sapPIPO","sapPSME")]
names(df)[c(1,3)] <- c("date","doy")
df <- left_join(df, daylength, by = "doy")
df <- df[df$doy >= 139 & df$doy <=236, ]

df <- df |>
  mutate(
    across(
      c(Tair_f, rH, VPD_f, Rg_f, SWC, daylength),
      ~ as.numeric(scale(.x))
    )
  )

wilcox.test(df[df$species == "PSME",]$growth,
            df[df$species == "PIPO",]$growth) #fail to reject the null hypothesis. they are not proven to be separate populations


lmermod <- lmer(isNewGrowth ~ 1 + (1|species), data = df)
summary(lmermod)
library(lmerTest)
ranova(lmermod)  #the random effect of species is significant

PSMEmodel <- glmer(
  isNewGrowth ~ 
    (Tair_f+rH+VPD_f+Rg_f+SWC)*daylength + 
    VPD_f:SWC + 
    (1|series), 
  data=df[df$doy >= 139 & df$doy <= 236 & df$species == "PSME", ], 
  family=binomial(link="logit")
)

PIPOmodel <- glmer(
  isNewGrowth ~ 
    (Tair_f+rH+VPD_f+Rg_f+SWC)*daylength + 
    VPD_f:SWC + 
    (1|series), 
  data=df[df$doy >= 139 & df$doy <= 236 & df$species == "PIPO", ], 
  family=binomial(link="logit")
)
summary(PSMEmodel)
summary(PIPOmodel)
library(car)

vifPSME <- vif(PSMEmodel)
vifPIPO <- vif(PIPOmodel)
vifPSME
vifPIPO

cor_matrix <- cor(
  PSMEmodel@frame[, c(
    "Tair_f", "rH", "VPD_f",
    "Rg_f", "SWC", "daylength"
  )],
  use = "complete.obs"
)
cor_df <- as.data.frame(as.table(cor_matrix))
colnames(cor_df) <- c("Variable1", "Variable2", "Correlation")
cor_df <- cor_df[
  as.numeric(factor(cor_df$Variable1)) <
    as.numeric(factor(cor_df$Variable2)), ]
cor_df <- cor_df[
  order(abs(cor_df$Correlation), decreasing = TRUE), ]
cor_df
