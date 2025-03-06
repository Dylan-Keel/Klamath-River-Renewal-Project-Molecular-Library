# The following is the script for the main analysis described in the manuscript:
#A molecular specimen bank for contemporary and future study to capture landscape-scale 
#fish biodiversity baselines before Klamath River dam-removal
# 

# Dylan Jon Keel 
# dkeel@res.us
# Resource Environmental Solutions
# 
# Katie Karpenko
# Cramer Fish Sciences
# katie.karpenko@fishsciences.net
# 03/04/2025


#library required packages


# library(doSNOW)
# library(doParallel)
# library(doMPI)
library(tidyverse)
library(vegan)
library(tictoc)
library(fields)
library(ggplot2)
library(sp)
library(sf)


library(sn)
library(lattice)
install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"))
library(INLA)
library(plotly)
library(geometry)
library(viridis)
#library(devtools) #might need devtools to install some packaages
library(tictoc)
library(kableExtra)

library(rgdal)

#sometimes this happens...
#install.packages("gstat",force=TRUE)
library(gstat)
library(remotes)
#library(INLAOutputs)


remotes::install_github("timcdlucas/INLAutils")
remotes::install_github("inbo/inlatools")
remotes::install_github("gfalbery/ggregplot")

#library(inlamesh3d)
library(inlatools)
library(INLAutils)
library(ggregplot)

#getwd()

###############

# Download GitHub Repository and set as working directory
setwd("...My File Path Here... /Downloads/Klamath-River-Renewal-Project-Molecular-Library-main/Klamath-River-Renewal-Project-Molecular-Library-main")

# source scripts for analysis
source("R_code/HighstatLibV13.R")
source("R_code/INLA_plotting_functions.R")

#read in site and sample data
data1=read.csv("data/Klamath_ML_Sample.ID.Metadata.2024.csv")

#read in sequencing data
data2=read.csv("data/Klamath_ML_Raw.Reads.2023.csv")

#process and flow correct sequence data

# subset metadata to include only 2023
covar2023<-subset(data1, Year=="2023")

# subset fish only hits from final metabar df
d_fish <- subset(data2, Type=="Fish")

# clean up names of taxa for plotting
d_fish$FigureName <- ifelse(grepl("spp.", d_fish$ReportName), d_fish$ReportName, d_fish$CommonName)

d_fish <-
  d_fish %>%
  mutate(FigureName = recode(FigureName, 
                             "Klamath Large Scale/Klamath Small Scale, Shortnose Sucker, Lost River Sucker" = "Catostomidae spp.",
                             "Rainbow Trout/steelhead" = "Rainbow Trout"))


d_clean <- d_fish %>%
  select(Site.Name, Sample.ID, ESVsize, FigureName) %>%  
  drop_na(Site.Name)  # Remove NAs (likely NTCs)

#
d_clean_unique <- d_clean %>%
  distinct(Sample.ID, FigureName, .keep_all = TRUE)  

d_wide <- d_clean_unique %>%
  pivot_wider(names_from = "FigureName", values_from = "ESVsize", values_fill = 0) %>%
  as.data.frame()

d_wide <- d_wide %>%
  mutate(Sample.ID = str_replace(Sample.ID, "[AB]$", ""))


# Filter for Replicate Numbers 4, 5, and 6
filtered_covar <- covar2023 %>%
  filter(Replicate.Number %in% c(4, 5, 6))%>%
  select(-Site.Name)


d_merge=left_join(d_wide,filtered_covar,by="Sample.ID")

d_merge=d_merge%>%
  select(1:which(colnames(d_merge) == "Sample.Volume"))  # Select all columns up to Sample.Volume

# Identify species columns (all except first two and last two)
species_columns <- colnames(d_merge)[3:(ncol(d_merge)-2)]


# calculate the volume corrected average read number per site
d_merg_vol <- d_merge %>%
  group_by(Site.Name)%>%
  summarise(across(all_of(species_columns), ~ sum(.x *Sample.Volume)/sum(Sample.Volume)),Total_Volume=sum(Sample.Volume))



# Group by Site.Name and summarize Sample Volume and flow.CMS
filtered_covar <- covar2023 %>%
  filter(Replicate.Number %in% c(4, 5, 6))
covar_summary <- filtered_covar %>%
  group_by(Site.Name) %>%
  summarise(Total_Sample_Volume = sum(Sample.Volume, na.rm = TRUE),
            Total_Flow_CMS = mean(flow.CMS, na.rm = TRUE))




# Merge d_wide with covar_summary to get Sample Volume and Flow per Site
merged_data <- d_merg_vol %>%
  left_join(covar_summary, by = "Site.Name")

# Step 2: Calculate Reads per Liter (RPL)
# Identify species columns (all except first two and last two)
species_columns <- colnames(d_merge)[3:(ncol(d_merge)-3)]

#calculate mean flow adjusted reads (reads/second)
merged_data_RPS <- merged_data %>%
  mutate(across(all_of(species_columns), ~ .x *1000*1000*Total_Flow_CMS))

# Arrange alphabetically and remove unwanted columns
merged_data_RPS_clean <- merged_data_RPS %>%
  select(-Total_Volume,-Total_Sample_Volume, -Total_Flow_CMS) %>%  # Remove the specified columns
  arrange(Site.Name)  # Sort alphabetically by Site.Name


# Define native species
native_species <- c("Catostomidae spp.",
                    "Entosphenus spp.", "Rainbow Trout",
                    "Chinook Salmon", "Coho Salmon", "Cottus spp.", 
                    "Speckled Dace", "Tui Chub")

# Select only the native species columns
merged_data_RPS_native <- merged_data_RPS_clean %>%
  select(Site.Name, all_of(native_species))



flow.data=merged_data_RPS_native


# Process the data
meta_table <- data1 %>%
  select(-c(Sample.ID, ncol(data1)-1, ncol(data1))) %>%  # Remove Sample.ID and last two columns
  filter(Year == 2023) %>%  # Keep only Year 2023
  group_by(Year, Site.Name) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE),  # Summarize numerics with mean
            across(where(is.character), first),  # Retain first value for character columns
            .groups = "drop") %>%  # Drop grouping after summarization
  mutate(
    Reference = ifelse(Site.Character == "Reference, anadromous", "Reference", "Impact"),
    Control = ifelse(Site.Character == "Dam removal, above barrier", "Control", "Anadromous")
  )


#Step 2: Calculate Shannon's Diversity Index and Species Richness for all rows
diversity_data.flow <- flow.data %>%
  rowwise() %>%
  mutate(
    ShannonsDI = diversity(c_across(where(is.numeric)), index = "shannon"),  # Shannon's Diversity Index
    Richness = sum(c_across(where(is.numeric)) > 0)  # Species Richness (non-zero species count)
  ) %>%
  ungroup()  # Ensure row-wise operation doesn't affect later operations


# Rename empty columns in meta_table
names(meta_table)[names(meta_table) == ""] <- paste0("Unnamed_", seq_along(names(meta_table))[names(meta_table) == ""])
meta_table$Unnamed_1=NULL

flow.data <- dplyr::left_join(diversity_data.flow, meta_table, by = "Site.Name")


par(mar=c(2,2,2,2))
# hist(total.data$ShannonsDI,breaks=50)
# hist(total.data$Richness,breaks=10)
# 
# hist(vol.data$ShannonsDI,breaks=50)
# hist(vol.data$Richness,breaks=10)

hist(flow.data$ShannonsDI,breaks=50)
hist(flow.data$Richness,breaks=10)

## check for multicolinearity
colnames(flow.data)

temp_data <- data.frame(lapply(flow.data, function(x) {
  if (is.character(x) || is.factor(x)) {
    as.numeric(factor(x))  # Convert character or factor to numeric based on unique values
  } else {
    x  # Leave numeric columns unchanged
  }
}))
colnames(temp_data)
MyVar.temp<-c("Air.Temperature"   ,     "Water.Temperature"    ,  "Dissolved.Oxygen"   ,   
              "Specific.Conductance" ,  "pH"    ,       # "Q.CMS",              
              "Habitat"        ,        "Size"       ,            "Reference"     ,         "Control"   )


MyVar.cont <- c("Air.Temperature"  ,    "Water.Temperature" ,   "Dissolved.Oxygen"  ,   "Specific.Conductance", "pH"  # , "Total_Volume"    ,           
                #"Q.CMS"   
                )

Mypairs(temp_data[,MyVar.temp])

# procede when finding none


#scale the covariates

flow.data.scale <- dplyr::select(flow.data, all_of(MyVar.cont)) 
names(flow.data.scale) <- paste(MyVar.cont, '.z', sep='')

scale.par = do.call(data.frame, 
                    list(mean = apply(flow.data.scale, 2, mean),
                         sd = apply(flow.data.scale, 2, sd),
                         min = apply(flow.data.scale, 2, min),
                         max = apply(flow.data.scale, 2, max)))

flow.data.scale <- as.data.frame(scale(flow.data.scale))
flow.data <- cbind(flow.data, flow.data.scale)
#total.data <-cbind(total.data,flow.data.scale)
#flow.data<-cbind(flow.data,flow.data.scale)

#now check all the covariates for multicolinearity
colnames(flow.data)
#
MyVar <- c("Air.Temperature.z"    ,  # "Total_Volume.z", #"Q.CMS.z",
           "Water.Temperature.z"    ,  "Dissolved.Oxygen.z"  ,     "Specific.Conductance.z"  ,
           "pH.z"             ,           
           "Habitat"    ,            "Size"            ,       "Reference"       ,      
           "Control"    )


#Mypairs(flow.data[,MyVar])
# plot them against the response variable


# MyMultipanel.ggp2(Z = vol.data, 
#                   varx = MyVar, 
#                   vary = "ShannonsDI", 
#                   ylab = "ShannonsDI",
#                   addSmoother = T,
#                   addRegressionLine = F,
#                   addHorizontalLine = FALSE)
# 
# 
# MyMultipanel.ggp2(Z = vol.data, 
#                   varx = MyVar, 
#                   vary = "Richness", 
#                   ylab = "Richness",
#                   addSmoother = T,
#                   addRegressionLine = F,
#                   addHorizontalLine = FALSE)
# 
# 
# MyMultipanel.ggp2(Z = total.data, 
#                   varx = MyVar, 
#                   vary = "ShannonsDI", 
#                   ylab = "ShannonsDI",
#                   addSmoother = T,
#                   addRegressionLine = F,
#                   addHorizontalLine = FALSE)
# 
# 
# MyMultipanel.ggp2(Z = total.data, 
#                   varx = MyVar, 
#                   vary = "Richness", 
#                   ylab = "Richness",
#                   addSmoother = T,
#                   addRegressionLine = F,
#                   addHorizontalLine = FALSE)


MyMultipanel.ggp2(Z = flow.data, 
                  varx = MyVar, 
                  vary = "ShannonsDI", 
                  ylab = "ShannonsDI",
                  addSmoother = T,
                  addRegressionLine = F,
                  addHorizontalLine = FALSE)


MyMultipanel.ggp2(Z = flow.data, 
                  varx = MyVar, 
                  vary = "Richness", 
                  ylab = "Richness",
                  addSmoother = T,
                  addRegressionLine = F,
                  addHorizontalLine = FALSE)



nonzero_values <- flow.data$ShannonsDI[flow.data$ShannonsDI > 0]

# Plot histogram
hist(flow.data$ShannonsDI, breaks = 10, col = "red", border = "black", 
     main = "Histogram with Zero Values", 
     xlab = "Value", ylab = "Frequency")
hist(nonzero_values, breaks = 10, col = "blue", border = "black", 
     main = "Histogram of Non-Zero Values", 
     xlab = "Value", ylab = "Frequency")




#make sure data is in the right format
####FIX BEFORE RUNNING OTHER SITES####
flow.data$Habitat=as.factor(flow.data$Habitat)

#flow.data$Habitat<-recode_factor(All_Data$YearR, "2021"=2, "2022"=3,"2023"=4,"2024"=5)
#YearR=as.numeric(All_Data$YearR)
#All_Data=subset(All_Data,select=-YearR)

flow.data$Size=as.factor(flow.data$Size)
flow.data$Reference=as.factor(flow.data$Reference)
flow.data$Control=as.factor(flow.data$Control)
flow.data$Stream=gsub("[0-9]", "", flow.data$Site.Name)
flow.data$Stream=gsub("Klamath_Ranch", "Klamath", flow.data$Stream)

# Assign the new specific stream names
flow.data$Stream[flow.data$Site.Name %in% c( "Klamath34", "Klamath36")] <- "JCBoyle"
flow.data$Stream[flow.data$Site.Name %in% c("Klamath12", "Klamath10", "Beaver1")] <- "Copco"
flow.data$Stream[flow.data$Site.Name %in% c("Klamath6", "Klamath4", "Jenny1", "Scotch1", "Klamath2")] <- "IronGate"

flow.data$Stream=as.factor(flow.data$Stream)


flow.data$ShannonsDI <- flow.data$ShannonsDI + 1e-6

ShannonsDI=flow.data$ShannonsDI

# 
# #No quads yet


plot(flow.data$ShannonsDI~flow.data$Habitat)
# NumZones.z=as.numeric(All_Data$NumZones.z)
# Number.of.Days.Inundated.z=as.numeric(All_Data$Number.of.Days.Inundated.z)
# MaxD.cm.z=as.numeric(All_Data$MaxD.cm.z)
# Area.m.z=as.numeric(All_Data$Area.m.z)
# NeighborCount.z=as.numeric(All_Data$NeighborCount.z)
# PerimeterM.z=as.numeric(All_Data$PerimeterM.z)
# PerimAreaRatio.z=as.numeric(All_Data$PerimAreaRatio.z)
# MaxD.cm.z2=as.numeric(All_Data$MaxD.cm.z)*as.numeric(All_Data$MaxD.cm.z)
# Number.of.Days.Inundated.z2=as.numeric(All_Data$Number.of.Days.Inundated.z)*as.numeric(All_Data$Number.of.Days.Inundated.z)


#make base model forms
base.form.pcp <- ShannonsDI ~ Habitat + Size + Control+ Reference +
  Dissolved.Oxygen.z + Specific.Conductance.z+ f(Stream,model="iid")

#pcprior <- list(theta = list(prior = "pc.prec", param = c(1, 0.01)))

control.fixed <- list(
  mean = 0,  # Mean for all fixed effects
  prec = 1e-6  # Very low precision (high variance), effectively flat
)

# Set uninformative prior for the precision of the Gamma distribution
# Inverse-gamma prior for precision
control.family <- list(
  hyper = list(
    prec = list(
      prior = "loggamma",  # Use a log-gamma distribution for precision
      param = c(1, 0.00001),  # Uninformative prior: low shape (1) and high variance (small rate)
      initial = log(1),  # Set initial value for the precision
      fixed = FALSE  # Let INLA estimate the precision
    )
  )
)



######here we fit some models without spatial effects
####normal non-spatial with distance from shore
# Prepare data: Add a column to categorize based on ShannonsDI

# Update Habitat labels

flow.data <- flow.data %>%
  mutate(Habitat_Label = ifelse(as.character(Habitat) == "Reservoir", "Reservoir", "Stream"),
         Category = ifelse(ShannonsDI > 0.05, "Above 0.05", "0.05 or Below"))
#names(flow.data)

# Ensure Stream factor order aligns with Site.Name order
flow.data$Stream <- factor(flow.data$Stream, levels = unique(flow.data$Stream[order(flow.data$Site.Name)]))

# Generate the plot with Site.Name as the x-axis labels
ggplot(flow.data, aes(x = reorder(Site.Name, ShannonsDI), y = ShannonsDI, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Habitat_Label, scales = "free_x") +  # Facet by Habitat_Label (Reservoir/Stream)
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate labels for readability
    strip.text = element_text(size = 12),  # Enhance facet label size
    panel.grid.major.x = element_blank()  # Remove x-axis major grid lines for a cleaner look
  ) +
  labs(
    title = "Shannon's Diversity Index (Native Fishes) by Site and Habitat Type",
    x = "Site Name",
    y = "Shannon's DI"
  ) +
  scale_fill_manual(values = c("Above 0.05" = "#00577b", "0.05 or Below" = "#97a641"))


set.seed(1993)



flow.data <- flow.data %>%
  mutate(across(where(is.factor), ~ as.numeric(as.factor(.))))

# Step 1: Fit the zero-inflation model
zero_inflated <- flow.data  # replace flow.data with your actual data frame
zero_inflated$has_non_zero <- as.numeric(flow.data$ShannonsDI > 0.05)

# Logistic regression for zero-inflation
binary_formula <- has_non_zero ~ Habitat + Size + Control+ Reference +
  Dissolved.Oxygen.z + Specific.Conductance.z+  f(Stream,model="iid")

new_data=zero_inflated

binary_fit <- inla(binary_formula, 
                   family = "binomial", 
                   data = zero_inflated, 
                   control.compute = list(waic = TRUE, cpo = TRUE, config = TRUE,return.marginals.predictor=TRUE),
                   control.predictor = list(link = 1, compute = TRUE))

summary(binary_fit)


# Extract fitted values (mean posterior predictions)
fitted_values <-binary_fit$summary.fitted.values$mean

# Combine fitted values with original data
zero_inflated$fitted <- fitted_values

library(pROC)

# Fit your INLA model (assuming you already have a model)
# model <- inla(has_non_zero ~ Habitat, family = "binomial", data = zero_inflated)

# Extract the fitted probabilities from your INLA model
fitted_probs <- binary_fit$summary.fitted.values$mean

# Calculate the ROC curve and AUC
roc_curve <- roc(zero_inflated$has_non_zero, fitted_probs)

# Display AUC value
auc_value <- auc(roc_curve)
print(paste("AUC:", auc_value))

# Optionally plot the ROC curve
par(mar=c(3,3,3,3))
par(mfrow = c(1,1))
plot(roc_curve, col = "black", lwd = 2, main = "ROC Curve for INLA Binomial Model")
abline(a = 0, b = 1, col = "grey", lty = 2)  # Add a diagonal reference line



# Step 1: Generate posterior predictive replications
n_replications <- 10000  # Number of replicated datasets
predictions <- inla.posterior.sample(n_replications, binary_fit)

# Step 2: Compute discrepancy measure (e.g., mean residuals for simplicity)
obs_discrepancy <- abs(zero_inflated$has_non_zero - binary_fit$summary.fitted.values$mean)

rep_discrepancy <- sapply(predictions, function(rep) {
  replicated_data <- sapply(rep$latent, mean)  # Average over replicated values
  abs(zero_inflated$has_non_zero - replicated_data)  # Calculate discrepancy
})

# Step 3: Calculate Bayesian p-value
bayesian_p_value <- mean(rep_discrepancy > obs_discrepancy)

bayesian_p_value


summary(binary_fit)


exp(7.585)
exp(31.511)


1- exp(-7.198)
1 -exp(-.873)

plot(flow.data$ShannonsDI~flow.data$Size)
######

# After accounting for the random effect of waterbody...

# look for spatial autocorrelation
#average distance between sites


# check out the model
par(mar=c(3,3,3,3))

observed<-zero_inflated$has_non_zero

D<-INLAutils::plot_inla_residuals(binary_fit, observed)



# Let's make a simplified world to view our data in 2D to check for
# spatial-autocorrelation. Let's assume that the river is a straight
# line and that sinuosity is not important.


#average distance between sites

xy <- with(flow.data, data.frame(Site.Name, X.Coordinate, Y.Coordinate))

# Convert the lat/long dataframe into a spatial object
latlong_sf <- st_as_sf(xy, coords = c("X.Coordinate", "Y.Coordinate"), crs = 4326)  # EPSG:4326 is WGS84

# Transform the lat/long coordinates to UTM
# sf automatically selects the appropriate UTM zone based on your data location
utm_sf <- st_transform(latlong_sf, crs = 32610)  # Example: UTM zone 10N, adjust accordingly

# View the results
#print(utm_sf)

# Extract the UTM coordinates
utm_coords <- st_coordinates(utm_sf)

# Replace the lat/long in the xy object with UTM coordinates
xy$UTM_Easting <- utm_coords[,1]  # UTM Easting (X coordinate)
xy$UTM_Northing <- utm_coords[,2] # UTM Northing (Y coordinate)

# Remove old Latitude and Longitude columns if you want
xy$Latitude <- NULL
xy$Longitude <- NULL

# View the updated dataframe
#print(xy)


#nonspat w/out velocity adj
Pi1 <- binary_fit$summary.fitted.values[, "mean"]



D1  <- (observed - Pi1) / sqrt(Pi1)

summary(D1)

#D1
# e <- Cs.Copies-Pi1

# dev.res <- sign(e)*sqrt(-2*(Cs.Copies*log(Pi1) + (1 - Cs.Copies)*log(1 - Pi1)))
# #D1
#Pi1
### shrink x 2 orders of magnitude
MyData <- data.frame(D1  = as.numeric(D1),
                     X = as.numeric(xy$UTM_Easting),
                     Y = as.numeric(xy$UTM_Northing))



###### here we look for spatial patterns in the residuals
MyData$MySize <- 2 * abs(MyData$D1) / max(MyData$D1)
MyData$MyCol <- ifelse(MyData$D1> 0, 1, 2)
lattice::xyplot(Y ~ X,
                data = MyData,
                cex = MyData$MySize,
                col = MyData$MyCol,
                pch = 1)

## we can see that sites 3 and 6 have clear spatial patterns
# in residuals
#View(MyData)


V1a=NULL

## now we make a semivariogram

coordinates(MyData) <- c("X","Y")
V1a <- gstat::variogram(D1 ~ X + Y, 
                        data = MyData, 
                        cutoff = 6000,
                        cressie = TRUE)

hist(V1a$dist)

p <- ggplot()
p <- p + geom_point(data = V1a,
                    aes(x = dist,
                        y = gamma))
p <- p + geom_smooth(data = V1a,
                     span = 0.9,
                     se = FALSE,
                     aes(x = dist,
                         y = gamma))
p <- p + xlab("Distance(m)") + ylab("Sample variogram")
p1 <- p + theme(text = element_text(size=15))+ylim(0,.25)
p1 

#with a sill of ~0.625 and a nugget of ~.25 we find:

Sill.Nug.R=.075/.125
1-Sill.Nug.R

# Your plot
p <- ggplot() +
  geom_point(data = V1a, aes(x = dist, y = gamma)) +
  geom_smooth(data = V1a, span = 0.9, se = FALSE, aes(x = dist, y = gamma)) +
  xlab("Distance (m)") +
  ylab("Sample variogram") +
  theme(text = element_text(size = 15)) +
  ylim(0, .25) +
  # Add horizontal line for the Sill
  geom_hline(yintercept = 0.125, color = "green", linetype = "dashed") +
  annotate("text", x = max(V1a$dist, na.rm = TRUE) * 0.9, y = 0.125, label = "Sill", color = "green", hjust = 0) +
  # Add horizontal line for the Nugget
  geom_hline(yintercept = 0.075, color = "orange", linetype = "dashed") +
  annotate("text", x = max(V1a$dist, na.rm = TRUE) * 0.9, y = 0.075, label = "Nugget", color = "orange", hjust = 0) +
  # Add vertical line for the Range
  geom_vline(xintercept = 2500, color = "red", linetype = "dashed") +
  annotate("text", x = 2500, y = max(V1a$gamma, na.rm = TRUE) * 0.95, label = "Range", color = "red", angle = 90, hjust = -0.5)

# Display the plot
p


#we see that spatial autocorrelation increases from 0 to 2500 meters
# in our data

Loc<-NULL
Loc <- cbind(xy$UTM_Easting,xy$UTM_Northing)
#Loc
#what are the distances between the points?
D <- dist(Loc)
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(D, 
     freq = TRUE,
     main = "", 
     xlab = "Distance between samples (m)",
     ylab = "Frequency")
#View(Loc)


RangeGuess <-2500


# There is still a spatial correlation in the residuals of spat.inla.  The first thing that you should try is to
# use a mesh with more triangles. Use a smaller MaxEdge and also a smaller
# cutoff. By doing that we allow for smaller-scale correlation.
# Right now we have:

require(splancs)



sum(is.na(Loc))
Loc <- na.omit(Loc)

Hull <- inla.nonconvex.hull(Loc, convex = -0.2)
MaxEdge  <- RangeGuess/10
mesh2d     <- inla.mesh.2d(boundary = Hull,
                           max.edge = c(1,40) * MaxEdge,
                           cutoff = MaxEdge/20 ,
                           offset = c(5, 100) ,
                           max.n=20000)  #control max.n to be <3000
# 
mesh2d$n
# 

# #back to 2D
par(mfrow = c(1,1), mar=c(0,0,0,0))
plot(mesh2d, asp=1, main = "")
points(Loc, col = 2, pch = 1, cex = 1.5)



####fix na values in covariates for spatial model
zero_inflated$Specific.Conductance.z[is.na(zero_inflated$Specific.Conductance.z)] <- mean(zero_inflated$Specific.Conductance.z, na.rm = TRUE)



# Define projector matrices for the mesh.
A <- inla.spde.make.A(mesh2d, loc = Loc)


spde <- inla.spde2.pcmatern(mesh2d, 
                            prior.range = c(RangeGuess, 0.5), 
                            prior.sigma = c(1, .05))

w.index <- inla.spde.make.index('w', n.spde = spde$n.spde)
#w.index


stackform<-as.formula(binary_formula)


modrando<-c("Stream")

terms <- base.form.pcp[3]
terms
terms <- gsub(" ", "", as.character(terms))
terms <- unlist(strsplit(terms, '\\+'))
terms <- terms[-grep('iid', terms)]
if(length(grep(':', terms))>0) terms <- terms[-grep(':', terms)] #apparently we can only use the main effects
terms <- c(terms, modrando)
terms

## make the formula for the stack
stackform <- formula(paste('~', paste(terms, collapse='+')))
stackform

Xm <- model.matrix(stackform, data = zero_inflated,na.action=na.pass)


Xm <- as.data.frame(Xm)
N <- nrow(zero_inflated)
#w.index


head(zero_inflated)


has_non_zero=zero_inflated$has_non_zero

StackFit <- inla.stack(
  remove.unused = T,
  tag = "Fit",
  data = list(
    has_non_zero = has_non_zero
  ),
  A = list(1, 1, A),
  effects = list(
    Intercept = rep(1, N),
    Xm = Xm[,-1],  # Use covariates without the intercept (includes dummy variables)
    w = w.index    # Ensure w.index is correctly aligned with random effects
  )
)
# 
# StackFit <- inla.stack(
#   remove.unused = TRUE,
#   tag = "Fit",
#   data = list(Shannon = vol.data$ShannonsDI),  # Ensure this variable exists in vol.data
#   A = list(1,1, A),  # Double-check that A is correctly sized for vol.data
#   effects = list(
#     Intercept = rep(1, nrow(vol.data)),
#     Xm = model.matrix(~ Habitat + Size + Control + Reference + Dissolved.Oxygen.z + Specific.Conductance.z - 1, data = vol.data), 
#     w = w.index
#   )
# )
# 
# # Use only relevant terms from Xm and provide the response separately
# StackFit <- inla.stack(
#  # remove.unused = TRUE,
#   tag = "Fit",
#   data = list(Shannon = ShannonDI),  # Response variable as separate vector
#   A = list(1, 1, A),  # A matrix for projections (make sure this aligns with w.index)
#   effects = list(
#     Intercept = rep(1, nrow(Xm)),  # Intercept term
#     Xm = Xm[,-1],  # Use the covariates from Xm without the intercept
#     w = w.index    # Ensure w.index aligns correctly
#   )
# )

#inla.stack.data(StackFit)

#Fit the global model
# 
# base.form.pcp <- ShannonsDI ~ Habitat + Size + Control+ Reference +
#   Dissolved.Oxygen.z + Specific.Conductance.z+  f(Stream,model="iid")

base.form.pcp <- has_non_zero ~ Habitat + Size + Control+ Reference +
  Dissolved.Oxygen.z + Specific.Conductance.z +#+  f(Stream,model="iid")
  f(w, model = spde)#+f(Site,model="iid",hyper=pcprior)



tic()
spat.inla.diversity<- inla(formula = update(base.form.pcp, . ~ . +
                                            f(Stream,model="iid")),
                         family = "binomial",
                         data = inla.stack.data(StackFit),
                         control.compute = list(waic = TRUE, config=TRUE,cpo=TRUE, po=TRUE,return.marginals.predictor=TRUE),
                         # control.family = control.family,
                         # control.inla = list(strategy = "gaussian"), 
                         control.predictor = list(link = 1, compute = TRUE,
                                                  A = inla.stack.A(StackFit)))

toc()


summary(spat.inla.diversity)



# Extract fitted values (mean posterior predictions)
fitted_values <-spat.inla.diversity$summary.fitted.values$mean[1:42]

# Combine fitted values with original data
zero_inflated$fitted <- fitted_values

library(pROC)

# Fit your INLA model (assuming you already have a model)
# model <- inla(has_non_zero ~ Habitat, family = "binomial", data = zero_inflated)

# Extract the fitted probabilities from your INLA model
fitted_probs <- spat.inla.diversity$summary.fitted.values$mean[1:42]

# Calculate the ROC curve and AUC
roc_curve <- roc(zero_inflated$has_non_zero, fitted_probs)

# Display AUC value
auc_value <- auc(roc_curve)
print(paste("AUC:", auc_value))

# Optionally plot the ROC curve
par(mar=c(3,3,3,3))
par(mfrow = c(1,1))
plot(roc_curve, col = "black", lwd = 2, main = "ROC Curve for INLA Binomial Model")
abline(a = 0, b = 1, col = "grey", lty = 2)  # Add a diagonal reference line



# Step 1: Generate posterior predictive replications
n_replications <- 1000  # Number of replicated datasets
predictions <- inla.posterior.sample(n_replications, spat.inla.diversity)

# Step 2: Compute discrepancy measure (e.g., mean residuals for simplicity)
obs_discrepancy <- abs(zero_inflated$has_non_zero - spat.inla.diversity$summary.fitted.values$mean)

rep_discrepancy <- sapply(predictions, function(rep) {
  replicated_data <- sapply(rep$latent, mean)  # Average over replicated values
  abs(zero_inflated$has_non_zero - replicated_data)  # Calculate discrepancy
})

# Step 3: Calculate Bayesian p-value
bayesian_p_value <- mean(rep_discrepancy > obs_discrepancy)

bayesian_p_value


summary(spat.inla.diversity)

#habitat
exp(7.706)
exp(40.779)

#size
1- exp(-35.4)
1 -exp(-3.727)


#reference
1- exp(-11.076)
1 -exp(-0.547)


# After accounting for the random effect of waterbody and spatial position...

# The log of the odds that streams have non-low native fish diversity is 
# ~8-41x greater than the log of the odds that reservoirs have non-low
# fish diversity that is extremely low (Shannon's Diversity Index < 0.05).

# The odds that mainstem sites have non-low native fish diversity is a near certainty
# to be lower than the odds that tributaries have non-low
# fish diversity that is extremely low (Shannon's Diversity Index < 0.05).

# Although the model finds that the odds that eference sites have 
# non-low native fish diversity that is a near certainty
# to be lower than the odds that impact sites have non-low
# fish diversity that is extremely low (Shannon's Diversity Index < 0.05), 
# low sample size of both makes these parameter estimates difficult to interpret.

# Our model fit the data well, with an AUC of 0.96 indicating a 96% true positive rate. 
# Additionally, our bayesian p-value was between 0.1 and 0.9 (0.62), indicating 
# goodness of fit.

# Because our posterior probability distributions spanned zero for Control/Impact,
#  Dissolved Oxygen Concentration, and Specific Conductance, there
# is no evidence that these variables described a meaningful amount of variance in
#native fish species richness.

summary(flow.data$ShannonsDI)













#### Now we fit the richness models in the same way
# Logistic regression for zero-inflation
base.form.rich <- Richness ~ Habitat + Size + Control+ Reference +
  Dissolved.Oxygen.z + Specific.Conductance.z+  f(Stream,model="iid")


nonspat.rich <- inla(base.form.rich,
                     
                     control.compute = list(waic=TRUE, cpo=TRUE, po=TRUE,config=TRUE,return.marginals.predictor=TRUE),
                     family = "poisson", 
                     
                     control.predictor = list(link = 1, compute = TRUE),
                     data = flow.data)
summary(nonspat.rich)



nonspat.d2.nbn <- inla(base.form.rich,
                       
                       control.compute = list(waic=TRUE, cpo=TRUE, po=TRUE,config=TRUE),
                       family = "nbinomial", 
                       
                       control.predictor = list(link = 1, compute = TRUE),
                       data = flow.data)
summary(nonspat.d2.nbn)

#we'll check to see if zero inflation is warranted

nonspat.d2.zip <- inla(base.form.rich,
                       
                       control.compute = list(waic=TRUE, cpo=TRUE, po=TRUE,config=TRUE),
                       family = "zeroinflatedpoisson1", 
                       
                       control.predictor = list(link = 1, compute = TRUE),
                       data = flow.data)
summary(nonspat.d2.zip)



nonspat.d2.zinbn <- inla(base.form.rich,
                         
                         control.compute = list(waic=TRUE, cpo=TRUE, po=TRUE,config=TRUE),
                         family = "zeroinflatednbinomial1", 
                         
                         control.predictor = list(link = 1, compute = TRUE),
                         data = flow.data)
summary(nonspat.d2.zinbn)


###########################Check for Overdispersion

#source("Modified Distribution Check inlatools.R")
#source("Modified Dispersion Check inlatools.R")



nonspat.d2_dc <- dispersion_check(nonspat.rich)
nonspat.d2.nbn_dc <- dispersion_check(nonspat.d2.nbn)
nonspat.d2.zip_dc <- dispersion_check(nonspat.d2.zip)
nonspat.d2.zinbn_dc <- dispersion_check(nonspat.d2.zinbn)



# Here's what we want to test our models for overdisperssion
plot(nonspat.d2_dc)

mean(nonspat.d2_dc$data)
mean(nonspat.d2_dc$model)
#the poisson model is clearly underdispersed we want the 
#  0.1 < P(D|data>D|model) < 0.9

mean(nonspat.d2.nbn_dc$model)

# that's low, but underdispersion can arise from an autocorrelation
# structure. Let's not rule out the NB distribution until we check
# for spatial autocorrelation


# check out the model
par(mar=c(3,3,3,3))

observed<-flow.data$Richness

D<-INLAutils::plot_inla_residuals(nonspat.rich, observed)


#P<-autoplot(nonspat.d2.zinbn,which=(c(1:5)))
#P

#summary(nonspat.d2.zinbn)
summary(nonspat.rich)


# Let's make a simplified world to view our data in 2D to check for
# spatial-autocorrelation. Let's assume that the river is a straight
# line and that sinuosity is not important.


#average distance between sites

xy <- with(flow.data, data.frame(Site.Name, X.Coordinate, Y.Coordinate))

# Convert the lat/long dataframe into a spatial object
latlong_sf <- st_as_sf(xy, coords = c("X.Coordinate", "Y.Coordinate"), crs = 4326)  # EPSG:4326 is WGS84

# Transform the lat/long coordinates to UTM
# sf automatically selects the appropriate UTM zone based on your data location
utm_sf <- st_transform(latlong_sf, crs = 32610)  # Example: UTM zone 10N, adjust accordingly

# View the results
#print(utm_sf)

# Extract the UTM coordinates
utm_coords <- st_coordinates(utm_sf)

# Replace the lat/long in the xy object with UTM coordinates
xy$UTM_Easting <- utm_coords[,1]  # UTM Easting (X coordinate)
xy$UTM_Northing <- utm_coords[,2] # UTM Northing (Y coordinate)

# Remove old Latitude and Longitude columns if you want
xy$Latitude <- NULL
xy$Longitude <- NULL

# View the updated dataframe
#print(xy)


#nonspat w/out velocity adj
Pi1 <- nonspat.rich$summary.fitted.values[, "mean"]



D1  <- (observed - Pi1) / sqrt(Pi1)

summary(D1)

#D1
# e <- Cs.Copies-Pi1

# dev.res <- sign(e)*sqrt(-2*(Cs.Copies*log(Pi1) + (1 - Cs.Copies)*log(1 - Pi1)))
# #D1
#Pi1
### shrink x 2 orders of magnitude
MyData <- data.frame(D1  = as.numeric(D1),
                     X = as.numeric(xy$UTM_Easting),
                     Y = as.numeric(xy$UTM_Northing))



###### here we look for spatial patterns in the residuals
MyData$MySize <- 2 * abs(MyData$D1) / max(MyData$D1)
MyData$MyCol <- ifelse(MyData$D1> 0, 1, 2)
lattice::xyplot(Y ~ X,
                data = MyData,
                cex = MyData$MySize,
                col = MyData$MyCol,
                pch = 1)

## we can see that sites 3 and 6 have clear spatial patterns
# in residuals
#View(MyData)


V1a=NULL

## now we make a semivariogram

coordinates(MyData) <- c("X","Y")
V1a <- gstat::variogram(D1 ~ X + Y, 
                        data = MyData, 
                        cutoff = 6000,
                        cressie = TRUE)

hist(V1a$dist)

p <- ggplot()
p <- p + geom_point(data = V1a,
                    aes(x = dist,
                        y = gamma))
p <- p + geom_smooth(data = V1a,
                     span = 0.9,
                     se = FALSE,
                     aes(x = dist,
                         y = gamma))
p <- p + xlab("Distance(m)") + ylab("Sample variogram")
p1 <- p + theme(text = element_text(size=15))+ylim(0,1.2)
p1 

#with a sill of ~0.625 and a nugget of ~.25 we find:

Sill.Nug.R=.25/.625
1-Sill.Nug.R

# Your plot
p <- ggplot() +
  geom_point(data = V1a, aes(x = dist, y = gamma)) +
  geom_smooth(data = V1a, span = 0.9, se = FALSE, aes(x = dist, y = gamma)) +
  xlab("Distance (m)") +
  ylab("Sample variogram") +
  theme(text = element_text(size = 15)) +
  ylim(0, 1.2) +
  # Add horizontal line for the Sill
  geom_hline(yintercept = 0.625, color = "green", linetype = "dashed") +
  annotate("text", x = max(V1a$dist, na.rm = TRUE) * 0.9, y = 0.625, label = "Sill", color = "green", hjust = 0) +
  # Add horizontal line for the Nugget
  geom_hline(yintercept = 0.25, color = "orange", linetype = "dashed") +
  annotate("text", x = max(V1a$dist, na.rm = TRUE) * 0.9, y = 0.25, label = "Nugget", color = "orange", hjust = 0) +
  # Add vertical line for the Range
  geom_vline(xintercept = 2500, color = "red", linetype = "dashed") +
  annotate("text", x = 2500, y = max(V1a$gamma, na.rm = TRUE) * 0.95, label = "Range", color = "red", angle = 90, hjust = -0.5)

# Display the plot
p


#we see that spatial autocorrelation increases from 0 to 2500 meters
# in our data

Loc<-NULL
Loc <- cbind(xy$UTM_Easting,xy$UTM_Northing)
#Loc
#what are the distances between the points?
D <- dist(Loc)
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(D, 
     freq = TRUE,
     main = "", 
     xlab = "Distance between samples (m)",
     ylab = "Frequency")
#View(Loc)


RangeGuess <-2500


# There is still a spatial correlation in the residuals of spat.inla.  The first thing that you should try is to
# use a mesh with more triangles. Use a smaller MaxEdge and also a smaller
# cutoff. By doing that we allow for smaller-scale correlation.
# Right now we have:

require(splancs)



sum(is.na(Loc))
Loc <- na.omit(Loc)

Hull <- inla.nonconvex.hull(Loc, convex = -0.2)
MaxEdge  <- RangeGuess/10
mesh2d     <- inla.mesh.2d(boundary = Hull,
                           max.edge = c(1,40) * MaxEdge,
                           cutoff = MaxEdge/20 ,
                           offset = c(5, 100) ,
                           max.n=20000)  #control max.n to be <3000
# 
mesh2d$n
# 

# #back to 2D
par(mfrow = c(1,1), mar=c(0,0,0,0))
plot(mesh2d, asp=1, main = "")
points(Loc, col = 2, pch = 1, cex = 1.5)



####fix na values in covariates for spatial model
flow.data$Specific.Conductance.z[is.na(flow.data$Specific.Conductance.z)] <- mean(flow.data$Specific.Conductance.z, na.rm = TRUE)



# Define projector matrices for the mesh.
A <- inla.spde.make.A(mesh2d, loc = Loc)


spde <- inla.spde2.pcmatern(mesh2d, 
                            prior.range = c(RangeGuess, 0.5), 
                            prior.sigma = c(1, .05))

w.index <- inla.spde.make.index('w', n.spde = spde$n.spde)
#w.index

base.form.pcp <- Richness ~ Habitat + Size + Control+ Reference +
  Dissolved.Oxygen.z + Specific.Conductance.z+ f(Stream,model="iid")

stackform<-as.formula(base.form.pcp)


modrando<-c("Stream")

terms <- base.form.pcp[3]
terms
terms <- gsub(" ", "", as.character(terms))
terms <- unlist(strsplit(terms, '\\+'))
terms <- terms[-grep('iid', terms)]
if(length(grep(':', terms))>0) terms <- terms[-grep(':', terms)] #apparently we can only use the main effects
terms <- c(terms, modrando)
terms

## make the formula for the stack
stackform <- formula(paste('~', paste(terms, collapse='+')))
stackform

Xm <- model.matrix(stackform, data = flow.data,na.action=na.pass)


Xm <- as.data.frame(Xm)
N <- nrow(flow.data)
#w.index

Richness=flow.data$Richness

StackFit <- inla.stack(
  remove.unused = T,
  tag = "Fit",
  data = list(
    Richness = Richness
  ),
  A = list(1, 1, A),
  effects = list(
    Intercept = rep(1, N),
    Xm = Xm[,-1],  # Use covariates without the intercept (includes dummy variables)
    w = w.index    # Ensure w.index is correctly aligned with random effects
  )
)
# 
# StackFit <- inla.stack(
#   remove.unused = TRUE,
#   tag = "Fit",
#   data = list(Shannon = vol.data$ShannonsDI),  # Ensure this variable exists in vol.data
#   A = list(1,1, A),  # Double-check that A is correctly sized for vol.data
#   effects = list(
#     Intercept = rep(1, nrow(vol.data)),
#     Xm = model.matrix(~ Habitat + Size + Control + Reference + Dissolved.Oxygen.z + Specific.Conductance.z - 1, data = vol.data), 
#     w = w.index
#   )
# )
# 
# # Use only relevant terms from Xm and provide the response separately
# StackFit <- inla.stack(
#  # remove.unused = TRUE,
#   tag = "Fit",
#   data = list(Shannon = ShannonDI),  # Response variable as separate vector
#   A = list(1, 1, A),  # A matrix for projections (make sure this aligns with w.index)
#   effects = list(
#     Intercept = rep(1, nrow(Xm)),  # Intercept term
#     Xm = Xm[,-1],  # Use the covariates from Xm without the intercept
#     w = w.index    # Ensure w.index aligns correctly
#   )
# )

#inla.stack.data(StackFit)

#Fit the global model
# 
# base.form.pcp <- ShannonsDI ~ Habitat + Size + Control+ Reference +
#   Dissolved.Oxygen.z + Specific.Conductance.z+  f(Stream,model="iid")

base.form.pcp <- Richness ~ Habitat + Size + Control+ Reference +
  Dissolved.Oxygen.z + Specific.Conductance.z +#+  f(Stream,model="iid")
  f(w, model = spde)#+f(Site,model="iid",hyper=pcprior)



tic()
spat.inla.shannon<- inla(formula = update(base.form.pcp, . ~ . +
                                            f(Stream,model="iid")),
                         family = "poisson",
                         data = inla.stack.data(StackFit),
                         control.compute = list(waic = TRUE, config=TRUE,cpo=TRUE, po=TRUE,return.marginals.predictor=TRUE),
                         # control.family = control.family,
                         # control.inla = list(strategy = "gaussian"), 
                         control.predictor = list(link = 1, compute = TRUE,
                                                  A = inla.stack.A(StackFit)))

toc()

spat.inla.rich=spat.inla.shannon
summary(spat.inla.rich)



autoplot(spat.inla.rich)

#INLAutils::plot_inla_residuals(spat.inla.d2.d.m, observed)

par(mar=c(4,4,4,4))
D<-INLAutils::plot_inla_residuals(spat.inla.rich, Richness)



########Seems we still have trouble with the lowest values

SpatField.w <- inla.spde2.result(inla = spat.inla.rich,
                                 name = "w",
                                 spde = spde,
                                 do.transfer = TRUE)
Kappa <- inla.emarginal(function(x) x,
                        SpatField.w$marginals.kappa[[1]] )

Sigma_u <- inla.emarginal(function(x) sqrt(x),
                          SpatField.w$marginals.variance.nominal[[1]] )

Range <- inla.emarginal(function(x) x,
                        SpatField.w$marginals.range.nominal[[1]] )


result = inla.spde2.result(spat.inla.rich, "w", spde)

Range

par(mar=c(4,4,4,4))
plot(result[["marginals.range.nominal"]][[1]], type = "l",
     main = "Nominal range, posterior density")


LocMesh <- mesh2d$loc[,1:2]
D<-dist(LocMesh)
#D <- as.matrix(D)
#D
# Using the estimated parameters from the model (see above)
# we can calculate the imposed Matern correlation values.
d.vec <- seq(0, #max(D)
             144392, length = 10000)
Cor.M <- (Kappa * d.vec) * besselK(Kappa * d.vec, 1)
Cor.M[1] <- 1

Sigma_u
Range
max(D)

cor.plot.data=data.frame(d.vec,Cor.M)
cor.plot.data

cor.plot=ggplot()+geom_line(data=cor.plot.data,aes(x=d.vec,y=Cor.M))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlim(0,5000)+
  ylim(0,1)+
  geom_abline(aes(intercept = 0.05, slope = 0),linetype=3)+
  xlab("Distance (m)")+
  ylab("Matern Correlation Values")+
  theme(axis.title = element_text(face="bold"))+
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))+
  geom_text(x = 1000, y = 0.1, aes(label = "5% autocorrelation"))

cor.plot
# 
# par(mfrow=c(1,1), mar = c(5,5,2,2))
# plot(x = d.vec,
#      y = Cor.M,
#      pch = 16,
#      type = "l",
#      cex.lab = 1.5,
#      xlab = "Distance",
#      ylab = "Correlation",
#      xlim = c(0, 10000))
# abline(h = 0.1, lty = 2)


cor.plot.data
Kappa

fitIndex <- inla.stack.index(StackFit, tag='Fit')$data


fitted.b<- spat.inla.rich$summary.fitted.values[fitIndex,]

#fitted.nonspat<-Chin.inla.null$summary.fitted.values[fitIndex,]



Pred   <- fitted.b[, "mean"]
VarY <- Pred
resid   <- (Richness - Pred) / sqrt(VarY)

MyData3 <- data.frame(resid = resid,
                      Xkm = utm_coords[,1],
                      Ykm = utm_coords[,2])


INLAutils::plot_inla_residuals(spat.inla.rich, observed)

# let's look at the residuals again
resid.plot=ggplot_inla_residuals(spat.inla.rich,Richness,CI = TRUE,binwidth = 0.1)

ggplot_inla_residuals2(spat.inla.rich,Richness, CI=TRUE,method = NA)


MyData3$MySize <- 2 * abs(MyData3$resid) / max(MyData3$resid)
MyData3$MyCol <- ifelse(MyData3$resid> 0, 1, 2)

#View(MyData3)
lattice::xyplot(MyData3$Ykm ~ MyData3$Xkm,
                data = MyData3,
                cex = MyData3$MySize,
                col = MyData3$MyCol,
                pch = 1)
par(mfrow=c(1,1))

hist(MyData3$resid, breaks = 20)

summary(spat.inla.rich)
#After accounting for the random effects of Site and Space 
#for no tested effects were included in our best model.


#fitted.d=nospat.inla.t.m$summary.fitted.values[fitIndex,]
rmse <- sqrt(mean((Richness - fitted.b$mean)^2))
rmse

#



fitIndex <- inla.stack.index(StackFit, tag='Fit')$data
fitted <- spat.inla.rich$summary.fitted.values[fitIndex,]

Pi1=fitted$mean

# Pi1
# 
# observed
# 

D1  <- (Richness - Pi1) / sqrt(Pi1)

summary(D1)



fitted.d<- spat.inla.rich$summary.fitted.values
RRtot=sum((observed-mean(observed))^2)
RRes=sum((observed-fitted.b$mean)^2)
pseudo_r2_val.d=1-RRes/RRtot
pseudo_r2_val.d
txt="R^{2} == 0.705"


CSD2=data.frame(flow.data,fitted.b)
# 49%% of the variance is accounted for by random effects
max(CSD2$mean)

plot=ggplot(CSD2,aes(x=mean,y=Richness))+
  geom_point()+
  theme_bw()+
  xlim(c(0,8))+
  ylim(c(0,8))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  
  geom_abline(aes(intercept = 0, slope = 1))+
  xlab("Fitted native fish species richness")+
  ylab("Observed native fish species richness")+
  theme(axis.title = element_text(face="bold"))

plot=plot+geom_text(x =2, y = 6, aes(label = txt),parse=T)+
  theme(axis.text = element_text(size = 16),axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 20))
plot



# Generate posterior predictive replications
n_replications <- 1000  # Number of replicated datasets
predictions <- inla.posterior.sample(n_replications, spat.inla.rich)

# Calculate observed discrepancies
obs_discrepancy <- abs(flow.data$Richness - spat.inla.rich$summary.fitted.values$mean)

# Calculate discrepancies for replicated data
rep_discrepancies <- sapply(predictions, function(rep) {
  replicated_data <- rep$latent[flow.data$Richness]  # Use appropriate index for observed data points
  abs(flow.data$Richness - replicated_data)
})

# Calculate Bayesian p-value
bayesian_p_value <- mean(rowMeans(rep_discrepancies) > obs_discrepancy)
print(bayesian_p_value)

# 
# # Create the plot (for a subset of replications, if necessary)
# discrepancy_data <- data.frame(
#   Observed = obs_discrepancy[1:42],
#   Replicated = rep_discrepancies[, sample(1:n_replications, length(obs_discrepancy))]
# )
# 
# ggplot(discrepancy_data, aes(x = Observed, y = Replicated.205)) +
#   geom_point(color = 'blue') +
#   geom_abline(intercept = 0, slope = 1, color = 'red', linetype = "dashed") +
#   labs(
#     title = "Bayesian p-value Comparison",
#     x = "Observed Discrepancy",
#     y = "Replicated Discrepancy"
#   ) +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 


summary(spat.inla.rich)



exp(0.828)

exp(2.368)


1-exp(-0.828)
1-exp(-0.149)



# Generate the plot with Site.Name as the x-axis labels
ggplot(flow.data, aes(x = reorder(Site.Name, Richness), y = Richness, fill = Habitat_Label)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Habitat_Label, scales = "free_x") +  # Facet by Habitat_Label (Reservoir/Stream)
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate labels for readability
    strip.text = element_text(size = 12),  # Enhance facet label size
    panel.grid.major.x = element_blank()  # Remove x-axis major grid lines for a cleaner look
  ) +
  labs(
    title = "Species Richness by Site and Habitat Type",
    x = "Site Name",
    y = "Richness"
  ) +
  scale_fill_manual(values = c("#00577b", "#97a641"))  # Adjust colors as needed





# After accounting for the random effects of waterbody and spatial position....

# Streams are expected to have 2.3 - 7.9 times greater native fish richness than 
#reservoirs.

# Tributaries are expected to have between 15 and 60% less native fish richness 
#than mainstem river sites.

# Because our posterior probability distributions spanned zero for Control/Impact,
# Reference/Impact, Dissolved Oxygen Concentration, and Specific Conductance, there
# is no evidence that these variables described a meaningful amount of variance in
#native fish species richness.


# Our model described approximately 64% of the variance in native fish species
# richness, had normally distributed residuals, and had a bayesian p-value between
# 0.1 and 0.9 (0.45), indicating that the model fit the data well.





library(ggplot2)
library(patchwork)
library(cowplot)
# AUC Plot B (ROC Curve using ggplot2 with step function, flipped axes)
roc_coords <- coords(roc_curve, ret = c("specificity", "sensitivity"), transpose = FALSE)

roc_df <- data.frame(
  specificity = roc_coords$specificity,
  sensitivity = roc_coords$sensitivity
)

# Extract fixed effect summary statistics
fixed_effects <- spat.inla.diversity$summary.fixed

# Extract the posterior mean, 2.5% and 97.5% credible intervals
forest_data <- data.frame(
  Variable = rownames(fixed_effects),
  OR = exp(fixed_effects$mean),  # Convert log-odds to odds ratio
  CI_Lower = exp(fixed_effects$`0.025quant`),  # 2.5% credible interval
  CI_Upper = exp(fixed_effects$`0.975quant`)  # 97.5% credible interval
)

print(forest_data)  # View extracted data

# Remove the intercept and reformat the variable names without reentering numeric values
forest_data <- forest_data[-1, ]  # Remove the intercept row

# Reverse the order of Variable if not already reversed
forest_data$Variable <- factor(forest_data$Variable, levels = rev(forest_data$Variable))

# Reorder and rename the Poisson model variables as per your requirements
forest_data$Variable <- recode(forest_data$Variable,
                                       "Habitat" = "Habitat (reservoir vs. stream)",
                                       "Size" = "Stream size (mainstem vs. tributary)",
                                       "Reference" = "Reference sites vs. impact sites",
                                       "Control" = "Control sites vs. impact sites",
                                       "Dissolved.Oxygen.z" = "Scaled dissolved oxygen (mg/L)",
                                       "Specific.Conductance.z" = "Scaled specific conductance (\u03BCS/cm)")

# Add asterisks to rows where the credible interval does not span 1 (Poisson model)
forest_data$Label <- ifelse(forest_data$CI_Lower > 1 | forest_data$CI_Upper < 1, "*", "")

# Extract fixed effect summary statistics for the Poisson model
fixed_effects_poisson <- spat.inla.rich$summary.fixed

# Extract the posterior mean, 2.5% and 97.5% credible intervals for the Poisson model
forest_data_poisson <- data.frame(
  Variable = rownames(fixed_effects_poisson),
  OR = exp(fixed_effects_poisson$mean),  # Convert log-odds to odds ratio (Poisson case)
  CI_Lower = exp(fixed_effects_poisson$`0.025quant`),  # 2.5% credible interval
  CI_Upper = exp(fixed_effects_poisson$`0.975quant`)  # 97.5% credible interval
)

# Remove the intercept and reformat the variable names without reentering numeric values
forest_data_poisson <- forest_data_poisson[-1, ]  # Remove the intercept row

# Reverse the order of Variable if not already reversed
forest_data_poisson$Variable <- factor(forest_data_poisson$Variable, levels = rev(forest_data_poisson$Variable))

# Reorder and rename the Poisson model variables as per your requirements
forest_data_poisson$Variable <- recode(forest_data_poisson$Variable,
                                       "Habitat" = "Habitat (reservoir vs. stream)",
                                       "Size" = "Stream size (mainstem vs. tributary)",
                                       "Reference" = "Reference sites vs. impact sites",
                                       "Control" = "Control sites vs. impact sites",
                                       "Dissolved.Oxygen.z" = "Scaled dissolved oxygen (mg/L)",
                                       "Specific.Conductance.z" = "Scaled specific conductance (\u03BCS/cm)")

# Add asterisks to rows where the credible interval does not span 1 (Poisson model)
forest_data_poisson$Label <- ifelse(forest_data_poisson$CI_Lower > 1 | forest_data_poisson$CI_Upper < 1, "*", "")

# Add a column to distinguish between models (Binomial and Poisson)
forest_data$Model <- "Non-Extreme Low Diversity (Binomial)"
forest_data_poisson$Model <- "Species Richness (Poisson)"

# Combine both datasets
combined_forest_data <- rbind(forest_data, forest_data_poisson)

# Reorder Variable to maintain the current order
combined_forest_data$Variable <- factor(combined_forest_data$Variable, levels = rev(unique(combined_forest_data$Variable)))

# Reorder the variables in both forest plots (A and C)
forest_data$Variable <- factor(forest_data$Variable, levels = c("Habitat (reservoir vs. stream)",
                                                               "Stream size (mainstem vs. tributary)",
                                                                "Reference sites vs. impact sites",
                                                                "Control sites vs. impact sites",
                                                                "Scaled dissolved oxygen (mg/L)",
                                                                "Scaled specific conductance (\u03BCS/cm)"))  # Specify the correct levels
forest_data_poisson$Variable <- factor(forest_data_poisson$Variable, levels = c("Habitat (reservoir vs. stream)",
                                                                                "Stream size (mainstem vs. tributary)",
                                                                                "Reference sites vs. impact sites",
                                                                                "Control sites vs. impact sites",
                                                                                "Scaled dissolved oxygen (mg/L)",
                                                                                "Scaled specific conductance (\u03BCS/cm)"))  # Match levels

# Update the AUC Plot B (ROC Curve with consistent font and size for text)
auc_plot <- ggplot(roc_df, aes(x = specificity, y = sensitivity)) +
  geom_step(color = "black", size = 1) +
  geom_abline(slope = -1, intercept = 1, linetype = "dashed", color = "grey") +
  annotate("text", x = 0.00, y = 0.05, label = "b", size = 8, fontface = "bold", hjust = 0) +
  annotate("text", x = 0.6, y = 0.2, label = "AUC = 0.97", size = 6, fontface = "plain") +  # Standardized font and size
  theme_classic() +
  labs(x = "Specificity", y = "Sensitivity") +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA)
  )

# Update Forest Plot A (Non-Extreme Low Diversity)
plot_binomial <- ggplot(forest_data, aes(x = Variable, y = OR)) +
  geom_pointrange(aes(ymin = CI_Lower, ymax = CI_Upper), size = 1, color = "black") +
  geom_text(aes(label = Label), nudge_x = 0.2, color = "black", size = 8) +
  coord_flip() +
  theme_classic() +
  labs(y = "Odds ratio (non-extreme low diversity)",
       x = NULL) +
  annotate("text", x = 1, y = 1e30, label = "a", size = 8, fontface = "bold", hjust = 1) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  scale_y_log10() +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black")

# Update Forest Plot C (Species Richness)
plot_poisson <- ggplot(forest_data_poisson, aes(x = Variable, y = OR)) +
  geom_pointrange(aes(ymin = CI_Lower, ymax = CI_Upper), size = 1, color = "black") +
  geom_text(aes(label = Label), nudge_x = 0.2, color = "black", size = 8) +
  coord_flip() +
  theme_classic() +
  labs(x = "Variable", y = "Rate ratio (species richness)") +
  annotate("text", x = 6, y = 8, label = "c", size = 8, fontface = "bold", hjust = 1) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(hjust = 1.2),
    axis.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  scale_y_log10() +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black")

# Update Richness Plot D (Observed vs Fitted Richness)
plot_richness <- ggplot(CSD2, aes(x = mean, y = Richness)) +
  geom_point() +
  theme_bw() +
  xlim(c(0, 8)) +
  ylim(c(0, 8)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  xlab("Fitted native fish species richness") +
  ylab("Observed native fish species richness") +
  annotate("text", x = 0.1, y = 7.6, label = "d", size = 8, fontface = "bold", hjust = 0) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  geom_text(x = 2, y = 6, aes(label = txt), parse = TRUE) +
  theme(
    strip.text.x = element_text(size = 20)
  )

# Combine the plots using patchwork with consistent spacing
combined_plot <- (plot_binomial + auc_plot) / (plot_poisson + plot_richness) +
  plot_layout(guides = "collect") & theme(legend.position = 'none')

# Display the final plot
print(combined_plot)






###---Heatmap of fish detections, shannon diversity, and species richness---###

### Required libraries
library(vegan)
library(tidyverse)
library(ggrepel)
library(readxl)
library(ggtext)
library(ggpubr)
library(patchwork)
library(ggh4x)
library(stringr)

## Import data

# import metadata and metabarcoding data
m <- data1
d <- data2

###---Data Prep---###

## Prepare Fish metabarcoding data set
# subset fish only reads from metabarcoding data
d_fish <- subset(d, Type=="Fish")
# select 2023 only for metadata and pull habitat type and site ID columns to merge to final dataframe later
m_2023 <- 
  m %>%
  subset(Year == 2023) %>%
  select(Site.Name, Habitat) %>%
  unique()

# clean up taxa CommonName for plotting
d_fish <-
  d_fish %>%
  mutate(CommonName = recode(CommonName, 
                             "Klamath Large Scale/Klamath Small Scale, Shortnose Sucker, Lost River Sucker" = "Catostomidae spp.",
                             "Rainbow Trout/steelhead" = "Rainbow Trout",
                             "Black/Brown Bullhead" = "Ameiurus spp.",
                             "Carps, Minnows"="Cyprinidae spp.",
                             "Lamprey"="Entosphenus spp.",
                             "Sculpin"="Cottus spp."))

# create a dataframe of species status
# grab unique common names
CommonName <- unique(d_fish$CommonName)
CommonName

# assign status in order of appearance in the list generated in previous step
status = c("Exotic", #Black Crappie
           "Exotic", #Ameiurus spp.
           "Exotic", #Bluegill Sunfish
           "Exotic", #Brown Trout
           "Exotic", #Cyprinidae spp.
           "Native", #Chinook Salmon
           "Native", #Coho Salmon
           "Exotic", #Goldfish
           "Exotic", #Fathead Minnow
           "Exotic", #Golden Shiner
           "Exotic", #Green Sunfish
           "Native", #Catostomidae spp.
           "Exotic", #Largemouth Bass
           "Native", #Entosphenus spp.
           "Native", #Rainbow Trout
           "Native", #Cottus spp.
           "Native", #Speckled Dace
           "Exotic", #Pumpkinseed Sunfish
           "Native", #Tui Chub
           "Exotic", #Yellow Bullhead
           "Exotic" #Yellow Perch
)

# merge common name list and status list into dataframe
# this will be added to final dataframe at a later step
status_df <- data.frame(CommonName, status)
# check that it assigned correctly
view(status_df)

## Create dataframe to add levels of site in a downstream to upstream order (for plotting)
# this step groups sites based on sampling reach and assigned a grouping level based on downstream to upstream position of Site.Name
# site names
unique(m_2023$Site.Name) 
# create new dataframe that places sites in correct order and assign levels, groupings (within reach) and groupin levels for nested axis
levels <- data.frame("Site.Name" = c("Klamath2", "Scotch1", "Klamath4", 
                                     "Jenny1", "Klamath6","Klamath10", 
                                     "Beaver1","Klamath12","Klamath34", 
                                     "Klamath36", "Klamath_Ranch","Klamath1", 
                                     "Bogus1", "Scotch2", "Scotch3", 
                                     "Jenny2", "Jenny3","Jenny4", 
                                     "Klamath8", "Fall1", "Fall2", 
                                     "Beaver2", "Klamath14", "Klamath16", 
                                     "Klamath18", "Klamath20","Klamath22", 
                                     "Klamath24", "Klamath26", "Klamath28", 
                                     "Klamath30", "Klamath32" ,"Shovel1" ,
                                     "Shovel2", "Hayden1",  "Spencer1" ,
                                     "Spencer2", "Spencer3", "Scott1", 
                                     "Scott2", "Scott3" ,"Kelsey1" ,"Kelsey2", "Canyon1"),
                     "level"=1:44,
                     "grouping"=rep(c("Iron Gate","Copco","JC Boyle", "Klamath (Below IGR)",
                                      "Bogus Creek","Scotch Creek", "Jenny Creek","Klamath (Below CR)",
                                      "Fall Creek","Beaver Creek",
                                      "Klamath (Below JCB)","Shovel Creek","Hayden Creek","Spencer Creek",
                                      "Scott River (Control)"),
                                    times=c(5,3,2,2,1,2,3,1,2,1,10,2,1,3,6)),
                     "grouping_level"=rep(c(1:15), times=c(5,3,2,2,1,2,3,1,2,1,10,2,1,3,6)))

# create vector of unique values for grouping level
grouping_level <- 
  levels %>%
  select(grouping_level, grouping) %>%
  unique()

# sum up reads per detected taxa per site
# not sure if this is done before or after flow corrected reads caluclated..
d_clean <-
  d_fish %>%
  group_by(Site.Name, CommonName) %>%        #combine reads for same spp w/in siteName
  mutate(sum_ESVsize = sum(ESVsize)) %>%          #sum reads
  select(CommonName, Site.Name, sum_ESVsize) %>%   # select desired columns 
  unique()  # drop duplicate rows 

head(merged_data_RPS_clean)
head(d_clean)

## Replace with flow corrected reads
# Reshape merged_data_RPS_clean to match d_clean format
d_clean <- merged_data_RPS_clean %>%
  pivot_longer(cols = -Site.Name,  # Pivot all columns except Site.Name
               names_to = "CommonName", 
               values_to = "sum_ESVsize") %>%
  filter(sum_ESVsize > 0)  # Remove zero counts if needed

# add option to use giga-reads
d_clean$giga_reads<-d_clean$sum_ESVsize/10^9


## Merge dataframes to create heatmap dataframe
# Merge various files for flow  corrected reads heat map
hm_df <-
  merge(d_clean, status_df, by = "CommonName") %>%
  merge(levels, by="Site.Name", all.y=TRUE)%>%
  merge(m_2023, by="Site.Name") 

# convert 0 reads to NA
# important for making 0 values plot transparent in heatmap
hm_df <- 
  hm_df %>%
  mutate(across(c("sum_ESVsize", "giga_reads"), ~ na_if(.,0))) 

# set levels of site.name to order x axis correctly. Based on levels_new order from levels dataframe
hm_df$Site.Name <- factor(hm_df$Site.Name, levels=levels$Site.Name[order(hm_df$level)], ordered=TRUE)
hm_df$grouping <- factor(hm_df$grouping, levels=grouping_level$grouping[order(hm_df$grouping_level)], ordered=TRUE)

#----------Heat map of reads ----------#

# labels for facet_grid
hm_labels_habitat<-c("Reservoir"="Reservoir", "Stream"="Stream", "Native"="Native", "Exotic"="Exotic")

# edit spp names to have scientific names display as italic in plot
spplist <- c("Goldfish"="Goldfish", 
             "Cyprinidae spp."="<i>Cyprinidae</i> spp.", 
             "Yellow Perch"="Yellow Perch", 
             "Speckled Dace"="Speckled Dace",
             "Golden Shiner"="Golden Shiner",   
             "Yellow Bullhead"="Yellow Bullhead",
             "Largemouth Bass"="Largemouth Bass", 
             "Pumpkinseed Sunfish"="Pumpkinseed Sunfish", 
             "Fathead Minnow"="Fathead Minnow", 
             "Brown Trout"="Brown Trout", 
             "Bluegill Sunfish"="Bluegill Sunfish", 
             "Catostomidae spp."="<i>Catostomidae</i> spp.",
             "Cottus spp."="<i>Cottus</i> spp.", 
             "Rainbow  Trout"="Rainbow Trout", 
             "Coho Salmon"="Coho Salmon", 
             "Black Crappie"="Black Crappie", 
             "Tui Chub"="Tui Chub",
             "Ameiurus spp."="<i>Ameiurus</i> spp.", 
             "Chinook Salmon"="Chinook Salmon", 
             "Entosphenus spp."="<i>Entosphenus</i> spp.", 
             "Green Sunfish"="Green Sunfish")

# no fish were detected at Scotch Creek but we want it to show up in the figure. Adding this placeholder data to make sure it is not excluded. Read count = 0 so nothing will be plotted for this site.
hm_df <- hm_df %>%
  mutate(CommonName = ifelse(grouping == "Scotch Creek", "Yellow Perch", CommonName), 
         status = ifelse(grouping == "Scotch Creek", "Exotic", status))

# plot with site names only
klamath.hm <- 
  hm_df %>%
  ggplot(., aes(x=interaction(level,grouping_level), y=CommonName, fill= sum_ESVsize)) + # change to fill=giga_reads after flow corrected calc
  geom_tile(colour="white", linewidth=0.25)+
  scale_y_discrete(limits=rev, label=spplist)+
  scale_x_discrete(guide = "axis_nested")+
  scale_fill_gradient(low = "#d8e1cf", high = "#438484", na.value = "white", 
                      name="10<sup>9</sup> DNA reads/sec.")+
  theme_bw()+
  theme(panel.spacing = unit(0, "lines"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.title = element_markdown(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_markdown()) +
  labs(x="Location", y="Taxa Detected") +
  facet_grid(cols=vars(Habitat), rows= vars(status), scales = "free", 
             space= "free", labeller=as_labeller(hm_labels_habitat))


# remove tick marks and labels from X axis
klamath.hm.1_notic <- 
  hm_df %>%
  ggplot(., aes(x=interaction(level, grouping_level), y=CommonName, fill= sum_ESVsize)) +
  # change to fill=giga_reads after flow corrected calc
  geom_tile(colour="white", linewidth=0.25)+
  scale_y_discrete(limits=rev, label=spplist)+
  #scale_x_discrete(guide = "axis_nested")+
  scale_fill_gradient(low = "#d8e1cf", high = "#438484", na.value = "white",
                      name="10<sup>9</sup> DNA reads/sec.", breaks=c(500,1000,1500,2000))+ 
  theme_bw()+
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        strip.placement = "outside",
        strip.text = element_text(size = 11, face = "bold"),
        legend.position = "bottom",
        legend.justification = "center", 
        legend.title = element_markdown(size=10, vjust= 0.75),
        legend.text = element_text(size=10),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_markdown())+ 
  labs(y="Taxa Detected")+
  facet_grid(cols=vars(Habitat), rows= vars(status), scales = "free", 
             space= "free", drop=TRUE, labeller=as_labeller(hm_labels_habitat))

# export heatmap 
# aspect_ratio=2.5
# ggsave("Analysis/R/figs/klamath.hm.1.numeric_notic.png", 
#     plot=klamath.hm.1_notic, width=12, height=aspect_ratio*2.25, dpi=600)


#---------- Richness and Shannon Diversity ----------#

# subset necessary columns from larger dataframe
hm_df_subset <- 
  hm_df %>%
  select(Site.Name,grouping_level,grouping, Habitat, level)

# native and exotic richness and diversity per site separately
shan.rich <- 
  hm_df %>%                              #grab original df
  filter(status=="Native") %>%
  group_by(Site.Name) %>%            #group by site name and status
  summarise(richness = specnumber(sum_ESVsize),
            shannon = diversity(sum_ESVsize, , index="shannon")) %>%# calculate diversity and richness
  mutate(shan.rounded = round(shannon, 2)) %>%
  merge(hm_df_subset, by="Site.Name", all.y=TRUE) %>% #add some additional data for plotting 
  mutate(across(c(richness,shannon,shan.rounded), ~ na_if(.,0))) %>%
  unique() %>%
  merge(data.frame(shan.name = "Native Species Diversity",
                   rich.name = "Native Species Richness"))

# set levels of site.name to order x axis correctly. Based on levels_new order from levels table
shan.rich$grouping = factor(shan.rich$grouping, levels=grouping_level$grouping[order(shan.rich$grouping_level)], ordered=TRUE)

### heatmap of shannon diversity
shan.plot <- 
  shan.rich %>%
  ggplot(., aes(x=interaction(level, grouping_level), y=shan.name, fill= shannon)) + 
  geom_tile(colour="black", linewidth=0.2)+
  #  geom_text(aes(label = shan.rounded), color ="black", size=1.75) +  
  scale_y_discrete(limits=rev)+
  #  scale_x_discrete(guide = "axis_nested")+
  scale_fill_gradient(low = "#C5FFC2", high = "#196B24", na.value = "white", 
                      name = "Native Species Diversity", limits=c(0,2), breaks=c(0.5,1,1.5,2))+
  theme_bw()+
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "none",
        strip.text = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=10),
        legend.title = element_text(size=10, vjust= 0.75),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),) + 
  facet_grid(cols=vars(Habitat),  scales = "free", 
             space= "free")


### heatmap of richness
rich.plot <- 
  shan.rich %>%
  ggplot(., aes(x=interaction(level, grouping_level), y=rich.name, fill= richness)) + 
  geom_tile(colour="black", linewidth=0.2)+
  # geom_text(aes(label = richness), color="black", size=2) +  
  scale_y_discrete(limits=rev)+
  scale_x_discrete(guide = "axis_nested")+
  scale_fill_gradient(low = "#cdcbcb", high = "#696969", na.value = "white", 
                      name="Native Species Richness", limits=c(0,8), breaks=c(2,4,6,8))+
  theme_bw()+
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "none",
        legend.position = "bottom",
        legend.text = element_text(size=10),
        legend.title = element_text(size=10, vjust= 0.75),
        axis.ticks.x=element_blank(),
        strip.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text( vjust = 1)) + 
  labs(x="Location") +
  facet_grid(cols=vars(Habitat),  scales = "free", 
             space= "free")

# Combine heatmap of reads, shannon diversity index and richness into one plot
hm.rich.shan <-
  (klamath.hm.1_notic / plot_spacer() /shan.plot / plot_spacer()/rich.plot) +  
  plot_layout(widths = c(12,12, 12,12, 12), heights = unit(c(3.75,-0.3, .2, -0.3, .2), c('in', 'in')))  + 
  plot_layout(guides = "collect") & theme(legend.position = 'bottom', 
                                          legend.box = 'vertical')

# Export combined plot. Remaining edits made in graphics editor.
aspect_ratio<-2.5
ggsave("figures/klamath.heatmap.richshan_notext3.png", 
       plot=hm.rich.shan, width=9.5, height=aspect_ratio*2.5, units="in", dpi=600)

# Edits made in graphics editor include:
# 1. Removing x axis labels and moving nested axis labels up to sit below x axis
# 2. Aligning legends vertically a bit better and moving to left of plot
# 3. Adding list of locations that correspond to nested axis:
hm_df %>%
  select(grouping_level, grouping) %>%
  unique() %>%
  arrange(grouping_level)

#---------END-------------#







# Klamath molecular library Manuscript
# Cramer Fish Sciences
# Katie Karpenko
# katie.karpenko@fishsciences.net
# Last Updated: Feb/27/2025

### NMDS plot ###

# Libraries
library(vegan)
library(tidyverse)
library(ggrepel)
library(readxl)
library(ggtext)
library(ggpubr)

#----------NMDS ----------#

# import final metabarcoding data set and metadata file
d <- data2 #update path!
m <-data1 #update path!

# subset metadata to include only 2023
covar2023<-subset(m, Year=="2023")

# subset fish only hits from final metabar df
d_fish <- subset(d, Type=="Fish")

# clean up names of taxa for plotting
d_fish$FigureName <- ifelse(grepl("spp.", d_fish$ReportName), d_fish$ReportName, d_fish$CommonName)

d_fish <-
  d_fish %>%
  mutate(FigureName = recode(FigureName, 
                             "Klamath Large Scale/Klamath Small Scale, Shortnose Sucker, Lost River Sucker" = "Catostomidae spp.",
                             "Rainbow Trout/steelhead" = "Rainbow Trout"))

# clean up data and select desired columns/rows
d_clean <- 
  d_fish %>%
  select(Site.Name, ESVsize, FigureName) %>%      #select columns needed, lumping data by siteName
  group_by(Site.Name, FigureName) %>%             #combine reads for same spp w/in siteName
  summarize(sum_ESV = sum(ESVsize)) %>%          #sum reads
  drop_na(Site.Name)                              #remove NA (these are NTCs)

## Select rarefaction depth

# sampling coverage
sampling_coverage<-
  d_fish %>%
  group_by(Sample.ID)%>%
  summarize(n_seqs=sum(ESVsize))

# density plot of sequencing coverage
sampling_coverage %>%
  ggplot(aes(x=n_seqs)) +
  geom_density() 

# histogram plot of sequencing coverage
sampling_coverage %>%
  ggplot(aes(x=n_seqs)) +
  geom_histogram(binwidth = 10000)+
  coord_cartesian(xlim=c(0,200000))

# scatter plot
sampling_coverage %>%
  ggplot(aes(x=1, y=n_seqs)) +
  geom_jitter()+
  scale_y_log10()

sampling_coverage %>%
  arrange(n_seqs) %>%
  ggplot(aes(x=1:nrow(.), y=n_seqs))+
  geom_line()

sampling_coverage %>%
  arrange(n_seqs) %>%
  print(n=20)

coverage_stats<-
  d_clean %>%
  group_by(Site.Name) %>%
  summarise(n_seqs = sum(sum_ESV),
            n_sings=sum(sum_ESV==100),
            goods = 100*(1-n_sings/n_seqs)) %>%
  filter(n_seqs>7364)

# at 7364, coverage is good
coverage_stats %>%
  ggplot(aes(x=n_seqs, y=goods))+
  geom_point()

coverage_stats %>%
  arrange(goods)

# setting min sequencing depth at the min amount of sequences observed per site. 
# this will be sampling depth for rarefaction.
min_seqs = 7364

# Pivot data wider to make matrix
d_wide <-
  d_clean %>%
  pivot_wider(names_from="FigureName", values_from = "sum_ESV", values_fill=0) %>%
  as.data.frame()

rownames(d_wide) <- d_wide$Site.Name 
d_wide <- d_wide[,-1]
d_wide <- as.matrix(d_wide)

### NMDS ###
set.seed(5523)
dist <- avgdist(d_wide, dmethod="bray", sample = min_seqs) #Takes "min_seqs" # of sequences from each sample and calculates distances, the repeats.  Number of iterations can be specified. Default is iterations = 100

### TO RERUN NMDS, remove # symbol from code below ###

# set.seed(26343)
# nmds_k2 = metaMDS(dist,distance = "bray", k = 2, maxit = 999, trymax = 999, wascores = TRUE)
# write_rds(nmds_k2, "Analysis/R/data_clean/klamath.nmds.rds") # update directory

# import NMDS result to keep it consistent with figure. There may be slight differences in position of points if it is rerun.
nmds_k2<- read_rds("results/klamath.nmds.rds")

# Check how well ordination plot represents real data (goodness function and Shepards plot)
goodness(nmds_k2) # vegan version
plot(nmds_k2, type="t", main = "goodness of fit")

#Shepards diagram
stressplot(nmds_k2) 

# extract stress value
stress_value <- nmds_k2$stress

# Plotting points in ordination space
plot(nmds_k2, "sites")   # Produces distance 
orditorp(nmds_k2, "sites")   # Gives points labels


### Visualization of NMDS with ggplot ###

# NMDS by habitat type
nmds.plot <-
  scores(nmds_k2) %>%                                    
  as_tibble(rownames = "Site.Name") %>%        #tibbles don't have row names. Mirrors row names into a samples column
  inner_join(., covar2023, by="Site.Name") %>%   #bring in metadata 
  ggplot (aes(x=NMDS1, y = NMDS2)) +                   
  geom_point(aes(fill = Habitat), size = 3, alpha=0.9, shape=21) +
  stat_ellipse(aes(fill = Habitat, color = Habitat), geom = "polygon", alpha = 0.2) +
  scale_fill_manual(name="Habitat", values = c(
    "grey",    #Reservoir
    "#438484"))+    #Stream
  scale_color_manual(name="Habitat", values = c(
    "grey",    #Reservoir
    "#438484"))+    #Stream
  theme( 
    panel.background = element_blank(),  #the background of the plot itself
    panel.border = element_rect(color = "gray", fill=NA), #border of plot panel
    axis.text = element_text(size = 12),               # Axis text size
    axis.title = element_text(size = 12),              # Axis titles
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  # made same color as background so as to disappear
    strip.background = element_rect(fill="#d1e5f0"),
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),      #color of background behind is facet title
    legend.position="bottom",                 #put "none" here to remove legend
    legend.key = element_blank(),
    legend.title = element_text(size = 14),      
    legend.text = element_text(size = 12)   ) +
  # geom_text(aes(label = siteName, 
  #              y = NMDS2 + 0.05),    # label points with site names (for QC), include if needed
  #          size = 3.5, color = "#304242", vjust = 0) +  
  annotate("text", x = 0.5, y = 0.6, label = "Stress: 0.11", size = 4.5, color = "#304242", hjust = 0.5)


###------Significant spp vectors-------#

### permutational analyses may result in slight variation from figure results
## import results from figure to avoid changes
sig.spp.scrs <- read.csv("results/sig.spp.scrs.csv") # update path

## Remove # to run code below if new significant species vectors need to be generated.

## Add sig spp vectors

# dist.spp.fit <- envfit(nmds_k2, d_wide, permutations = 999)                # this fits species vectors (i.e., intrinsic variables)
#spp.scrs <- as.data.frame(scores(dist.spp.fit, display = "vectors"))       #save species intrinsic values into dataframe
#spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))                  #add species names to dataframe
#spp.scrs <- cbind(spp.scrs, pval = dist.spp.fit$vectors$pvals)             #add pvalues to dataframe so you can select species which are significant
#spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
#sig.spp.scrs <- subset(spp.scrs, pval<=0.001) #subset data to show species significant at 0.001
#head(spp.scrs)

# Export sig spp scores
# write.csv(sig.spp.scrs, "Analysis/R/data_clean/sig.spp.scrs.csv")

# Create status field for significant species dataframe
# grab unique common names
CommonName.sigspp <- unique(d_fish$FigureName)

# set status for each species
status <- c("Exotic", #Black Crappie
            "Exotic", #Ameiurus spp.
            "Exotic", #Bluegill Sunfish
            "Exotic", #Brown Trout
            "Exotic", #Cyprinidae spp.
            "Native", #Chinook Salmon
            "Native", #Coho Salmon
            "Exotic", #Goldfish
            "Exotic", #Fathead Minnow
            "Exotic", #Golden Shiner
            "Exotic", #Green Sunfish
            "Native", #Catostomidae spp.
            "Exotic", #Largemouth Bass
            "Native", #Entosphenus spp.
            "Native", #Rainbow Trout
            "Native", #Cottus spp.
            "Native", #Speckled Dace
            "Exotic", #Pumpkinseed Sunfish
            "Native", #Tui Chub
            "Exotic", #Yellow Bullhead
            "Exotic" #Yellow Perch
)
# create df of common names and status
status_df_sigspp <- data.frame(CommonName.sigspp, status)
# check that it assigned correctly
view(status_df_sigspp)

# Add status to sig.spp.scrs df to define linetype by status for visualizing verctors
sig.spp.scrs<-merge(sig.spp.scrs, status_df_sigspp, by.x="X", by.y="CommonName.sigspp")

# adding species to nmds plot from above
nmds.vectors.names<-
  nmds.plot +
  geom_segment(data = sig.spp.scrs, 
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2, 
                   linetype = status), 
               arrow = arrow(length = unit(0.20, "cm")), colour = "#304242", lwd = 0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs, 
                           aes(x = NMDS1, y = NMDS2, label = Species), 
                           cex = 3, direction = "both", segment.size = 0.1, 
                           color = "#304242") + #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  scale_linetype_manual(values = c("Native" = "solid", "Exotic" = "dashed"), 
                        labels = c("Exotic", "Native"))+
  guides(linetype = guide_legend(title = "Status"))

# create object of NMDS plot with vectors without spp names. Spp names added in graphics editor after export.
nmds.vectors<-
  nmds.plot +
  geom_segment(data = sig.spp.scrs, 
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2, 
                   linetype = status), 
               arrow = arrow(length = unit(0.20, "cm")), colour = "#304242", lwd = 0.3) + #add vector arrows of significant species
  scale_linetype_manual(values = c("Native" = "solid", "Exotic" = "dashed"), 
                        labels = c("Exotic", "Native"))+
  guides(linetype = guide_legend(title = "Status"))

#export NMDS without vectors
aspect_ratio<-2.5
ggsave("figures/klamath.NMDS2.png", 
       plot=nmds.plot, width=8.5, height=aspect_ratio*2.5, dpi=600) # update path

#export NMDS with vectors and species names
ggsave("figures/klamath.NMDS.vectors.spp2.png", 
       plot=nmds.vectors.names, width=8.5, height=aspect_ratio*2.5, dpi=600) # update path

#export NMDS with vectors only (no spp. names, added in text editor to manually control placement)
ggsave("figures/klamath.NMDS.vectors2.png", 
       plot=nmds.vectors, width=8.5, height=aspect_ratio*2.5, dpi=600) # update path

#---------END-------------#

