library(tidyr)
library(dplyr)


##This script build the dataframe for estimating unknown costs of IAS in EU countries
#Need to have the combination of 
#species/sector/type of cost

# Determine what combination of species/type of cost etc ------------------

#Open InvaCost database used to build the model
eu.df.full <- read.csv("/Users/morganehenry/Documents/GitHub/EU_extrapolation/data/EU_extrapol_processed_data.csv")

#Only include data at the country level

eu.df.country <- eu.df.full[eu.df.full$Spatial_scale == "Country",]

# Country only data + combination of species, sector and type of cost ----------------------------------------

comb4 <- paste(eu.df.country$Species, ",", eu.df.country$Impacted_sector, ",", eu.df.country$Type_of_cost_merged) #Add comma so we can easily separate
print(paste("Number of unique combination when considering species and type of cost:", length(unique(comb4))))

comb.unique4 <- data.frame(comb.unique= unique(comb4)) #need as dataframe to use the function separate

##Convert this vector to a dataframe
df.comb4 <- separate(comb.unique4, col= "comb.unique", c("Species", "Impacted_sector", "Type_of_cost_merged"), sep=" , ")

# Prepare data for extrapolation from species.csv -------------------------

##Open dataframe containing the presence of invasive species in Europe

sp <- read.csv("/Users/morganehenry/Documents/GitHub/EU_extrapolation/data/species_EU_countries.csv", header=T)

# How many specie are present in both the InvaCost data and the extrapolation df --------
nb.sp4 <- intersect(df.comb4$Species, sp$Species) ; print(length(unique(nb.sp4))) #Country only data, combination of species, impacted sector and type of cost: 46

# Add country data to the dataframe species.csv ---------------------------
country.df <- read.csv("/Users/morganehenry/Documents/GitHub/EU_extrapolation/data/Country_new.csv")
country.df <- country.df[,-1]
country.df <- country.df %>% rename(Official_country = Country)

##Merge sp dataframe with the country information (GDP, population, etc)
sp <- merge(sp, country.df, by="Official_country")

##Merge sp with % gdp per sector 
df.sect <- read.csv("/Users/morganehenry/Documents/GitHub/EU_extrapolation/data/sector_gdp.csv")

sp <- merge(sp, df.sect, by= "Official_country")

# Build dataframe for extrapolation ---------------------------------------
sp4 <- sp[sp$Species %in% df.comb4$Species,]
df.extrapolation4 <- merge(df.comb4, sp4, by= "Species")
write.csv(df.extrapolation4, "./data/df_extrapolation_Countrydata_SpeciesSectorTypeCost.csv", row.names = F)
