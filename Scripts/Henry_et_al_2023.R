############ Spatial and sectoral projections ############

library(doBy)
library(tidyr)
library(dplyr)

### Functions necessary in the script "funcs_all.R"

source("Henry_2023_funcs.R")

eu.df.full <- read.csv("EU_extrapol_processed_data.csv") #change file directory
eu.df.full$comb2 <- paste(eu.df.full$Species, eu.df.full$Impacted_sector, eu.df.full$Type_of_cost_merged) #create a column with the combination of species name, sector impacted, and type of costs (damage, management, mixed)
eu.df.full <- eu.df.full[order(eu.df.full$comb2),] #reorder so everything matches later
dfsp <- read.csv("df_extrapolation_Countrydata_SpeciesSectorTypeCost.csv") #data for extrapolation

##Give all the combinations of predictors (in order to tell scale_fit_fun what it should run)
pos.pred <- c("GDP_2019", "Population_2019", "Surface_area_km2", "Impacted_sector")
full_gamma <-c(0,0,0,0,NA) #named "full" gamma vector to make the optimization function (scale_fit_fun) faster
names(full_gamma)=c(pos.pred, "sd")
pred.l <- list(m0= NULL, m1= pos.pred[1], m2= pos.pred[2], m3= pos.pred[3], m4= pos.pred[4], m5= pos.pred[c(1,2)], m6= pos.pred[c(1,4)], m7= pos.pred[c(1,3)], 
               m8= pos.pred[c(4,2)], m9= pos.pred[c(4,3)], m10= pos.pred[c(2,3)], m11= pos.pred[c(1,3,4)], m12= pos.pred[c(1,2,3)], m13= pos.pred[c(1,2,4)],
               m14= pos.pred[c(2,3,4)], m15= pos.pred) #create all combinations of predictors

#We don't need to keep all variables in the df
pred <- c("Cost_ID", "Official_country", "Species", "Spatial_scale", "Cost_estimate_per_year_2017_USD_exchange_rate", "Impacted_sector", "Type_of_cost_merged", "cost_mil", "Population_2019", "Surface_area_km2", "GDP_2019",
          "Agri_forest_fish_perc_gdp", "Health_expend_perc_gdp", "comb2")
eu.df <- eu.df.full[eu.df.full$Spatial_scale == "Country",pred]


for(m in 1:length(pred.l)){ #loop over each model
  pred.list <- pred.l[[m]] #to obtain the correct predictors
  
  par1 <- c(rep(0.1, length(pred.list)), sd(eu.df$cost_mil, na.rm=T)) #starting values for optimization function
  
  names(par1) <- c(pred.list, "sd")
  
  df.extrap <- do_df_extrap(df.sp= dfsp, df.fit= eu.df) #this contains ONLY the extrapolation data, no invacost data
  
  df.fit1 <- eu.df
  eu.df$type <- "invacost"
  df.fit1$type <- "extrapol"
  df.fit.full <- rbind(eu.df, df.fit1)
  
  if(length(par1) == 1){ #cannot use optim with a single precictor
    par1 <- c(0,1000) #reassign par1 beacause optim need to know the min and max to look for the optimized value
    op.gamma <- optimize(scale_fit_fun, par1, maximum = F, df= df.fit.full, pred.list= pred.list, extrap= F, fit.test=F) #optimisation function
    gamma1 <- op.gamma$minimum #extract gamma1
    
  } else { #when more than a single parameter to fit
    op.gamma <- optim(par1, scale_fit_fun, df= df.fit.full, pred.list= pred.list, extrap= F, fit.test=F) ##optimisation function
    gamma1 <- op.gamma$par #extract gamma1
  }
  names(gamma1) <- c(pred.list, "sd")
  
  ##Metrics of fit --> calculate AIC for the given set of parameters gamma
  fit.t1 <- scale_fit_fun(gamma1 = gamma1, df.fit.full, pred.list= pred.list, extrap = F, fit.test=T) #extract df to calculate measure of fit
  v.AIC <- aic.fit(fit.t1, gamma1= gamma1)
  
  
  ##Estimation of projected costs
  ratio.fit <- get_ratio(df.fit.full, pred.list= pred.list) #obtain all the ratios
  df.extrap <- df.extrap[names(eu.df)] #make sure the columns names match
  df.full.extrap <- rbind(eu.df, df.extrap)
  
  #We can't estimate the cost for a certain comb2 that is NOT present in the invacost df
  #Need to remove comb2 that are present in extrapol but absent in invacost and vice versa
  absent <- setdiff(unique(df.full.extrap$comb2[df.full.extrap$type == "invacost"]), unique(df.full.extrap$comb2[df.full.extrap$type == "extrapol"])) #those are in the invacost dataset but NOT in the extrapol
  absent2 <- setdiff(unique(df.full.extrap$comb2[df.full.extrap$type == "extrapol"]), unique(df.full.extrap$comb2[df.full.extrap$type == "invacost"]))
  df.full.extrap <- df.full.extrap[!(df.full.extrap$comb2 %in% absent),]
  df.full.extrap <- df.full.extrap[!(df.full.extrap$comb2 %in% absent2),]
  
  #Run the function scale_fit_fun to obtain predictde costs
  extrap <- scale_fit_fun(gamma1 = gamma1, df.full.extrap, pred.list=pred.list, extrap = T, fit.test= F)
  
  sum_pred_cost <- sum(extrap$cost_pred, na.rm= T)
  
  gamma.df2 <- c(gamma1, v.AIC, sum_pred_cost)
  names(gamma.df2) <- c(pred.list, "sd", "v.AIC", "sum_pred_cost")
  
  write.csv(gamma.df2, file=paste0("./","EU_cost", "_", names(pred.l[m]),".csv"))

}


############################################################################################################

############ Temporal projection ############

#delete last 3 years to delete lag time 
extrapolation<-expanded_obs2 %>% filter(Impact_year <= "2017") 

#Robust Regression extrapolation 
#costs 
df <- extrapolation %>% group_by(Impact_year) %>% summarise(Cost= sum(cost_mil)) 
df2<-df[!(df$Impact_year<1980),]
robust.linear <- robustbase::lmrob(Cost~ Impact_year, maxit.scale=1000, data = df2) 

confidence.interval = 0.95 
prediction.years <- data.frame(Impact_year = 1980:2040) 
pred.robust.linear <- try(predict(robust.linear, 
                                  prediction.years, 
                                  interval = "confidence", 
                                  level = confidence.interval), 
                          silent = TRUE)

rownames(pred.robust.linear) <- prediction.years[, 1]
# Robust regression - quadratic effect 
#incomplete.year.threshold = NULL 
#incomplete.year.threshold <- maximum.year + 1 
#yearly.cost$transf.cost <- yearly.cost$Annual.cost 
#yearly.cost.calibration <- yearly.cost[-which(yearly.cost[, "Year"] >= incomplete.year.threshold), ] 
robust.quadratic <- robustbase::lmrob(Cost ~ Impact_year + I(Impact_year^2), maxit.scale=900, data = df2, 
                                      cov = ".vcov.w") # Covariance matrix estimated using asymptotic normality of the coefficients 
# See ?lmrob and Koller & Stahel 2011 
pred.robust.quadratic <- try(predict(robust.quadratic, 
                                     prediction.years, 
                                     interval = "confidence", 
                                     level = confidence.interval), 
                             silent = TRUE)
rownames(pred.robust.quadratic) <- prediction.years[, 1]
#model.RMSE["robust.quadratic", "RMSE.calibration"] <- sqrt(mean(residuals(robust.quadratic)^2)) 
#model.RMSE["robust.quadratic", "RMSE.alldata"] <- sqrt( 
# mean((pred.robust.quadratic[match(yearly.cost$Year, rownames(pred.robust.quadratic)), 
# "fit"] - yearly.cost$transf.cost)^2)) 
#plotting 
new_data_linear <- bind_cols(prediction.years, pred.robust.linear) 
colnames(new_data_linear) <- c("Impact_year", "fit_linear", "lwr_linear", "upr_linear") 
new_data_quadratic <- bind_cols(prediction.years, pred.robust.quadratic) 
colnames(new_data_quadratic) <- c("Impact_year", "fit_quadratic", "lwr_quadratic", "upr_quadratic")
new_data<-left_join(new_data_linear,new_data_quadratic)
myarrow=arrow(angle = 15, type = "closed")
EU_total<-ggplot(mapping = aes(x = Impact_year, y = Cost)) + 
  geom_point(data = df, aes())+ 
  geom_line(data = new_data, aes(y = fit_linear, x = Impact_year), arrow=myarrow, colour = "blue", 
            size = 1, alpha = 0.95, show.legend = T)+ 
  geom_ribbon(data = new_data, 
              mapping = aes(ymin = lwr_linear, ymax = upr_linear, x = Impact_year), fill = "blue", 
              inherit.aes = FALSE, alpha = 0.2)+ 
  geom_line(data = new_data, aes(y = fit_quadratic, x = Impact_year), linetype = "dashed", arrow=myarrow, colour = "orange", 
            size = 1, alpha = 0.95, show.legend = T)+ 
  geom_ribbon(data = new_data, 
              mapping = aes(ymin = lwr_quadratic, ymax = upr_quadratic, x = Impact_year), fill = "orange", 
              inherit.aes = FALSE, alpha = 0.2)+ 
  xlab("Year") + ylab("Total cost in the EU (in US$ 2017 value)")+ 
  scale_fill_manual()+ 
  coord_cartesian(ylim = c(0, 200)) + 
  theme(legend.position = "top")+ 
  theme_ggeffects()+ 
  theme_pubclean()
EU_total<-EU_total + ggtitle("InvaCost 2.0 - EU extrapolation") 
EU_total 
svg("EU_RR_HR_management_without_outliers_cost.svg") 
EU_total 
dev.off()