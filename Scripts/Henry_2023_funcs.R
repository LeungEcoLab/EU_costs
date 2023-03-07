
##Script to group all the necessary functions to run the scripton the supercomputer

# do.df.extrap creates the dataframe that will be used for estimating unknown costs in EU. We need to make sure that the combination of species/sector/type of cost match between this df and the fitting df
# get_ratio gives the min and max ration for each predictors
# bound_ratio truncates the ratios values if it's outside the range of the data
# scale_fit_fun is the main function to optimised the parameters gamma1, calculates metrics of fit, and estimate unknown costs
# aic.fit is a custom-made AIC function to calculate the metrics of fit for each model (to later do model averaging)


do_df_extrap <- function(df.sp, df.fit){
  #Keep only species/sector/type of cost that are NOT present in the invacost database
  df.sp$comb2 <-  paste(df.sp$Species, df.sp$Impacted_sector, df.sp$Type_of_cost_merged)
  df.sp$comb <- paste(df.sp$Official_country, df.sp$Species, df.sp$Impacted_sector, df.sp$Type_of_cost_merged)
  df.fit$comb <- paste(df.fit$Official_country, df.fit$Species, df.fit$Impacted_sector, df.fit$Type_of_cost_merged)
  
  same <- which(df.sp$comb %in% df.fit$comb)
  df.extrap <- df.sp[-same,] #remove all rows of data that are present in df
  df.extrap <- df.extrap[!df.extrap$comb == "San Marino Myocastor coypus Agriculture Damage",] #remove this line because GDP of agriculture in san marino is 0
  
  #Remove countries not in the EU
  extra.c <- c("Andorra", "San Marino", "Monaco", "Liechtenstein") #NOT IN EU
  df.extrap <- df.extrap[!(df.extrap$Official_country %in% extra.c),]
  
  df.fit$type <- "invacost" #to keep track of waht is what
  df.extrap$type <- "extrapol" #to keep track of what is what
  
  df.extrap$Cost_estimate_per_year_2017_USD_exchange_rate <- NA
  df.extrap$cost_mil <- NA
  df.extrap$Cost_ID <- paste("extrap", "_", 1:nrow(df.extrap), sep="")
  df.extrap$Spatial_scale <- unique(df.fit$Spatial_scale)
  
  #reorder columns to use rbind
  df.extrap <- df.extrap[names(df.fit)]
  
  ##We want to keep only the extrapolation part of the dataframe, because the "invacost" df will be the bootstrap
  return(df.extrap)
}


get_ratio <- function(dff, pred.list){
  
  ln <- pred.list #replace pred.list by ln (list.name) so that we don't overwrite it
  dff$nosector <- 1
  
  pos <- which(dff$type == "invacost")
  x <- dff[pos,] #invacost data
  y <- dff[-pos,] #extrapolation data
  
  sector <- c(Health = "Health_expend_perc_gdp", `Authorities-Stakeholders` = "nosector",
              Agriculture = "Agri_forest_fish_perc_gdp", Fishery= "Agri_forest_fish_perc_gdp",
              `Public and social welfare` = "nosector", Environment= "nosector",
              Diverse= "nosector", Forestry= "Agri_forest_fish_perc_gdp")
  
  mxr=matrix(NA, nrow=nrow(y), ncol=length(ln)*2)
  
  for(i in 1:nrow(y)){
    if(is.null(ln)){
      return(mxr)
      
    } else {
      if("Impacted_sector" %in% pred.list == T){
        ln[ln == "Impacted_sector"] <- sector[y$Impacted_sector[i]] #if Impacted_sector is present, change it to the specific sector of interest
      }
      
      for(j in 1:length(ln)){
        mxr[i, j] <- min(y[i, ln[j]]/x[,ln[j]])
        mxr[i, j+length(ln)] <- max(y[i, ln[j]]/x[,ln[j]])
      }
    }
  }
  mxr2 <- as.data.frame(mxr)
  colnames(mxr2) <- c(paste(rep("min", length(pred.list)), pred.list, sep="_"), paste(rep("max", length(pred.list)), pred.list, sep="_"))
  return(mxr2)
}

bound_ratio <- function(v, min.r, max.r){
  v[v < min.r] <- min.r
  v[v > max.r] <- max.r
  return(v)
}


scale_fit_fun <- function(gamma1, df, pred.list, extrap = F, fit.test =F){ #gamma1= vector of parameters to optimize
  
  ln <- pred.list #replace pred.list by ln (list.name) so that we don't overwrite it
  
  if(extrap == F){ #default (extrap = F) is when we want to estimate the gamma parameters (with optim)
    if(length(gamma1) == 1){
      names(gamma1) <- "sd"
    }
    full_gamma[names(gamma1)]=gamma1
    
    scale_fun <- function(dff, fit.test){
      #Add a column "nosector" when we don't have specific sector values (e.g. "authorities-stakeholders", "Public and social welfare", ...)
      dff$nosector <- 1
      
      pos <- which(dff$type == "invacost")
      x <- dff[pos,] #invacost data
      y <- dff[-pos,] #extrapolation data
      
      sector <- c(Health = "Health_expend_perc_gdp", `Authorities-Stakeholders` = "nosector", 
                  Agriculture = "Agri_forest_fish_perc_gdp", Fishery= "Agri_forest_fish_perc_gdp",
                  `Public and social welfare` = "nosector", Environment= "nosector", 
                  Diverse= "nosector", Forestry= "Agri_forest_fish_perc_gdp")
      
      if(nrow(y) == 1){
        y$cost_pred <- NA 
        
        if(fit.test ==T){ #return the full dataset to calculate metrics of fit
          return(y)
        } else {
          ll <- sum(log(dnorm(y$cost_mil - y$cost_pred, sd= full_gamma["sd"]))) #calculate log-likelihood sd= gamma1[length(gamma1)]
          return(-ll) #return -ll because optim finds the minimum
        }
        
      } else {
  
        for(i in 1:nrow(y)){
          if("Impacted_sector" %in% pred.list == T){
            ln[ln == "Impacted_sector"] <- sector[y$Impacted_sector[i]] #assign the proper impacted_sector
          }
          pred1 <- x$cost_mil[-i] * (y$GDP_2019[i] / x$GDP_2019[-i])^full_gamma["GDP_2019"] * 
            (y$Population_2019[i] / x$Population_2019[-i])^full_gamma["Population_2019"] * 
            (y$Surface_area_km2[i] / x$Surface_area_km2[-i])^full_gamma["Surface_area_km2"] * 
            (y[i, sector[y$Impacted_sector[i]]] / x[,sector[x$Impacted_sector[1]]][-i])^full_gamma["Impacted_sector"] #parameter fit to sector ##Not sure why [1]??? 
          y$cost_pred[i] <- mean(pred1)
        }
        
        if(fit.test ==T){
          return(y)
          
        } else {
          ll <- sum(log(dnorm(y$cost_mil - y$cost_pred, sd= full_gamma["sd"])))
          return(-ll) 
        }
      }
    }
    
    df.by <- by(df, df$comb2, scale_fun, fit.test= fit.test) #apply the function "scale_fun" to each "comb2" column. comb2 is the combination of species/sector/type of cost
    
    if(fit.test == T){
      df.full <- do.call(rbind, df.by) 
      return(df.full)
      
    } else {
      measure.fit <- sum(df.by, na.rm = T) #sum of the log-likelihood
      return(measure.fit)
    }
    
  } else { #When extrap = T we want to obtain an estimation of costs
    
    scale_fun <- function(dff){
      
      ln <- pred.list
      #Add a column "nosector" when we don't have specific sector values (e.g. "authorities-stakeholders", "Public and social welfare", ...)
      dff$nosector <- 1
      
      pos <- which(dff$type == "invacost")
      x <- dff[pos,] #invacost data
      y <- dff[-pos,] #extrapolation data
      
      sector <- c(Health = "Health_expend_perc_gdp", `Authorities-Stakeholders` = "nosector", 
                  Agriculture = "Agri_forest_fish_perc_gdp", Fishery= "Agri_forest_fish_perc_gdp",
                  `Public and social welfare` = "nosector", Environment= "nosector", 
                  Diverse= "nosector", Forestry= "Agri_forest_fish_perc_gdp")
      
      v <- list() #v will contains the ratios
      for(i in 1:nrow(y)){
        if(is.null(ln)){
          pred1 <- x$cost_mil
          y$cost_pred[i] <- mean(pred1)
        } else {
          
          if("Impacted_sector" %in% pred.list == T){
            ln[ln == "Impacted_sector"] <- sector[y$Impacted_sector[i]]
          }
        
        for(j in 1:length(ln)){
          v[[j]] <- bound_ratio(y[i, ln[j]]/x[,ln[j]], min.r= min(ratio.fit[,paste("min", pred.list[j], sep="_")]), max.r= max(ratio.fit[,paste("max", pred.list[j], sep="_")])) #use pred.list instead of ln because the name for sector is "sector" not the specific one
        }
          
          vm <- matrix(unlist(v), ncol=length(ln)) #create a matrix where each column is a predictor (so that we can multiply across rows)
          colnames(vm) <- ln
  
          prod.vm <-apply(vm, MARGIN= 1, function(x){ prod(x^gamma1[-length(gamma1)])}) #multiply each row and apply the gamma1 exp
          
          pred1 <- x$cost_mil * prod.vm
          y$cost_pred[i] <- mean(pred1)
          
        }
      }
      
      return(y) #Here we want to keep the entire dataframe to compare/analyse the results
    }
    
    df.by <- by(df, df$comb2, scale_fun)
    df.full <- do.call(rbind, df.by) 
    
    return(df.full)
  }
}

aic.fit <- function(x, gamma1){
  LL <- sum(log(dnorm(x$cost_mil - x$cost_pred, sd= gamma1[length(gamma1)])), na.rm=T)
  aic <- -2*(LL) + 2*(length(gamma1)-1) #-1 because one of the parameters is the standard deviation
  return(aic)
}


