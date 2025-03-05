inlaDat <- function(mdl){
  mdlForm <- mdl$.args$formula
  mdlterms <- attr(terms(mdlForm), 'term.labels')
  resp <- gsub('\\(|\\)',"",mdlForm[2])
  
  dat <- mdl$.args$data
  temp <- NULL
  for(i in 1:length(dat)){
    if(names(dat)[[i]]==resp) {
      nr = length(dat[[i]])
      next
    }
    temp <- cbind(temp, dat[[i]])
  }
  temp <- data.frame(temp)
  names(temp) <- names(dat)[2:length(dat)]
  
  temp <- temp[,-which(colSums(is.na(temp)) > nr)]
  temp <- subset(temp, select = names(temp) %in% mdlterms)
  temp <- temp[complete.cases(temp),]
  #temp[,resp] <- dat[[resp]]
  return(temp)
}

## mdl.results - output from the model selection code run above
## orig.dat - the data frame with all the unsmoothed data (can be scaled or unscaled)
## data.points - the number of data points that you want the response plot to have
## intChar - the character used to define an interaction between two covariates. 
##            The default is an underscore ('_') because a colon (':') was reserved
##            in formulas. The underscore was what I used in the inlaGAMFunc, but I
##            think I also made that an option in the function. If there is an 
##            interaction it will be smoothed as a tensor by the GAM.
## newData - a new data frame that has values for all the covariates in the model for 
##           which you want model predictions of the response. 
## newNames - names of the rows for the the newData. For example, if there are names
##            for each cell that you want to assign to each row that you want to refer
##            in the output, you can enter that here.
## coordCols - the columns that contain any data that you would like to be included in
##             the newData output from the predict function. I'm assuming it will be
##             latitude longitude coordinates, which is why I have named it 'coordCols'.
## scaleNewData - Do you want to z-transform the new data?
## plotQuantiles - What quantiles of the original data do you want to plot?
## unscale - is the model data scale and do you want to backtransform it. Note: If data is
##           scaled for analysis, we typically don't backtransform it until we plot
##            it.
## scale.par - the mean and standard deviation for the covariates from the original 
##             that will be used to scale and unscale the data 
#mdl=spat.inla.thal.d
newDat.lincomb <- function(mdl, orig.dat, data.points, intChar="_",
                           newData=NULL, newNames=NULL, 
                           coordCols=NULL, scaleNewData=F,
                           plotQuantiles=c(0.01, 0.99),
                           unscale=F, scale.par=NULL){
  
  #extract the formula from the model results
  mdlForm <- mdl$.args$formula
  mdlDat <- inlaDat(mdl)
  terms.obj <- terms(mdlForm)
  mdl.terms <- attr(terms.obj, 'term.labels')
  mdl.terms <- mdl.terms[-grep('\\(|\\)', mdl.terms)]
  mdl.terms<-c(mdl.terms, modrando)
    
  #paran.offset <- attr(terms.obj, 'variables')[[attr(terms.obj, 'offset')+1]]
 # offsetName <- gsub(".*\\(|\\).*", "", paran.offset)[2]
  
  #determine if there are any smoothed terms and the number of knots for those terms
  smth.terms <- substr(mdl.terms, 1, nchar(mdl.terms)-2)
  nKnots <- data.frame(table(smth.terms), stringsAsFactors = F)
  nKnots <- nKnots[which(nKnots$Freq>1),]
  
  # unsmoothed terms
  nonsmth <- mdl.terms[mdl.terms %in% names(orig.dat)]
  
  #factor terms
  fact <- mdl.terms[-which(mdl.terms %in% names(orig.dat))]
  fact <- fact[-which(fact=="Intercept")]
  
  #combine smoothed and unsmoothed 
  dat.terms <- c(nonsmth, as.character(nKnots$smth.terms), fact)
  
  #deal with interaction terms
  intTerms <- dat.terms[grep(intChar, dat.terms)]
  if(length(intTerms)>0) {
    for(i in 1:length(intTerms)){
      intAct <- intTerms[i]
      temp <- unlist(strsplit(intAct, intChar))
      dat.terms <- dat.terms[-which(dat.terms==intAct)]
      dat.terms <- c(dat.terms, temp)
    }
  }
  
  lcDat <- NULL
  lcnames <- NULL
  
  #this is the loop to make the new dataset for the non-smoothed terms if 
  # no new dataset is provided
  if(is.null(newData)){
    for (i in 1:length(dat.terms)){
      
      plotvar <- dat.terms[i]
      
      # get the range of values for the term that you will want to plot eventually
      if (dat.terms[i] %in% fact){
        var.range <- unique(mdlDat[,plotvar])
      } else {  
        var.range <- seq(from = quantile(mdlDat[, plotvar], plotQuantiles[1]), 
                         to = quantile(mdlDat[, plotvar], plotQuantiles[2]),
                         length=data.points)
        if(unscale==T) var.range <- (var.range * scale.par[plotvar, 'sd']) + scale.par[plotvar, 'mean']
      }
      
      templs <- list() 
      
      #get the median value for every other term
      for(j in 1:length(dat.terms)){
        templs[[j]] <- var.range
        if (dat.terms[i] != dat.terms[j]){
          if (is.factor(mdlDat[,dat.terms[j]])){
            templs[[j]]<-myMode(mdlDat[,dat.terms[j]])
          } else {
            templs[[j]]<-quantile(mdlDat[,dat.terms[j]], probs = 0.5)
            #templs[[j]]<-mean(mdlDat[,dat.terms[j]])
          }
        }  
      }
      tempdf <- expand.grid(templs)
      lcDat <- rbind(lcDat, tempdf)
      
      #save names so you can extract these terms later in the plotting function
      lcnames <- c(lcnames, 
                   paste(dat.terms[i], 
                         sprintf('%0.2d', 1:length(var.range)), 
                         sep=''))
    }
    names(lcDat) <- dat.terms
    
    if(length(intTerms) > 0) {
      for(i in 1:length(intTerms)){
        intAct <- intTerms[i]
        temp <- unlist(strsplit(intAct, intChar))
        intDF <- expand.grid(lcDat[grep(temp[1], lcnames), temp[1]],
                             lcDat[grep(temp[2], lcnames), temp[2]])
        names(intDF) <- temp
        
        for(j in 1:length(dat.terms)){
          if (!(dat.terms[j] %in% temp)){
            if (is.factor(mdlDat[,dat.terms[j]])){
              val <- myMode(mdlDat[,dat.terms[j]])
              intDF[,dat.terms[j]] <- rep(val, nrow(intDF))
            } else {
              val <- quantile(mdlDat[,dat.terms[j]], probs = 0.5)
              intDF[,dat.terms[j]] <- rep(val, nrow(intDF))
            }
          }  
        }
        lcDat <- bind_rows(lcDat, intDF)
        
        #save names so you can extract these terms later in the plotting function
        lcnames <- c(lcnames, 
                     paste(intAct, sprintf('%0.2d', 1:nrow(intDF)), sep=''))
      }
    }
    
    new <- lcDat
    ## This is what to do if the user inputs a new data set
  } else {
    
    lcDat <- newData[names(newData) %in% dat.terms]
    if(scaleNewData==T){
      for(j in 1:ncol(lcDat)){
        varj <- names(lcDat)[j]
        lcDat[,j] <- (lcDat[,j] - scale.par[varj, 'mean'])/scale.par[varj, 'sd']
      }
    }
    
    coords <- newData[names(newData) %in% coordCols]
    new <- cbind(coords, lcDat)
    
    if(is.null(newNames)) {
      lcnames <- paste('lc', sprintf('%0.2d', 1:nrow(newData)), sep='')
    } else {
      lcnames <- newNames
    }
  }
  
  #now initialize the terms for the smoothing loop
  smth.terms <- as.character(nKnots$smth.terms)
  
  Xcr <- list()
  lcs <- list()
  smth.names <- NULL
  drop <- NULL
  
  #this is the loop to smoothed the terms based on the model
  if(length(smth.terms)>0){
    for(i in 1:length(smth.terms)){
      
      if(length(grep(intChar, smth.terms[i])) > 0) {
        nk <- sqrt(nKnots$Freq[i]+1) #the number of knots for each smoothed term
        
        smthvar <- unlist(strsplit(smth.terms[i], intChar))
        smth <- te(smthvar[1], smthvar[2], bs = "cr", k = nk, fx = TRUE)
        smther <- 'te'
        drop <- c(drop, smthvar)
        
        for(j in 1:length(smthvar)){
          smth$margin[[j]]$term <- smthvar[j]
          smth$margin[[j]]$label <- paste('s(', smthvar[j], ')', sep='')
        }
        
      } else {
        nk <- nKnots$Freq[i] + 1 #the number of knots for each smoothed term
        
        smthvar <- smth.terms[i]
        smth <- s(smthvar, bs = "cr", k = nk, fx = TRUE)
        smther <- 's'
        drop <- c(drop, smth.terms[i])
      }
      
      smth$term <- smthvar
      smth$label <- paste(smther,'(', paste(smthvar, collapse=','), ')', sep='')
      
      #This is an mgcv function that actually does the smoothing based on the 
      # number of knots that we're using
      sm <- smoothCon(smth, 
                      data = orig.dat,  
                      absorb.cons = TRUE)[[1]]
      
      #another mgcv function to predict smooths with new data
      predmat <- data.frame(PredictMat(sm, new))
      rowind <- grep(smth.terms[i], lcnames)
      Xcr[[i]] <- predmat[rowind,]
      
      nc = ncol(Xcr[[i]])
      nr = nrow(Xcr[[i]])
      xnames <- paste(paste(smthvar, collapse=intChar), sprintf('%0.2d', 1:nc), sep='')
      names(Xcr[[i]]) <- xnames
      
      ## this is actually for the linear combinations data set
      predmat <- data.frame(PredictMat(sm, lcDat))
      tempnames <- names(lcDat)
      lcDat <- cbind(lcDat, predmat)
      names(lcDat) <- c(tempnames, xnames)
    }
  }
  
  if(!is.null(drop)) lcDat <- lcDat[,-which(names(lcDat) %in% drop)]
  ## make the linear combinations for the non-smoothed terms
  lcs <- inla.make.lincombs(lcDat)
  
  ## add the variable names to the linear combinations list so that the data can
  ## be refered to later in the plotting function
  names(lcs) <- lcnames
  
  ## output the results
  output <- list(lcs, lcDat, lcnames, new)
  names(output) <- c('lincomb', 'lincombDat', 'lincombRowNames', 'newDat')
  
  return(output)
}


### this is the predict function to re-fit the inla model but this will call the
### newDat.lincomb function to make all the necessary linear combinations to make
### the response plots

## species - what species are you making predictions for. If you have output for multiple
##             species and you only want to do a few, you can specify that here.
## mdl - what is the mdl that you want use to make predictions for.
## mdl.stack  - this is the inla stack that was used to fit the inla model 
## spatialMod - Logical (TRUE/FALSE) - Was a spatial autocorrelation model
##                used or not. 
## orig.dat   -  this is the original (unsmoothed) data
## data.pts   - how many data.pts do you want in you response plot?
## newData - a new data frame that has values for all the covariates in the model for 
##           which you want model predictions of the response. 
## newNames - names of the rows for the the newData. For example, if there are names
##            for each cell that you want to assign to each row that you want to refer
##            in the output, you can enter that here.
## coordCols - the columns that contain any data that you would like to be included in
##             the newData output from the predict function. I'm assuming it will be
##             latitude longitude coordinates, which is why I have named it 'coordCols'.
## scaleNewData - Do you want to z-transform the new data?
## plotQuantiles - What quantiles of the original data do you want to plot?
## unscale - is the model data scale and do you want to backtransform it. Note: If data is
##           scaled for analysis, we typically don't backtransform it until we plot
##            it.
## scale.par - the mean and standard deviation for the covariates from the original 
##             that will be used to scale and unscale the data 


inlaPredict <- function(mdl, mdl.stack, spatialMod,
                        orig.dat, data.pts, plotQuantiles=c(0.01, 0.99),
                        newData=NULL, newNames=NULL, coordCols=NULL,
                        scaleNewData=NULL, unscale=F, scale.par=NULL) {
  
  
  fam <- mdl$.args$family
  #extract the model formula and determine the covariates that are included in the model
  mdlForm <- mdl$.args$formula 
  
  # run the newDat.lincomb function to make the linear combinations that will
  # be used to make predictions for each covariate in the model
  new <- newDat.lincomb(mdl=mdl, orig.dat=orig.dat, data.points=data.pts,
                        newData=newData, newNames=newNames, 
                        scaleNewData = scaleNewData, coordCols=coordCols,
                        scale.par=scale.par, unscale=unscale,
                        plotQuantiles=plotQuantiles)
  
  lcs <- new$lincomb
  newDat <- new$newDat
  newRowNames <- new$lincombRowNames
  
  # run the model based on whether it includes a spatial covariance structure or not
  if(spatialMod==T){
    mod <- inla(mdlForm,
                family = fam,
                data=inla.stack.data(mdl.stack),
                lincomb = lcs,
                #verbose=T,
                #control.inla = list(lincomb.derived.only=F),
                control.compute = list(config=T),
                control.predictor = list(A = inla.stack.A(mdl.stack),
                                         compute=TRUE, 
                                         quantiles = c(0.025, 0.975)))
  } else {
    mod <-  inla(mdlForm,
                 family = fam,
                 data=mdlDat,
                 lincomb = lcs,
                 #control.inla=list(lincomb.derived.only=F),
                 control.compute = list(config=T),
                 control.predictor = list(A = inla.stack.A(mdl.stack),
                                          compute=TRUE, 
                                          quantiles = c(0.025, 0.975)))
  }
  
  output <- list(newDat, newRowNames, mod)
  names(output) <- c('newData', 'newRowNames', 'model')
  return(output)
}



## Function to the plot the results from the inlaPredict function as a 
## line plot (no interaction) or a contour plot (interaction)

## inlaPrediction - results from inlaPrediction
## plotVar - the covariate that you want to plot
## orig.dat   -  this is the original (unsmoothed) data
## unscale - is the model data scaled and do you want to backtransform it. 
##           Note: If data is scaled for analysis, we typically don't 
##           backtransform it until we plot it.
## scale.par - the mean and standard deviation for the covariates from the original 
##             that will be used to scale and unscale the data 
## colPalette - the color palette to use to make the plot
## legend - logical (TRUE/FALSE) - do you want to include the legend with the plot
## xlabl - what do you want the text of the x-axis label to be?
## ylabl - what do you want the text of the y-axis label to be?
## xlimits - Do you want to limit the range of the xaxis? The default is the range
##            of the data.
## ylimits - Do you want to limit the range of the yaxis? The default is the range
##            of the data.
## pritPlotobject - logical (TRUE/FALSE) - Do you want to print the object 
##                    directly after making it or do you want to store it?

# inlaPrediction = sharksPred
# plotVar = 'Fishing'
# orig.dat = Sharks06
# unscale=F
# scale.par=NULL
# xlabl='Fishing'
# ylabl='probability of presence'
# intChar='_'
# colPalette = 'Set2'
# legend=F
# xlimits=NULL
# ylimits=NULL
# printplotobject=FALSE


inlaPlot <- function(inlaPrediction, plotVar, orig.dat, intChar='_', 
                     unscale=F, scale.par=NULL,
                     colPalette = 'Set2', legend=TRUE,
                     xlabl, ylabl, xlimits=NULL, ylimits=NULL, 
                     printplotobject=FALSE){
  
  #only use one color if there is only one species but use multiple colors if there
  # are multiple species
  plotColors <- 'black'
  
  # extract the model output from the inlaPrediction function
  mod <- inlaPrediction$model
  modlink <- mod$misc$linkfunctions$names
  newDat <- inlaPrediction$newData
  newRowNames <- inlaPrediction$newRowNames
  
  ### now extract the linear combination results, which are the predictions that we
  ### will use to make the response plots
  keep <- c('mean', '0.025quant', '0.975quant')
  
  # make the plot data frame based on the predictions from the linear combinations
  plotDat.temp <- mod$summary.lincomb.derived[,keep]
  # we will only keep the rows that include the names of our variable (this is why
  # I named every row of the linear combinations list based on the covariate)
  #rowind <- grep(plotVar, row.names(plotDat.temp))
  newRowInd <- grep(plotVar, newRowNames)
  plotDat.temp <- plotDat.temp[newRowInd,]
  
  if(length(grep(intChar, plotVar))>0) {
    # need to split and extract values for both variables 
    # if it is an interaction
    intVars <- unlist(strsplit(plotVar, intChar))
    
    xvals <- newDat[newRowInd, intVars]
    
    intPlot <- T
  } else {
    newColInd <- grep(plotVar, names(newDat))
    xvals <- newDat[newRowInd, newColInd]
    intPlot <- F
  }
  
  plotDat.temp <- cbind(xvals, plotDat.temp)
  
  if(intPlot==T){
    names(plotDat.temp) <- c(intVars, 'modPred', '0.025quant', '0.975quant')
    
    for(j in 1:length(intVars)){
      ## backtransform values for both axes so they are on the original 
      ## scale of the data
      if(unscale == T) {
        plotDat.temp[,intVars[j]] <- (plotDat.temp[,intVars[j]] * scale.par[intVars[j], 'sd']) + scale.par[intVars[j], 'mean']
      } else {
        plotDat.temp$xaxis <- plotDat.temp[,intVars[j]]
      }
    }
  } else {
    names(plotDat.temp) <- c(plotVar, 'modPred', '0.025quant', '0.975quant')
    
    ## backtransform the x-axis so that it is on the original scale of the data
    if(unscale == T) {
      plotDat.temp$xaxis <- (plotDat.temp[,plotVar] * scale.par[plotVar, 'sd']) + scale.par[plotVar, 'mean']
    } else {
      plotDat.temp$xaxis <- plotDat.temp[,plotVar]
    }
  }
  
  ## transform the predictions to the correct scale 
  if(modlink=='logit'){
    plotDat.temp <- within(plotDat.temp, {
      predictedProb <- plogis(modPred)
      lowCredInt <- plogis(`0.025quant`)
      upCredInt <- plogis(`0.975quant`)
    })
  } else {
    if(modlink=='log'){
      plotDat.temp <- within(plotDat.temp, {
        predictedProb <- exp(modPred)
        lowCredInt <- exp(`0.025quant`)
        upCredInt <- exp(`0.975quant`)
      })
    }
  }
  
  
  ## reduce the data to just the columns that we need so that we can combine all 
  ## species together into a single data frame

  if(intPlot==F) {
    keep <- c('xaxis', 'predictedProb', 'lowCredInt', 'upCredInt')
  } else {
    keep <- c(intVars, 'predictedProb', 'lowCredInt', 'upCredInt')
  }
  plotDat.red <- plotDat.temp[,keep]
  
 ## combine this species with all other species
 plotDat <- plotDat.red 
 
 # determine if plotvar is a factor or not
 vartype <- ifelse(is.factor(orig.dat[,plotVar]), 'fact', 'cont')
  
  if(intPlot==T){
    require(scales)
    
    minscale <- 0
    maxscale <- 1
    medscale <- 0.5
    
    p <- ggplot(data=plotDat, aes_string(x=intVars[1], y=intVars[2])) +
      geom_tile(aes(fill=predictedProb)) + labs(x=xlabl, y=ylabl) +
      scale_fill_gradientn(colours = c("cornflowerblue", "white","orange"),
                           values = rescale(c(minscale, medscale, maxscale)), 
                           breaks=c(round(minscale,3), round(medscale,3), round(maxscale,3)),
                           name="Predicted Probability", guide = "colorbar", limits= c(minscale, round(maxscale,3))) +   
      scale_color_gradientn(colours = c("cornflowerblue", "white","orange"),
                            values = rescale(c(minscale, medscale, maxscale)), 
                            breaks=c(round(minscale, 3), round(medscale, 3), round(maxscale, 3)),
                            name="Predicted Probility", guide = "colorbar", limits= c(minscale, round(maxscale,3))) +
      scale_x_continuous(expand = c(0,0))+
      scale_y_continuous(expand = c(0,0)) +
      theme(legend.key.size = unit(0.25, "cm"), axis.title = element_text(size=12), axis.text = element_text(size=11)) 
    
  } else {
    
    if(vartype=='cont'){
      # set the plot limits
      if (is.null(xlimits)) xlimits=c(min(plotDat$xaxis), max(plotDat$xaxis))
      if (is.null(ylimits)) ylimits=c(min(plotDat$lowCredInt), max(plotDat$upCredInt))
      
      # set the axis labels
      if (is.null(xlabl)) xlabl=plotVar
      
      ## make the plot
      p <- ggplot(plotDat, aes(x=xaxis, y=predictedProb)) +
        geom_line(size=1) +
        geom_ribbon(aes(ymin = lowCredInt, ymax = upCredInt), alpha = .2) +
        #scale_colour_manual("Taxon", values=plotColors) +
        #scale_fill_manual('Taxon', values=plotColors) +
        #ylim(0, ymax) +
        scale_x_continuous(limits=xlimits) + scale_y_continuous(limits=ylimits) +
        theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        xlab(xlabl)
    } else {
      # set the plot limits
      plotDat$xaxis <- factor(plotDat$xaxis)
      #if (is.null(xlimits)) xlimits=factor()
      if (is.null(ylimits)) ylimits=c(min(plotDat$lowCredInt), max(plotDat$upCredInt))
      
      # set the axis labels
      if (is.null(xlabl)) xlabl=plotVar
      
      ## make the plot
      p <- ggplot(plotDat, aes(x=xaxis, y=predictedProb)) +
        geom_point(size=1) +
        geom_errorbar(aes(ymin=lowCredInt, ymax=upCredInt, width=0.2)) +
        scale_y_continuous(limits=ylimits) +
        theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        xlab(xlabl)
    }
    
    if(!is.null(ylabl)) p <- p + ylab(ylabl)
  }
  
  ## keep or remove the legend depending on what you want
  if(legend==FALSE) p <- p + theme(legend.position = 'none')
  
  ## print or don't print, it's up to you
  if (printplotobject == TRUE) {
    print(p)
  } else {
    return(p)
  }
}
