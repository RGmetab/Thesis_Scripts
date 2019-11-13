#Title:  Metabolite Set Enrichment Analysis (MSEA)
#Author: Rudi Grosman - UoL NMR Centre
#Usage: pathwayFisher(Metabolites, PathTable)


pathwayFisher <- function(Metabolites, PathTable, ConvTable=NULL,withEASE=T, PvalOnly=F, ShowMatches=F){
  
  if (is.data.frame(PathTable)){
    PathTable <- PathTable
  } else if (is.character(PathTable)){
    PathTable <- read.csv(PathTable, stringsAsFactors = F, sep = '\t')
  } else if (is.null(PathTable)){
    return(stop('No Pathway Table supplied. Aborting script.'))
  }
  
  BGlist <- apply(PathTable,2,function (x) length(grep('C[0-9]{5}',x)))
  BGtotal <- sum(apply(PathTable,2,function (x) length(grep('C[0-9]{5}',x))))
  MetNo <- length(Metabolites)
  MetMatch <- apply(PathTable, 2, function (x) length(intersect(Metabolites, x)))
  BGlist <- BGlist[names(MetMatch)]
  MetMatch <- MetMatch[MetMatch != 0]
  ResultDF <- data.frame(PathCode=numeric(0),
                         Low.Conf.Int=numeric(0),
                         High.Conf.Int=numeric(0),
                         P.val=numeric(0))
  
    for (i in 1:length(MetMatch)){
    if(withEASE){
      ContTable <- matrix(c(MetMatch[i]-1,BGlist[i],
                            MetNo-MetMatch[i],BGtotal-BGlist[i]),
                          byrow=T, nrow=2)
    } else {
      ContTable <- matrix(c(MetMatch[i],BGlist[i],
                            MetNo-MetMatch[i],BGtotal-BGlist[i]),
                          byrow=T, nrow=2)
    }
    
    rownames(ContTable) <- c("In Pathway","Not in Pathway")
    colnames(ContTable) <- c("Query Metabolites","Pathway Metabolites")
    FisherResult <- fisher.test(ContTable, alternative='greater')
    tmpOut <- c(names(MetMatch)[i], FisherResult$conf.int, FisherResult$p.value)
    ResultDF[i,] <- tmpOut
  }
  
  ResultDF$P.val <- as.numeric(ResultDF$P.val)
  ResultDF$adjPval <- p.adjust(ResultDF$P.val, method='BH')
  
  ResultDF$Significance <- 'No'
  ResultDF$Significance[ResultDF$adjPval<0.05] <- '*'
  ResultDF$Significance[ResultDF$adjPval<0.01] <- '**'
  ResultDF$Significance[ResultDF$adjPval<0.001] <- '***'
  

  
  if (is.null(ConvTable)){
    return(ResultDF)
  } else if (is.character(ConvTable)){
    ConvTable <- read.csv(ConvTable, sep='\t', stringsAsFactors=F, header=F)
    ResultDF$PathName <- ConvTable[match(ResultDF$PathCode, ConvTable[,1]),2]

  } else if (is.data.frame(ConvTable)){
    ResultDF$PathName <- ConvTable[match(ResultDF$PathCode, ConvTable[,1]),2]
  }
  
  if(ShowMatches){
    MatchList <- sapply(ResultDF$PathCode, function(x) paste(sort(intersect(PathTable[,x], Metabolites)),collapse = ', '))
    ResultDF$Matches <- MatchList
    ResultDF$Hits <- MetMatch
    ResultDF$Total <- sapply(ResultDF$PathCode, function(x) length(grep('C[0-9]{5}',(PathTable[,x]))))
    ResultDF$Percentage <- round((ResultDF$Hits/ResultDF$Total)*100,2)
    #This is work in progress was applied as it was explained in DAVID, not
    #sure if it is good and exactly what kind of information it contributes.
    ResultDF$FoldEnrich <- round((ResultDF$Hits/MetNo)/
                                   (ResultDF$Total/BGtotal),4)
  }
  
  return(ResultDF)
  
}
