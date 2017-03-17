###################################################################
#   Complexity masures for regression problemas                   #
#   Proposed by Ana Carolina Lorena and Ivan Costa                #
#   Implemented by Aron Ifanger Maciel and Ana Carolina Lorena    #
###################################################################

library(igraph)
library(FNN)

###################################################################
#   CALLS                                                         #
#       Normalize(dataset)                                        #
#       FormatDataset(dataset, output)                            #
#       MaxPosition(array)                                        #
#       MinPosition(array)                                        #
#       spearman_from_rank(array)                                 #
#       ExamplesRemovedNumber(x,y,minCorrelation)                 #
#                                                                 #
#       C1(dataset)                                               #
#       C2(dataset)                                               #
#       C3(dataset)                                               #
#       C4(dataset)                                               #
#                                                                 #
#       L1(dataset)                                               #
#       L2(dataset)                                               #
#                                                                 #
#       S1(dataset)                                               #
#       S2(dataset)                                               #
#       S3(dataset)                                               #
#                                                                 #
#       L3(dataset)                                               #
#       S4(dataset)                                               #
#       T2(dataset)                                               #
###################################################################

###################################################################
#   First group of functions                                      #
#   These are auxiliary functions                                 #
###################################################################

#-----------------------------------------------------------------#
#   01 - Normalize                                                #
#       Normalize the columns of a dataset (within [0,1])         #
#-----------------------------------------------------------------#
Normalize = function(dataset) {
    dataset = as.matrix(dataset)
    numberColumn = ncol(dataset)

    for (column in 1:numberColumn)
        dataset[,column] = (dataset[,column] - min(dataset[,column])) /
                           (max(dataset[,column]) - min(dataset[,column]))

    dataset
}

#-----------------------------------------------------------------#
#   02 - FormatDataset                                            #
#       Format the entry dataset                                  #
#-----------------------------------------------------------------#
FormatDataset = function(dataset, output){
    
    dataset = as.matrix(dataset)
    numberColumn = ncol(dataset)
  
    if(!is.null(output)){
        input = dataset
    } else {
        input  = as.matrix(dataset[,-numberColumn])
        output = as.matrix(dataset[,numberColumn])
        numberColumn = ncol(input)
    }

    list(input = Normalize(input), output = Normalize(output), 
        numberColumn = numberColumn, numberRows = nrow(input))
}

#-----------------------------------------------------------------#
#   03 - MaxPosition                                              #
#       The larger element index                                  #
#-----------------------------------------------------------------#
MaxPosition = function(array) order(-array)[1]

#-----------------------------------------------------------------#
#   04 - MinPosition                                              #
#       The smaller element index                                  #
#-----------------------------------------------------------------#
MinPosition = function(array) order(array)[1]

#-----------------------------------------------------------------#
#   05 - spearman_from_rank                                       #
#       Computes the Spearman correlation between the differences #
#       of two ranks (with no ties)                               #
#-----------------------------------------------------------------#
spearman_from_rank = function(rank){
  size=length(rank)
  results=1-6*sum(rank^2)/(size^3-size)
  results
}
    
#-----------------------------------------------------------------#
#   06 - ExamplesRemovedNumber                                    #
#       Remove examples from a dataset until a specific           #
#       correlation to the output is achived                      #
#-----------------------------------------------------------------#
ExamplesRemovedNumber = function(x,y,minCorrelation)
{
  
  numberRows = length(x)
  if(numberRows == length(y))
  {
    remainingRows = numberRows
    maxPosition = 0
    xorder=rank(x)
    yorder=rank(y)
    diff=xorder-yorder
    
    correlation = spearman_from_rank(diff)
    
    if(correlation < 0){
      yorder=rank(-y)
      diff=xorder-yorder
      correlation = spearman_from_rank(diff)
    }
    
    while(abs(correlation) < minCorrelation && !is.na(correlation))
    {
      
      maxPosition = which.max(abs(diff))
      
      diff=diff +
        ((yorder>yorder[maxPosition]) -
           (xorder>xorder[maxPosition]))
      
      yorder=yorder[-maxPosition]
      xorder=xorder[-maxPosition] 
      diff=diff[-maxPosition]
      remainingRows = remainingRows - 1
      correlation = spearman_from_rank(diff)
      
      if(is.na(correlation))
        correlation
    }
    
    (numberRows-remainingRows)/numberRows
  } else NA
}

###################################################################
#   Second group of functions                                     #
#   Complexity Measures                                           #
###################################################################

#-----------------------------------------------------------------#
#   07 - C1 - Maximum Feature Correlation to the Output           #
#       First one has to calculate the absolute value of the Spe  #
#       arman correlation between each feature and the outputs.   #
#       The absolute value is taken because both extremes of the  #
#       correlation measure, which range between [-1; 1], repre   #
#       sent a strong correlation, either inverse or direct. C1   #
#       is given by the maximum correlation value obtained among  #
#       all the features.                                         #
#-----------------------------------------------------------------#
C1 = function(dataset, output = NULL) {   
  
  formatedDataset = FormatDataset(dataset, output)
  input = formatedDataset$input
  output = formatedDataset$output
  numberColumn = formatedDataset$numberColumn
  
  correlations = array(0,numberColumn)
  
  for (column in 1:numberColumn){
    correlation = cor(output, input[,column], method = "spearman")
    
    if(!is.na(correlation))
      correlations[column] = abs(correlation)
  }
  
  max(correlations)
}

#-----------------------------------------------------------------#
#   08 - C2 Average Feature Correlation to the Output             #
#       Similar to C1, but computes the average of all the corre  #
#       lations, as opposed to just taking the maximum value      #
#       among them. Therefore, it is a measure of the relation    #
#       ship of all the features to the outputs, even though the  #
#       correlations are calculated individually for each feature.#
#       Higher values of C2 indicate simpler problems.            #
#-----------------------------------------------------------------#
C2 = function(dataset, output = NULL)
{   
    formatedDataset = FormatDataset(dataset, output)
    input = formatedDataset$input
    output = formatedDataset$output
    numberColumn = formatedDataset$numberColumn

    correlations = array(0,numberColumn)
    
    for (column in 1:numberColumn)
        correlations[column] = abs(cor(output, input[,column], 
            method = "spearman"))

    naRemove = !is.na(correlations)

    mean(correlations[naRemove])
}

#-----------------------------------------------------------------#
#   09 - C3 Individual Feature Efficiency                         #
#       Calculates the number of examples that must be removed    #
#       from the dataset until a high correlation value to the    #
#       output is achieved, divided by the total number of exam   #
#       ples. This computation is done for each individual featu  #
#       re and C3 returns the minimum of the values found, corres #
#       ponding to the feature more related to the output. Lower  #
#       values of C3 indicate simpler problems.                   #
#-----------------------------------------------------------------#
C3 = function(dataset, output = NULL, minCorrelation = 0.9)
{
    formatedDataset = FormatDataset(dataset, output)
    input = formatedDataset$input
    output = formatedDataset$output
    numberColumn = formatedDataset$numberColumn

    correlations = array(1,numberColumn)
    
    for (column in 1:numberColumn)
        correlations[column] = ExamplesRemovedNumber(output, input[,column], minCorrelation)
    
    naRemove = !is.na(correlations)
    
    min(correlations[naRemove])
}

#-----------------------------------------------------------------#
#   10 - C4 Collective Feature Efficiency                         #
#       The feature most correlated to the output is identified   #
#       and all examples with a small residual value after a      #
#       linear fit are removed. The next most correlated          #
#       feature to the remaining data is found and the previous   #
#       process is repeated until all features have been analyzed #
#       or no example remains (different from C3, after analyzing #
#       one feature, the examples that were removed are disregar  #
#       ded). The ratio of examples for which a linear fit with   #
#       small residual was not achieved is returned. Lower values #
#       indicate simpler problems.                                #
#-----------------------------------------------------------------#
C4 = function(dataset, output = NULL, minResidual = 0.1){
    
    formatedDataset = FormatDataset(dataset, output)
    input = formatedDataset$input
    output = formatedDataset$output
    numberColumn = formatedDataset$numberColumn
    numberRows = formatedDataset$numberRows

    stop = FALSE
    looked = 0
    
    while(!stop){
        
        correlations = cor(output,input,
                             method="spearman")
        indexMostCorrelated = MaxPosition(abs(correlations))
        if(!is.na(indexMostCorrelated)){
            
            looked = looked + 1
            linearModel = lm(output~input[,indexMostCorrelated])
            indexRemove = abs(linearModel$residuals) > minResidual
            input = input[indexRemove,]
            output = output[indexRemove]
            if(sum(!indexRemove) == length(indexRemove) | length(output)==1 |  looked == numberColumn)
                stop = TRUE
        }
    }
    if(length(output) == 1)
       0
    else
       length(output)/numberRows
}

#-----------------------------------------------------------------#
#   11 - L1 Distance of the Data Itens to the Linear Function     #
#       The sum of the absolute values of the residues of a multi #
#       ple linear regressor. Lower values indicate simpler pro   #
#       blems, which can be fit by a linear function.             #
#-----------------------------------------------------------------#
L1 = function(dataset, output = NULL){
    
    formatedDataset = FormatDataset(dataset, output)
    input = formatedDataset$input
    output = formatedDataset$output
    naRemove = !is.na(cor(output, input, method = "spearman"))
    numberColumns = formatedDataset$numberColumn
    if(numberColumns > 1)
        linearModel = lm(output~input[,naRemove])
    else 
        linearModel = lm(output~input[])
    mean(abs(linearModel$residuals))
}

#-----------------------------------------------------------------#
#   12 - L2 Average Error of Linear Regressor                     #
#       L2 sums the square of the residuals from a multiple       #
#       linear regression. Smaller values indicate simpler        #
#       (linear) problems.                                        #
#-----------------------------------------------------------------#
L2 = function(dataset, output = NULL){
    
    formatedDataset = FormatDataset(dataset, output)
    input = formatedDataset$input
    output = formatedDataset$output

    naRemove = !is.na(cor(output, input, method = "spearman"))

    numberColumns = formatedDataset$numberColumn
    if(numberColumns > 1)
      linearModel = lm(output~input[,naRemove])
    else 
      linearModel = lm(output~input[])
        
    mean(linearModel$residuals^2)
}

#-----------------------------------------------------------------#
#   13 - S1 Output Distribution                                   #
#       First a Minimum Spanning Tree (MST) is generated from     #
#       data. Herewith, each data item corresponds to a vertex of #
#       the graph. The edges are weighted according to the        #
#       Euclidean distance between the examples. The MST will     #
#       greedily connect examples closer to each other. Next S1   #
#       monitors whether the examples joined in the MST have      #
#       close output values. Lower values indicate simpler        #
#       problems, where the outputs of close examples in the      #
#       input space are also next to each other.                  #
#-----------------------------------------------------------------#
S1 = function(dataset, output = NULL){
    
    formatedDataset = FormatDataset(dataset, output)
    input = formatedDataset$input
    output = formatedDataset$output
    numberRows = formatedDataset$numberRows
    
    fullGraph = graph.full(numberRows, directed = FALSE, loops = FALSE)
    E(fullGraph)$weight= dist(input, method = "euclidian")
    mst = minimum.spanning.tree(fullGraph, algorithm = "prim")
    edgelist = get.edgelist(mst);
    mean(abs(output[edgelist[,1]] - output[edgelist[,2]]))
}

#-----------------------------------------------------------------#
#   14 - S2 Input Distribution                                    #
#       S2 first orders the data points according to their output #
#       values yi and then estimates the distance between pairs of# 
#       examples that are neighbors in the obtained ordering. S2  #
#       complements S1 by measuring how similar in the input space# 
#       are data items with close outputs. Lower values indicate  # 
#       simpler problems.                                         #
#-----------------------------------------------------------------#
S2  = function(dataset, output = NULL){
    
    getDistances = function(dataset, numberRows, numberColumns){

        distances = array(0,numberRows-1)
        for(line in 2:numberRows){
            if(numberColumns > 1)
               distances[line-1] = dist(dataset[(line-1):line,])
            else
               distances[line-1] = dist(dataset[(line-1):line])
          } 
        distances
    }

    formatedDataset = FormatDataset(dataset, output)
    input = formatedDataset$input
    output = formatedDataset$output
    numberRows = formatedDataset$numberRows
    numberColumns = formatedDataset$numberColumn
   
    order = order(output)
    if(numberColumns > 1)
        distances = getDistances(input[order,], numberRows,numberColumns)
    else
        distances = getDistances(input[order], numberRows,numberColumns)
    
    mean(distances)
}

#-----------------------------------------------------------------#
#   15 - S3 Error of a Nearest Neighbor regressor                 #
#       S3 calculates the mean squared error of a nearest neighbor#
#       regressor, using leave-one-out. It measures how the       #
#       examples are close together and high values imply that    #
#       there are many gaps in the input space. Lower values      #
#       indicate simpler problems.                                #
#-----------------------------------------------------------------#
S3  = function(dataset, output = NULL){

    formatedDataset = FormatDataset(dataset, output)
    input = formatedDataset$input
    output = formatedDataset$output
    numberRows = formatedDataset$numberRows
    predictions = matrix(0,numberRows)
    distances = as.matrix(dist(input, method = "euclidian",diag=TRUE,upper=TRUE))

    diag(distances) <- Inf
    for(line in 1:numberRows){
        predictions[line] = output[MinPosition(distances[line,,drop=FALSE])]
    }
    mean((predictions-output)^2)
}

#-----------------------------------------------------------------#
#   16 - L3 Non-linearity of a Linear regressor                   #
#       L3 selects pairs of examples with close outputs (with low #
#       |yi − yj|) and creates a new random test point by randomly#
#       interpolating them. It is based on the Non-linearity of a #
#       Linear Classifier (L3) measure from [2]. A linear         #
#       regressor using the original data is trained and has its  #
#       mean square error measured in the new data. L3 measures   #
#       how sensitive the regressor is to the new points. If the  #
#       original training points are distributed smoothly, their  #
#       interpolated variants will be close to the original data  #
#       items. Lower values indicate simpler problems.            #
#-----------------------------------------------------------------#
L3 = function(dataset, output = NULL){

    formatedDataset = FormatDataset(dataset, output)
    input = formatedDataset$input
    output = formatedDataset$output
    numberRows = formatedDataset$numberRows
    naRemove = !is.na(cor(output, input, method = "spearman"))
    numberColumn = sum(naRemove)

    order = order(output)
      
    output = output[order] 
    numberColumns = formatedDataset$numberColumn
    randomUniform = runif(numberRows - 1)
    
    if(numberColumns > 1){
      input = input[order,]  
      linearModel = lm(output~input[,naRemove])
      newInput = randomUniform*input[2:numberRows-1,naRemove] + 
        (1-randomUniform)*input[2:numberRows,naRemove]
      
    }
    else{
      input = input[order]
      linearModel = lm(output~input[])
      newInput = randomUniform*input[2:numberRows-1] + 
        (1-randomUniform)*input[2:numberRows]
      
    }
    newOutput = randomUniform*output[2:numberRows-1] + 
               (1-randomUniform)*output[2:numberRows]
    newPredict = predict.lm(linearModel,newdata = as.data.frame(input<-newInput))
    
    mean((newPredict - newOutput)^2)
}

#-----------------------------------------------------------------#
#   17 - S4 Non-linearity of Nearest Neighbor regressor           #
#       Employs the same procedure as before, but using a nearest #
#       neighbor regressor instead in the output predictions.     #
#-----------------------------------------------------------------#
S4 = function(dataset, output = NULL){

    formatedDataset = FormatDataset(dataset, output)
    input = formatedDataset$input
    output = formatedDataset$output
    numberRows = formatedDataset$numberRows
    naRemove = !is.na(cor(output, input, method = "spearman"))
    numberColumn = sum(naRemove)

    order = order(output)
    output = output[order,] 
    numberColumns = formatedDataset$numberColumn
    randomUniform = runif(numberRows - 1)
    
    if(numberColumns > 1){
      input = input[order,]  
      newInput = randomUniform*input[2:numberRows-1,] + 
        (1-randomUniform)*input[2:numberRows,]
      
          }
    else{
      input = input[order]
      newInput = randomUniform*input[2:numberRows-1] + 
        (1-randomUniform)*input[2:numberRows]
      
      }
    newOutput = randomUniform*output[2:numberRows-1] + 
               (1-randomUniform)*output[2:numberRows]
    newPredict = knn.reg(as.data.frame(input), as.data.frame(newInput), output, k = 1)$pred
    mean((newPredict - newOutput)^2)
}

#-----------------------------------------------------------------#
#   18 - T2 Average Number of Examples per dimension              #
#       T2 is the logarithm of the average number of examples per #
#       dimension. It gives an indicative on data sparsity.       #
#-----------------------------------------------------------------#
T2 = function(dataset, output = NULL){
    formatedDataset = FormatDataset(dataset, output)
    (formatedDataset$numberRows / formatedDataset$numberColumn)
}
