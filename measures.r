###################################################################
#   Complexity masures for regression problemas                   #
#   Proposed by Ana Carolina Lorena and Ivan G. Costa             #
#   Implemented by Aron Ifanger Maciel                            #
###################################################################

library(igraph)
library(FNN)

###################################################################
#   CALLS                                                         #
#       Normalize(dataset)                                        #
#       FormatDataset(dataset, output)                            #
#       MaxPosition(array)                                        #
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
#       H1(dataset)                                               #
###################################################################

###################################################################
#   First group of functions                                      #
#   These functions are used in measures evaluation               #
###################################################################

#-----------------------------------------------------------------#
#   01 - Normalize                                                #
#       Normalize the columns                                     #
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
#   04 - ExamplesRemovedNumber                                    #
#       Remove the examples from the dataset until a high correla #
#       tion value to the output is achieved                      #
#-----------------------------------------------------------------#

ExamplesRemovedNumber = function(x,y,minCorrelation)
{
    numberRows = length(x)   
    if(numberRows == length(y))
    {
        remainingRows = numberRows
        maxPosition = 0
        x = Normalize(x)
        y = Normalize(y)
        
        correlation = cor(x,y,method = "spearman")
        if(is.na(correlation))
            NA
        if(correlation < 0)
            y = 1-y
                    
        indexDistance = abs(x - y) 
        while(abs(correlation) < minCorrelation && !is.na(correlation))
        {
            maxPosition = MaxPosition(indexDistance)
            x = x[1:remainingRows != maxPosition]
            y = y[1:remainingRows != maxPosition]
            
            indexDistance = abs(x - y)
            
            remainingRows = length(x)
            correlation = cor(x,y,method = "spearman")
            
            if(is.na(correlation))
                correlation
        }
        
        (numberRows-remainingRows)/numberRows
    } else NA
}


###################################################################
#   Second group of functions                                     #
#   Measures evaluation                                           #
###################################################################

#-----------------------------------------------------------------#
#   05 - C1 - Maximum Feature Correlation to the Output           #
#       First one has to calculate the absolute value of the Spe  #
#       arman correlation between each feature and the outputs.   #
#       The absolute value is taken because both extremes of the  #
#       correlation measure, which range between [􀀀1; 1], repre   #
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
#   06 - C2 Average Feature Correlation to the Output             #
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
        correlations[column] = cor(output, input[,column], 
            method = "spearman")

    naRemove = !is.na(correlations)

    mean(correlations[naRemove])
}

#-----------------------------------------------------------------#
#   07 - C3 Individual Feature Efficiency                         #
#       Calculates the number of examples that must be removed    #
#       from the dataset until a high correlation value to the    #
#       output is achieved, divided by the total number of exam   #
#       ples. This computation is done for each individual featu  #
#       re and C3 returns the minimum of the values found, corres #
#       ponding to the feature more related to the output. Lower  #
#       values of C3 indicate simpler problems.                   #
#-----------------------------------------------------------------#

C3 = function(dataset, output = NULL, minCorrelation = 0.7)
{
    formatedDataset = FormatDataset(dataset, output)
    input = formatedDataset$input
    output = formatedDataset$output
    numberColumn = formatedDataset$numberColumn

    correlations = array(1,numberColumn)
    
    for (column in 1:numberColumn)
        correlations[column] = ExamplesRemovedNumber(output, 
                                input[,column], minCorrelation)        
    
    naRemove = !is.na(correlations)
    
    min(correlations[naRemove])
}

#-----------------------------------------------------------------#
#   08 - C4 Collective Feature Efficiency                         #
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

    remove = FALSE
    looked = array(FALSE, numberColumn)
    
    while(sum(!remove) > 0 & sum(looked) != numberColumn){
        
        indexMostCorrelated = MaxPosition(abs(cor(output,input,
                                                method="spearman")))
        
        if(!is.na(indexMostCorrelated)){
            
            looked[indexMostCorrelated] = TRUE
            linearModel = lm(output~input[,indexMostCorrelated])
            indexRemove = abs(linearModel$residuals) > minResidual
        
            input = input[indexRemove,]
            output = output[indexRemove]
        
            if(sum(!indexRemove) == length(indexRemove) | length(output)==1)
                remove = TRUE
            else
                remove = TRUE
        }
    }
    
    length(output)/numberRows
}


#-----------------------------------------------------------------#
#   09 - L1 Distance of the Data Itens to the Linear Function     #
#       The sum of the absolute values of the residues of a multi #
#       ple linear regressor. Lower values indicate simpler pro   #
#       blems, which can be fit by a linear function.             #
#-----------------------------------------------------------------#

L1 = function(dataset, output = NULL){
    
    formatedDataset = FormatDataset(dataset, output)
    input = formatedDataset$input
    output = formatedDataset$output

    naRemove = !is.na(cor(output, input, method = "spearman"))
    linearModel = lm(output~input[,naRemove])

    mean(abs(modelo$residuals))
}

#-----------------------------------------------------------------#
#   10 - L2 Average Error of Linear Regressor                     #
#       L2 sums the square of the residuals from a multiple       #
#       linear regression. Smaller values indicate simpler        #
#       (linear) problems.                                        #
#-----------------------------------------------------------------#

L2 = function(dataset, output = NULL){
    
    formatedDataset = FormatDataset(dataset, output)
    input = formatedDataset$input
    output = formatedDataset$output

    naRemove = !is.na(cor(output, input, method = "spearman"))
    linearModel = lm(output~input[,naRemove])

    mean(modelo$residuals^2)
}

#-----------------------------------------------------------------#
#   11 - S1 Output Distribution                                   #
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
#   12 - S2 Input Distribution                                    #
#       S2 first orders the data points according to their output #
#       values yi and then estimates the distance between pairs of# 
#       examples that are neighbors in the obtained ordering. S2  #
#       complements S1 by measuring how similar in the input space# 
#       are data items with close outputs. Lower values indicate  # 
#       simpler problems.                                         #
#       Ordena-se os dados pela sua saída y i e em seguida estima-#
#       se a distância entre pares de exemplos que sejam vizinhos.#
#       Com isso é possível medir o quão dados similares possuem  #
#       saídas diferentes. Valores menores indicam problemas mais #
#       simples.                                                  # 
#-----------------------------------------------------------------#

S2  = function(dataset, output = NULL){
    
    getDistaces = function(dataset, numberRows){

        distances = array(0,numberRows-1)
        
        for(line in 2:numberRows)
                distances[line] = dist(dataset[(line-1):line,])

        distances
    }

    formatedDataset = FormatDataset(dataset, output)
    input = formatedDataset$input
    output = formatedDataset$output
    numberRows = formatedDataset$numberRows
    
    naRemove = !is.na(input,output)
    input = input[,naRemove]
    order = order(output)
    
    distances = getDistaces(output[order,], numberRows)
    
    mean(distances)
}

#-----------------------------------------------------------------#
#   13 - S3 Error of a Nearest Neighbor regressor                 #
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

    for(line in 1:numberRows)
        predictions[line] = knn.reg(input[-line,], 
                                    input[line,], 
                                    output[-line], 
                                    k = 1)$pred

    mean(abs(predictions-output))
}

#-----------------------------------------------------------------#
#   14 - L3 Non-linearity of a Linear regressor                   #
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
    input = input[order,]    
    output = output[order,] 
    
    linearModel = lm(output~input[,naRemove])

    randomUniform = runif(numberRows - 1)
    
    newInput = randomUniform*input[2:numberRows-1,] + 
               (1-randomUniform)*input[2:numberRows,]
    newOutput = randomUniform*output[2:numberRows-1] + 
               (1-randomUniform)*output[2:numberRows]
    
    newPredict = predict.lm(linearModel,as.data.frame(newInput))

    mean(abs(newPredict - newOutput))
}

#-----------------------------------------------------------------#
#   15 - S4 Non-linearity of Nearest Neighbor regressor           #
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
    input = input[order,]    
    output = output[order,] 
    
    linearModel = lm(output~input[,naRemove])

    randomUniform = runif(numberRows - 1)
    
    newInput = randomUniform*input[2:numberRows-1,] + 
               (1-randomUniform)*input[2:numberRows,]
    newOutput = randomUniform*output[2:numberRows-1] + 
               (1-randomUniform)*output[2:numberRows]
    
    newPredict = knn.reg(input, newInput, output, k = 1)$pred

    mean(abs(newPredict - newOutput))
}

#-----------------------------------------------------------------#
#   16 - T2 Average Number of Examples per dimension              #
#       T2 is the average number of examples per dimension. It    #
#       gives an indicative on data sparsity                      #
#-----------------------------------------------------------------#

T2 = function(dataset, output = NULL){
    formatedDataset = FormatDataset(dataset, output)
    formatedDataset$numberRows / formatedDataset$numberColumn
}
