
test_gra_displayHeatmap <- function() {

    ## Setup data
    dataset <- replicate(8, rnorm(n=500, mean=11, sd=2.5))
    
    ## Run function... as test?
    print(displayHeatmap(dataset, "euclidean", title="Test Heatmap Euclidean"))
    print(displayHeatmap(dataset, "manhattan", title="Test Heatmap Manhattan"))
    print(displayHeatmap(dataset, "euclidean"))
    
}

test_gra_displayLinearReg <- function() {
    
    ## Setup data
    dataset <- data.frame(replicate(2, rnorm(n=500, mean=11, sd=2.5)))
    
    ## Run function... as test?
    print(displayLinearReg(dataset))
    print(displayLinearReg(dataset, 1, 2))
    print(displayLinearReg(dataset, 1, 2, 4, 16))
    print(displayLinearReg(dataset, 1, 2, 2, 3, TRUE))
    
}
