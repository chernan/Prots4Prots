
test_gra_displayHeatmap <- function() {

    ## Setup data
    dataset <- replicate(8, rnorm(n=500, mean=11, sd=2.5))
    
    ## Run function... as test?
    print(displayHeatmap(dataset, "euclidean", title="Test Heatmap Euclidean"))
    print(displayHeatmap(dataset, "manhattan", title="Test Heatmap Manhattan"))
    print(displayHeatmap(dataset, "euclidean"))
    
}
