


## calculate effect decompositions
calc.effects <- function(theta, beta, alpha){

    ## Direct effect
    DE <- theta

    ## Mediated (indirect) effect
    IE <- beta * alpha
    
    ## total effect
    TE <- IE + DE

    df <- data.frame(Indirect=IE, Direct=DE, Total=TE)
    return(df)
}



##  calculates summary statistics of posterior distribution
calc.stats <- function(x){
    means <- colMeans(x)
    cis <- HPDinterval(mcmc(x))
    return(data.frame(means, cis, dose=factor(c(0, 1.6, 3.5, 7, 14))))
}
