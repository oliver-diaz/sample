library(gmp)
library(Rmpfr)

# is Prime?

isPrime0 <- function(x){
  if(x==2){
    return(TRUE)
  }
  for(d in 2:(x-1)){
    if(x %% d ==0){ 
      return(FALSE) }
  } 
  return(TRUE)
}

isPrime0 <- Vectorize(isPrime0, vectorize.args = "x")

isPrime1 <- function(x){
  if(x==2){
    return(TRUE)
  }
  xb <- ceiling(sqrt(x))
  for(d in 2:xb){
    if(x %% d ==0){ 
      return(FALSE)
    }
    } 
  return(TRUE)
}
isPrime1 <- Vectorize(isPrime1, vectorize.args = "x")

# Prime factorization

primeFactors0 <- function(x){
  pfactors <- c()
  d <- 2
  while(x >1){
    if( isPrime1(d)){
      if(x %% d == 0){
        pfactors <- c(pfactors,d)
        x <- x/d
        } else{
          d <- d+1 
        }
      } else{
        d <- d+1
      }
  }
  return(pfactors)
}

primeFactors1 <- function(x){
  pfactors <- c()
  d <- 2
  while(x >1){
    if( d <= sqrt(x) +1){
      if(isPrime1(d)){
        if(x %% d == 0){
          pfactors <- c(pfactors,d)
          x <- x/d
          } else{
            d <- d+1
            }
        } else{
          d <- d+1
          }
    } else{
      return(c(pfactors,x))
    }
  }
  return(pfactors)
}


### timing process
rt <- Sys.time() 
isPrime0(4102317001) # 4102317001
Sys.time() - rt
 
rt <- Sys.time() 
isPrime1(4102317001)
Sys.time() - rt

rt <- Sys.time() 
(p <- primeFactors0(4102317001))
Sys.time() - rt

ptm <- proc.time() 
(p <- primeFactors1(4102317001))
proc.time() - ptm

