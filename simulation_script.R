## script to model the deterioaraty of a rna through time

## itinial size of strand
size0=10000
## degradation rate = average number of break per unit of time per bp
mu=0.01

## function that returns a simulated strand,
## input: size0: initial size of the strand, t time since it has been deteriorating, mu degradation rate
size <- function(size0,t,mu){
    n=rpois(1,t*mu*size0) # number of breaks = poisson distributed
    if(n>1)
    {
        if(n<size0)
        {
            k=c(0,sort(sample(1:size0,n)),size0) # loci of break = random}
            return(k[2:(n+2)]-k[1:(n+1)])
        }
        else
            return(0)
        
    }
    else
        return(size0)
}

## model the input of new rna through time
## if a = 0 constant flow
## otherwise linear model
concentration <- function(t,a,b)
{ ## population size model
    return(a*t+b)
}

##simulation
## buffering
tot=c() ## save all sizes produces
for(tt in 1:1000) ## for strand produced between now and 1000 time unites ago
{
    ##    print("--")
    ##    print(tt)
    ##   print(length(tot))
    for (j in 1:concentration(tt,0,1000)) ## for each strand introcued j time unit ago
    {# simulate its break and add it to the lot
        tot=c(tot,size(size0,tt,mu))
    }
    tot=tot[tot>30] ## eliminate very small bits to speed up simulations
}


plot(as.data.frame(table(tot)),xlab="RNA size")




mydist=as.data.frame(table(tot))

plot(mydist[,1],(mydist[,2]))


plot(mydist[,1],(mydist[,2])/sum(mydist[,2]),ylab="frequency", xlab="RNA size",col=rainbow(200)[mydist[,1]],type="l")


## what is the relationship between a strand 
u=1:500
plot(mydist[u,2],mydist[2*u,2], xlab="freq of RNA size x", ylab="freq of RNA size 2x",col=rainbow(500)[mydist[u,1]])



summary(tot)

plot(tot)

