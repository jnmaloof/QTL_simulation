# QTL simulations

### Intro

The goal is to determine the best design for QTL mapping in an F2 population.  Given that a certain number of plants will be grown, is it better for all of those plants to be F2s or, instead, to use a smaller number of F2s and clonally replicate them?

### Load Libraries

```r
library(qtl)
library(ggplot2)
library(reshape2)
```

### Set up parameters.

```r
chromosomes <- 5
individuals <- 600 # number of individual plants that would be planted
repl <- 3 # number of prelications per F2
```

### Simulate a genetic map.
100 markers per chromosome

```r
map <- sim.map(len=rep(100,chromosomes),n.mar=100,include.x = FALSE)
```

### Simuate a QTL Cross
One QTL per chromosome, in the middle.  
The QTL effect size varies, increasing as we move from chromosome to chromosome.

```r
cross <- sim.cross(map=map,
                   model = matrix(c(1:chromosomes, #QTL location (chromosome)
                                    rep(50,chromosomes), #QTL position on chromosome
                                    1:chromosomes, # additive effect sizes, ranging from 1 to # of chromosomes
                                    rep(0,chromosomes)), # dominance deviation
                                  nrow = chromosomes, byrow = FALSE), #1 QTL per chromosome at 50cM, increasing effect size 1:12
                   n.ind=individuals,
                   type="f2")
cross <- sim.geno(cross)
```

### Add noise to phenotypes
Add some noise.  We create several different phenotypes, each with more noise added.

```r
phenotypes <- lapply(0:chromosomes, function(noise) {
  #create the plant replicates
  pheno.matrix <- t(apply(cross$pheno,1,rep,repl)) 
  
  #add the nosise
  pheno.matrix <- pheno.matrix + rnorm(n=prod(dim(pheno.matrix)), 
                                       mean=0,
                                       sd=noise)
}) #each element of the list is a phenotype matrix of 3 replicates, with increasing amounts of noise added.

#now create separate cross objects for each of the two designs
cross.no.rep <- cross

#for the no rep case take the first replicate
cross.no.rep$pheno <- as.data.frame(sapply(phenotypes, FUN = function(x) x[,1])) 

#as a reality check also make an object with F2s but with the 1/replication number of individiuals
cross.no.rep.small <- cross.no.rep[,1:round(individuals/3)]

cross.with.rep <- cross

#for the case with replicates take the average of the replicates
cross.with.rep$pheno <- as.data.frame(sapply(phenotypes, rowMeans))
cross.with.rep <- cross.with.rep[,1:round(individuals/3)] # so that the total number of individuals grown is the same 
```

### Check it
correlation should be better in the cross with replication

```r
cor(cross.no.rep$pheno)
```

```
##           V1        V2        V3        V4        V5        V6
## V1 1.0000000 0.9820973 0.9374766 0.8758875 0.7725978 0.7259740
## V2 0.9820973 1.0000000 0.9208799 0.8567686 0.7674631 0.7175219
## V3 0.9374766 0.9208799 1.0000000 0.8145670 0.7206982 0.6642565
## V4 0.8758875 0.8567686 0.8145670 1.0000000 0.6797688 0.6114985
## V5 0.7725978 0.7674631 0.7206982 0.6797688 1.0000000 0.5544198
## V6 0.7259740 0.7175219 0.6642565 0.6114985 0.5544198 1.0000000
```

```r
cor(cross.no.rep.small$pheno)
```

```
##           V1        V2        V3        V4        V5        V6
## V1 1.0000000 0.9836894 0.9388166 0.8841488 0.7872309 0.7514146
## V2 0.9836894 1.0000000 0.9222580 0.8652147 0.7750367 0.7225872
## V3 0.9388166 0.9222580 1.0000000 0.8149353 0.7210697 0.6718041
## V4 0.8841488 0.8652147 0.8149353 1.0000000 0.6897966 0.6633860
## V5 0.7872309 0.7750367 0.7210697 0.6897966 1.0000000 0.6038300
## V6 0.7514146 0.7225872 0.6718041 0.6633860 0.6038300 1.0000000
```

```r
cor(cross.with.rep$pheno)
```

```
##           V1        V2        V3        V4        V5        V6
## V1 1.0000000 0.9936731 0.9775825 0.9600793 0.9317744 0.9095759
## V2 0.9936731 1.0000000 0.9713712 0.9514925 0.9257797 0.8935442
## V3 0.9775825 0.9713712 1.0000000 0.9408865 0.9106251 0.8916907
## V4 0.9600793 0.9514925 0.9408865 1.0000000 0.8909694 0.8747427
## V5 0.9317744 0.9257797 0.9106251 0.8909694 1.0000000 0.8389174
## V6 0.9095759 0.8935442 0.8916907 0.8747427 0.8389174 1.0000000
```


### Do the QTL mapping
Use scanone, imputation method.

```r
no.rep.scanone <- scanone(cross.no.rep,pheno.col=1:(chromosomes+1),method = "imp")
no.rep.scanone.small <- scanone(cross.no.rep.small,pheno.col=1:(chromosomes+1),method = "imp")
with.rep.scanone <- scanone(cross.with.rep,pheno.col=1:(chromosomes+1),method = "imp")
```

### Plot the results

```r
no.rep.scanone.melt <- melt(no.rep.scanone,id.vars = c("chr","pos"))
no.rep.scanone.melt$type <- paste("no rep.",individuals,"F2s")

no.rep.scanone.small.melt <- melt(no.rep.scanone.small,id.vars = c("chr","pos"))
no.rep.scanone.small.melt$type <- paste("no rep.",round(individuals/3),"F2s")

with.rep.scanone.melt <- melt(with.rep.scanone,id.vars = c("chr","pos"))
with.rep.scanone.melt$type <- paste(repl,"replicates.",individuals, "total plants")
combined.scanone <- rbind(no.rep.scanone.melt,no.rep.scanone.small.melt,with.rep.scanone.melt)
pl <- ggplot(combined.scanone,aes(x=pos,y=value,color=variable))
pl <- pl + facet_grid(type ~ chr) 
pl <- pl + geom_line()
pl + ylab("LOD score")
```

![](simulations_files/figure-html/map2-1.png)

With no noise (phenotype "V1") replicaton (first row) and no replication--small (middle row) are the same.  Neither is as good as no replication--full (bottom row). With increasing noise we see the benefits of replication if the total number of F2s is the same (see V6). But if the total number of _plants_ is the same then the no replication design is much better.


