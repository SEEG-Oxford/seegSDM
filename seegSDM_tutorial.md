Example workflow / tutorial for the seegSDM package
========================================================

This document runs through a typical workflow for distribution modelling using the ```seegSDM``` package.

It will include importing and checking data, running BRT ensembles in parallel and examining the fitted models.
This is a work in progress, so sorry if it stops halfway through or doesn't make any sense yet! Please report any issues via the [issues tracking system](https://github.com/SEEG-Oxford/seegSDM/issues).

The structure is:

##### [Installing the package](#install)
##### [Loading data](#load)
##### [Quality control](#quality)
##### [Generating pseudo-absences](#pseudo)
##### [Extracting covariate data](#extract)
##### [Running a single BRT model](#BRT)
##### [Running a BRT ensemble in parallel](#ensemble)
##### [Visualising the BRT ensemble](#vis)
##### [Outputting the results](#output)



### <a id="install"></a>Installing the package

To install seegSDM straight from github we use the ```install_github``` function in the ```devtools``` package.

```r
# if it isn't already installed, install devtools from CRAN
# install.packages('devtools')

# and load it
library(devtools)

# use install_github to install seegSDM, giving the name of repo & owner
# and installing all the packages it depends on

# install_github('seegSDM', 'SEEG-Oxford', dependencies = 'Depends')

# seegSDM should now be installed, so we just need to load it
library(seegSDM)
```



### <a id="load"></a>Loading data

Next we load in some occurrence data. Here we'll use fake occurrence data provided with the package, though you can import your own using e.g. ```read.csv```. The occurrence object has a number of different columns, giving coodinates of observations of a fake disease (both polygon and point data) as well as information needed for model fitting.


```r
# load the data
data(occurrence)

# look at the first 6 lines
head(occurrence)
```

```
##   UniqueID Admin Year     x     y Area
## 1        1 -9999 1990 -6.15 -4.75   NA
## 2        2 -9999 2013 -4.65  5.05   NA
## 3        3 -9999 2009  2.65 -2.95   NA
## 4        4 -9999 1992 -5.35  5.05   NA
## 5        5 -9999 2008 -6.95 13.35   NA
## 6        6 -9999 1999 -9.75  3.75   NA
```


Most of the ```seegSDM``` functions use ```SpatialPoints*``` objects (from the ```sp``` package) so we convert ```occurrence``` into one of these. We can do this using the function ```occurrence2SPDF``` which makes some assumptions about whats in ```occurrence``` (see the helpfile for details).


```r
# convert to a SpatialPoints object
occ <- occurrence2SPDF(occurrence)
```


Next we load a bunch of raster files containing covariates for the model. Again, we use some fake data rasters which are provided with the package. You can import your own using ```raster``` from the ```raster``` package or maybe using the ````seegSDM``` function ```importRasters``` to make things a little easier.


```r
# load the covariate rasters
data(covariates)

# see a summary
covariates
```

```
## class       : RasterBrick 
## dimensions  : 200, 200, 40000, 3  (nrow, ncol, ncell, nlayers)
## resolution  : 0.1, 0.1  (x, y)
## extent      : -10, 10, -5, 15  (xmin, xmax, ymin, ymax)
## coord. ref. : +init=epsg:3395 
## data source : in memory
## names       :   cov_a,   cov_b,   cov_c 
## min values  : -1.2989, -0.6801,  1.0000 
## max values  : -0.2735,  0.6020,  6.0000
```

```r

# and plot them
plot(covariates)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 



### <a id="quality"></a>Quality control

There are currently two functions to check the quality of incoming data: ```checkRasters``` which checks that rasters match up with an expected template raster and ```checkOccurrence``` which runs a number of checks on the incoming occurrence data to make sure there aren't any serious errors, and to tidy up some more minor ones.

First we'll run ```checkRasters``` to make sure the ```covariates``` lines up with a template raster that we know is correct.


```r
# load a template raster to check covariates against
data(template)

# run checkRasters
checkRasters(covariates, template)
```

```
## class       : RasterBrick 
## dimensions  : 200, 200, 40000, 3  (nrow, ncol, ncell, nlayers)
## resolution  : 0.1, 0.1  (x, y)
## extent      : -10, 10, -5, 15  (xmin, xmax, ymin, ymax)
## coord. ref. : +init=epsg:3395 
## data source : in memory
## names       :   cov_a,   cov_b,   cov_c 
## min values  : -1.2989, -0.6801,  1.0000 
## max values  : -0.2735,  0.6020,  6.0000
```


If everything is fine the original object is returned (so here R prints a summary), otherwise an error is thrown. See ```?checkRasters``` for more details of the checks that are done.

Next we use ```checkOccurrence``` run a whole bunch of checks (```?checkOccurrence``` for details) on the incoming occurrence data. We need a few other raster layers to do this, including an evidence consensus layer and a raster brick of admin unit maps (again, these are fake data).

We load and plot the evidence consensus layer:

```r
data(consensus)
plot(consensus)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 


and a raster brick with (fake) GAUL codes at 4 different levels.

```r
# load and plot the four layers
data(admin)
plot(admin)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 


Now we can run ```checkOccurrence```

```r
# overwriting the occ object we defined before
occ <- checkOccurrence(occurrence, consensus, admin)
```

```
## 25 polygons had areas greater than the area threshold of 1 and will be removed.
```


some of the occurrence points were in polygons with areas geater than the default allowed maximum. The maximum allowed area can be altered. Other than this all of the checks were passed. If any major problems had been detected, ```checkOccurrence``` would have thrown an error message.

### <a id="pseudo"></a>Generating pseudo-absences & extracting the data

There are various different schools of thought on how to select pseudo-absences for presence-only species distribution modelling.

The recent Dengue distribution paper by [Bhatt et al.](http://dx.doi.org/10.1038/nature12060) used a distance threshold and an evidence consensus layer to select both pseudo-absence and pseudo-presence points. The function ```extractBhatt``` takes a vector of three required parameters and applies this procedure. It also extracts the values of the covariates for all of these points. For polygon occurrence records the multiple values are efficiently extracted and summarised using the function ```extractAdmin```. For point data the ```raster``` package function ```extract``` is used.

Other pseduo-data generation methods can be carried out using the more general functions ```bgSample``` and ```bgDistance```. Here we apply ```extractBhatt``` to our data using some arbitrary parameter settings.


```r
# to make sure this tutorial is reproducible, we set the seed for the
# random number generator
set.seed(1)

# run extractBhatt, defining the last covariate as a factor
lis <- extractBhatt(c(2, 1, 5), occ, covariates, consensus, admin, factor = c(FALSE, 
    FALSE, TRUE), return_points = TRUE)
```


With ```return_points = TRUE```, ```extractBhatt``` produces a list with three elements: a dataframe for modelling and the locations of the pseudo-presence and pseudo-absence records. We have a look at the object and plot the points.


```r
# look what's in the list returned
names(lis)
```

```
## [1] "data"            "pseudo_absence"  "pseudo_presence"
```

```r

# summarise the dataframe which we'll use for modelling
summary(lis$data)
```

```
##        PA          cov_a            cov_b        cov_c  
##  Min.   :0.0   Min.   :-1.265   Min.   :-0.583   1:  3  
##  1st Qu.:0.0   1st Qu.:-0.859   1st Qu.:-0.305   2: 70  
##  Median :0.5   Median :-0.708   Median :-0.149   3:102  
##  Mean   :0.5   Mean   :-0.739   Mean   :-0.126   4: 93  
##  3rd Qu.:1.0   3rd Qu.:-0.592   3rd Qu.: 0.026   5: 30  
##  Max.   :1.0   Max.   :-0.361   Max.   : 0.457   6:  2
```

```r

# evidence consensus layer as a background
plot(consensus)

# add the pseudo-absences in light blue (note they aren't in high scoring
# consensus regions)
points(lis$pseudo_absence, pch = 16, col = "light blue")

# and pseudo-presences in purple (not in the very low scoring regions)
points(lis$pseudo_presence, pch = 16, col = "purple")

# and add the occurrence points
points(occ, pch = 16)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10.png) 



### <a id="BRT"></a>Running a single BRT model 

We're now ready to run a BRT model. The ```gbm.step``` function in the ```dismo``` package (which ```seegSDM``` loads) runs a cross-validation procedure to pick the best number of trees (an important parameter in BRT) and runs the final model. ```seegSDM``` provides a wrapper function ```runBRT``` for ```gbm.step``` with a set of SEEG-preferred default settings. Only four arguments need to be provided: the dataframe, the indicies for the presence/pseudo-absence and covariate columns and a ```RasterBrick``` object to predict to.


```r
brt <- runBRT(lis$data, 2:4, 1, covariates)
```

```
## 
##  
##  GBM STEP - version 2.9 
##  
## Performing cross-validation optimisation of a boosted regression tree model 
## for PA with dataframe data and using a family of bernoulli 
## Using 300 observations and 3 predictors 
## creating 10 initial models of 10 trees 
## 
##  folds are stratified by prevalence 
## total mean deviance =  1.386 
## tolerance is fixed at  0.0014 
## now adding trees... 
## fitting final gbm model with a fixed number of  430  trees for  PA
```


```runBRT``` returns a list giving the model, a raster of the predicted probability of presence and data to plot the relative influence and covariate effects. The last two are used in the BRT ensemble modelling, but we can visualise the single model using functions from the ```gbm``` package.

We can plot the individual marginal effect curves for each covariate...


```r
par(mfrow = c(1, nlayers(covariates)))
for (i in 1:nlayers(covariates)) plot(brt$model, i)
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12.png) 


...the 2-dimensional interaction between the first two covariates...

```r
plot(brt$model, 1:2)
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13.png) 


...the relative influence of each covariate...

```r
summary(brt$model)
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14.png) 

```
##         var rel.inf
## cov_a cov_a  87.688
## cov_b cov_b  10.255
## cov_c cov_c   2.057
```


...and the map of predicted habitat suitability produced by ```runBRT```.

```r
plot(brt$pred, zlim = c(0, 1))
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15.png) 



### <a id="ensemble"></a>Running a BRT ensemble in parallel

The plots above show one of the drawbacks of single BRT models: they fit very jerky effects of environmental covariates. We can get smoother and more realistic effect curves through fitting and averaging an ensemble of multiple BRT models. Ensembling also allows us to reduce reliance of the model on the arbitrary selection of parameters for pseudo-data generation and enables us to produce estimates of uncertainty in both the model and its predictions.

Next we'll set up an ensemble of models fitted using pseudo-data generated in different ways by altering the parameters passed to ```extractBhatt```. We first define the ranges of parameters to use then get all the different combinations using R's ```expand.grid``` function.


```r
# pseudo-absences per occurrence:
na <- c(1, 5, 10)
# pseudo-presences per occurrence:
np <- c(0.1, 0.05, 0.01)
# distance (in decimal degrees) within which to sample these
mu <- c(1, 3, 5)

pars <- expand.grid(na = na, np = np, mu = mu)

# each row contains a different configuration od parameters
head(pars)
```

```
##   na   np mu
## 1  1 0.10  1
## 2  5 0.10  1
## 3 10 0.10  1
## 4  1 0.05  1
## 5  5 0.05  1
## 6 10 0.05  1
```

```r

# now there are 3 * 3 * 3 = 27 different combinations and models to run
nrow(pars)
```

```
## [1] 27
```


To use all of the ```seegSDM``` functions for running ensembles we need to change these into a list, which we can do using ```lapply``` and a small function to subset ```pars```.


```r
sub <- function(i, pars) pars[i, ]
par_list <- lapply(1:nrow(pars), sub, pars)
```


Unfortunately running BRT ensembles can be time consuming and it's preferable to run them in parallel across multiple processors. There are a number of different R packages to help with this, but here we use the ```snowfall``` package.

Loading seegSDM has already loaded snowfall, so the first thing we need to do is set up a parallel cluster using the function ```sfInit```, which is as easy as this:


```r
# set up a cluster of two cpus in with parallel execution.  you may want
# to run a different number depending on your computer!
sfInit(cpus = 4, parallel = TRUE)
```

```
## R Version:  R version 3.0.1 (2013-05-16)
```

```
## snowfall 1.84-4 initialized (using snow 0.3-12): parallel execution on 4
## CPUs.
```


We need to export any functions we want to run in parallel out to each processor. Everything we'll use is in ```seegSDM```, so we export the whole package using ```sfLibrary```


```r
sfLibrary(seegSDM)
```

```
## Library seegSDM loaded.
```

```
## Library seegSDM loaded in cluster.
```


Now we're ready to run some code in parallel. We use the function ```sfLapply``` which acts like ```lapply```, except that each element of the list is processed in parallel. We run ```extractBhatt``` over these different parameter settings.




### <a id="vis"></a>Visualising the BRT ensemble


### <a id="output"></a>Outputting the results



