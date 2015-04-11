setwd("~/Desktop/data_incubator_project")

library(Matrix)
library(rworldmap)

## functions may be used for the project

markov.clustering = function ( trans.mat, inflate=2, multiple=2, n.iterations=100, tolerate=1e-15 ){
  # trans.mat = the transition matrix with row sums ones
  # inflate   = power of inflation per iteration
  # multiple  = power of multiplication per iteration
  # n.iterations = number of iterations
  # tolerate = set the values less than tolerate to be zeros
  # output: the list 'clusters', each element contains the indices of points in one cluster
  n = dim(trans.mat)[1]
  for (i in 1:n.iterations){
    print(i)
    for (k in 1:multiple){
      trans.mat = trans.mat %*% trans.mat  
    }
    trans.mat = trans.mat ^ inflate
    trans.mat = trans.mat * 1e9
    s = rowSums(trans.mat)
    s[s==0]=1
    sm = matrix(rep(s,n),n,n)
    trans.mat = trans.mat/sm
  }
  c.mat = as.matrix(unique(data.frame(t(as.matrix(trans.mat)))))
  ncl = dim(c.mat)[1]
  clusters=list()
  for (j in 1:ncl){
    clusters[[j]] = which(c.mat[j,]!=0)
  }
  return (clusters)
}

## read datasets
routes = read.table('./routes.dat',sep=',',header=FALSE,
                    col.names=c('airline','airline.id','depart.airport','depart.airport.id',
                                'dest.airport','dest.airport.id','code.share','stops','equipment'))
airports = read.table('./airports.dat',sep=',',header=FALSE,
                      col.names=c('id','name','city','country','code3','code4','latitude','longitude',
                                  'altitude','time.zone','dst','tz.name'))

## pre-process the datasets
# remove the redundant airports 
airports.recorded = as.character(airports$code3)
dept.airports = as.character(unique(routes$depart.airport))
dest.airports = as.character(unique(routes$dest.airport))
dept.airports.recorded = intersect(dept.airports,airports.recorded)
dest.airports.recorded = intersect(dest.airports,airports.recorded)
airports.reduced = union(dept.airports,dest.airports)

# update airports and routes to reduced airports
airports = subset(airports,airports$code3 %in% as.factor(airports.reduced))
n.airports = dim(airports)[1] # number of airports I care about
airports$id = 1:n.airports # reset the id's of the airports
routes = subset(routes,routes$depart.airport %in% dept.airports.recorded)
routes = subset(routes,routes$dest.airport %in% dest.airports.recorded)                
n.routes = dim(routes)[1] # number of routes

# add airport.id information into 'routes' according to 'airports'
routes$dept.id = rep(0,n.routes)
routes$dest.id = rep(0,n.routes)
for (i in 1:n.routes){
  routes$dept.id[i] = which(as.character(airports$code3)==as.character(routes$depart.airport[i]))
  routes$dest.id[i] = which(as.character(airports$code3)==as.character(routes$dest.airport[i]))
}

# create the pairwise matrix corresponding to the flight network
mat = Diagonal(n.airports,1)
for (i in 1:n.routes){
  v1 = routes$dept.id[i]
  v2 = routes$dest.id[i]
  tmp = mat[v1,v2]
  mat[v1,v2] = tmp + 1
}
r.sums = rowSums(mat)
div = matrix(rep(r.sums,n.airports),n.airports,n.airports)
mat = mat/div # make the row sums of mat ones by dividing the row sums


## conduct the markov chain 
clusters = markov.clustering(mat)
l = length(clusters)
cl= rep(0,n.airports)
for (i in 1:l){
  cl[ clusters[[i]] ] = i  
}
ncl = sapply(clusters,length)
ncl = ncl[ncl>=10]
ncl = sort(ncl,decreasing=TRUE)


## figures
png(filename='clustering_map.png',width=1000,height=600)
newmap <- getMap(resolution = "hight")
plot(newmap)
points(airports$longitude,airports$latitude,col=cl%%25,cex=0.5,pch=16)
dev.off()

png(filename='ncl.png',width=1000,height=600)
barplot(ncl,ylab='number of airports',main='Cluster Scales')
dev.off()

## save results
save.image("./results.RData")






