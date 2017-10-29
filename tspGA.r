set.seed(50)
n_cities <- 10
distance <- matrix((sample(n_cities^2)),n_cities,n_cities)
for (i in 1:n_cities){
        for (j in i:n_cities){
        distance[i,j]=distance[j,i]
        }
}
diag(distance) = 0

pathDistance <- function(path,distance) {
        d <- 0
        for (i in 1:(length(path)-1))
                d <- d + distance[ path[i] , path[i+1] ]
        return(d)
}


Ackeley <- function(x) {
        return(-20*exp(-0.2*sqrt(0.5*(x[1]^2+0^2))) - exp(0.5*(cos(2*pi*x[1])+cos(2*pi*0))) + exp(1) + 20)
}

binToDecMapping <- function(bitArray , range = c(-5,5)) {
	range[1] + (  ( abs(diff(range)) / ( 2^(length(bitArray)) - 1 ) ) ) * sum(2^((length(bitArray)-1):0) * bitArray )
}

randomTrials <- function(count = 100, n_cities = 10) {
	trialStacks <- NULL
	for (i in 1:count)
		trialStacks <- rbind( trialStacks ,  sample(n_cities))
	trialStacks
}

getDecStack <- function(trialStacks , range = c(-5,5) ) {
	decStack <- NULL
	for (i in 1:nrow(trialStacks))
		decStack <- c(decStack, binToDecMapping(trialStacks[i,], range = range) )
	decStack
}

getFitenssStack <- function(trialStacks, func = pathDistance , distanceMatrix=distance) {
	fitnessStack <- NULL
	for (i in 1:nrow(trialStacks))
		fitnessStack <- c(fitnessStack, func ( trialStacks[i,], distanceMatrix )  )
	fitnessStack
}

naturalSelection <- function( trialStacks , type = "tournament", func = pathDistance, distanceMatrix = distance ) {
	if ( type=="tournament" ) {
		fitnessStack <- getFitenssStack (trialStacks, func , distanceMatrix)
		len <- length(fitnessStack)
		winners <- NULL
		for (round in 1:len) {
			pair <- sample(len,2)
			if ( fitnessStack[ pair[1] ] <= fitnessStack[ pair[2] ] )
				winners <- rbind(winners,trialStacks[pair[1],])
			else
				winners <- rbind(winners,trialStacks[pair[2],])
		}
		return (winners)
	}
}


crossOver <- function(trialStacks, crossProb = 0.8 ,type = "one_point") {
	if (type=="one_point") {
		len <- nrow(trialStacks)
		digits <- ncol(trialStacks)
		for (round in 1:len) {
			pair <- sample(len,2)
			crossCut <- sample(digits, 1)
			if (runif(1) < crossProb) {	
				t1 <- trialStacks[ pair[1], ] [1:crossCut]
				indices <- NULL
				for (i in t1)
					indices <- c(indices, which(trialStacks[pair[2] , ] == i))
				t1 <- c(t1, trialStacks[ pair[2] , ][-indices])

				t2 <- trialStacks[ pair[2], ] [1:crossCut]
				indices <- NULL
				for (i in t2)
					indices <- c(indices, which(trialStacks[pair[1] , ] == i))
				t2 <- c(t2, trialStacks[ pair[1] , ][-indices])

				trialStacks[ pair[1], ] <- t1
				trialStacks[ pair[2], ] <- t2
			}
		}
		return(trialStacks)	
	}
}

mutation <- function(trialStacks, mutate_prob=(1/nrow(trialStacks)) , type = "bit_flip" ) {
	if (type=="bit_flip") {
		len <- nrow(trialStacks)
		digits <- ncol(trialStacks)
		for (i in 1:len) {
			for (j in 1:digits) {
				if (runif(1) < mutate_prob) {
					swapIndices <- sample(digits,1)
					temp <- trialStacks[i,j]
					trialStacks[i,j] <- trialStacks[i,swapIndices]
					trialStacks[i,swapIndices]  <- temp
				}
			}
		}
		return(trialStacks)
	}
}

GeneOpt <- function( func = pathDistance, distanceMatrix=distance, iterations=50, count=100, n_cities=10, crossProb=0.8, mutate_prob=0.001 , display=TRUE , cordList) {
	trialStacks <- randomTrials(count , n_cities )
	for ( i in 1:iterations) {
		trialStacks <-naturalSelection ( trialStacks , type = "tournament", func , distanceMatrix)
		trialStacks <- crossOver (trialStacks, crossProb = crossProb ,type = "one_point")
		trialStacks <- mutation (trialStacks, mutate_prob=(1/nrow(trialStacks)) , type = "bit_flip" )
		if (display==TRUE) {
			fitnessStack <- getFitenssStack (trialStacks, func , distanceMatrix)
			min_f = min(fitnessStack) 
			arg_min = trialStacks[which(fitnessStack ==  min_f)[1] ,  ]
			plot(c( cordList[arg_min[1],1], cordList[arg_min[2],1] ) , c( cordList[arg_min[1],2], cordList[arg_min[2],2] ), type='b', xlim = c(10,100) , ylim = c(90,100) )
			for (i in 2:(n_cities)) 
				lines(c ( cordList[arg_min[i] , 1], cordList[arg_min[i+1] , 1] ) , c( cordList[arg_min[i] , 2], cordList[arg_min[i+1] , 2] ) , type='b')
		}
	}
	fitnessStack <- getFitenssStack (trialStacks, func , distanceMatrix)
	min_f = min(fitnessStack) 
	arg_min_f = trialStacks[which(fitnessStack ==  min_f)[1] ,  ]
	list( min = min_f , arg_min = arg_min_f)
}

a <- read.csv('f14.csv')
a <- as.matrix(a)
b <- runif(14,25,50)
c <- runif(14,25,50)
a <- rbind(a,cbind(b,c))
b <- runif(14,60,90)
c <- runif(14,60,90)
a <- rbind(a,cbind(b,c))
distance <- as.matrix(dist(a))
n_cities <- nrow(distance)
Answer <- GeneOpt ( func = pathDistance, distanceMatrix=distance, iterations=200, count=500, n_cities, crossProb=0.8, mutate_prob=0.001, display=TRUE, cordList=a )
print(Answer)
