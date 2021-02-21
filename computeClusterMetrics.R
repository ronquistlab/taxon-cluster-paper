# Set the p-distance threshold for PI clusters; 1.5% of 313 bp = 4.7 bp
p_threshold <- 0.015

# Set proportion of specimens threshold for major haplotype clusters
major_threshold <- 0.20

# Set distance threshold for checking (if this distance or larger from checked haplotypes, check)
checking_threshold <- 2

# Read in functions for deriving specimen checking values
source("specimenSelection.R")

# Read in the cluster table
E <- read.table("cluster_table.tsv", header=TRUE, sep="\t")

# Get the number of specimens per haplotype (1-based, in order)
# The table function simply counts the number of entries (specimens)
# for each unique haplotype label in the E$OC.0 column
# The output is a table with two columns, one with the label
# of the haplotype and one with the count, labeled $Freq.
D <- as.data.frame( table( E$OC.0 ) )
specimens <- D$Freq
num_haplotypes <- length(specimens)

# Get the morpho species for each haplotype (1-offset)
# Note that haplotype labels are not consecutive integers
# as one might have assumed. So we need to match on the
# unique haplotype labels rather than on a consecutive
# series of integers. To avoid excessive matching throughout
# the script, a costly operation, we sort the labels and
# assume that the info associated with haplotypes uses the
# same sort order always. We assume labeling of morphospecies
# in consecutive 0-offset indices for simplicity. We check
# that this assumption holds and throw an error if not.
unique_haplotypes <- sort(unique( E$OC.0 ))
x <- 1:length(unique(E$MORPHO))
x <- x-1
if ( sum( x %in% E$MORPHO ) != length(x) )
    stop( "Unexpected labeling of morphological species" )
morphospecies <- E$MORPHO[ match( unique_haplotypes, E$OC.0 ) ] + 1

# Get the column names for the clustering algorithms
# The first two columns are specimens and morphospecies; we skip these
colNames <- colnames(E)
colNames <- colNames[ 3:length(colNames) ]

# Now loop over each clustering column, accumulating the statistics we need
num_morpho <- length( unique( E$MORPHO ) )
specimensToCheck <- numeric()
speciesMissed <- numeric()
numClusters <- numeric()
numActualClusters <- numeric()
numSplitClusters <- numeric()
numMatchedClusters <- numeric()
numMergedClusters <- numeric()
matchScore <- numeric()
altMatchScore <- numeric()
for ( clusterMethod in colNames )
{
    cat( "Now processing clusterMethod ", clusterMethod, "\n" )
    cat( "==========================================\n" )

    # Get the cluster info
    clusters <- E[[clusterMethod]][ match( unique_haplotypes, E$OC.0 ) ]

    # Find out how many specimens we need to check of each species
    num_checks <- compute_checks( specimens, clusters, morphospecies, p_threshold, major_threshold, checking_threshold )

    # Append the data to the result vectors
    specimensToCheck[length(specimensToCheck)+1] <- sum(num_checks) 
    speciesMissed[length(speciesMissed)+1] <- sum( num_checks == 0 )
    numClusters[length(numClusters)+1] <- length( unique(clusters) )
    
    # Statistics from cluster table
    numActualClusters[length(numActualClusters)+1] <- length( unique(E[[clusterMethod]]) )
    X <- as.data.frame( table( E$MORPHO, E[[clusterMethod]] ) )
    X <- X[ X$Freq != 0, ]
    morpho_splits      <- unique( X$Var1[ duplicated(X$Var1) ] )
    morpho_singletons  <- unique(X$Var1)[ !(unique(X$Var1) %in% morpho_splits) ]
    cluster_splits     <- unique( X$Var2[ duplicated(X$Var2) ] )
    cluster_singletons <- unique(X$Var2)[ !(unique(X$Var2) %in% cluster_splits) ]

    num_split    <- length( morpho_splits)
    num_matched  <- sum( X$Var2[ match(morpho_singletons, X$Var1) ] %in% cluster_singletons )
    num_merged   <- num_morpho - num_matched - num_split
    alt_num_matched <- sum( X$Var2[ match( morpho_singletons, X$Var1 ) ] %in% cluster_singletons )

    numSplitClusters[length(numSplitClusters)+1] <- num_split
    numMatchedClusters[length(numMatchedClusters)+1] <- num_matched 
    numMergedClusters[length(numMergedClusters)+1] <- num_merged
    matchScore[length(matchScore)+1] <- (2*num_matched) / (num_morpho + length(unique(clusters)))
    altMatchScore[length(altMatchScore)+1] <- (2*alt_num_matched) / (length(unique(X$Var1)) + length(unique(X$Var2)))
}

clusterPerformance <- data.frame( colNames, specimensToCheck, speciesMissed, numClusters, numActualClusters, numSplitClusters, numMatchedClusters, numMergedClusters, matchScore, altMatchScore )
save(clusterPerformance, file="clusterMethodPerformance.RData")

