# Help functions to compute specimen selection statistics

load("distmat.RData")

# Read in the cluster table
E <- read.table("cluster_table.tsv", header=TRUE, sep="\t")

# Key from haplotype (1-offset) to 1-offset index into distance matrix
# Note that we need match haplotype ID to index in sorted array of unique
# haplotype IDs. The haplotype IDs are given in the column
# corresponding to OC at 0%. They are unique integers but not in a
# consecutive series.
unique_haplotypes <- sort( unique( E$OC.0 ) )
index2haplotype <- match( E$OC.0[ match( haploTypeSpecimens, E$Specimens ) ], unique_haplotypes )
haplotype2index <- match( 1:length(unique_haplotypes), index2haplotype )

# Convenient function to convert row and col (1-based) to
# element index into compact distance vector (1-based)
distvec_index <- function( i, j )
{
    row <- max(i,j)
    col <- min(i,j)
    
    ((row-1)*(row-2))/2+col
}


# Compute p-distance between haplotype (col index)
compute_pdistance <- function( i, j )
{
    i1 <- haplotype2index[i]
    j1 <- haplotype2index[j]
    pdistances[ distvec_index(i1,j1) ]
}


# Compute bp-distance between haplotypes (col index)
compute_bpdistance <- function( i, j )
{
    i1 <- haplotype2index[i]
    j1 <- haplotype2index[j]
    distances[ distvec_index(i1,j1) ]
}


# Compute maximum distance within a cluster
max_p_dist <- function( cluster_haplos )
{
    max_p <- 0
    for ( i in 1:(length(cluster_haplos)-1) )
    {
        i1 <- cluster_haplos[i]
        for ( j in (i+1):length(cluster_haplos) )
        {
            j1 <- cluster_haplos[j]
            d <- compute_pdistance(i1,j1)
            if ( d > max_p )
                max_p <- d
        }
    }
    return( max_p )
}


# Find and return the indices of the most distant haplotypes in a cluster
max_dist_haplotypes <- function( cluster_haplos )
{
    max_p <- 0
    max_i <- 0
    max_j <- 0

    for ( i in 1:(length(cluster_haplos)-1) )
    {
        i1 <- cluster_haplos[i]
        for ( j in (i+1):length(cluster_haplos) )
        {
            j1 <- cluster_haplos[j]
            if ( compute_pdistance(i1,j1) > max_p )
            {
                max_i <- i1
                max_j <- j1
            }
        }
    }

    c(max_i, max_j)
}


# Find and return all the major haplotypes in a cluster 
major_haplotypes <- function( cluster_indices, specimens, threshold_prop )
{
    num_specimens <- sum( specimens[cluster_indices] )
    cluster_indices[ specimens[ cluster_indices ] >= threshold_prop*num_specimens ]
}


#####################
# compute_checks: function for computing the specimens to be checked for morphology
#                 for each clustering
#
# parameters:
#   specimens           - a vector of length #haplotypes with number of specimens for each haplotype (1-offset)
#   clusters            - a vector of length #haplotypes with cluster index (potentially 0-offset) for each haplotype (1-offset)
#   morphospecies       - a vector of length #haplotypes with morphospecies index for each haplotype (1-offset)
#   p_threshold         - max-p-distance threshold for finding problematic clusters
#   major_threshold     - proportion threshold for recognizing major clusters
#   checking_threshold  - bp-distance threshold to nearest neighbor; if bp>= threshold, then checking is needed
#
# return value      - vector of number of specimens to check for each morphospecies; 0 if none (we missed this species)

compute_checks <- function( specimens, clusters, morphospecies, p_treshold, major_threshold, checking_threshold )
{
    # Initialize haplotype integer vector (1-offset)
    haplotypes <- 1:length( specimens )

    # Initialize return vector
    num_checks <- rep( 0, length( unique( morphospecies ) ) )

    for ( cluster_index in unique(clusters) )
    {
        # Get the indices of the haplotypes in the cluster
        cluster_haplos <- haplotypes[ clusters == cluster_index ]
        cat( "cluster_index=", cluster_index, "\n")
        cat( cluster_haplos, "\n" )

        # If this is a single-haplotype cluster, we are done; it validates
        # and we check a single specimen in the cluster to make sure it does not match some other cluster
        if ( length( cluster_haplos ) == 1 )
            check_haplos <- cluster_haplos
        else
        {
            # More than one haplotype 
            # Compute max distance within cluster
            max_p <- max_p_dist( cluster_haplos )
            cat("max_p=", max_p, "\n")

            # Find indices of the haplotype pair in the cluster with the maximum distance
            max_haplos <- max_dist_haplotypes( cluster_haplos )
            cat("max_haplos=", max_haplos, "\n")

            # Find indices of major haplotypes
            major_haplos <- major_haplotypes( cluster_haplos, specimens, major_threshold )
            cat("major_haplos=", major_haplos, "\n")

            # Assume that we check the most distant haplotypes
            check_haplos <- max_haplos
            
            # Now check if we should replace a distant haplotype with a major haplotype
            # We first remove any potential duplicates in major_haplos that are already in check_haplos
            # If we replace one, then we remove that from remaining major_haplos because we do not
            # want to substitute in the same major haplotype twice
            # Question for Emily: Should it be < checking_threshold, or <= checking_threshold?
            replacement_haplos <- major_haplos[ !(major_haplos %in% check_haplos) ]
            if ( length(replacement_haplos) > 0 )
            {
                for ( j in 1:length(check_haplos) )
                {
                    for ( k in 1:length(replacement_haplos) )
                    {
                        if ( compute_bpdistance(check_haplos[j],replacement_haplos[k]) <= checking_threshold )
                        {
                            check_haplos[j] <- replacement_haplos[k]
                            replacement_haplos <- replacement_haplos[ -k ]  # Remove this major haplotype; we do not want it twice
                            break;
                        }
                    }
                    if ( length(replacement_haplos) == 0 )
                        break;
                }
            }

            # Check whether it is a PI cluster; if so, add any remaining major haplotypes
            # Here we only use max_p as this can easily be applied across clustering algorithms
            # To add stability, we need a definition that applies across clustering algorithms
            if ( max_p >= p_threshold )
                check_haplos <- unique( c( check_haplos, major_haplos ) )

            if ( length(check_haplos) < 2 )
                cat( "ERROR: too few specimens to check...\n")

            # Check if we have more than one species so we need to do more thorough checking
            if ( length( unique( morphospecies[check_haplos] ) ) > 1 || max_p >= p_threshold )
            {
                # We need to do thorough checking
                cat( "Thorough checking needed\n" )
                for ( n in 1:length(check_haplos) )
                    cat( n, " -- haplotype=", check_haplos[n], " -- species=", morphospecies[check_haplos[n]], "\n" )

                # Add major clusters if we thought first that this was simple (max_p < p_threshold )
                if ( max_p < p_threshold && length(major_haplos) > 0 )
                    check_haplos <- unique( c( check_haplos, major_haplos ) )

                # Cycle through all haplotypes in the cluster that are not checked
                # If it differs from a checked haplotype by less than the checking_threshold
                # we do not have to check it; otherwise it needs to be checked
                for ( j in cluster_haplos )
                {
                    if ( is.na( match( j, check_haplos ) ) )
                    {
                        needs_checking <- TRUE
                        for ( k in check_haplos )
                        {
                            if ( compute_bpdistance(j,k) < checking_threshold )
                            {
                                needs_checking <- FALSE
                                break
                            }
                        }
                        if ( needs_checking )
                            check_haplos <- c( check_haplos, j )
                    }
                }

                # If species boundaries are unclear, we check even more specimens
                # We know that if specimens have not been checked, they are at most 1 bp away from the closest species
                # If they are within 2 bp of more than one species, then we add them in to clarify boundaries
                for ( j in cluster_haplos )
                {
                    if ( is.na( match( j, check_haplos ) ) )
                    {
                        neighbours <- numeric();
                        for ( k in check_haplos )
                        {
                            if ( compute_bpdistance(j,k) <= checking_threshold )
                                neighbours <- c( neighbours, morphospecies[k] )
                        }
                        if ( length( unique( neighbours ) ) > 2 )    
                            check_haplos <- c( check_haplos, j )
                    }
                }
            }
        }

        # Now we can summarize how many checks and for which species
        numPrevious = sum(num_checks)
        for ( i in check_haplos )
            num_checks[ morphospecies[i] ] <- num_checks[ morphospecies[i] ] + 1

        cat( "haplotypes=", length(cluster_haplos), "-- haplotypes to check=", length(check_haplos), " -- sum_delta=", sum(num_checks) - numPrevious, "\n" ) 
        if ( sum(num_checks) - numPrevious != length(check_haplos) )
        {
            cat("ERROR\n")
            for ( i in check_haplos ) cat( i, " -- ", morphospecies[i], "\n" )
            break
        }
    }

    return( num_checks )
}


