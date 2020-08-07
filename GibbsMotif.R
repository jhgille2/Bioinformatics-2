# Gibbs sampler for motif detection - 
# Unsupervised detection of a sequence motif of length 14 from a set of 10 
# 300 bp DNA sequences
#
# Author: Jay Gillenwater
# Last Edited: 2020-04-30

library(doParallel)
library(foreach)
library(parallel)
library(snow)
library(doSNOW)
library(pbapply)


# The 10 sequences 
Seq0 <- c("CGGACTGTAAGATTAAGTCCAGCCGCGAAATACCTCAGTTAGGGTGGCAACATTCTTTTAGCGATTCACGCACTCCGCGTAGGGCGCATGGCGAATATATCAAATGACGAAGTCCACGCTGACCGGTCACGCTGTTCTGGTAGGCTATGGGAAAATTGAGTACAGATTACAGAGATTGCTGGTCAGGGTAAGAATACTCTCTGATCGCCGAAGGCAGGCACCTCCCGGTTCTGGCACACGGAAGCTGTAAAAAGTCCTGCCCGTCCTTACTCTGGCCACGCTCATCGCTCTACCACAGAA")
Seq1 <- c("TGAGGGCCACCTGGCTGGAGACATGTAGAGCTACGATCGGCGCCCCCCACCGTAGTGCCGTCATCTGGATTTCACACTGCTGCATGAACCTTCAAAAAATGATGCAGTATATACGATTACAGATTACAGCCTGATTGGAGAGAAAGATTCTTAAAGGGACTTCAAACGAGCCAAGCCACAGCATAGTGAGTACCGCGTTATGCGTGTTGCTTGCGCCTCGTATATCAGGTTAAGTCGGGACCGACCTGGAACTTGGGTTGTGGACCCACTGATGAAGTCGCCTAGGTGGTAGCATAGAGA")
Seq2 <- c("GAATCCGGGAGAAACGCAACAGGCCGCGTTAGTGCTACACGGAACGGCCCCTTGGTTTGCAGATAGGTAGTCTGTATGATTTAATCAGTCTACTGTGTCGCCAGAAACTGAACTATGTACGTCGGTACGGCTATTCATGCGTGGCTGCAAATATTAGTCTATAATTTCTCCACTTGCTCCGAGCGGATCAATGCCTATGCTGACTGACATTGGTGTACAGATTACAACGTTCAGTATGAGTTTCCCTAATCCTTCTTGTCTCTACACGCCCGCTCGTCCTGGTGATACCGGACATTCCTT")
Seq3 <- c("CGACGGGGATTCCGCGCAGAGTATCTAGTTGACGCTCGCAATGCTGACTAAGGCTCAAAAACGTGTCGAATGATTTGTCAGGCACGCCTTCATTGCAAACGTCTAATAATCTTGCCTGACCCGAACTTAGATGGGCTGATATGCAAGTACGGCAAACATGCAGTTCTTGAGGCAGCTCCGGATCCGGTGAGTCTGATTAACCCACTTGCAGTCAATCCCGAGTTATTTGTTAACCAGGTTCCATAGTGCACCAGATGCGGACTCCATGTATAGCGCCTATGACTACAGATTACAAGACAT")
Seq4 <- c("GACCTCTGGGAATTTCGCCGGATCTTTGAGACCGAAAGATTACAGATAACAACACAGGTCCCGTAGCCTTTACCCCAGGTACAGTTGTCCCTAACAGGCGCCTCGTTTACTTGGACATGAGGGCCATCGAGTATCAGACCCTGGGCTATTAGATGTGAAGATTCGTGACTCCCTCCAGGTCGATTTAGTCTGTTAATACTAAGCTTAAGGTCTTTAACCGCGGGGTTGTTTGCTCCTATATGTGATAACATGGTCTTGCGTTGATCGCTATACTCTACTATCACGAGGATAAAACCCAAT")
Seq5 <- c("CCCCCTCGACTTCGTGAAATGTGTGTTGTAAGTTGTTGTAGAAGCGATGATTTAGATTGGTCAATGTGCTAGGGTCGGTGAAAATGCTTCAGAGGGCTGGAGCGCTACCCTTTAATGGTTCCGAACATACCCCTATATAAGGGGCTTTTACGACTTTAGTGAGCCCCCCCCTCTTCGAAGCTGAATACTTTGAACTGGCCATGGACTTCCTGCGCGTGGGGCTTGGCATATATGCATATAAAAAGTCTCTGGGTACCGCTCTGAGTGATCAAAGATGACAGCGTCGTCTTTCCAACCCAG")
Seq6 <- c("TTACATACTCTTCGCGTAGAGCATACAGATTATATCATTTCGCTACCACCCTATAGTACTTCCAAATCGCGCAGTACCGTCGAGGATAACGGTCCGGCAAGCACGTAGACCTACTTGCCTGTTAGCTAACCGTCCATCCGGAGCACTTGCGACTACCCGAAACGTGCCCTCCCACCTTCCCTAGCCTTTGTCCGATGGAGGTTGCAACCTCTTGGAGTCGCAGGGATCTTTAACTCACGTCTAGACCAAGTACCCAGGAGCGACGGGATCGGGAGACTTTATTTGGCGTTCAGGGCGAGT")
Seq7 <- c("TAGAGTAGACGAGAGCAAAGACACAGTGGGGTGGCTACCACGCTGTTCGTATATGCACTTGGCGCCGCGTTGATCCTAGTGGGTTCAGTTTACAGAGTTAGATCACAGATAACATCGATTTCGGTCTGGTGATACGGAAAAATGAATTGAAGATGATCCTAGCTATTCCGGGTACAAATTTAATATACATAGGTGCGCCGGCTCATGACATCCATACTTGCCCAGTACAATAGTTCTATTTTTAGGAGTTCAGTATAGTGTCTCTTCATTAGCGTACTCCTCGGCCCTTGGTCCAGAACT")
Seq8 <- c("CCGATGACCTTGCGCCTCAACAACCGAGGCCACCTACCATGTCAGGTGTGACGTACCTCGCCATTACGATCACAGATTACAAATGCAGGCGCCTAATATTTTCACGCTAAATAGAATAGTTTTTTCAGATGAGAAAGACGTATGACAACCAAGGTTAATTCGTAATAGACAGGCGCTGGGGCTAGGGATGCCTTTCAGTACTTTACAACACACCTTCGGCCTGCGTGGGCGGAGCGGAGTAATAGCTGTTGTCCCATTGGCCATAGAGGAGTTTGGTACCCTGAATCGCGGTCCTGGGCT")
Seq9 <- c("CCGATACCATGAAACTTTCGCCATCCAACGTGATTCAGCGTTTATTAGCAGTGCTGGAGACCAAACTCAGATCGCATATTACAACAGGTCGCTCAGTGGAGGCCCATTAGAACACCGAGCACTATTTCATCGTATTGCGGTCTAGTTGCAATATGAACGATCCGCATCGTTCGGGTGGCACTGTGAAGATACGCATCGTGAGGAAACCCTGATGAATAGTCCTGTTCATAATGAACAGAGTCTGCTGTGTGCTGGATTCGCGGGTGATATACGAAGTGACACCGAAATGCATTGGCTCGC")


# Combine these sequences into a list and then seperate each string into a matrix of characters
AllSeqs    <- list(Seq0, Seq1, Seq2, Seq3, Seq4, Seq5, Seq6, Seq7, Seq8, Seq9)
AllSeqs    <- lapply(AllSeqs, function(x) t(as.matrix(unlist(strsplit(x, ""))))) # Split strings into characters (nucleotides)

# The length of the motif to search for
MotifLen <- 14

# A function to implement a sum of pairs scoring scheme, Match is the score for matches, MisMatch is the penalty for mismatches
SumOfPairs <- function(MotifMatrix, Match = 5, MisMatch = -2){
  
  SeqLen <- ncol(MotifMatrix)
  
  # Dataframe to reference every row combination
  SeqPairs <- expand.grid(c(1:nrow(MotifMatrix)), c(1:nrow(MotifMatrix)))
  SeqPairs <- SeqPairs[!c(SeqPairs[, 1] == SeqPairs[, 2]), ]
  
  SeqPairs$PairScore <- NA
  
  # A function to score two sequences
  ScorePair <- function(Seq1, Seq2){
    
    MatchNum      <- sum(Seq1 == Seq2) # Number of matches
    MatchScore    <- Match*MatchNum    # Score of matches
    MisMatchScore <- MisMatch*(SeqLen - MatchNum) # Score of mismatches
    
    MatchScore + MisMatchScore # Final score
  }
  
  SeqPairs$PairScore <- apply(SeqPairs, 1, function(x) ScorePair(MotifMatrix[x[1], ], MotifMatrix[x[2], ]))
  
  return(sum(SeqPairs$PairScore))
}


# The Gibbs sampler function
# Seqs: A list of DNA sequences
# MotifLen: Legnth of the motif to detect
# TrackPerformance: Record the score for the motif for each iteration? (for plotting)
GibbsMotif2 <- function(Seqs = AllSeqs, MotifLen = 14, TrackPerformance = FALSE){

  # A function that fnds the frequency of each nucleotide base in a vector
  CountBase <- function(x, offset = 1){
    nA <- sum(x == "A")
    nT <- sum(x == "T")
    nG <- sum(x == "G")
    nC <- sum(x == "C")
    
    # An "offset" to prevent 0 probabilities in calculations
    (c("A" = nA, "T" = nT, "G" = nG, "C" = nC) + (offset/4))/(length(x) + offset)
    
  }
  
  # The background probability of observing each nucleotide in the full data
  BackgroundProbs <- CountBase(unlist(Seqs))
  
  # Sample random starting places for each motif
  StartingPos <- sample(1:(300-MotifLen), length(Seqs), replace = TRUE)
  
  # The ending position of each substring
  EndingPos <- StartingPos + MotifLen - 1
  
  # A vector to compare starting positions of k-mers
  CheckStarts <- rep(FALSE, length(Seqs))
 
  # Probability associated with each kmer in the current set of motifs
  SeqProbs <- rep(0, length(Seqs))
  
  # Initialize a matrix to hold the proposed motifs
  ProposalMat <- matrix(nrow = length(AllSeqs), ncol = MotifLen)
  
  # Fill this matrix with the initial substrings
  for(i in 1:length(Seqs)){
    ProposalMat[i,] <- Seqs[[i]][StartingPos[[i]]:EndingPos[[i]]]
  }
  
  # Vector to hold a score for each iteration
  IterationScore <- c()
  
  while(!all(CheckStarts)){

    # Sample a motif to exclude from the current set
    i <- sample(c(1:nrow(ProposalMat)), 1)
    
    # Exclude the current kmer from the proposal matrix
    CurrentProposalMat <- ProposalMat[-i, ] 
    
    # A matrix which holds the probability of observing each base at each position in the proposal matrix
    ProposalProbs <- apply(CurrentProposalMat, 2, CountBase)/BackgroundProbs
    ProposalProbs <- ProposalProbs/colSums(ProposalProbs)
    
    # A list o hold the probability of each tested k-mer in a moving window
    AllKmerprobs <- vector("list", length = 300-14)
    
    # Calculate weights for each 14-mer of the sampled sequence
    for(j in 1:(300 - MotifLen)){
      
      # Pull a 14-mer from the current sequence
      CurrentKMer <- Seqs[[i]][j:(j + (MotifLen - 1))] 
      
      # I'll find the probability of observing the sequence in the background by summing the background probabilities
      # of each base in the sequence
      BackgroundProb_seq <- prod(BackgroundProbs[CurrentKMer])
      
      # For the sequence probability in the current subset, I need to lookup the probability for each base according to the 
      # observed probability at it's position
      CurrentProb_seq <- prod(diag(ProposalProbs[CurrentKMer,]))
      
      # Update the list holding the probabilities for each k-mer
      AllKmerprobs[[j]] <- CurrentProb_seq/BackgroundProb_seq
    }
    
    # Use probabilities from the list to sample a new 14-mer
    StartWeights <- unlist(AllKmerprobs)/sum(unlist(AllKmerprobs))
     
    # Cumulative sum of probabities so that a starting site can be sampled by 
    # drawing a sample from the uiform distribution
    StartWeights_prob <- cumsum(StartWeights)
    
    # Sample a new starting place for this kmer weighted by the probabilities
    NewStart <- min(which(StartWeights_prob >= runif(1)))
    
    # The probability of the current sequence
    SeqProbs[[i]] <- StartWeights[NewStart]
    
    
    # Update Starting position check
    CheckStarts[[i]] <- ifelse(NewStart == StartingPos[[i]], TRUE, FALSE)
    StartingPos[[i]] <- NewStart

    # Place the new k-mer back in the proposal matrix
    ProposalMat[i, ] <- Seqs[[i]][NewStart:(NewStart + (MotifLen-1))]
    
    if(TrackPerformance){
      IterationScore <- c(IterationScore, SumOfPairs(ProposalMat))
    }
  }
  
  BaseProbs <- apply(ProposalMat, 2, function(x) CountBase(x, offset = 0))
  
  IterationScoredf <- data.frame(Iteration = 1:length(IterationScore), Score = IterationScore)
  
  FinalResults <- list(Motif = ProposalMat, StartingPositions = StartingPos, SequenceProbabilities = SeqProbs, BaseProbs = BaseProbs, ScoreTrace = IterationScoredf)
  FinalResults
}



# In Parallel, Repeat sampling many times to try to find better alignments
# Skip this and go to the next section to do the same task sequentially 
nCores <- detectCores() - 1 # Decrease this number if running in parallel makes the computer run too slow
if(nCores == 0){
  nCores <- 1
}
if(getDoParWorkers() != nCores){
  cl <- makeCluster(nCores)
  registerDoSNOW(cl)
}

# How many times to run the Gibbs sampler
iterations <- 50

# A progress bar
pb       <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts     <- list(progress = progress)

# Collect many samples
system.time({MultipleSamples <- foreach(i = 1:iterations, .options.snow = opts) %dopar%{
  GibbsSample <- GibbsMotif2()
}
})

# Remove the cluster
registerDoSEQ()

###################################################################################
# SEQUENTIAL IMPLEMENTATION 
# TAKES A LONG TIME, SKIP THIS SECTION IF ALREADY RUN IN PARALLEL

# The number of iterations can be reduced down to around 100-150 and 
# still have a very good chance of giving a good motif (See the end of the script)

iterations <- 50
MultipleSamples <- vector("list", length = iterations)

pb <- txtProgressBar(max = iterations, style = 3)

system.time({for(i in 1:iterations){
  MultipleSamples[[i]] <- GibbsMotif2()
  setTxtProgressBar(pb, i)
}})
###################################################################################



# Score each motif alignment
MaxScoreMotifs <- pbsapply(MultipleSamples, function(x) SumOfPairs(x[[1]], MisMatch = -2))

# Which alignments have the maximum score from the current set
MaxScoreMotifs <- which(MaxScoreMotifs == max(MaxScoreMotifs))

# Compare these max scoring motifs to see if there are different max scoring matrices
# Pull the starting positions for each of these
MaxStartingPositions <- matrix(ncol = 10, nrow = length(MaxScoreMotifs))
for(i in 1:nrow(MaxStartingPositions)){
  MaxStartingPositions[i, ] <- MultipleSamples[[MaxScoreMotifs[[i]]]]$StartingPositions
}

# Check if there are multiple "best" starting positions that scored the same
FirstStartingPos <- MaxStartingPositions[1, ]
CheckDiff        <- apply(MaxStartingPositions, 1, function(x) all(x == FirstStartingPos))

# Return the best motif(s), with starting positions and nucleotide frequencies
FinalMotif <- c()
if(all(CheckDiff)){
  FinalMotif <- MultipleSamples[[MaxScoreMotifs[[1]]]]
}else{
  nUnique <- c(1, which(!CheckDiff))
  FinalMotif <- vector("list", length = length(nUnique))
  for(i in 1:length(FinalMotif)){
    FinalMotif[[i]] <- MultipleSamples[[MaxScoreMotifs[[nUnique[[i]]]]]]
  }
}

# The alignment of the best motif, the starting position in each of the provided sequence
# The probability of each starting position, and the frequency of each base in the starting position
FinalMotif


# Extra section...... 
# How many iterations should you run this algorithm to find the "max" result?
# Model as a negative binomial where a "success" is if the motif produced is the
# maximum value for the whole set

# The scores for every iteration of the Gibbs sampler
MotifScores  <- pbsapply(MultipleSamples, function(x) SumOfPairs(x[[1]], MisMatch = -2))
MotifSuccess <- MotifScores == max(MotifScores) # Which produced a success (max score)?

p_success <- mean(MotifSuccess) # TRUE = 1

nbinomialfn <- function(q){
  qnbinom(q, 1, p_success)
}

# Plot of the quantile function
curve(nbinomialfn, from = 0.5, to = 1-0.0001)

# How many iterations for a 99.9999% chance of converging to the max motif?
qnbinom(0.999999, 1, p_success) # 147 (with 2000 samples)


# Can probably run the sampler fewer times (around 150 iterations), and still get the best possible result
# ......or at least the same "best" result as you'd get after running it 2000 times

