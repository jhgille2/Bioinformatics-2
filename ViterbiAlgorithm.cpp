/* Implementation of the forward algorithm for calculating the log-likelihood
   of a given sequence, with provided HMM parameters. The HMM is relatively simple
   with only two hidden states: AT-rich and GC-rich. Each region emits bases with 
   different probabilities (AT rich is more likely to emit A or T and vice versa) */

#include <iostream>
#include <string>
#include <fstream>
#include <math.h> 
#include <limits> 

/* For this algorithm, I will do calculations in log-space to prevent numerical underflow errors in the recursion portion for long sequences.
   The link to a review about these calculations is posted below.

   http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf.
              
   I will implement the pseudocode described in the review, and use the same notation in my calculations*/

   // Extended exponential
   double eexp(double x)
   {
       double result;
       if(isnan(x))
       {
           result = 0;
       }else{
           result = exp(x);
       }
       return(result);
   }

   // Extended logarithm
   double eln(double x)
   {
       double result;
       if(x == 0)
       {
            result = std::numeric_limits<double>::quiet_NaN();
       }else if(x > 0)
       {
           result = log(x);
       }else{
           std::cout << std::endl << "Negative Input" << std::endl;
           return(0);
       }
       return(result);
   }

   // Extended logarithm sum
   // For this to work properly, all inputs have to be logarithms already e.g. x = eln(x), y = eln(y)
   double elnsum(double x, double y)
   {
       double result;

       if(isnan(x) || isnan(y))
       {
           if(isnan(x))
           {
               result = y;
           }else{
               result = x;
           }
       }else{
           if(x > y)
           {
               result = x + eln(1 + exp(y - x));
           }else{
               result = y + eln(1 + exp(x - y));
           }
       }
    return(result);
   }

   // Extended logarithm product
   // Same as the sub, this function assumes that inputs are logarithms
   double elnproduct(double x, double y)
   {
       double result;
       
       if(isnan(x) || isnan(y))
       {
           result = std::numeric_limits<double>::quiet_NaN();
       }else{
           result = x + y;
       }
       return(result);
   }


int main ( int argc, char **argv )
{
    // Read in the sequence from a text file
    /*
    std::ifstream file(argv[1]);
    std::string str;
    std::string Seq;

    while(std::getline(file, str))
    {
        Seq += str;
        Seq.push_back('\n');
    }

    std::cout << "Sequence read: " << Seq << std::endl;
    */

   // Need to fix file read. for now hardcode the sequence
   std::string Seq="GCGAGTTCCCCCATGCGTGTCGGCCCCCGCGTCGCTTTTACTATATCCACCCCATTCCATTGAAGGTATTTACTTGGGATCAAATACCGACTAGAGTTATTAAATCTTAATGACCTATCCTAAGTTTAAGTAACCCAGGTACAGTCGTTAGCTTAGTTCCAAGCGTCCCTCACGTGCACTACGGTATTCCCCCTCGTCATGTGCAACCCCCCCGCTCAGCGGATTTCCAGCCGGCCAACAGTTCGACGAGAAGGCATACCCGCCAGGCAACCTACCGGCCACTCCCTGCGCCCGGACCTTTGACTTGAAACTTTTCTTATGTTGTGGGCCCCTGCGCCTCTGTAGGTATGATGGGTTTAAGAGCTTTAGGCCGCCACTATCGGCACCACCATAGTGTAAAACATCGCATAGGCCTGGGCTAAGGGGGGTTTGATGCGGGTGGAGGCGGGTGAGGATAACGGAGCATGGCATACCGTATGTTTTTCCTACGGCGGCTGGGCCCGTAGTCGAAATCTCTAAATATCAACTATAAGGTGAGCGCCAGGGTCGCTGGGTCGAAGCGGCGACAATTTGTGACGGGCTACGCGACGGGTCCCTTATGCGCGAGGCTGTGCGCTGGATTCAACCTGGCGCCTCCATGCATGCCGCGAATCTGTAGTTAATTGCACGGATAATACGCTGGGGGCACCAAGACCTGGTTAACGTTTCATTCCTGCTGTTAAATCAGGTCGTCAGGGACCTAGTGACACCTGCCCGACATCCCGCGAACCTGCCCGGAGTTCACGGATACTATTGGTAATACGTAAAATGCGGGTCGGGACCGTATGGGATTCTCCTAATTAGTCATGTCTCGTAATTCCATATGTGGCTGTTTTACAATACGATTTTAATAACTCTTACACTGGTGTCTCCGGGGGCGTTTCTTCTTACACAATTATACGAAAACAATCTGACGCGATAAATTGTCATCCAGTGTAATATGGCTACGGCGAACGACACCCATTTGGGGAAGTGGTTAAACTGGGTGTATGATAACCCGCTTAGTCTAGCAGGTACACCGTTCGTCCTGCGACTTCTACAGGAGTACCCCTTATTCTTGTGATCTCACAGCGGAAGGCATAGTCCATCGTGCCACCGCTTATCGTGACTACTTCTGAGTAAGTTGTGTTGGTAAAAACCAGCGTATCATGTCGTCTGAAG";
   
   // Number of observations and hidden states
   int nObs = Seq.length();
   int nStates = 2;

   // Transition probabilities 2x2 matrix 
   /*       AT    GC
        AT 0.98  0.02
        GC 0.05  0.95
   */
   double TransitionProbs[2][2] = {{0.51, 0.49}, {0.49, 0.51}};

   // Emission probabilities
   /*    G    C    A    T
     AT 0.2  0.2  0.3 0.3
     GC 0.3  0.3  0.2  0.2
   */
   double EmissionProbs[2][4] = {{0.20, 0.20, 0.30, 0.30}, {0.30, 0.30, 0.20, 0.20}};

   // A string that will be used in the recursion step to pull elements from the emission probability array by matching indexes
   // I am almost positive there must be a better way to do this, but I believe this should work
   std::string EmissionHelper = "GCAT";

   // The recursion Array
   /*
        1   2   3......N
    AT
    GC
   */
   double RecursionArray[2][Seq.length()];

   // Initialize a traceback array of identical dimensions
   int TracebackArray[2][Seq.length()];
   TracebackArray[0][0] = 0;
   TracebackArray[0][1] = 0;

   // Array to hold probabilities for the recusrion step
   double CurrentProbs[2];

   // Initialization step
   RecursionArray[0][0] = eln(1) + eln(EmissionProbs[0][EmissionHelper.find(Seq.at(0))]); // P(y1 = AT) = 1 (as given in the homework)
   RecursionArray[1][0] = eln(0); // P(y1 = GC) = 1-P(y1 = AT) = 0

   // Recursion step
   // For each observed base
   for(unsigned i = 1; i < Seq.length(); ++i)
   {
       // For each state
       for(unsigned j = 0; j < nStates; ++j)
       {
           std::fill_n(CurrentProbs, 2, std::numeric_limits<double>::quiet_NaN());
           int EmissionIndex = EmissionHelper.find(Seq.at(i));

           for(unsigned k=0; k < nStates; ++k)
           {
               // For each previous state
               CurrentProbs[k] = RecursionArray[k][i-1] + eln(TransitionProbs[k][j]) +  eln(EmissionProbs[j][EmissionIndex]);
               /*std::cout << "Current Probability:" <<CurrentProbs[k] << std::endl;
               std::cout << "K: " << k << std::endl;
               std::cout << "J: " << j << std::endl;
               std::cout << "I: " << i << std::endl;
               std::cout << "Emission Prob: " << eln(EmissionProbs[j][EmissionIndex]) << std::endl;
               std::cout << "Transition Prob: " << eln(TransitionProbs[k][j]) << std::endl;
               std::cout << "Gamma: " << RecursionArray[k][i-1] <<std::endl << std::endl;*/
           }

           double MaxProb = 0;
           int MaxIndex = 0;
           for(unsigned m = 0; m < nStates; m++)
           {
               if(m==0)
               {
                   MaxProb = CurrentProbs[m];
               }
               if(CurrentProbs[m] > MaxProb)
               {
                   MaxProb = CurrentProbs[m];
                   MaxIndex = m;
               }
           }

           RecursionArray[j][i] = CurrentProbs[MaxIndex];
           TracebackArray[j][i] = MaxIndex;
       }
   }

   // Print the recursion array
   /*for(unsigned i = 0; i < nStates; ++i)
   {
       for(unsigned j = 0; j < Seq.length(); ++j)
       {
           std::cout << RecursionArray[i][j] << '\t';
       }
        std::cout << std::endl;
   }*/

   // Print the traceback array
   /*std::cout << "Traceback Array" <<std::endl;
   for(unsigned i = 0; i < nStates; ++i)
   {
       for(unsigned j = 0; j < Seq.length(); ++j)
       {
           std::cout << TracebackArray[i][j];
       }
        std::cout << std::endl;
   }*/

   // Final log probability of the most likely sequence is the maximum of the last column of the recursion array
   double FinalProb;
   int FinalIndex = 0;
   for(unsigned T=0; T < nStates; T++)
   {
       if(T == 0)
       {
           FinalProb = RecursionArray[T][Seq.length() - 1];
       }
       if(RecursionArray[T][Seq.length() - 1] > FinalProb)
       {
           FinalProb = RecursionArray[T][Seq.length() - 1];
           FinalIndex = T;
       }
   }

   // Return Final Sequence log-likelihood
   std::cout << "State sequence log-likelihood: " << FinalProb << std::endl;

   // Traceback
   int SequenceTrace[Seq.length() - 1];
   int NextIndex = 0;
   for(unsigned t=Seq.length(); t--;)
   {
       if(t == Seq.length() - 1)
       {
           NextIndex = TracebackArray[FinalIndex][(Seq.length() - 1)];
       }

       SequenceTrace[t] = TracebackArray[NextIndex][t];
       NextIndex = TracebackArray[NextIndex][t];
   }

   // Return the Traceback
   std::cout << "Traceback (0 indicates AT rich state, 1 indicates a GC rich state)" << std::endl;
   for(unsigned i = 0; i < Seq.length(); i++)
   {
       std::cout << SequenceTrace[i];
   }
   std::cout << std::endl;

    // Print Transition Matrix
    std::cout << "Transition Matrix" << std::endl;
    std::cout << '\t' << "AT" << '\t' << "GC" << std::endl;
    for(unsigned i= 0; i < 2; ++i)
    {
        for(unsigned j = 0; j < 2; ++j)
        {
            if(i == 0 & j == 0)
            {
                std::cout << "AT" << '\t' << TransitionProbs[i][j] << '\t';
            }else if(i == 1 & j == 0)
            {
                std::cout << "GC" << '\t' << TransitionProbs[i][j] << '\t';
            }else{
                std::cout << TransitionProbs[i][j] << '\t';
            }    
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // Print emission matrix
    std::cout << "Transition Matrix" << std::endl;
    std::cout << '\t' << "G" << '\t' << "C" << '\t' << "A" << '\t' << "T" << std::endl;
    for(unsigned i= 0; i < 2; ++i)
    {
        for(unsigned j = 0; j < 4; ++j)
        {
            if(i == 0 & j == 0)
            {
                std::cout << "AT" << '\t' << EmissionProbs[i][j] << '\t';
            }else if(i == 1 & j == 0)
            {
                std::cout << "GC" << '\t' << EmissionProbs[i][j] << '\t';
            }else{
                std::cout << EmissionProbs[i][j] << '\t';
            }    
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

}
