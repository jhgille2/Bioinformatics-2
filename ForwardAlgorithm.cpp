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

     

   /* Various function tests to see if I'm getting the output I expect to
   double TestX = 0.3;
   double TestY = 0.2;
   double TestNaN = std::numeric_limits<double>::quiet_NaN();

   std::cout << "Test of extended exponential " << eexp(TestX) << std::endl; // Should be 1.34985....
   std::cout << "Test of extended logarithm " << eln(TestX) << std::endl;    // Should be -1.20397.....
   std::cout << "Test of extended logarithm sum " << elnsum(TestX, TestY) << std::endl; // Should be eln(0.3) + eln(1+exp(eln(y) - eln(x))) = -1.20397 + eln(1 + exp(-1.6094 + 1.20397)) = -0.69313033284...
   std::cout << "Test of extended logarithm product " << elnproduct(TestX, TestY) << std::endl; // Should be eln(x) + eln(y) = -1.20397 - 1.6094 = -2.81337

   std::cout << "Test of extended sum including a NaN " << elnsum(TestX, TestNaN) << std::endl; // Should be -1.20397....
   std::cout << "Test of extended sum including a NaN " << elnproduct(TestX, TestNaN) << std::endl; // Should be NaN

   std::cout << "Tests of isnan (should return TRUE) " << isnan(TestNaN) << std::endl;

   std::cout << "Multiplication Test " << elnproduct(-1.20397, -0.0202027) << std::endl;
   std::cout << "Sum Test: " << elnsum(std::numeric_limits<double>::quiet_NaN(), elnproduct(0.3, 0.98)) << std::endl; 

   if(isnan(TestX) || isnan(TestNaN))
   {
       std::cout << "At least one value of TestX OR TestNaN is NaN " << std::endl;
       if(isnan(TestX))
       {
           std::cout << "TestX is NaN" << std::endl;
       }else if(isnan(TestNaN))
       {
           std::cout << "TestNaN is NaN" << std::endl;
       }
   } */
   


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
   double EmissionProbs[2][4] = {{0.2, 0.2, 0.3, 0.3}, {0.3, 0.3, 0.2, 0.2}};

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

   // Initialization step
   RecursionArray[0][0] = elnproduct(eln(1), eln(0.2)); // P(y1 = AT) = 1 (as given in the homework)
   RecursionArray[1][0] = eln(0); // P(y1 = GC) = 1-P(y1 = AT) = 0

   // Recursion step
   // For each observed base
   for(unsigned i = 1; i < Seq.length(); ++i)
   {
       // For each state
       for(unsigned j = 0; j < nStates; ++j)
       {
           double logalpha = std::numeric_limits<double>::quiet_NaN();
           int EmissionIndex = EmissionHelper.find(Seq.at(i));

           for(unsigned k=0; k < nStates; ++k)
           {
               // For each previous state
               logalpha = elnsum(logalpha, elnproduct(RecursionArray[k][i-1], eln(TransitionProbs[k][j])));
           }
           RecursionArray[j][i] = elnproduct(logalpha, eln(EmissionProbs[j][EmissionIndex]));
       }
   }

   /*for(unsigned i = 0; i < nStates; ++i)
   {
       for(unsigned j = 0; j < Seq.length(); ++j)
       {
           std::cout << RecursionArray[i][j] << " ";
       }
        std::cout << std::endl;
   }*/

   // Sum last column of the recusion array to calculate the final probability of the sequence
    double FinalLikelihood = 0;
    for(unsigned i = 0; i < nStates; ++i)
    {
        FinalLikelihood += RecursionArray[i][Seq.length() - 1];
    }

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


   // Print the final likelihood
   std::cout << std::endl << "Sequence Log-likelihood: " << FinalLikelihood << std::endl;
}
