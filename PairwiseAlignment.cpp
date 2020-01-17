/* Implementation of the Needleman-Wunsch algoithm for pairwise sequence alignment in C++.
   The function takes a FASTA file as it's only argument. The FASTA file should only contain two sequences, otherwise
   the alignment between the last two sequences of the file will be returned. */


#include <iostream>
#include <string>
#include <fstream>

int main ( int argc, char **argv )
{
    /* Get the two sequences to align from a FASTA file. 
    I borrowed heavily from this code http://rosettacode.org/wiki/FASTA_format#C.2B.2B for this portion. 
    In fact, the only alteration I made was to send the two sequences to the seq1 and seq2 strings instead
    of printing them to the console as in the original code.*/
        if( argc <= 1 ){
        std::cerr << "Usage: "<<argv[0]<<" [FASTA file]" << std::endl;
        return -1;
    }
 
    std::ifstream input(argv[1]);
    if(!input.good()){
        std::cerr << "Error opening '"<<argv[1]<<"'. Bailing out." << std::endl;
        return -1;
    }
 
    std::string line, name, content, seq1, seq2;
    while( std::getline( input, line ).good() ){
        if( line.empty() || line[0] == '>' ){ // Identifier marker
            if( !name.empty() ){ // Read the first sequence into seq1
                seq1 = content;
                name.clear();
            }
            if( !line.empty() ){
                name = line.substr(1);
            }
            content.clear();
        } else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){ // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            } else {
                content += line;
            }
        }
    }
    if( !name.empty() ){ // Read the second sequence into seq2
        seq2 = content;
    }
 
    std::cout << "-Sequences read from FASTA file-" << '\n' << "Seq1: " << seq1 << std::endl << "Seq2: " << seq2 << std::endl << '\n';

  // Begin implementation of the Needleman-Wunsch algorithm to align the two sequences from the FASTA file

  // Penalties
  int match = 0;
  int fiveprimegap = 8;
  int threeprimegap = 7;
  int gap = 11;
  int mismatch = 4;

  // An array to hold the scores of dimensions m * n where m is the number of bases in seq1, and n is the number of bases in seq2
  int ScoreArray[seq1.length() + 1][seq2.length() + 1];

  // An array to store the direction which was taken to arrive at a particular cell
  char DirectionArray[seq1.length() + 1][seq2.length() + 1];
  
  /* Initialize the score array so that the first row and column correspond to the blank scores. The first row is initialized with 
     values corresponding to the 5' gap penalty*/
  for (unsigned i=0; i<seq1.length()+1; ++i)
  {
      ScoreArray[i][0] = i*gap;
      DirectionArray[i][0] = 'U';
  }

  for (unsigned i = 0; i<seq2.length()+1; ++i)
  {
      ScoreArray[0][i] = i*fiveprimegap;
      DirectionArray[0][i] = 'L';
  }


  // Fill in the rest of the score array
  for (unsigned row=1; row<seq1.length() + 1; ++row)
  {
      for (unsigned column=1; column<seq2.length() + 1; ++column)
      {
          int DiagScore;
          int LeftScore;
          int UpScore;
          int MinScore;

          // The score of moving diagonally to a position (match or mismatch)
          DiagScore = seq1.at(row - 1)==seq2.at(column - 1) ? ScoreArray[row - 1][column - 1] + match : ScoreArray[row - 1][column - 1] + mismatch;

          // The score of moving from the cell to the left
          LeftScore = ScoreArray[row][column - 1] + gap;

          // The score of moving down from the top cell
          UpScore = ScoreArray[row - 1][column] + gap;

          // Select the minimum score and input it into the score array
          MinScore = std::min(DiagScore, std::min(LeftScore, UpScore));
          ScoreArray[row][column] = MinScore;

          // Fill in the appropriate direction in the direction array
          if(DiagScore == MinScore){
              DirectionArray[row][column] = 'D';
          } else if(LeftScore == MinScore)
          {
              DirectionArray[row][column] = 'L';
          }else{
              DirectionArray[row][column] = 'U';
          }
      }
  }

    // Print the score array
/*   for(unsigned i=0; i<seq1.length()+1; ++i)
  {
      for(unsigned j=0; j<seq2.length()+1; ++j)
      {
          std::cout << ScoreArray[i][j] << ' ';
      }
      std::cout << "\n";
  } */

  /* I believe this is a pretty inelegant way to go about doing this next part, but anyway.... what I want to do now is calculate the 
  scores for the last row of the array as if they had incurred a three-prime gap penalty. If any cell in the last row has a
  "corrected" score which takes this penalty into account that is still less than the score of the bottom right corner, I will
  use that cell as the starting point for the traceback instead.*/
  int CorrectedScores[seq2.length() + 1];
  for(unsigned i=0; i<seq2.length() + 1; ++i)
  {
      CorrectedScores[i] = ScoreArray[seq1.length()][i] + ((seq2.length() - i)*threeprimegap);
  }

  // The "uncorrected" optimal score
  int OptimalScore_start = ScoreArray[seq1.length()][seq2.length()];
  
  // Find the smallest element of the corrected scores. This is the column where the traceback will begin
  int SmallestIndex = seq2.length();
  int SmallestVal = OptimalScore_start;
  for(unsigned i=0; i<seq2.length(); ++i)
  {
      if(CorrectedScores[i] < SmallestVal)
      {
          SmallestIndex = i;
          SmallestVal = CorrectedScores[i];
      }
  }

// Print the direction array
/*   DirectionArray[0][0] = 'O';

  for(unsigned i=0; i<seq1.length()+1; ++i)
  {
      for(unsigned j=0; j<seq2.length()+1; ++j)
      {
          std::cout << DirectionArray[i][j] << ' ';
      }
      std::cout << "\n";
  }

  std::cout << '\n'; */

  // Trace back using the direction matrix to produce the alignment
  int CurrentRow = seq1.length();
  int CurrentCol = SmallestIndex;

  std::string Align1 = "";
  std::string Align2 = "";
  
  // Number of 3' gaps
  int GapNum = seq2.length() - SmallestIndex;

  // Check if there is a 3' gap and add it if there is
  if(SmallestIndex < seq2.length())
  {
      for(unsigned i=0; i < GapNum; ++i)
      {
          Align1.insert(0, 1, '-');
          Align2.insert(0, 1, seq2.at(seq2.length() - i - 1));
      }
  }

  while(CurrentRow >= 1 | CurrentCol >= 1)
  {
      if(DirectionArray[CurrentRow][CurrentCol] == 'D') // Moving diagonally
      {
          Align1.insert(0, 1, seq1.at(CurrentRow - 1));
          Align2.insert(0, 1, seq2.at(CurrentCol - 1));
          --CurrentRow;
          --CurrentCol;
      }
      else if(DirectionArray[CurrentRow][CurrentCol] == 'U') // Moving up
      {
          Align1.insert(0, 1, seq1.at(CurrentRow - 1));
          Align2.insert(0, 1, '-');
          --CurrentRow;
      }
      else // Moving left
      {
          Align1.insert(0, 1, '-');
          Align2.insert(0, 1, seq2.at(CurrentCol - 1));
          --CurrentCol;
      }
  }

  // Print the optimal alignment and its score
  std::cout << "-Optimal Alignment-"<< '\n' << Align1 << '\n' << Align2 << '\n';
  std::cout << "Alignment score: " << SmallestVal << '\n';

  return 0;
}
