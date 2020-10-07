# Smith-Waterman Algorithm
Program in Python3 that implements Smith-Waterman algorithm based on the blosum50 scoring matrix.
Made for IST Computational Biology subject 2020/2021.

## Requirements
requires the following libraries:
- os
- [Bio.SubsMat](https://biopython.org/wiki/Download)
- time

## Usage
1º insert the first amino acid sequence in all caps lock; 
2º insert the second amino acid sequence in all caps lock; 
3º insert the gap penalty (the program will always consider the integer as negative); 
4º outputs the solution(s) with the sequences where the best local combination is highlighted, the value of the best combination and the time it took to execute.

Example:

First amino acid sequence: WPIWPC 
Second amino acid sequence: IIWPI 
Gap penalty: 4 

Solution n1: 
WP|IWP|C 
 I|IWP|I 
 
Solution n2: 
  |WPI|WPC 
II|WPI| 

Value of alignment is: 30! 
--- 0.0026280879974365234 seconds ---

## Credits
Program made by Francisco Esteves.

Tested by Diogo Gonçalves, Clémentine Abel and Gonçalo Damas.
