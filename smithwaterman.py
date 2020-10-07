###############################################################################
#   Author: Francisco Esteves       Group: 3
#   Description: Smith-Waterman script in Python3.
#   Input: amino acid 1 (str); amino acid 2 (str); gap penalty (int).
#   Output: score of the optimal alignment(s) (int);
#           combination of the optimal alignment(s) (str).
###############################################################################

# imports
from os import system, name
from Bio.SubsMat import MatrixInfo
import time


# class that represents the complete Smith-Waterman matrix
# smith_matrix has the score of each point
# trace_back_matrix has the "arrows" that point to the previous combination
class Matrices:
    def __init__(self, s1, s2):
        self.smith_matrix = [[0 for i in range(len(s1)+1)] for j in range(len(s2)+1)]
        self.trace_back_matrix = [['X' for i in range(len(s1)+1)] for j in range(len(s2)+1)]


# clear the screen function
def clear():
    # windows
    if name == 'nt':
        system('cls')
    # mac/linux
    else:
        system('clear')


# score of given combination in the blosum50 matrix
# input: (char) A, (char) B
# returns: (int) the value of combination
def combination_score(a, b):
    for i, j in MatrixInfo.blosum50:
        if i == a and j == b:
            return(MatrixInfo.blosum50[(a, b)])
        else:
            if i == b and j == a:
                return(MatrixInfo.blosum50[(b, a)])
    # Error
    print("Error. Amino acid combination not in blosum50 matrix."
            + "\nMake sure you used valid capital letters\n.")
    exit(-1)


# fill the smith_matrix according to the Smith-Waterman and blosum50 matrix;
# fill the trace_back_matrix as follows:
#   'D' -> points in the diagonal
#   'U' -> points up
#   'L' -> points left
#    0  -> doesn't point anywhere
# input: (str) s1, (str) s2, (int) gap_pen
# returns: (object Matrix) matrices
def smith_algorithm(s1, s2, gap_pen):

    start_time = time.time()
    # create a Matrix object
    matrices = Matrices(s1, s2)
    # fills the matrices row by row
    for i in range(len(s2)):
        for j in range(len(s1)):
            combine_vals = combination_score(s1[j], s2[i]) + matrices.smith_matrix[i][j]
            gap_in_s1 = matrices.smith_matrix[i+1][j] + gap_pen
            gap_in_s2 = matrices.smith_matrix[i][j+1] + gap_pen
            # evaluates the best combination and where to point
            if (combine_vals >= gap_in_s2) and (combine_vals >= gap_in_s1):
                matrices.smith_matrix[i+1][j+1] = combine_vals
                matrices.trace_back_matrix[i+1][j+1] = 'D'
                if combine_vals < 0:
                    matrices.smith_matrix[i+1][j+1] = 0
                    matrices.trace_back_matrix[i+1][j+1] = 'X'
            elif (gap_in_s1 >= combine_vals) and (gap_in_s1 >= gap_in_s2):
                matrices.smith_matrix[i+1][j+1] = gap_in_s1
                matrices.trace_back_matrix[i+1][j+1] = 'L'
                if gap_in_s1 < 0:
                    matrices.smith_matrix[i+1][j+1] = 0
                    matrices.trace_back_matrix[i+1][j+1] = 'X'
            else:
                matrices.smith_matrix[i+1][j+1] = gap_in_s2
                matrices.trace_back_matrix[i+1][j+1] = 'U'
                if gap_in_s2 < 0:
                    matrices.smith_matrix[i+1][j+1] = 0
                    matrices.trace_back_matrix[i+1][j+1] = 'X'

    # use to print the matrices
    #for i in range(len(matrices.smith_matrix)):
    #    print(matrices.smith_matrix[i])
    #print("")
    #for i in range(len(matrices.trace_back_matrix)):
    #    print(matrices.trace_back_matrix[i])

    # finds the highest value in the matrix & how many times it repeasts itself
    highest = 0
    position = [[0,0]]
    for i in range(len(s2)):
        for j in range(len(s1)):
            if matrices.smith_matrix[i+1][j+1] == highest:
                position.append([i+1,j+1])
            elif matrices.smith_matrix[i+1][j+1] > highest:
                highest = matrices.smith_matrix[i+1][j+1]
                position = [[i+1,j+1]]

    # runs through the trace_back_matrix while writting the combination
    k = 0
    for i, j in position:
        k += 1
        comb_s1 = []
        comb_s2 = []
        while True:
            if matrices.trace_back_matrix[i][j] == 'D':
                comb_s1.insert(0, s1[j-1])
                comb_s2.insert(0, s2[i-1])
                i = i-1
                j = j-1
            elif matrices.trace_back_matrix[i][j] == 'U':
                comb_s1.insert(0, '-')
                comb_s2.insert(0, s2[i-1])
                i = i-1
            elif matrices.trace_back_matrix[i][j] == 'L':
                comb_s1.insert(0, s1[j-1])
                comb_s2.insert(0, '-')
                j = j-1
            else:
                break

        print("\n\nSolution n" + str(k) + ":\n")

        # i and j are the positions in the matrix where the alignment starts
        # if s1 starts it's comb before s2
        for a in range(i-j):
            print(" ", end='')
        # print the first letters that do not belong in the comb
        for a in range(j):
            print(s1[a], end='')
        # print the comb
        print("|", end="")
        for a in range(len(comb_s1)):
            print(comb_s1[a], end='')
        print("|", end="")
        # remove the '-' from the comb for the next if
        for hifen in comb_s1:
	           if(hifen == '-'):
		                 comb_s1.remove('-')
        # if there still are letters to print after the comb, print them
        if j-1+len(comb_s1) < len(s1):
            print(s1[j+len(comb_s1):], end='')
        print("")

        # if s1 starts it's comb before s2
        for a in range(j-i):
            print(" ", end='')
        # print the first letters that do not belong in the comb
        for a in range(i):
            print(s2[a], end='')
        # print the comb
        print("|", end="")
        for a in range(len(comb_s2)):
            print(comb_s2[a], end='')
        print("|", end="")
        # remove the '-' from the comb for the next if
        for hifen in comb_s2:
	           if(hifen == '-'):
		                 comb_s2.remove('-')
        # if there still are letters to print after the comb, print them
        if i-1+len(comb_s2) < len(s2):
            print(s2[i+len(comb_s2):], end='')
        print("")

    print("\nValue of alignment is: " + str(highest) + "!\n")

    print("--- %s seconds ---" % (time.time() - start_time))

    return matrices


def main():

    # clears the screen
    clear()
    # input of variables
    s1 = str(input("First amino acid sequence: "))
    s2 = str(input("Second amino acid sequence: "))
    gap_pen = input("Gap penalty: ")
    try:
        gap_pen = int(gap_pen)
        # always negative
        if gap_pen > 0:
            gap_pen = -gap_pen
    except ValueError:
        print("Error. Gap penalty not in integer format.")
        exit(-1)

    # compute the WHOLE smith-waterman algorithm
    matrices = smith_algorithm(s1, s2, gap_pen)

    return(0)


if __name__ == '__main__':
    main()
