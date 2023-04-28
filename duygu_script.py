# global alignment - (Needleman-Wunch)
# local alignment - (Smith-Waterman)

from numpy import *
import numpy as np
import sys

all_aas_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', "*"]


# create a zero matrix with given (sequence) lengths
def make_0matrix(len_seq1, len_seq2):
    return [[0] * len_seq1 for i in range(len_seq2)]


# gets 2 sequences from the file
def get_strings(file_name):
    with open(file_name, "r") as file:
        lines = file.readlines()
    first_seq, second_seq = lines[0].strip("\n "), lines[1].strip("")
    return first_seq, second_seq


# read matrix into a dictionary
def read_matrix(filename):
    # reads the substitution matrix with the filename into a dictionary for ex. BLOSUM62
    sm = {}
    f = open(filename, "r")
    line = f.readline()
    line = line.strip()
    tokens = line.split("  ")
    size = len(tokens)
    alphabet = []
    for i in range(0, size):
        alphabet.append(tokens[i][0])
    for i in range(0, size):
        line = f.readline()
        line = line.strip()
        tokens = line.split(" ")
        while '' in tokens:
            tokens.remove('')
        tokens = tokens[1:]
        for j in range(0, len(tokens)):
            k = alphabet[i] + alphabet[j]
            sm[k] = int(tokens[j])
    return sm


# get the scores of two nucleotides from the dictionary that takes the scores from input score file
def get_score_fom_matrix(a, b, matrix_dict):
    assert len(a) == len(b) == 1  # check if they are both characters
    if a not in all_aas_list:
        a = "*"
    if b not in all_aas_list:
        b = "*"
    concat_str = a + b

    return matrix_dict[concat_str]


def return_str(str2, str1, table, param, index1, index2):
    a1 = ""
    a2 = ""
    i = index1
    j = index2
    while table[i][j] != param:  # while i > 0 or j > 0:

        if table[i][j] == "d":  # diagonal move

            a1 = str1[i - 1] + a1
            a2 = str2[j - 1] + a2
            i = i - 1
            j = j - 1
        elif table[i][j] == "s":

            a1 = str1[i - 1] + a1
            a2 = "-" + a2
            i = i - 1
        elif table[i][j] == "e":

            a1 = "-" + a1
            a2 = str2[j - 1] + a2
            j = j - 1


    return a2, a1


def glob_return_str(str2, str1, table):

    a1 = ""
    a2 = ""
    i = len(str1)
    j = len(str2)
    while table[i][j] != 0:  # while i > 0 or j > 0:

        if table[i][j] == "d":  # diagonal move
            a1 = str1[i - 1] + a1
            a2 = str2[j - 1] + a2
            i = i - 1
            j = j - 1
        elif table[i][j] == "s":
            a1 = str1[i - 1] + a1
            a2 = "-" + a2
            i = i - 1
        elif table[i][j] == "e":
            a1 = "-" + a1
            a2 = str2[j - 1] + a2
            j = j - 1


    return a2, a1


def take_args():

    arg_list = sys.argv
    if len(arg_list) < 6:
        print("You gave incomplete set of arguments, please rerun the program from terminal with all arguments.")
        return 0
    else:
        input_txt = arg_list[1]
        align_algo = arg_list[2]  # It can be "local" or "global"
        score_matrix = arg_list[3]
        gap_op = int(arg_list[4])
        gap_ex = int(arg_list[5])
        # initialising substitution matrix
        subs_dict = read_matrix(score_matrix)
        # initialising sequences
        seq1, seq2 = get_strings(input_txt)

        # len_seq1, len_seq2 = len(seq1), len(seq2)

        if align_algo == "global":
            matrix_m = make_0matrix(len(seq1) + 1, len(seq2) + 1)  # no gap
            matrix_x = make_0matrix(len(seq1) + 1, len(seq2) + 1)  # gap in x
            matrix_y = make_0matrix(len(seq1) + 1, len(seq2) + 1)  # gap in y
            backtrack = make_0matrix(len(seq1) + 1, len(seq2) + 1)  # gap in y


            # initialising 3 +1 matrices to keep track of no gap or gap in x,y in affine alignment

            matrix_m[0][1], matrix_m[1][0] = gap_op, gap_op
            matrix_x[0][1], matrix_y[1][0] = gap_op, gap_op


            for i in range(2, len(matrix_m)):
                matrix_m[i][0] = matrix_m[i-1][0] + gap_ex
                matrix_y[i][0] = matrix_y[i-1][0] + gap_ex  # gaps in seq2

            for j in range(2, len(matrix_m[0])):

                matrix_m[0][j] = matrix_m[0][j-1] + gap_ex
                matrix_x[0][j] = matrix_x[0][j-1] + gap_ex

            for j in range(len(matrix_m)):
                matrix_x[j][0] = -Infinity

            for j in range(len(matrix_m[0])):
                matrix_y[0][j] = -Infinity

            for i in range(1, len(backtrack)):
                backtrack[i][0] = "s"  # gap_op + i * gap_ex
            for j in range(1, len(backtrack[0])):
                backtrack[0][j] = "e"  # gap_op + j * gap_ex


            for i in range(1, len(matrix_m)):  # len(seq1) + 1):
                for j in range(1, len(matrix_m[i])):  # len(seq2) + 1)
                    matrix_m[i][j] = get_score_fom_matrix(seq1[j - 1], seq2[i - 1], subs_dict) + \
                                     max(matrix_m[i - 1][j - 1], matrix_x[i - 1][j - 1], matrix_y[i - 1][j - 1]
                                         )
                    # upper table
                    matrix_x[i][j] = max(
                        gap_op + matrix_m[i][j - 1],
                        gap_ex + matrix_x[i][j - 1]
                        # gap_op + gap_ex + matrix_y[i][j - 1]
                    )
                    # lower table
                    matrix_y[i][j] = max(
                        gap_op + matrix_m[i - 1][j],  # gap_op + gap_ex + matrix_m[i - 1][j]
                        # gap_op + gap_ex + matrix_x[i - 1][j],
                        gap_ex + matrix_y[i - 1][j]
                    )
                    optimum_one = max(matrix_m[i][j], matrix_x[i][j],
                                      matrix_y[i][j])
                    if optimum_one == matrix_m[i][j]:
                        backtrack[i][j] = "d"
                    elif optimum_one == matrix_y[i][j]:
                        backtrack[i][j] = "s"
                    elif optimum_one == matrix_x[i][j]:
                        backtrack[i][j] = "e"

            return_seq_1, return_seq2 = glob_return_str(seq1, seq2, backtrack)

        elif align_algo == "local":
            matrix_m = make_0matrix(len(seq1) + 1, len(seq2) + 1)  # no gap
            matrix_x = make_0matrix(len(seq1) + 1, len(seq2) + 1)  # gap in x
            matrix_y = make_0matrix(len(seq1) + 1, len(seq2) + 1)  # gap in y
            backtrack = make_0matrix(len(seq1) + 1, len(seq2) + 1)

            # initialising 3 +1 matrices to keep track of no gap or gap in x,y in affine alignment
            matrix_m[0][1], matrix_m[1][0] = 0, 0
            matrix_x[0][1], matrix_y[1][0] = 0, 0
            backtrack[0][0] = 0

            for i in range(2, len(matrix_m)):
                matrix_m[i][0] = 0
                matrix_y[i][0] = 0 # gaps in seq2

            for j in range(2, len(matrix_m[0])):
                matrix_m[0][j] = 0
                matrix_x[0][j] = 0

            for j in range(len(matrix_m)):
                matrix_x[j][0] = 0

            for j in range(len(matrix_m[0])):
                matrix_y[0][j] = 0

            for i in range(1, len(backtrack)):
                backtrack[i][0] = 0  # gap_op + i * gap_ex
            for j in range(1, len(backtrack[0])):
                backtrack[0][j] = 0 # gap_op + j * gap_ex

            for i in range(1, len(matrix_m)):  # len(seq1) + 1):
                for j in range(1, len(matrix_m[i])):  # len(seq2) + 1)
                    matrix_m[i][j] = get_score_fom_matrix(seq1[j - 1], seq2[i - 1], subs_dict) + \
                                     max(matrix_m[i - 1][j - 1], matrix_x[i - 1][j - 1], matrix_y[i - 1][j - 1]
                                         )
                    # upper table
                    matrix_x[i][j] = max(
                        gap_op + matrix_m[i][j - 1],
                        gap_ex + matrix_x[i][j - 1]
                        # gap_op + gap_ex + matrix_y[i][j - 1]
                    )
                    # lower table
                    matrix_y[i][j] = max(
                        gap_op + matrix_m[i - 1][j],  # gap_op + gap_ex + matrix_m[i - 1][j]
                        # gap_op + gap_ex + matrix_x[i - 1][j],
                        gap_ex + matrix_y[i - 1][j]
                    )
                    optimum_one = 0
                    opt = max(matrix_m[i][j], matrix_x[i][j],
                                      matrix_y[i][j])
                    if opt > optimum_one:
                        optimum_one = opt

                    if optimum_one == matrix_m[i][j]:
                        backtrack[i][j] = "d"
                    elif optimum_one == matrix_y[i][j]:
                        backtrack[i][j] = "s"
                    elif optimum_one == matrix_x[i][j]:
                        backtrack[i][j] = "e"


            list1 = np.array(matrix_x)
            list2 = np.array(matrix_y)
            list3 = np.array(matrix_m)
            max1, max2, max3 = list1.max(), list2.max(), list3.max()
            max_numb = max(max1, max2, max3)
            if max_numb == max1:
                index1, index2 = np.argwhere(matrix_x == max1)[0]
            elif max_numb == max2:
                index1, index2 = np.argwhere(matrix_y == max2)[0]
            elif max_numb == max3:
                index1, index2 = np.argwhere(matrix_m == max3)[0]

            return_seq_1, return_seq2 = return_str(seq1, seq2, backtrack, 0, index1, index2)

        a = len(return_seq_1)
        str3 = ""

        total_score = 0
        matches = 0
        for i in (range(a)):
            if i < 1:
                if return_seq_1[i] == "-":
                    str3 += " "
                    total_score += gap_op
                elif return_seq2[i] == "-":
                    str3 += " "
                    total_score += gap_op
                elif return_seq_1[i] == return_seq2[i]:
                    matches +=1
                    str3 += "|"
                    total_score += get_score_fom_matrix(return_seq_1[i],return_seq2[i],subs_dict)
                else:
                    total_score += get_score_fom_matrix(return_seq_1[i], return_seq2[i], subs_dict)
                    str3 += " "
            elif i >= 1:
                if return_seq_1[i] == "-":
                    if return_seq_1[i-1] == "-":
                        str3 += " "
                        total_score += gap_ex
                    else:
                        str3 += " "
                        total_score += gap_op
                elif return_seq2[i] == "-":
                    if return_seq2[i - 1] == "-":
                        str3 += " "
                        total_score += gap_ex
                    else:
                        str3 += " "
                        total_score += gap_op
                elif return_seq_1[i] == return_seq2[i]:
                    matches += 1
                    str3 += "|"
                    total_score += get_score_fom_matrix(return_seq_1[i], return_seq2[i], subs_dict)
                else:
                    total_score += get_score_fom_matrix(return_seq_1[i], return_seq2[i], subs_dict)
                    str3 += " "
        print(return_seq_1)
        print(str3)
        print(return_seq2)
        print("Alignment score: " + str(float(total_score)))
        print("Identity value: " + str(matches) + "/" + str(len(return_seq_1)) +" "+ str(100*(float(matches)/float(len(return_seq_1)))) +"%")



take_args()
