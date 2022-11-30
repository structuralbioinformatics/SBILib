'''
Needleman-Wunsch global alignment algorithm

Based upon:
http://www.codesofmylife.com/
as defined in:
http://en.wikipedia.org/wiki/Needleman-Wunsch_algorithm

And adapted to the SBI library
A variation has been added so that initial movement of the alignment can
be penalized differently
'''
import numpy as np
import re

from SimilarityMatrix import SimilarityMatrix as SM
from SeqAli import SeqAli


def needleman_wunsch(sequence1, sequence2, similarity_matrix,
                     gap_init = None, gap_penalty = None):

    rows = len(sequence1) + 1
    cols = len(sequence2) + 1
    sequence1 = re.sub('[xX-]', '*', sequence1)
    sequence2 = re.sub('[xX-]', '*', sequence2)
    fMatrix = np.zeros((rows, cols), float)
    similarityMatrixMap = SM.get_matrix(similarity_matrix)

    if gap_init is None:
        gap_init = similarityMatrixMap['*']['A']
    if gap_penalty is None:
        gap_penalty = similarityMatrixMap['*']['*']

    for i in range(0, rows):
        fMatrix[i][0] = i * gap_init
    for j in range(0, cols):
        fMatrix[0][j] = j * gap_init

    for i in range(1, rows):
        for j in range(1, cols):
            mtch   = fMatrix[i - 1][j - 1]
            mtch  += similarityMatrixMap[sequence1[i - 1]][sequence2[j - 1]]
            delete = fMatrix[i - 1][j] + gap_penalty
            insert = fMatrix[i][j - 1] + gap_penalty
            fMatrix[i][j] = max(mtch, delete, insert)

    return trackBack(fMatrix, sequence1, sequence2,
                     gap_penalty, similarityMatrixMap)


def trackBack(fMatrix, seq1, seq2, gap, similarity_map):
    '''
    Tracks back to create the aligned sequence pair
    '''
    alignedSeq1 = ''
    alignedSeq2 = ''
    i = len(seq1)
    j = len(seq2)
    while i > 0 and j > 0:
        score = fMatrix[i][j]
        diagScore = fMatrix[i - 1][j - 1]
        upScore = fMatrix[i][j - 1]
        leftScore = fMatrix[i - 1][j]
        if score == diagScore + similarity_map[seq1[i - 1]][seq2[j - 1]]:
            alignedSeq1 = seq1[i - 1] + alignedSeq1
            alignedSeq2 = seq2[j - 1] + alignedSeq2
            i -= 1
            j -= 1
        elif score == leftScore + gap:
            alignedSeq1 = seq1[i - 1] + alignedSeq1
            alignedSeq2 = '-' + alignedSeq2
            i -= 1
        elif score == upScore + gap:
            alignedSeq1 = '-' + alignedSeq1
            alignedSeq2 = seq2[j - 1] + alignedSeq2
            j -= 1
        else:
            raise 'Not Possible'
    while i > 0:
        alignedSeq1 = seq1[i - 1] + alignedSeq1
        alignedSeq2 = '-' + alignedSeq2
        i -= 1
    while j > 0:
        alignedSeq1 = '-' + alignedSeq1
        alignedSeq2 = seq2[j - 1] + alignedSeq2
        j -= 1

    stats = alignment_stats(alignedSeq1, alignedSeq2, similarity_map)
    return SeqAli(sequences      = [alignedSeq1, alignedSeq2],
                  sequence_inits = [1, 1],  identities = stats['idt'],
                  positives = stats['pos'], gaps = stats['gaps'])


def alignment_stats(seq1, seq2, similarity_map):
    stats = {'idt': 0, 'pos': 0, 'gaps': 0}

    for i in range(len(seq1)):
        if seq1[i] == "-" or seq2[i] == "-":       stats['gaps'] += 1
        elif seq1[i] == seq2[i]:                   stats['idt']  += 1
        elif similarity_map[seq1[i]][seq2[i]] > 0: stats['pos']  += 1
    return stats
