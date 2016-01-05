#define -NPY_1_7_API_VERSION

#compile with: gcc -fpic -shared -I /usr/include/python2.7/ -o nw_align.so nw_align.c

__author__ = 'Fernando Puente-Sánchez'
__email__ = 'fpusan@gmail.com'
__version__ = '1.0'
__date__ = '05-Jan-2016'
__license__ = 'BSD-3'
__copyright__ = 'Copyright 2013-2015 Fernando Puente-Sánchez'

BSD3_LICENSE = """

    Copyright (c) 2015, Fernando Puente Sánchez
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this
      list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    * Neither the name of moira, Poisson binomial filtering and its Poisson
      approximation, nor the names of its contributors may be used to endorse or
      promote products derived from this software without specific prior
      written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    
"""

cimport cython
@cython.boundscheck(False)
cpdef nw_align(seq_1, seq_2, int match, int mismatch, int gap, refine_overlap = True, verbose = False):
    """
    Globally align two sequences using the Needleman-Wunsch algorithm.
    Note that there is no gap_extend penalty in NW. Mothur charges the gap_open penalty for every gapped position.
    If refine_overlap=True, mothur's solution for aligning sequences with big non-overlapping regions will be applied.
    - First row and column will be initialized to 0 instead of applying the gap penalty.
    - Last row or column will be modified with by the nw_overlap function to insert the correct number of gaps at 3'. 
    """

    #Populate matrices.
    seq_1 =' %s'%seq_1 #We need an extra space at the beginning for building the matrix.
    seq_2 =' %s'%seq_2
    cdef list score_matrix = []
    cdef list pointer_matrix = [] #Will contain tuples with the displacement in rows and columns (i.e. (0, -1) for going left).
    pointers = {(0, 0) : 'x', (-1, -1) : 'd', (0, -1) : 'l', (-1, 0) : 'u'}

    cdef unsigned int i, j
    cdef unsigned int L1 = len(seq_1)
    cdef unsigned int L2 = len(seq_2)
    #Start with first row and column.
    for i from 0 <= i < L1:
        if refine_overlap:
            score_matrix.append([0])
        else:
            score_matrix.append([<int>i * gap])
        pointer_matrix.append([(-1, 0)])
    for j from 1 <= j < L2:
        if refine_overlap:
            score_matrix[0].append(0)
        else:
            score_matrix[0].append(<int>j * gap)
        pointer_matrix[0].append((0, -1))

    cdef int diag_score, up_score, left_score
    #Then the rest.
    for i from 1 <= i < L1:
        for j from 1 <= j < L2:
            if seq_1[i] == seq_2[j]:
                diag_score = score_matrix[i-1][j-1] + match
            else:
                diag_score = score_matrix[i-1][j-1] + mismatch

            up_score = score_matrix[i-1][j] + gap
            left_score = score_matrix[i][j-1] + gap

            if diag_score >= up_score:
                if diag_score >= left_score:
                    score_matrix[i].append(diag_score)
                    pointer_matrix[i].append((-1, -1))
                else:
                    score_matrix[i].append(left_score)
                    pointer_matrix[i].append((0, -1))
            else:
                if up_score >= left_score:
                    score_matrix[i].append(up_score)
                    pointer_matrix[i].append((-1, 0))
                else:
                    score_matrix[i].append(left_score)
                    pointer_matrix[i].append((0, -1))

    if refine_overlap:
        nw_overlap(score_matrix, pointer_matrix)

    if verbose:
        print 'Needleman-Wunsch matrix.'
        if refine_overlap:
            print '(Mothur overlap corrections were applied)'
        print '\t' + '\t\t'.join(seq_2)
        for i in range(L1):
            print seq_1[i], '\t'.join([str((score_matrix[i][j], pointers[pointer_matrix[i][j]])) for j in range(L2)])
        print

    #Traceback:
    cdef int score = 0
    seq_1_aligned = []
    seq_2_aligned = []
    i, j = len(seq_1) - 1, len(seq_2) - 1

    while i > 0 or j > 0:

        score += score_matrix[i][j]
        pointer = pointer_matrix[i][j]
        if pointer[0] == -1:
            seq_1_aligned.append(seq_1[i])
        else:
            seq_1_aligned.append('-')
        if pointer[1] == -1:
            seq_2_aligned.append(seq_2[j])
        else:
            seq_2_aligned.append('-')

        i += pointer[0] #Note that these values are 0 or negative.
        j += pointer[1]

    seq_1_aligned = ''.join(reversed(seq_1_aligned)) #Since we filled it backwards.
    seq_2_aligned = ''.join(reversed(seq_2_aligned))

    return seq_1_aligned, seq_2_aligned, score

@cython.boundscheck(False)
cdef nw_overlap(list score_matrix, list pointer_matrix):

    """
    Clean up the 3' part of the alignment as mothur does, in order to avoid having a lot of scattered
    bases in the 3' region.
    
    Their solution was to look for the cells with the highest scores in the last row and column, and force gaps
    by making previous cells in those line/column point towards them.

    According to mothur's documentation, issues in the 5' region are solved by themselves during traceback.
    
    Note that they also initialize the first row and column of the NW matrix to 0 instead of applying gap penalties. 
    """
    
    cdef unsigned int i, j
    cdef unsigned int L1 = len(score_matrix)
    cdef unsigned int L2 = len(score_matrix[0])
    cdef int best_in_last_column_score = -10000
    cdef unsigned int best_in_last_column_index = 0
    for i from 0 <= i < L1:
        cell_score = score_matrix[i][-1]
        if cell_score >= best_in_last_column_score:
            best_in_last_column_score = cell_score
            best_in_last_column_index = i

    cdef int best_in_last_row_score = -10000
    cdef unsigned int best_in_last_row_index = 0
    for j from 0 <= j < L2:
        cell_score = score_matrix[-1][j]
        if cell_score >= best_in_last_row_score:
            best_in_last_row_score = cell_score
            best_in_last_row_index = j

    #Decide which sequence will get the gaps:
    if (best_in_last_column_index, best_in_last_row_index) == ((len(score_matrix) - 1), (len(score_matrix[0]) - 1)):
        pass

    elif best_in_last_column_score > best_in_last_row_score:
        for i in range(len(score_matrix) -1, best_in_last_column_index, -1):
            pointer_matrix[i][-1] = -1, 0

    else:
        for j in range(len(score_matrix[0]) - 1, best_in_last_row_index, -1):
            pointer_matrix[-1][j] = 0, -1

