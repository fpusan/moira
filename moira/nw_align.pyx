#define -NPY_1_7_API_VERSION
#cython: boundscheck=False, wraparound=False

#compile with: cython nw_align.pyx && gcc -fpic -shared -I /usr/include/python2.7/ -o nw_align.so nw_align.c

__author__ = 'Fernando Puente-Sánchez'
__email__ = 'fpusan@gmail.com'
__version__ = '1.1'
__date__ = '01-Apr-2016'
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
from cython.view cimport array as cvarray

cpdef nw_align(seq_1, seq_2, int match, int mismatch, int gap, refine_overlap = True, verbose = False):
    """
    Globally align two sequences using the Needleman-Wunsch algorithm.
    Note that there is no gap_extend penalty in NW. Mothur charges the gap_open penalty for every gapped position.
    If refine_overlap=True, mothur's solution for aligning sequences with big non-overlapping regions will be applied.
    - First row and column will be initialized to 0 instead of applying the gap penalty.
    - Last row or column will be modified with by the nw_overlap function to insert the correct number of gaps at 3'. 
    """
    cdef Py_ssize_t i, j, rd, cd
    #Populate matrices.
    seq_1 =' %s'%seq_1 #We need an extra space at the beginning for building the matrix.
    seq_2 =' %s'%seq_2
    cdef Py_ssize_t L1 = len(seq_1)
    cdef Py_ssize_t L2 = len(seq_2)
    cdef int [:,:] score_matrix = cvarray(shape=(L1,L2), itemsize=sizeof(int), format="i")
    cdef int [:,:] row_displacement = cvarray(shape=(L1,L2), itemsize=sizeof(int), format="i")
    cdef int [:,:] column_displacement = cvarray(shape=(L1,L2), itemsize=sizeof(int), format="i")
    pointers = {(0, 0) : 'x', (-1, -1) : 'd', (0, -1) : 'l', (-1, 0) : 'u'}

    #Start with first row and column.
    for i from 0 <= i < L1:
        if refine_overlap:
            score_matrix[i, 0] = 0
        else:
            score_matrix[i, 0] = i * gap
        row_displacement[i, 0] = -1
        column_displacement[i, 0] = 0
    for j from 1 <= j < L2:
        if refine_overlap:
            score_matrix[0, j] = 0
        else:
            score_matrix[0, j] = j * gap
        row_displacement[0, j] = 0
        column_displacement[0, j] = -1

    cdef int diag_score, up_score, left_score
    #Then the rest.
    for i from 1 <= i < L1:
        for j from 1 <= j < L2:
            if seq_1[i] == seq_2[j]:
                diag_score = score_matrix[i-1, j-1] + match
            else:
                diag_score = score_matrix[i-1, j-1] + mismatch

            up_score = score_matrix[i-1, j] + gap
            left_score = score_matrix[i, j-1] + gap

            if diag_score >= up_score:
                if diag_score >= left_score:
                    score_matrix[i, j] = diag_score
                    row_displacement[i, j] = -1
                    column_displacement[i, j] = -1
                else:
                    score_matrix[i, j] = left_score
                    row_displacement[i, j] = 0
                    column_displacement[i, j] = -1
            else:
                if up_score >= left_score:
                    score_matrix[i,j] = up_score
                    row_displacement[i, j] = -1
                    column_displacement[i, j] = 0
                else:
                    score_matrix[i,j] = left_score
                    row_displacement[i, j] = 0
                    column_displacement[i, j] = -1

    if refine_overlap:
        nw_overlap(score_matrix, row_displacement, column_displacement)

    if verbose:
        print 'Needleman-Wunsch matrix.'
        if refine_overlap:
            print '(Mothur overlap corrections were applied)'
        print '\t' + '\t\t'.join(seq_2)
        for i in range(L1):
            print seq_1[i], '\t'.join([str((score_matrix[i,j], pointers[(row_displacement[i,j],column_displacement[i,j])])) for j in range(L2)])
        print

    #Traceback:
    cdef int score = 0
    seq_1_aligned = []
    seq_2_aligned = []
    i, j = L1 - 1, L2 - 1

    while i > 0 or j > 0:

        rd, cd = row_displacement[i, j], column_displacement[i, j]
        score += score_matrix[i, j]
        if rd == -1:
            seq_1_aligned.append(seq_1[i])
        else:
            seq_1_aligned.append('-')
        if cd == -1:
            seq_2_aligned.append(seq_2[j])
        else:
            seq_2_aligned.append('-')

        i += rd #Note that these values are 0 or negative.
        j += cd

    seq_1_aligned = ''.join(reversed(seq_1_aligned)) #Since we filled it backwards.
    seq_2_aligned = ''.join(reversed(seq_2_aligned))

    return seq_1_aligned, seq_2_aligned, score


cdef nw_overlap(int [:,:] score_matrix, int [:,:] row_displacement, int [:,:] column_displacement):

    """
    Clean up the 3' part of the alignment as mothur does, in order to avoid having a lot of scattered
    bases in the 3' region.
    
    Their solution was to look for the cells with the highest scores in the last row and column, and force gaps
    by making previous cells in those line/column point towards them.

    According to mothur's documentation, issues in the 5' region are solved by themselves during traceback.
    
    Note that they also initialize the first row and column of the NW matrix to 0 instead of applying gap penalties. 
    """
    
    cdef unsigned int i, j
    cdef Py_ssize_t L1 = score_matrix.shape[0]
    cdef Py_ssize_t L2 = score_matrix.shape[1]
    cdef int best_in_last_column_score = -10000
    cdef Py_ssize_t best_in_last_column_index = 0
    cdef int cell_score
    for i from 0 <= i < L1:
        cell_score = score_matrix[i, L2-1]
        if cell_score >= best_in_last_column_score:
            best_in_last_column_score = cell_score
            best_in_last_column_index = i

    cdef int best_in_last_row_score = -10000
    cdef Py_ssize_t best_in_last_row_index = 0
    for j from 0 <= j < L2:
        cell_score = score_matrix[L1-1,j]
        if cell_score >= best_in_last_row_score:
            best_in_last_row_score = cell_score
            best_in_last_row_index = j

    #Decide which sequence will get the gaps:
    if (best_in_last_column_index, best_in_last_row_index) == ((len(score_matrix) - 1), (len(score_matrix[0]) - 1)):
        pass

    elif best_in_last_column_score > best_in_last_row_score:
        for i in range(len(score_matrix) -1, best_in_last_column_index, -1):
            row_displacement[i, L2-1] = -1
            column_displacement[i, L2-1] = 0

    else:
        for j in range(len(score_matrix[0]) - 1, best_in_last_row_index, -1):
            row_displacement[L1-1, j] = 0
            column_displacement[L1-1, j] = -1

