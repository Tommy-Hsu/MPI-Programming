#pragma once

#include <cstdlib>
#include <ctime>
#include <mpi.h>
#include <fstream>
#include <iostream>

void construct_matrices(int *n_ptr, int *m_ptr, int *l_ptr, int **a_mat_ptr, int **b_mat_ptr);
void matrix_multiply(const int n, const int m, const int l, const int *a_mat, const int *b_mat);
void destruct_matrices(int *a_mat, int *b_mat);
