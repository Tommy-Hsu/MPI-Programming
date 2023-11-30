#include "matrix_operations.h"

#define MASTER 0
#define FROM_MASTER 1
#define FROM_WORKER 2

void construct_matrices(int *n_ptr, int *m_ptr, int *l_ptr, int **a_mat_ptr, int **b_mat_ptr) {

    int world_rank, world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int numworkers = world_size - 1;
    int mtype;
    int offset = 0;
    int rows = 0;

    if(world_rank == MASTER){

        /* construct matrix from stdin */
        std::cin >> *n_ptr >> *m_ptr >> *l_ptr;
        *a_mat_ptr = new int[*n_ptr * *m_ptr];
        *b_mat_ptr = new int[*m_ptr * *l_ptr];
        for (int i = 0; i < *n_ptr * *m_ptr; i++) {
            std::cin >> (*a_mat_ptr)[i];
        }
        for (int i = 0; i < *m_ptr * *l_ptr; i++) {
            std::cin >> (*b_mat_ptr)[i];
        }
    }

    MPI_Bcast(n_ptr, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(m_ptr, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(l_ptr, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(world_rank == MASTER){
        /* Send matrix data to the worker tasks */
        int averow = *n_ptr/numworkers;
        int extra = *n_ptr%numworkers;
        mtype = FROM_MASTER;
        MPI_Request MPI_request[4 * numworkers];

        for (int dest=1; dest<=numworkers; dest++)
        {
            rows = (dest <= extra) ? averow+1 : averow;
            MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
            MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
            MPI_Send(&(*a_mat_ptr)[offset], rows*(*m_ptr), MPI_INT, dest, mtype, MPI_COMM_WORLD);
            MPI_Send(&(*b_mat_ptr), (*m_ptr)*(*l_ptr), MPI_INT, dest, mtype, MPI_COMM_WORLD);
            printf("Sent %d a_mat and %d b_mat to task %d \n", rows*(*m_ptr), (*m_ptr)*(*l_ptr), dest);
            offset = offset + rows;
        }
    }
    else{
        mtype = FROM_MASTER;
        MPI_Status MPI_status[4];
        printf("%d", world_rank);
        MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &MPI_status[0]);
        MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &MPI_status[1]);
        *a_mat_ptr = new int[rows  * *m_ptr];
        *b_mat_ptr = new int[*m_ptr * *l_ptr];
        MPI_Recv(&(*a_mat_ptr)[offset], rows * (*m_ptr), MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &MPI_status[2]);
        MPI_Recv(&(*b_mat_ptr), (*m_ptr) * (*l_ptr), MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &MPI_status[3]);
        printf("Received %d a_mat and %d b_mat from task %d\n", rows*(*m_ptr), (*m_ptr)*(*l_ptr), MASTER);
    }
}

void matrix_multiply(const int n, const int m, const int l, const int *a_mat, const int *b_mat) {
    
    int world_rank, world_size;
    MPI_Status MPI_status[4];
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int numworkers = world_size - 1;
    int offset;
    int source;
    int rows;
    int mtype;

    if(world_rank == MASTER){

        /* Receive results from worker tasks */
        int* c_mat = new int[n * l];
        mtype = FROM_WORKER;
        for (int i=1; i<=numworkers; i++)
        {
            source = i;
            // MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
            // MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
            // MPI_Recv(&c[offset], rows*l, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
            // printf("Received results from task %d\n",source);
        }

        /* Print results */
        // for(int i = 0; i < n; i++){
        //     for(int j = 0; j < l; j++){
        //         std::cout << c_mat[i * l + j] << " ";
        //     }
        //     std::cout << std::endl;
        // }

        delete[] c_mat;
    }
    else{

        /* display a_mat and b_mat */
        for(int i = 0; i < n; i++){
            for(int j = 0; j < m; j++){
                std::cout << a_mat[i * m + j] << " ";
            }
            std::cout << std::endl;
        }

        /* Perform matrix multiplication */
        // int *c_mat = new int[rows * l];
        // for(int i = 0; i < rows; i++){
        //     for(int j = 0; j < l; j++){
        //         c_mat[i * l + j] = 0;
        //         for(int k = 0; k < m; k++){
        //             c_mat[i * l + j] += a[i * m + k] * b[k * l + j];
        //         }
        //     }
        // }

        // mtype = FROM_WORKER;
        // MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
        // MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
        // MPI_Send(&c_mat, rows*l, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);

        // delete[] c_mat;
    }
}

void destruct_matrices(int *a_mat, int *b_mat) {

    if(a_mat == NULL || b_mat == NULL){
        return;
    }
    delete[] a_mat;
    delete[] b_mat;
}