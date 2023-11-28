// #include <mpi.h>
// #include <stdio.h>
// #include <stdlib.h>
// #include <time.h>
// #include <sys/types.h>
// #include <unistd.h>

// int main(int argc, char **argv)
// {
//     // --- DON'T TOUCH ---
//     MPI_Init(&argc, &argv);
//     double start_time = MPI_Wtime();
//     double pi_result;
//     long long int tosses = atoi(argv[1]);
//     int world_rank, world_size;
//     // ---

//     MPI_Win win;
//     long long int tosses_per_rank = tosses / (world_size);
//     long long int tosses_per_rank_rest = tosses % (world_size);
//     // long long int hits_per_rank_array[world_size] = {0};
//     long long int* nHits;

//     // TODO: MPI init
//     MPI_Comm_size(MPI_COMM_WORLD, &world_size);
//     MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
//     char processor_name[MPI_MAX_PROCESSOR_NAME];
//     int name_len;
//     MPI_Get_processor_name(processor_name, &name_len);


//     if (world_rank == 0)
//     {
//         // Master
//         MPI_Alloc_mem(sizeof(long long int), MPI_INFO_NULL, &nHits);
//         *nHits = 0;

//         MPI_Win_create(nHits, sizeof(long long int), sizeof(long long int), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
//         unsigned int seed = (unsigned int)time(NULL) + (unsigned long)pthread_self();
//         for( long long int i = 0; i < tosses_per_rank+tosses_per_rank_rest; i++ ) {
//             double x = (double)rand_r(&seed) / (double)RAND_MAX;
//             double y = (double)rand_r(&seed) / (double)RAND_MAX;
//             if ( x*x + y*y <= 1.f ) {
//                 *nHits++;
//             }
//         }
//         printf("rank %d, nHits %lld\n", world_rank, *nHits);
//     }
//     else
//     {
//         // Workers
//         MPI_Win_create(NULL, 0, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &win);
//         unsigned int seed = (unsigned int)time(NULL) + (unsigned long)pthread_self();
//         *nHits = 0;
//         for( long long int i = 0; i < tosses_per_rank; i++ ) {
//             double x = (double)rand_r(&seed) / (double)RAND_MAX;
//             double y = (double)rand_r(&seed) / (double)RAND_MAX;
//             if ( x*x + y*y <= 1.f ) {
//                 (*nHits)++;
//             }
//         }
//         printf("rank %d, nhits %lld\n", world_rank, *nHits);
//     }
//     MPI_Win_free(&win);

//     if (world_rank == 0)
//     {
//         // TODO: handle PI result
//         pi_result = 4.0 * static_cast<double>(*nHits) /static_cast<double>(tosses);

//         // --- DON'T TOUCH ---
//         double end_time = MPI_Wtime();
//         printf("%lf\n", pi_result);
//         printf("MPI running time: %lf Seconds\n", end_time - start_time);
//         // ---
//     }
    
//     MPI_Finalize();
//     return 0;
// }

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>

int main(int argc, char **argv)
{
    // --- DON'T TOUCH ---
    MPI_Init(&argc, &argv);
    double start_time = MPI_Wtime();
    double pi_result;
    long long int tosses = atoi(argv[1]);
    int world_rank, world_size;
    // ---

    // TODO: MPI init
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    MPI_Win win;
    long long int tosses_per_rank = tosses / world_size;
    long long int tosses_per_rank_rest = tosses % world_size;
    long long int nHits_local = 0, *nHits_global;

    unsigned int seed = (unsigned int)time(NULL) + (unsigned long)pthread_self();
    for (long long int i = 0; i < tosses_per_rank + (world_rank == 0 ? tosses_per_rank_rest : 0); i++) {
        double x = (double)rand_r(&seed) / (double)RAND_MAX;
        double y = (double)rand_r(&seed) / (double)RAND_MAX;
        if (x * x + y * y <= 1.0) {
            nHits_local++;
        }
    }
    printf("rank %d, nHits %lld\n", world_rank, nHits_local);

    if (world_rank == 0)
    {
        // Master
        MPI_Alloc_mem(world_size * sizeof(long long int), MPI_INFO_NULL, &nHits_global);
        *nHits_global = {0};
        MPI_Win_create(nHits_global, world_size * sizeof(long long int), sizeof(long long int), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);
        MPI_Put(&nHits_local, 1, MPI_LONG_LONG_INT, 0, world_rank, 1, MPI_LONG_LONG_INT, win);
        MPI_Win_unlock(0, win);

        int ready = 0;
        while (!ready)
        {
            MPI_Win_lock(MPI_LOCK_SHARED, 0, 0, win);
            ready = 1;
            for (int i = 1; i < world_size; i++)
            {
                if (nHits_global[i] == 0)
                {
                    ready = 0;
                    break;
                }
            }
            MPI_Win_unlock(0, win);
        }
    }
    else
    {
        // Workers
        MPI_Win_create(NULL, 0, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &win);
        MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);
        MPI_Put(&nHits_local, 1, MPI_LONG_LONG_INT, 0, world_rank, 1, MPI_LONG_LONG_INT, win);
        MPI_Win_unlock(0, win);
    }

    if (world_rank == 0)
    {
        // TODO: handle PI result
        pi_result = 0.0;
        double temp = 4.0 / static_cast<double>(tosses);
        for(size_t i = 0; i < world_size; i++)
        {
            pi_result += static_cast<double>(nHits_global[i]);
        }
        pi_result *= temp;

        // --- DON'T TOUCH ---
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        MPI_Free_mem(nHits_global);
    }

    MPI_Win_free(&win);
    MPI_Finalize();
    return 0;
}
