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
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // every rank has
    long long int tosses_per_rank = tosses / (world_size);
    long long int tosses_per_rank_rest = tosses % (world_size);
    long long int* hits_per_rank_array = (long long int*)malloc(sizeof(long long int) * world_size);
    long long int nHits = 0;

    // TODO: use MPI_Gather
    unsigned int seed = (unsigned int)time(NULL) + (unsigned long)pthread_self();
    for( long long int i = 0; i < tosses_per_rank; i++ ) {
        double x = (double)rand_r(&seed) / (double)RAND_MAX;
        double y = (double)rand_r(&seed) / (double)RAND_MAX;
        if ( x*x + y*y <= 1.f ) {
            nHits++;
        }
    }
    MPI_Gather(&nHits, 1, MPI_LONG_LONG_INT, hits_per_rank_array, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

    if (world_rank == 0)
    {
        // TODO: PI result
        for (int i = 1; i < world_size; i++) {
            nHits += hits_per_rank_array[i];
        }
        pi_result = 4.f * static_cast<double>(nHits) / static_cast<double>(tosses);

        // --- DON'T TOUCH ---
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        // ---
    }
    
    MPI_Finalize();
    return 0;
}
