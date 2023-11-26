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

    // TODO: binary tree redunction
    unsigned int seed = (unsigned int)time(NULL) + (unsigned long)pthread_self();
    long long int nHits = 0;
    for( long long int i = 0; i < tosses_per_rank; i++ ) {
        double x = static_cast<double>(rand_r(&seed)) / static_cast<double>(RAND_MAX);
        double y = static_cast<double>(rand_r(&seed)) / static_cast<double>(RAND_MAX);
        if ( x*x + y*y <= 1.f ) {
            nHits++;
        }
    }
    printf("rank %d, nhits %lld\n", world_rank, nHits);

    for(size_t step = 1; step < world_size; step *= 2)
    {
        if (world_rank % (2*step) == 0)
        {
            // Receive
            long long int recv_nHits;
            MPI_Recv(&recv_nHits, 1, MPI_LONG_LONG_INT, world_rank+step, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            nHits += recv_nHits;
        }
        else
        {
            // Send
            MPI_Send(&nHits, 1, MPI_LONG_LONG_INT, world_rank-step, 0, MPI_COMM_WORLD);
            break;
        }
    }

    if (world_rank == 0)
    {  
        unsigned int seed = (unsigned int)time(NULL) + (unsigned long)pthread_self();
        for( long long int i = 0; i < tosses_per_rank_rest; i++ ) {
            double x = static_cast<double>(rand_r(&seed)) / static_cast<double>(RAND_MAX);
            double y = static_cast<double>(rand_r(&seed)) / static_cast<double>(RAND_MAX);
            if ( x*x + y*y <= 1.f ) {
                nHits++;
            }
        }

        // TODO: PI result
        pi_result = 4.f * (double)nHits / (double)tosses;

        // --- DON'T TOUCH ---
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        // ---
    }

    MPI_Finalize();
    return 0;
}
