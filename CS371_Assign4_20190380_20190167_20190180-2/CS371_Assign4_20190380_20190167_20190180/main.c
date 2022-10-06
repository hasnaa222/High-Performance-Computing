#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "math.h"
int main(int argc, char * argv[])
{
    int my_rank;		/* rank of process	*/
    int process;			/* number of process	*/
    int portionSize,remSize=0;
    int i,Number, s=3;
    double globalSum = 0.0, localSum = 0.0;
    long double expected = 1.202056903159594;
    double time1,time2, duration =0.0, TP;

    MPI_Status status;
    MPI_Init( &argc, &argv );
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &process);
    MPI_Barrier(MPI_COMM_WORLD);

    if( my_rank == 0)
    {
        printf("please enter the number of N : ");
        scanf("%d",&Number);
        time1 = MPI_Wtime();
        portionSize = Number/(process);
        remSize = Number - (process)*portionSize;
        if(remSize >0)
        {
//            time1 = MPI_Wtime();
            for(i =Number-remSize+1 ; i<=Number; i++ )
            {
                localSum += 1.0 / pow(i,s);
            }
//            time2 = MPI_Wtime();
//            duration += time2 - time1;
        }

    }
    MPI_Bcast(&portionSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Number, 1, MPI_INT, 0, MPI_COMM_WORLD);

//    time1 = MPI_Wtime();
    for(i=1+(portionSize*my_rank); i<=portionSize+(portionSize*my_rank) ; i++)
    {
        localSum += 1.0 / pow(i, s);
    }
    MPI_Reduce(&localSum,&globalSum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
//    time2 = MPI_Wtime();
//    duration += time2 - time1;
//    MPI_Reduce(&duration,&TP,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

    if (my_rank == 0 )
    {
        printf("Total Number reslut %lf\n",globalSum );
        printf("Total error = expected - calculated : %Lf\n",expected - globalSum );
//        printf("Calculate the runtime : %lf\n", TP );
        time2 = MPI_Wtime();
        duration += time2 - time1;
        printf("Calculate the runtime : %lf\n", duration );

    }

    MPI_Finalize();
    return 0;
}
