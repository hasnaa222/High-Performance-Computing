#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "omp.h"


int main(int argc, char * argv[])
{
    FILE *fp;
    fp = fopen("CS371_Assign3_20190380_20190167_20190180.txt", "r");
    int flag=0, min=0, max=0,remSize=0;
    int pointsNum, barsNum,range,portionSize;
    int i,j=0,my_rank,process;
    char buffer[10];
    int * points;
    MPI_Status status;
    MPI_Init( &argc, &argv );
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &process);
    int Thread = omp_get_max_threads();
    if( my_rank == 0)
    {
        printf("Please enter the number of points : ");
        scanf("%d", &pointsNum);
        printf("Please enter the number of bars : ");
        scanf("%d", &barsNum);
        points = (int*) malloc(pointsNum * sizeof(int));

        for (i = 0; i < pointsNum; i++)
        {
            fscanf(fp,"%d\n", &points[i]);
            if(min>points[i])
                {
                    min = points[i];
                }
                if(max<points[i])
                {
                    max = points[i];
                }
        }
        range = ( max - min )/ barsNum;
        portionSize = pointsNum/(process);
        remSize = pointsNum - (process)*portionSize;
    }
    MPI_Bcast(&range, 1, MPI_INT, 0,MPI_COMM_WORLD);
    MPI_Bcast(&barsNum, 1, MPI_INT, 0,MPI_COMM_WORLD);
    MPI_Bcast(&portionSize, 1, MPI_INT, 0,MPI_COMM_WORLD);
    MPI_Bcast(&min, 1, MPI_INT, 0,MPI_COMM_WORLD);
    MPI_Bcast(&pointsNum, 1, MPI_INT, 0,MPI_COMM_WORLD);
    MPI_Bcast(&remSize, 1, MPI_INT, 0,MPI_COMM_WORLD);
    int* bars = (int*) malloc((process*barsNum) * sizeof(int));
    int* localBars = (int*) malloc(barsNum * sizeof(int));
    int* globalSum = (int*) malloc(barsNum * sizeof(int));
    for(i=0; i<barsNum; i++)
    {
        localBars[i]=0;
        globalSum[i]=0;
    }
    for (j=0; j<(barsNum*process); j++)
    {
        bars[j]=0;
    }


    int** localBarSum = (int **) malloc(Thread * sizeof(int *));
    for (i = 0; i < Thread; i++)
    {
        localBarSum[i] = (int *) malloc(barsNum * sizeof(int));
        for(j =0; j< barsNum; j++)
        {
            localBarSum [i][j]=0;
        }
    }
    if(my_rank==0 && remSize!=0)
    {
        #pragma omp parallel shared(pointsNum,remSize, barsNum, min, range,localBarSum,points) private(i, j)
        {

            #pragma omp for schedule(static)
            for(i=pointsNum-remSize ; i<pointsNum; i++)
            {
                #pragma omp parallel
                #pragma omp for schedule(static)
                for (j=1; j<=barsNum; j++)
                {
                    if(points[i] <= min +(range*j) && points[i] > min +(range*(j-1)))
                    {
                        localBarSum[omp_get_thread_num()][j-1]=localBarSum[omp_get_thread_num()][j-1]+1;
                    }
                    else if( points[i] == min)
                    {
                        localBarSum[omp_get_thread_num()][0]=1+localBarSum[omp_get_thread_num()][0];
                    }
                }
            }
        }
    }
    int* portionPoints = (int*) malloc(pointsNum * sizeof(int));
    MPI_Scatter(points, portionSize, MPI_INT, portionPoints, portionSize, MPI_INT, 0, MPI_COMM_WORLD);


    #pragma omp parallel shared(portionSize, barsNum, min, range,localBarSum) private(i, j)
    {

        #pragma omp for schedule(static)
        for(i = 0; i<portionSize; i++)
        {
            #pragma omp parallel
            #pragma omp for schedule(static)
            for (j=1; j<=barsNum; j++)
            {
                if(portionPoints[i] <= min +(range*j) && portionPoints[i] > min +(range*(j-1)))
                {
                    localBarSum[omp_get_thread_num()][j-1]=localBarSum[omp_get_thread_num()][j-1]+1;
                }
                else if( portionPoints[i] == min)
                {
                    localBarSum[omp_get_thread_num()][0]=1+localBarSum[omp_get_thread_num()][0];
                }
            }
        }
    }

    for(j=0; j< Thread; j++)
    {
        for(i =0; i< barsNum; i++)
        {
            localBars[i]+=localBarSum [j][i];
        }
    }
    MPI_Gather (localBars, barsNum, MPI_INT, bars, barsNum, MPI_INT, 0, MPI_COMM_WORLD );

    if(my_rank==0)
    {
        #pragma omp parallel shared(process, barsNum,globalSum ,bars,flag) private(j)
        {
            #pragma omp for schedule(dynamic)
            for (j=0; j<(barsNum*process); j++)
            {
                globalSum[flag] += bars[j];
                flag++;
                if(flag == barsNum)
                {
                    flag =0;
                }
            }
        }
    }
    if ( my_rank == 0)
    {
        for (j=0; j<barsNum; j++)
        {
            printf("The range start with %d, end with %d, with count %d.\n", min +(range*j),min +(range*(j+1)),globalSum[j]);
        }
    }

    fclose(fp);
    MPI_Finalize();
    free(globalSum);
    free(bars);
    free(localBars);
    free(portionPoints);
    for (i = 0; i < Thread; i++)
    {
        free(localBarSum[i]);
    }
    free(localBarSum);
    return 0;
}
