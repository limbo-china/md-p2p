#include "timer.h"
#include "mympi.h"
#include "parameter.h"
#include "info.h"
#include "atom.h"
#include "potential.h"
#include "system.h"

#include <stdio.h>
#include <unistd.h>
#include <mpi.h>

void updateMomenta(System* sys, Parameter* para); 
void updatePosition(System* sys, Parameter* para);

int main(int argc, char** argv){
	
	MPI_Init(&argc, &argv);
	initRank();

	// char processor_name[20];
	// int name_len;
	// MPI_Get_processor_name(processor_name, &name_len);

	// fprintf(stdout, "%s\n",processor_name);
	//fprintf(stdout, "rankNums: %d\n", getRankNums());
	//fprintf(stdout, "myRank: %d\n", getMyRank());


	beginTimer(total);

	Parameter* para = readParameter();
	printPara(stdout,para);

	
	//sleep(5);
	System* sys = initSystem(para);

	beginTimer(loop);
	for(int i=1;i<=para->stepNums;i++){
    	updateMomenta(sys, para); 

    	updatePosition(sys, para);

    	//beginTimer(adjustatom);
    	adjustAtoms(sys);
    	//endTimer(adjustatom);

    	beginTimer(force);
    	computeForce(sys);
    	endTimer(force);

    	updateMomenta(sys, para); 
    	if(i%para->printNums == 0){

    	//MPI_Allreduce(&sys->atoms->myNum, &sys->atoms->totalNum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    	//printTotalAtom(stdout,sys->atoms);
    		if(ifZeroRank())
			{
				printf("当前步数: %d 	已模拟时间: %g fs\t",i,i*para->stepTime);
			}		
    		computeTotalKinetic(sys);
    		printTemper(stdout,sys->energy,sys->atoms->totalNum);
    	}
    }
    printTemper2(stdout,sys->energy,sys->atoms->totalNum);
	endTimer(loop);

	//beginTimer(test);
    //MPI_Win_free(win);
    //endTimer(test);
	

	endTimer(total);

	double looptime = getGlobalTime(loop);
	double forcetime =  getGlobalTime(force);
	double commtime = getGlobalTime(communication);
	double globalloop=0.0,globalforce=0.0,globalcomm=0.0;

	MPI_Allreduce(&looptime, &globalloop, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&forcetime, &globalforce, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&commtime, &globalcomm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	if(ifZeroRank())
	{
		//fprintf(stdout, "total time: %g\n",getGlobalTime(total));
		fprintf(stdout, "\n------\n总执行时间: %g s\n"/*,globalloop*/,globalloop/getRankNums());
		//fprintf(stdout, "adjust time: %g\n",getGlobalTime(adjustatom));
		fprintf(stdout, "通信消耗时间: %g s\n"/*,globalcomm*/,globalcomm/getRankNums());
		fprintf(stdout, "计算力时间: %g s\n------\n"/*,globalforce*/,globalforce/getRankNums());
		//fprintf(stdout, "test time: %g\n",getGlobalTime(test));

	}


	MPI_Finalize();
	return 0;
}

void updateMomenta(System* sys, Parameter* para){

	double t = 0.5*para->stepTime;

	for (int nCell=0; nCell<sys->cells->myCellNum; nCell++)
      	for (int n=MAXPERCELL*nCell,count=0; count<sys->cells->atomNum[nCell]; count++,n++)
      		for(int i=0;i<3;i++)
         		sys->atoms->momenta[n][i] += t*sys->atoms->force[n][i];
}
void updatePosition(System* sys, Parameter* para){

	double t = para->stepTime;
	double m = sys->lattice->atomM;

	for (int nCell=0; nCell<sys->cells->myCellNum; nCell++)
      	for (int n=MAXPERCELL*nCell,count=0; count<sys->cells->atomNum[nCell]; count++,n++)
      		for(int i=0;i<3;i++)
         		sys->atoms->pos[n][i] += t*sys->atoms->momenta[n][i]/m;
}