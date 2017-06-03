#include "energy.h"
#include "system.h"

#include <mpi.h>
#include <stdlib.h>

// 初始化结构体
void initEnergy(Energy** ener){

	*ener = (Energy*)malloc(sizeof(Energy));
    Energy* energy = *ener;

    energy->kineticEnergy = 0.0;
    energy->potentialEnergy = 0.0;
}

// 计算体系的总动能
void computeTotalKinetic(struct SystemStr* sys){

printf("1\n");
  // double ener[2];
  // double res_ener[2];
  // ener[0] = 0.0;
  // ener[1] = sys->energy->potentialEnergy;
	 double myKineticEnergy = 0.0;
	 double globalKineticEnergy = 0.0;
	double atomM = sys->lattice->atomM;

	// 计算本空间的原子总动能
   	for (int nCell=0; nCell<sys->cells->myCellNum; nCell++)
      	for (int n=MAXPERCELL*nCell,count=0; count<sys->cells->atomNum[nCell]; count++,n++)
      		for(int i=0; i<3; i++)
         		myKineticEnergy += sys->atoms->momenta[n][i]*sys->atoms->momenta[n][i]
         			*0.5/atomM;
printf("2\n");
    // AllReduce, 得到整个体系的总动能和势能
    MPI_Allreduce(&myKineticEnergy, &globalKineticEnergy, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    printf("3\n");

   	sys->energy->kineticEnergy = globalKineticEnergy;
    //sys->energy->potentialEnergy = ener[1];
}