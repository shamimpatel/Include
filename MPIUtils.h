//
//  MPIUtils.h
//  DiffractionMain
//
//  Created by Shamim Patel on 11/03/2013.
//  Copyright (c) 2013 Shamim Patel. All rights reserved.
//

#ifndef _MPIUTILS_H_
#define _MPIUTILS_H_
#include <string.h>
#include <sstream>
#include "mpi.h"
using namespace std;

std::string CreateProcessorUniqueFilename( std::string Filename, std::string Extension)
{
    std::stringstream S;
    int ProcNum;
    MPI_Comm_rank(MPI_COMM_WORLD,&ProcNum);
    S << Filename << "_P" << ProcNum << Extension;
    return S.str();
}


#endif
