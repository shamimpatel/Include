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
#include "MPICout.h"
using namespace std;

std::string CreateProcessorUniqueFilename( const char* Filename, const char* Extension)
{
    std::stringstream S;
    int ProcNum;
    MPI_Comm_rank(MPI_COMM_WORLD,&ProcNum);
    S << Filename << "_P" << ProcNum << Extension;
    return S.str();
}


std::string CreateConcatCommand( const char* Filename, const char* Extension)
{
    int NumProcessors;
    MPI_Comm_size(MPI_COMM_WORLD,&NumProcessors);
    
    stringstream S;
    
    S << "cat ";
    for(int i=0;i<NumProcessors;i++)
    {
        S << Filename << "_P" << i << Extension << " ";
    }
    
    S << "> " << Filename << Extension;
    
    
    return S.str();
    
}

#endif
