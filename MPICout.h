//
//  MPICout.h
//  TestCode
//
//  Created by Shamim Patel on 22/03/2013.
//  Copyright (c) 2013 Shamim Patel. All rights reserved.
//
#ifndef _MPICOUT_H_
#define _MPICOUT_H_


#include <iostream>
#include <mpi.h>
using namespace std;

typedef basic_ostream<char,char_traits<char>> cout_type;
class __MPI_cout
{
    //Processor that will be allowed to output;
    int Processor;    
public:
        
    __MPI_cout()
    {
        this->Processor = 0;
    }
    
    __MPI_cout(unsigned int Processor)
    {
        this->Processor = Processor;
    }
    
    void SetProcessor(unsigned int Processor)
    {
        this->Processor = Processor;
    }
   
    __MPI_cout& operator<<(cout_type& (*__pf)(cout_type&))
    {
        int ProcNum;
        MPI_Comm_rank(MPI_COMM_WORLD,&ProcNum);
        if(ProcNum == this->Processor)
        {
            cout << __pf;
            return *(this);
        }
        else
        {
            return *(this);
        }
    }
            
    __MPI_cout& operator<<(basic_ios<char,char_traits<char>>&
                          (*__pf)(basic_ios<char,char_traits<char>>&))
    {
        int ProcNum;
        MPI_Comm_rank(MPI_COMM_WORLD,&ProcNum);
        if(ProcNum == this->Processor)
        {
            cout << __pf;
            return *(this);
        }
        else
        {
            return *(this);
        }
    }
            
    __MPI_cout& operator<<(ios_base& (*__pf)(ios_base&))
    {
        int ProcNum;
        MPI_Comm_rank(MPI_COMM_WORLD,&ProcNum);
        if(ProcNum == this->Processor)
        {
            cout << __pf;
            return *(this);
        }
        else
        {
            return *(this);
        }
    }
    __MPI_cout& operator<<(bool __n)
    {
        int ProcNum;
        MPI_Comm_rank(MPI_COMM_WORLD,&ProcNum);
        if(ProcNum == this->Processor)
        {
            cout << __n;
            return *(this);
        }
        else
        {
            return *(this);
        }
    }
    __MPI_cout& operator<<(short __n)
    {
        int ProcNum;
        MPI_Comm_rank(MPI_COMM_WORLD,&ProcNum);
        if(ProcNum == this->Processor)
        {
            cout << __n;
            return *(this);
        }
        else
        {
            return *(this);
        }
    }
    __MPI_cout& operator<<(unsigned short __n)
    {
        int ProcNum;
        MPI_Comm_rank(MPI_COMM_WORLD,&ProcNum);
        if(ProcNum == this->Processor)
        {
            cout << __n;
            return *(this);
        }
        else
        {
            return *(this);
        }
    }
    __MPI_cout& operator<<(int __n)
    {
        int ProcNum;
        MPI_Comm_rank(MPI_COMM_WORLD,&ProcNum);
        if(ProcNum == this->Processor)
        {
            cout << __n;
            return *(this);
        }
        else
        {
            return *(this);
        }
    }
    __MPI_cout& operator<<(unsigned int __n)
    {
        int ProcNum;
        MPI_Comm_rank(MPI_COMM_WORLD,&ProcNum);
        if(ProcNum == this->Processor)
        {
            cout << __n;
            return *(this);
        }
        else
        {
            return *(this);
        }
    }
    __MPI_cout& operator<<(long __n)
    {
        int ProcNum;
        MPI_Comm_rank(MPI_COMM_WORLD,&ProcNum);
        if(ProcNum == this->Processor)
        {
            cout << __n;
            return *(this);
        }
        else
        {
            return *(this);
        }
    }
    __MPI_cout& operator<<(unsigned long __n)
    {
        int ProcNum;
        MPI_Comm_rank(MPI_COMM_WORLD,&ProcNum);
        if(ProcNum == this->Processor)
        {
            cout << __n;
            return *(this);
        }
        else
        {
            return *(this);
        }
    }
    __MPI_cout& operator<<(long long __n)
    {
        int ProcNum;
        MPI_Comm_rank(MPI_COMM_WORLD,&ProcNum);
        if(ProcNum == this->Processor)
        {
            cout << __n;
            return *(this);
        }
        else
        {
            return *(this);
        }
    }
    __MPI_cout& operator<<(unsigned long long __n)
    {
        int ProcNum;
        MPI_Comm_rank(MPI_COMM_WORLD,&ProcNum);
        if(ProcNum == this->Processor)
        {
            cout << __n;
            return *(this);
        }
        else
        {
            return *(this);
        }
    }
    __MPI_cout& operator<<(float __f)
    {
        int ProcNum;
        MPI_Comm_rank(MPI_COMM_WORLD,&ProcNum);
        if(ProcNum == this->Processor)
        {
            cout << __f;
            return *(this);
        }
        else
        {
            return *(this);
        }
    }
    __MPI_cout& operator<<(double __f)
    {
        int ProcNum;
        MPI_Comm_rank(MPI_COMM_WORLD,&ProcNum);
        if(ProcNum == this->Processor)
        {
            cout << __f;
            return *(this);
        }
        else
        {
            return *(this);
        }
    }
    __MPI_cout& operator<<(long double __f)
    {
        int ProcNum;
        MPI_Comm_rank(MPI_COMM_WORLD,&ProcNum);
        if(ProcNum == this->Processor)
        {
            cout << __f;
            return *(this);
        }
        else
        {
            return *(this);
        }
    }
    __MPI_cout& operator<<(const void* __p)
    {
        int ProcNum;
        MPI_Comm_rank(MPI_COMM_WORLD,&ProcNum);
        if(ProcNum == this->Processor)
        {
            cout << __p;
            return *(this);
        }
        else
        {
            return *(this);
        }
    }
    __MPI_cout& operator<<(basic_streambuf<char,char_traits<char>>* __sb)
    {        
        int ProcNum;
        MPI_Comm_rank(MPI_COMM_WORLD,&ProcNum);
        if(ProcNum == this->Processor)
        {
            cout << __sb;
            return *(this);
        }
        else
        {
            return *(this);
        }
    }            
    __MPI_cout& operator<<(const char* String)
    {
        int ProcNum;
        MPI_Comm_rank(MPI_COMM_WORLD,&ProcNum);
        if(ProcNum == this->Processor)
        {
            cout << String;
            return *(this);
        }
        else
        {
            return *(this);
        }
    }
};
            
            
__MPI_cout MPI_cout;
#define _USING_MPI_


#endif
