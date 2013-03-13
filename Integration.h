//
//  Integration.h
//  TestCode
//
//  Created by Shamim Patel on 13/03/2013.
//  Copyright (c) 2013 Shamim Patel. All rights reserved.
//

#ifndef _INTEGRATION_H_
#define _INTEGRATION_H_
#include <cmath>


//use trapezium rule to estimate the value of a 1D integral
double Compute1DIntegral( double (*F)(double), double lowerLimit, double upperLimit, double dx )
{
    double Sum = 0.0;
    int nSteps = (upperLimit - lowerLimit)/dx;
    
    double x0 = F(lowerLimit);
    double xN = F(upperLimit);
    
    double ExtraTerms = (x0 + xN)*dx*0.5;
    
    for( int Step = 1; Step <= nSteps; Step++)
    {
        double x = lowerLimit + dx*Step;
        Sum += F(x)*dx;
    }
    return Sum + ExtraTerms;
}




#endif
