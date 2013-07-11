//
//  Pinhole.h
//  CCDPinholeProcess
//
//  Created by Shamim Patel on 05/03/2013.
//  Copyright (c) 2013 Shamim Patel. All rights reserved.
//

#ifndef _PINHOLE_H_
#define _PINHOLE_H_

#include "Vector.h"

class PinholePlane
{
    Vector PinholeOrigin;
    Vector PinholeNormal;
    double DistToPlane;
    Vector InputPinholeNormal;
    double PinholeRadius;
    
public:
    PinholePlane( Vector InputOrigin, Vector InputNormal,
                  double InputPinholeRadius  )
    {
        this->InputPinholeNormal = InputNormal;
        this->PinholeOrigin = InputOrigin;
        this->PinholeRadius = InputPinholeRadius;
        
        this->PinholeNormal = InputNormal;
        
        if( fabs(PinholeNormal.Magnitude() - 1.0) > 1E-8 )
        {
            cout << "Pinhole Input Normal not normalized (Diff: " << PinholeNormal.Magnitude() - 1.0 <<  ")" << endl;
            PinholeNormal = PinholeNormal.Normalized();
        }
        
        
        //This isn't that important for now. But later we might care about
        //the proper orientation of the pinhole's x,y axes
        if(PinholeOrigin.Dot(PinholeNormal) < 0.0)
        {
            cout << "Input Normal points towards Origin" << endl;
            //need to check like this or 0.0 --> -0.0
            //-0.0 with the trig functions goes to the other end of the output range.
            if(PinholeNormal.x != 0.0)
            {
                PinholeNormal.x *= -1.0;
            }
            if(PinholeNormal.y != 0.0)
            {
                PinholeNormal.y *= -1.0;
            }
            if(PinholeNormal.z != 0.0)
            {
                PinholeNormal.z *= -1.0;
            }
        }
    
        DistToPlane = PinholeNormal.Dot(PinholeOrigin);
    }


    //Return true for ray going through pinhole (or not hitting the pinhole plane at all).
    //false for ray hitting plane but failing pinhole test
    bool TestRayPinholeIntersect( Vector RaySource, Vector RayDirection,
                                      double &RayLength, Vector &IntersectPoint)
    {
        if( RayPlaneIntersect(PinholeNormal, DistToPlane,
                              RaySource, RayDirection, RayLength, IntersectPoint) )
        {
            if( (IntersectPoint - PinholeOrigin).Magnitude() <= PinholeRadius )
            {
                return true; //ray passes through pinhole
            }
            else
            {
                return false;
            }
        }
        else
        {
            return true; //ray does not hit pinhole plane
        }
    }


};




#endif
