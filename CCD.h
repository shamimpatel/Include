#ifndef _CCD_H_
#define _CCD_H_

#include <cstdlib>
#include "Vector.h"


class CCD
{
    Vector CCDNormal, CCDOrigin, CCDXAxis, CCDYAxis;
    double DistToPlane;

    double XPixelWidth, YPixelWidth;
    double InputCCDXMin, InputCCDXMax, InputCCDYMin, InputCCDYMax;

public:
    CCD(Vector InputCCDOrigin, Vector InputCCDNormal, double InputCCDAngle,
        double XPixelWidth, double YPixelWidth,
        double InputCCDXMin, double InputCCDXMax,
        double InputCCDYMin, double InputCCDYMax)
    {
        this->XPixelWidth = XPixelWidth;
        this->YPixelWidth = YPixelWidth;

        this->InputCCDXMin = InputCCDXMin;
        this->InputCCDYMin = InputCCDYMin;
        this->InputCCDXMax = InputCCDXMax;
        this->InputCCDYMax = InputCCDYMax;


        CCDNormal = InputCCDNormal;
        CCDOrigin = InputCCDOrigin;

        CCDXAxis = Vector(1,0,0);
        CCDYAxis = Vector(0,1,0);

        if( InputCCDNormal.Magnitude() != 1.0)
        {
            cout << "Input Normal not normalized" << endl;
            CCDNormal = InputCCDNormal.Normalized();
        }

        if(InputCCDOrigin.Dot(InputCCDNormal) < 0.0)
        {
            cout << "Input Normal points towards Origin" << endl;
            //need to check like this or 0.0 --> -0.0
            //-0.0 with the trig functions goes to the other end of the output range.
            if(CCDNormal.x != 0.0)
            {
                CCDNormal.x *= -1.0;
            }
            if(CCDNormal.y != 0.0)
            {
                CCDNormal.y *= -1.0;
            }
            if(CCDNormal.z != 0.0)
            {
                CCDNormal.z *= -1.0;
            }
            //CCDNormal = CCDNormal * -1.0;
        }

        cout << "CCDNormal:\t";
        CCDNormal.Print();

        DistToPlane = CCDNormal.Dot(CCDOrigin);

        double NormalTheta, NormalPhi;
        GetAlternateThetaPhi( CCDNormal, NormalTheta, NormalPhi );

        //Rotate the CCD Axes so that they're perpendicular to the normal.
        CCDXAxis = RotateY(CCDXAxis,NormalTheta);
        CCDXAxis = RotateZ(CCDXAxis,NormalPhi);

        CCDYAxis = RotateY(CCDYAxis,NormalTheta);
        CCDYAxis = RotateZ(CCDYAxis,NormalPhi);

        //Need small tolerance as it'll come out as 1E-8 instead of zero
        if(fabs(CCDXAxis.Dot(CCDYAxis.Cross(CCDNormal)) - 1.0) >= 0.0001)
        {
            cout << "Error: CCD Axes and CCD Normal not orthogonal!" << endl;
            exit(0);
        }

        //spin around input normal not the "real" normal which may be modified.
        //alternatively could have span it before rotating the ccd x/y vectors?
        CCDXAxis = RotationAxis(CCDXAxis,InputCCDNormal.Normalized(),InputCCDAngle);
        CCDYAxis = RotationAxis(CCDYAxis,InputCCDNormal.Normalized(),InputCCDAngle);

        cout << "CCDXAxis:\t";
        CCDXAxis.Print();
        cout << "CCDYAxis:\t";
        CCDYAxis.Print();

    }

    void GenerateCCDCorners( Vector &Corner1, Vector &Corner2, Vector &Corner3, Vector &Corner4)
    {
        Corner1 = CCDOrigin + CCDXAxis*InputCCDXMin + CCDYAxis*InputCCDYMin;
        Corner2 = CCDOrigin + CCDXAxis*InputCCDXMin + CCDYAxis*InputCCDYMax;
        Corner3 = CCDOrigin + CCDXAxis*InputCCDXMax + CCDYAxis*InputCCDYMin;
        Corner4 = CCDOrigin + CCDXAxis*InputCCDXMax + CCDYAxis*InputCCDYMax;
    }

    Vector GetCCDXAxis()
    {
        return CCDXAxis;
    }

    Vector GetCCDYAxis()
    {
        return CCDYAxis;
    }

    bool RayCCDIntersect(Vector RaySource, Vector RayDirection,
                         double &DistanceFromSource, Vector &Intersection)
    {
        return RayPlaneIntersect(CCDNormal,DistToPlane,RaySource,
                                 RayDirection,DistanceFromSource,Intersection);
    }

    void ProjectToCCD( Vector V, double &X, double &Y)
    {

        Vector Point = V - CCDOrigin;
        X = Point.Dot(CCDXAxis);
        Y = Point.Dot(CCDYAxis);
    }

    //this should be give the X/Y coords after being projected onto the CCD Plane.
    void GetPixelFromCCDXY( double inX, double inY,
                            int &outX,  int &outY)
    {
        outX = inX/XPixelWidth;
        outY = inY/YPixelWidth;

        if(inX < 0)
        {
            outX = outX - 1;
        }
        if(inY < 0)
        {
            outY = outY - 1;
        }
    }

    bool GetPixelRayCCDIntersect(Vector RaySource, Vector RayDirection,
                                 double &RayLength, Vector &IntersectPoint,
                                 double &CCDIntersectX, double &CCDIntersectY,
                                 int &XPixel, int &YPixel)
    {
        if( RayCCDIntersect(RaySource, RayDirection,
                            RayLength, IntersectPoint) == false)
        {
            return false;
        }

        ProjectToCCD(IntersectPoint,CCDIntersectX,CCDIntersectY);
        GetPixelFromCCDXY(CCDIntersectX,CCDIntersectY,XPixel,YPixel);

        return true;
    }
};

#endif // _CCD_H_
