#ifndef _CCD_H_
#define _CCD_H_

#include <cstdlib>
#include "Vector.h"
#include "FileReading.h"


class CCD
{
public:
    Vector CCDNormal, CCDOrigin, CCDXAxis, CCDYAxis;
    double DistToPlane;

    double XPixelWidth, YPixelWidth;
    double InputCCDXMin, InputCCDXMax, InputCCDYMin, InputCCDYMax, InputCCDXSize, InputCCDYSize;
    
    int nXPixels;
    int nYPixels;

public:
    CCD(Vector InputCCDOrigin, Vector InputCCDNormal, double InputCCDAngle,
        int numXPixels, int numYPixels,
        double InputCCDXSize, double InputCCDYSize)
    {
        this->XPixelWidth = InputCCDXSize/double(numXPixels);
        this->YPixelWidth = InputCCDYSize/double(numYPixels);

        this->InputCCDXSize = InputCCDXSize;
        this->InputCCDYSize = InputCCDYSize;
        
        this->InputCCDXMin = -InputCCDXSize/2.0;
        this->InputCCDYMin = -InputCCDYSize/2.0;
        this->InputCCDXMax = InputCCDXSize/2.0;
        this->InputCCDYMax = InputCCDYSize/2.0;

        this->nXPixels = numXPixels;
        this->nYPixels = numYPixels;
        
        
        cout << "XPixelWidth:\t" << XPixelWidth << "\tYPixelWidth:\t" << YPixelWidth << endl;
        
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
            //CCDNormal = CCDNormal * -1.0; Don't do it this way!
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
            exit(1);
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
        
        outX += nXPixels/2;
        outY += nYPixels/2;
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


CCD GenerateCCDFromInputScript( std::string Filename )
{
    ifstream datafile(Filename.c_str());
    if(datafile.is_open() == false)
    {
        cout << "Error: Failed to open InputScript.txt" << endl;
        exit(1);
    }
    std::map<std::string,std::string> InputData;
    AddToMapFromFile(datafile, InputData);
    datafile.close();
    
    Vector InputCCDOrigin(110,-1,50); //origin of CCD
    Vector InputCCDNormal(0,0,1); //direction that CCD points in.
    double InputCCDAngle = 0;
    
    
    double InputCCDXSize = 0.0;
    double InputCCDYSize = 0.0;
    
    
    //double InputCCDXPixelSize = 0.05;
    //double InputCCDYPixelSize = 0.05;
    int NumXPixels, NumYPixels;
    
    VectorFromMap("CCDOrigin",InputData,InputCCDOrigin);
    VectorFromMap("CCDNormal",InputData,InputCCDNormal);
    DoubleFromMap("CCDAngle", InputData,InputCCDAngle);

    //DoubleFromMap("CCDXPixelSize", InputData, InputCCDXPixelSize);
    //DoubleFromMap("CCDYPixelSize", InputData, InputCCDYPixelSize);
    
    IntFromMap("CCDNumXPixels", InputData, NumXPixels);
    IntFromMap("CCDNumYPixels", InputData, NumYPixels);
    
    
    DoubleFromMap("CCDXSize", InputData, InputCCDXSize);
    DoubleFromMap("CCDYSize", InputData, InputCCDYSize);
    
    
    CCD CCDCamera(InputCCDOrigin, InputCCDNormal, InputCCDAngle,
                  NumXPixels, NumYPixels,
                  InputCCDXSize, InputCCDYSize);
    
    return CCDCamera;

}

#endif // _CCD_H_
