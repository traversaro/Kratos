//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Albert Puigferrat Perez
//                   Ignasi de Pouplana
//

#if !defined(KRATOS_ELEMENT_UTILITIES )
#define  KRATOS_ELEMENT_UTILITIES

// Project includes
#include "includes/define.h"

// System includes
//#include <cmath>

// Project includes
//#include "utilities/math_utils.h"
#include "includes/element.h"

// Application includes

namespace Kratos
{

class ElementUtilities
{


public:

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void InterpolateVariableWithComponents(array_1d<double,2>& rVector,const Matrix& Ncontainer,
                                        const array_1d<array_1d<double,3>, 3>& VariableWithComponents,const unsigned int& GPoint)
    {
        noalias(rVector) = ZeroVector(2);

        for(unsigned int i=0; i<3; i++)
        {
            rVector[0] += Ncontainer(GPoint,i)*VariableWithComponents[i][0];
            rVector[1] += Ncontainer(GPoint,i)*VariableWithComponents[i][1];
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void InterpolateVariableWithComponents(array_1d<double,2>& rVector,const Matrix& Ncontainer,
                                        const array_1d<array_1d<double,3>, 4>& VariableWithComponents,const unsigned int& GPoint)
    {
        noalias(rVector) = ZeroVector(2);

        for(unsigned int i=0; i<4; i++)
        {
            rVector[0] += Ncontainer(GPoint,i)*VariableWithComponents[i][0];
            rVector[1] += Ncontainer(GPoint,i)*VariableWithComponents[i][1];
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    static inline void InterpolateVariableWithComponents(array_1d<double,3>& rVector,const Matrix& Ncontainer,
                                        const array_1d<array_1d<double,3>, 4>& VariableWithComponents,const unsigned int& GPoint)
    {
        noalias(rVector) = ZeroVector(3);

        for(unsigned int i=0; i<4; i++)
        {
            rVector[0] += Ncontainer(GPoint,i)*VariableWithComponents[i][0];
            rVector[1] += Ncontainer(GPoint,i)*VariableWithComponents[i][1];
            rVector[2] += Ncontainer(GPoint,i)*VariableWithComponents[i][2];
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void InterpolateVariableWithComponents(array_1d<double,3>& rVector,const Matrix& Ncontainer,
                                        const array_1d<array_1d<double,3>, 8>& VariableWithComponents,const unsigned int& GPoint)
    {
        noalias(rVector) = ZeroVector(3);

        for(unsigned int i=0; i<8; i++)
        {
            rVector[0] += Ncontainer(GPoint,i)*VariableWithComponents[i][0];
            rVector[1] += Ncontainer(GPoint,i)*VariableWithComponents[i][1];
            rVector[2] += Ncontainer(GPoint,i)*VariableWithComponents[i][2];
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // TODO
    // Triangle2D3 version.
    static inline double ProjectedLength(const Geometry<Node<3> >& rGeometry,
                                                            const array_1d<double,3>& rVelocity)
    {
        double Hvel = 0.0;
        double CosAngle = 0.0;
        double EdgeNorm = 0.0;
        double NormVel = norm_2 (rVelocity);

        const unsigned int NumNodes = 3;

        if (NormVel > 1e-12)
        {
            // Loop over edges looking for maximum 'projected' length
            array_1d<double,3> Edge(3,0.0);
            for(unsigned int i = 0; i < NumNodes; ++i)
            {
                unsigned int j = (i+1) % NumNodes;
                Edge = rGeometry[j] - rGeometry[i];
                double InnerProd = inner_prod (rVelocity, Edge);
                double EdgeNormAux = norm_2 (Edge);

                double CosAngleAux = InnerProd / (NormVel * EdgeNormAux);

                if (EdgeNormAux > EdgeNorm)
                {
                    EdgeNorm = EdgeNormAux;
                }

                if (CosAngleAux > CosAngle)
                {
                    CosAngle = CosAngleAux;
                }

                double HvelAux = EdgeNorm * CosAngle;

                if (HvelAux > Hvel)
                {
                    Hvel = HvelAux;
                }

            }
        }

        return Hvel;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    static inline void FillArray1dOutput(array_1d<double,3>& rOutputValue, const array_1d<double,2>& ComputedValue)
    {
        rOutputValue[0] = ComputedValue[0];
        rOutputValue[1] = ComputedValue[1];
        rOutputValue[2] = 0.0;
    }

    //----------------------------------------------------------------------------------------

    static inline void FillArray1dOutput(array_1d<double,3>& rOutputValue, const array_1d<double,3>& ComputedValue)
    {
        rOutputValue[0] = ComputedValue[0];
        rOutputValue[1] = ComputedValue[1];
        rOutputValue[2] = ComputedValue[2];
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

}; /* Class ElementUtilities*/
} /* namespace Kratos.*/

#endif /* KRATOS_ELEMENT_UTILITIES defined */