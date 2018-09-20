//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


#ifndef KRATOS_MPM_STRESS_PRINCIPAL_INVARIANTS_UTILITY
#define KRATOS_MPM_STRESS_PRINCIPAL_INVARIANTS_UTILITY

// System includes
#include <cmath>
#include <set>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"

namespace Kratos
{

   class MPMStressPrincipalInvariantsUtility
   {

      public:

            typedef Matrix MatrixType;

            typedef Vector VectorType;

            typedef unsigned int IndexType;

            typedef unsigned int SizeType;

            static inline void SortPrincipalStress(Vector& rPrincipalStress, Vector& rMainStrain, Matrix& rMainDirections)
            {
                  // Create Copy
                  Matrix PrincipalDirection1 = ZeroMatrix(3,1);
                  Matrix PrincipalDirection2 = ZeroMatrix(3,1);
                  Matrix PrincipalDirection3 = ZeroMatrix(3,1);

                  for(unsigned int i=0; i<3; i++)
                  {
                        PrincipalDirection1(i,0) = rMainDirections(0,i);
                  }
                  for(unsigned int i=0; i<3; i++)
                  {
                        PrincipalDirection2(i,0) = rMainDirections(1,i);
                  }
                  for(unsigned int i=0; i<3; i++)
                  {
                        PrincipalDirection3(i,0) = rMainDirections(2,i);
                  }

                  // Reorder and swap
                  if(rPrincipalStress[0]<rPrincipalStress[1])
                  {
                        std::swap(rPrincipalStress[0],rPrincipalStress[1]);
                        std::swap(rMainStrain[0],rMainStrain[1]);
                        Matrix TempMatrix = PrincipalDirection1;
                        PrincipalDirection1 = PrincipalDirection2;
                        PrincipalDirection2 = TempMatrix;
                  }

                  if(rPrincipalStress[1]<rPrincipalStress[2])
                  {
                        std::swap(rPrincipalStress[1],rPrincipalStress[2]);
                        std::swap(rMainStrain[1],rMainStrain[2]);
                        Matrix TempMatrix = PrincipalDirection2;
                        PrincipalDirection2 = PrincipalDirection3;
                        PrincipalDirection3 = TempMatrix;
                  }

                  if(rPrincipalStress[0]<rPrincipalStress[1])
                  {
                        std::swap(rPrincipalStress[0],rPrincipalStress[1]);
                        std::swap(rMainStrain[0],rMainStrain[1]);
                        Matrix TempMatrix = PrincipalDirection1;
                        PrincipalDirection1 = PrincipalDirection2;
                        PrincipalDirection2 = TempMatrix;
                  }

                  // Copy back to original matrix
                  for(unsigned int i=0; i<3; i++)
                  {
                        rMainDirections(i,0) = PrincipalDirection1(i,0);
                  }
                  for(unsigned int i=0; i<3; i++)
                  {
                        rMainDirections(i,1) = PrincipalDirection2(i,0);
                  }
                  for(unsigned int i=0; i<3; i++)
                  {
                        rMainDirections(i,2) = PrincipalDirection3(i,0);
                  }
            }

            static inline void CalculateStressInvariants( const Vector& rStress, double& I1, double& J2)
            {
                  // Volumetric Equivalent
                  I1 = 0;
                  for (unsigned int i = 0; i < 3; ++i)
                        I1 += rStress[i];
                  I1 /= 3.0;

                  // Deviatoric Equivalent
                  J2 = 0;
                  for (unsigned int i = 0; i < 3; i++)
                        J2 += std::pow( rStress[i] - I1, 2);

                  if ( rStress.size() == 6 )
                  {
                        for (unsigned int i = 3; i < 6; i++)
                              J2 += 2.0 * std::pow( rStress[i], 2);
                  }

                  J2 = std::sqrt( J2/2.0 );
            }

            static inline void CalculateStressInvariants( const Vector& rStress, double& I1, double& J2, double& rLodeAngle)
            {
                  // CalculateStressInvariants - for I1 and J2
                  CalculateStressInvariants( rStress, I1, J2);
                  
                  // Compute deviatoric stress
                  Matrix deviatoric_stress_tensor = MathUtils<double>::StressVectorToTensor( rStress); 
                  for (unsigned int i = 0; i < 3; ++i)
                        deviatoric_stress_tensor(i,i) -= I1;

                  // Compute Lode Angle
                  rLodeAngle = MathUtils<double>::Det(deviatoric_stress_tensor); // J_3 = det(deviatoric_stress_tensor)
                  rLodeAngle = -1.0 * rLodeAngle / 2.0 * std::pow(3.0/J2, 1.5);

                  double epsilon = 1.0e-9;
                  if ( std::abs(J2) < epsilon ) {                                               // if J2 is 0
                        rLodeAngle = GetPI() / 6.0;
                  } 
                  else if ( std::abs( rLodeAngle ) > (1.0 - epsilon) ) {                        // if current rLodeAngle magnitude is larger than 1.0
                        rLodeAngle = ( GetPI() / 6.0 ) * rLodeAngle / std::abs(rLodeAngle);
                  }                 
                  else {                                                                        // otherwise
                        rLodeAngle = std::asin( rLodeAngle ) / 3.0;
                  }

            }

            static double GetPI()
            {
                  return std::atan(1.0)*4.0;
            }

   }; // end Class MPMStressPrincipalInvariantsUtility

} // end namespace Kratos

#endif // KRATOS_MPM_STRESS_PRINCIPAL_INVARIANTS_UTILITY
