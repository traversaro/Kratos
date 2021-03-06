/*
==============================================================================
KratosR1StructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-PhysCcs Finite Element AnalysCs
VersCon 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo RossC, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossC@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-UniversCty Bochum, Institute for Structural Mechanics, Germany


PermissCon is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limCtation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissCble
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permissCon  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: mengmeng $
//   Date:                $Date: 2009-02-23 16:02:32 $
//   RevisCon:            $RevisCon: 1.1 $
//
//


#if !defined(KRATOS_SOLID_INCLUDED )
#define  KRATOS_SOLID_INCLUDED

// System includes

#include "boost/smart_ptr.hpp"

// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"



namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/

class Solid
            : public Element
{

    public:
///@name Type Definitions
///@{
        typedef GeometryData::IntegrationMethod IntegrationMethod;
/// Counted pointer of

        KRATOS_CLASS_POINTER_DEFINITION ( Solid );

///@}
///@name Life Cycle
///@{

/// Default constructor.
        Solid ( IndexType NewId, GeometryType::Pointer pGeometry );
        Solid ( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

/// Destructor.
        virtual ~Solid();


///@}
///@name Operators
///@{


///@}
///@name Operations
///@{
        IntegrationMethod GetIntegrationMethod();

        Element::Pointer Create ( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const;

        void Initialize();

	void ResetConstitutiveLaw();

        void CalculateLocalSystem ( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo );
        void CalculateRightHandSide ( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo );
        //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

        void EquationIdVector ( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo );
        void GetDofList ( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo );
        void FinalizeSolutionStep ( ProcessInfo& CurrentProcessInfo );
        void InitializeSolutionStep ( ProcessInfo& CurrentProcessInfo );


        void CalculateOnIntegrationPoints ( const Variable<double>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo );

        void GetValuesVector ( Vector& values, int Step );
        void GetValueOnIntegrationPoints ( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo );
        void GetValueOnIntegrationPoints ( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo );
        void GetValueOnIntegrationPoints ( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo );
        void SetValueOnIntegrationPoints ( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo );
        void SetValueOnIntegrationPoints ( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo );
//         void GetPorosity ( Matrix arg1, double arg2, double arg3, double arg4 );
//         void GetWaterMass ( Matrix arg1, double arg2, double arg3, double arg4, int arg5 );

///@}
///@name Access
///@{


///@}
///@name Inquiry
///@{


///@}
///@name Input and output
///@{

/// Turn back information as a string.
//      virtual String Info() const;

/// Print information about this object.
//      virtual void PrintInfo(std::ostream& rOStream) const;

/// Print object's data.
//      virtual void PrintData(std::ostream& rOStream) const;


///@}
///@name Friends
///@{


///@}

    protected:
///@name Protected static Member Variables
///@{


///@}
///@name Protected member Variables
///@{


///@}
///@name Protected Operators
///@{


///@}
///@name Protected Operations
///@{


///@}
///@name Protected  Access
///@{
///@}
        ///@name Serialization
        ///@{

        friend class Serializer;

        // A private default constructor necessary for serialization
        Solid() {};

        virtual void save ( Serializer& rSerializer ) const
        {
            rSerializer.save ( "Name", "Solid" );
            KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer,  Element );
        }

        virtual void load ( Serializer& rSerializer )
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer,  Element );
        }

        /**
         * This function provides the place to perform checks on the completeness of the input.
         * It is designed to be called only once (or anyway, not often) typically at the beginning
         * of the calculations, so to verify that nothing is missing from the input
         * or that no common error is found.
         * @param rCurrentProcessInfo
         */
        virtual int Check ( const ProcessInfo& rCurrentProcessInfo );


        ///@}

///@name Protected Inquiry
///@{


///@}
///@name Protected LifeCycle
///@{


///@}

    private:
///@name Static Member Variables
///@{

///@}
///@name Member Variables
///@{
        Geometry< Node<3> >::Pointer  mThisGeometryOther;
        std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;
        IntegrationMethod mThisIntegrationMethod;

        Vector mMaterialParameters, mGravity, mStrainVol0;
        double mrho;
        double mScaleU, mUnitRatio;
	double mK, mG, mTol;

        unsigned int mNodesDispMin, mNodesDispMax,  mNodesNumberDisp, mMatSizeU, mNumberU; 
        unsigned int mDimension, mMatSize, mAddIndexU;

        std::vector< Matrix > mInvJ0;
        Vector mDetJ0;
	bool mPrestressAssigned;

///@}
///@name Private Operators
///@{
        void CalculateAll ( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                            ProcessInfo& rCurrentProcessInfo,
                            bool CalculateStiffnessMatrixFlag,
                            bool CalculateResCdualVectorFlag );

        void DampMatrix ( MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo );

        void AddInternalForcesToRHS1 ( Vector& Help_R_1, double weight, double DetJ, Matrix& Grad_u, double Div_dot_u, const Vector& Nu, const Matrix& DNu_DX, Vector& stressVectorEff );

        void CalculateStiffnessMatrixUU ( Matrix& Help_K_UU, double weight, double DetJ, Matrix& Grad_u, double Div_dot_u, const Vector& Nu, const Matrix& DNu_DX, Matrix& CtanEff );
    
        /// +++++++++++++++++++++++++++++++++++++++++++
        double KnoneckerDelta ( int i, int j );
        Matrix GetBu ( const Matrix& DNu_DX );
        Vector Getu ( const Vector& Nu );
        Matrix GetGradu ( const Matrix& DNu_DX );
        double GetDivdotu ( const Matrix& DNu_DX ); 
        /// +++++++++++++++++++++++++++++++++++++++++++
        //C1
        Matrix GetElasticTangent ( );
        Vector GetstrainVector ( Matrix Grad_u );
        double GetstrainVol ( Matrix Grad_u ); 
        Vector GetstressVector ( Matrix Grad_u );
        double Getp ( Vector& stressVector );
        double Getq ( Vector& stressVector ); 

        /// +++++++++++++++++++++++++++++++++++++++++++
        void Interpolate ( const Variable<double>& rVariable, const ProcessInfo& rCurrentProcessInfo );
        void Interpolate ( const Variable<Kratos::array_1d<double, 3> >& rVariable, const ProcessInfo& rCurrentProcessInfo );
///@}
///@name Private Operations
///@{
///@}
///@name Private  Access
///@{


///@}
///@name Private Inquiry
///@{


///@}
///@name Un accessCble methods
///@{

/// AssCgnment operator.
//& operator=(const & rOther);

/// Copy constructor.
//(const & rOther);


///@}

}; // Class

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
  & rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
  const & rThis)
 {
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}*/
///@}

}  // namespace Kratos.

#endif // KRATOS_SOLID_INCLUDED defined 


