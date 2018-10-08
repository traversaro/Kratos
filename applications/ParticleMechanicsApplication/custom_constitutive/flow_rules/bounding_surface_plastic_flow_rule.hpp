//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


#if !defined(KRATOS_BOUNDING_SURFACE_PLASTIC_FLOW_RULE_H_INCLUDED )
#define      KRATOS_BOUNDING_SURFACE_PLASTIC_FLOW_RULE_H_INCLUDED

// System includes
#include <cmath>

// External includes

// Project includes
#include "custom_constitutive/flow_rules/MPM_flow_rule.hpp"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

//struct MCStressInvariants {

//double MeanStress;
//double J2InvSQ;
//double LodeAngle;

//};

//struct MCSmoothingConstants {

//double A;
//double B;

//};
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
class BoundingSurfacePlasticFlowRule
    :public MPMFlowRule
{



public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NonLinearAssociativePlasticFlowRule
    KRATOS_CLASS_POINTER_DEFINITION( BoundingSurfacePlasticFlowRule );

    struct MaterialParameters
    {
        double SpecificVolume;

    public:
        void PrintInfo()
        {
            KRATOS_INFO("MPMFlowRule.MaterialParameters") << "SpecificVolume = " <<  SpecificVolume  << std::endl;
        }

    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BoundingSurfacePlasticFlowRule();

    /// Initialization constructor.
    BoundingSurfacePlasticFlowRule(YieldCriterionPointer pYieldCriterion);

    /// Copy constructor.
    BoundingSurfacePlasticFlowRule(BoundingSurfacePlasticFlowRule const& rOther);

    /// Assignment operator.
    BoundingSurfacePlasticFlowRule& operator=(BoundingSurfacePlasticFlowRule const& rOther);

    // CLONE
    MPMFlowRule::Pointer Clone() const override;

    /// Destructor.
    ~BoundingSurfacePlasticFlowRule() override;

    bool CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, const Matrix& rIncrementalDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen) override;

    bool UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables ) override;

    Matrix GetElasticLeftCauchyGreen(RadialReturnVariables& rReturnMappingVariables) override;

    unsigned int GetPlasticRegion() override;

    void ComputeElastoPlasticTangentMatrix(const RadialReturnVariables& rReturnMappingVariables, const Matrix& rNewElasticLeftCauchyGreen, const double& alfa, Matrix& rConsistMatrix) override;
    
    void CalculatePrincipalStressTrial(const RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Matrix& rStressMatrix) override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    // /// Turn back information as a string.
    // virtual std::string Info() const;

    // /// Print information about this object.
    // virtual void PrintInfo(std::ostream& rOStream) const;

    // /// Print object's data.
    // virtual void PrintData(std::ostream& rOStream) const;


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    Vector mElasticPrincipalStrain;
    Vector mPlasticPrincipalStrain;

    Vector mPrincipalStressUpdated;
    Vector mPrincipalStressTrial;

    unsigned int mRegion;
    bool mLargeStrainBool;

    MaterialParameters mMaterialParameters;

    double mInitialVolumetricStrain;

    Vector mCenterOfHomologyStress;
    Vector mPreviousStress;
    Vector mImagePointStress;
    
    double mPreviousMeanStressP;
    double mPreviousDeviatoricStressQ;

    double mPlasticMultiplier;

    bool mImagePointComputedBool;
    bool mIsOnceUnloaded;

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
    void InitializeMaterial(YieldCriterionPointer& pYieldCriterion, HardeningLawPointer& pHardeningLaw, const Properties& rProp) override;

    void InitializeMaterialParameters();


    void CalculatePrincipalStressVector(const Vector& rPrincipalStrain, Vector& rPrincipalStress);

    void CalculatePrincipalStrainFromStrainInvariants(Vector& rPrincipalStrain, const double& rVolumetricStrain, const double& rDeviatoricStrain, const Vector& rDirectionVector);

    void CalculateStrainInvariantsFromPrincipalStrain(const Vector& rPrincipalStrain, double& rVolumetricStrain, double& rDeviatoricStrain, Vector& rDeviatoricStrainVector);

    void ReturnStressFromPrincipalAxis(const Matrix& rEigenVectors, const Vector& rPrincipalStress, Matrix& rStressMatrix);

    void CalculateTransformationMatrix(const Matrix& rMainDirection, Matrix& rA);

    
    bool CalculateConsistencyCondition(const RadialReturnVariables& rReturnMappingVariables, const Vector& rPrincipalStress, const Vector& rPrincipalStrain, unsigned int& region, Vector& rPrincipalStressUpdated);
 
    void CalculatePlasticMultiplier(const Vector& rDirectionN, const Vector& rDirectionM, const double& rHardening, const Matrix& rElasticMatrix, const Vector rPrincipalStrain, double& rPlasticStrainMultiplier);

    void CalculateImagePointStress(const Vector& rCenterOfHomologyStress, const Vector& rCurrentStress, Vector& rImagePointStress, double& rConstantB, const bool& rBIsKnown = false);
    
    void CalculateCenterOfHomologyStress(Vector& rCenterOfHomologyStress);


    void ComputeElasticMatrix(const double& rMeanStressP, Matrix& rElasticMatrix);

    void ComputeInverseElasticMatrix(const double& rMeanStressP, Matrix& rInverseElasticMatrix);

    void ComputePlasticMatrix(const Vector& rDirectionN, const Vector& rDirectionM, const double& rHardening, const Matrix& rElasticMatrix, Matrix& rPlasticMatrix);

    void CalculateModificationMatrix(Matrix& rModMatrixT);

    
    void CalculateLoadingDirection(const Vector& rPrincipalStressVector, Vector& rLoadingDirection);

    void CalculatePlasticFlowDirection(const Vector& rPrincipalStressVector, const Vector& rImagePointStressVector, Vector& rPlasticFlowDirection);

    void CalculateYieldSurfaceDerivatives(const Vector& rPrincipalStressVector, Vector& rFirstDerivative);

    void CalculatePlasticPotentialDerivatives(const Vector& rPrincipalStressVector, const Vector& rImagePointPrincipalStressVector, Vector& rFirstDerivative);

    void CalculatePlasticPotentialInvariantDerivatives(const Vector& rPrincipalStressVector, const Vector& rImagePointPrincipalStressVector, Vector& rFirstDerivative);

    void CalculatePlasticPotentialSecondDerivatives(const Vector& rPrincipalStressVector, const Vector& rImagePointPrincipalStressVector, Matrix& rSecondDerivative);

    void CalculatePlasticPotentialInvariantSecondDerivatives(const Vector& rPrincipalStressVector, const Vector& rImagePointPrincipalStressVector, Vector& rSecondDerivative);


    double CalculateCriticalStateLineSlope(const double& rLodeAngle);

    double GetAlphaParameter();

    double GetDirectionParameter(const Vector& rPrincipalStressVector, const Vector& rImagePointPrincipalStressVector);

    double GetDirectionAngle(const Vector& rPrincipalStressVector);

    double GetPI();
  
    //virtual void GetPrincipalStressAndStrain(Vector& PrincipalStresses, Vector& PrincipalStrains);
    ///@}
    ///@name Protected  Access
    ///@{


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


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class NonLinearAssociativePlasticFlowRule

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
// 				    NonLinearAssociativePlasticFlowRule& rThis);

// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
// 				    const NonLinearAssociativePlasticFlowRule& rThis)
// {
//   rThis.PrintInfo(rOStream);
//   rOStream << std::endl;
//   rThis.PrintData(rOStream);

//   return rOStream;
// }
///@}

///@} addtogroup block



///@}
///@ Template Operations
///@{


///@}


}  // namespace Kratos.

#endif // KRATOS_BOUNDING_SURFACE_PLASTIC_FLOW_RULE_H_INCLUDED  defined 
