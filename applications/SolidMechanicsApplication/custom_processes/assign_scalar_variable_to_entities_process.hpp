//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2016 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_ASSIGN_SCALAR_VARIABLES_TO_ENTITIES_PROCESS_H_INCLUDED)
#define  KRATOS_ASSIGN_SCALAR_VARIABLES_TO_ENTITIES_PROCESS_H_INCLUDED



// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// The base class for assigning a value to scalar variables or array_1d components processes in Kratos.
/** This function assigns a value to a variable belonging to all of the nodes in a given mesh
*/
class AssignScalarVariableToEntitiesProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AssignScalarVariableToEntitiesProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignScalarVariableToEntitiesProcess);

    enum EntityType { NODES, CONDITIONS, ELEMENTS };

    enum AssignmentType { DIRECT, ADDITION, SUBTRACTION, MULTIPLICATION, DIVISION };


    ///@}
    ///@name Life Cycle
    ///@{
    AssignScalarVariableToEntitiesProcess(ModelPart& rModelPart) : Process(Flags()), mrModelPart(rModelPart) {}

    AssignScalarVariableToEntitiesProcess(ModelPart& rModelPart,
                                          Parameters rParameters) : Process(Flags()) , mrModelPart(rModelPart)
    {
        KRATOS_TRY

        Parameters default_parameters( R"(
            {
                "model_part_name":"MODEL_PART_NAME",
                "variable_name": "VARIABLE_NAME",
                "entity_type": "NODES",
                "value" : 1.0,
                "compound_assignment": "direct"
            }  )" );


        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mvariable_name = rParameters["variable_name"].GetString();

        if( rParameters["entity_type"].GetString() == "NODES" ){
          mEntity = NODES;
        }
        else if(  rParameters["entity_type"].GetString() == "CONDITIONS" ){
          mEntity = CONDITIONS;
        }
        else if(  rParameters["entity_type"].GetString() == "ELEMENTS" ){
          mEntity = ELEMENTS;
        }
        else{
          KRATOS_ERROR <<" Entity type "<< rParameters["entity_type"].GetString() << " is not supported " << std::endl;
        }

        if( mEntity == NODES ){

          if( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(mvariable_name) ) //case of component variable
          {
            typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
            component_type var_component = KratosComponents< component_type >::Get(mvariable_name);

            if( rModelPart.GetNodalSolutionStepVariablesList().Has( var_component.GetSourceVariable() ) == false )
            {
              KRATOS_ERROR << "trying to set a variable that is not in the model_part - variable name is " << mvariable_name << std::endl;
            }

            mdouble_value = rParameters["value"].GetDouble();
          }
          else if( KratosComponents< Variable<double> >::Has( mvariable_name ) ) //case of double variable
          {
            if( rModelPart.GetNodalSolutionStepVariablesList().Has( KratosComponents< Variable<double> >::Get( mvariable_name ) ) == false )
            {
              KRATOS_ERROR << "trying to set a variable that is not in the model_part - variable name is " << mvariable_name << std::endl;
            }

            mdouble_value = rParameters["value"].GetDouble();
          }
          else if( KratosComponents< Variable<int> >::Has( mvariable_name ) ) //case of int variable
          {
            if( rModelPart.GetNodalSolutionStepVariablesList().Has( KratosComponents< Variable<int> >::Get( mvariable_name ) ) == false )
            {
              KRATOS_ERROR << "trying to set a variable that is not in the model_part - variable name is " << mvariable_name << std::endl;
            }

            mint_value = rParameters["value"].GetInt();

          }
          else if( KratosComponents< Variable<bool> >::Has( mvariable_name ) ) //case of bool variable
          {
	    if( rModelPart.GetNodalSolutionStepVariablesList().Has( KratosComponents< Variable<bool> >::Get( mvariable_name ) ) == false )
            {
              KRATOS_ERROR << "trying to set a variable that is not in the model_part - variable name is " << mvariable_name << std::endl;
            }

	    mbool_value = rParameters["value"].GetBool();
          }

        }
        else if( mEntity == CONDITIONS || mEntity == ELEMENTS ){

          if( KratosComponents< Variable<double> >::Has( mvariable_name ) ) //case of double variable
          {
            mdouble_value = rParameters["value"].GetDouble();
          }
          else if( KratosComponents< Variable<int> >::Has( mvariable_name ) ) //case of int variable
          {
            mint_value = rParameters["value"].GetInt();
          }
          else if( KratosComponents< Variable<bool> >::Has( mvariable_name ) ) //case of bool variable
          {
            mbool_value = rParameters["value"].GetBool();
          }
          else{
            KRATOS_ERROR << "trying to set a variable that is not in the model_part - variable name is " << mvariable_name << std::endl;
          }

        }

        this->SetAssignmentType(rParameters["compound_assignment"].GetString(), mAssignment);

        KRATOS_CATCH("")
    }


    /// Destructor.
    virtual ~AssignScalarVariableToEntitiesProcess() {}


    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{


    /// Execute method is used to execute the AssignScalarVariableToEntitiesProcess algorithms.
    void Execute()  override
    {

        KRATOS_TRY

        if( mEntity == NODES ){

          if( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(mvariable_name) ) //case of component variable
          {
            typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
            component_type var_component = KratosComponents< component_type >::Get(mvariable_name);
            AssignValueToNodes< component_type, double>(var_component , mdouble_value);
          }
          else if( KratosComponents< Variable<double> >::Has( mvariable_name ) ) //case of double variable
          {
            AssignValueToNodes<>(KratosComponents< Variable<double> >::Get(mvariable_name), mdouble_value);
          }
          else if( KratosComponents< Variable<int> >::Has( mvariable_name ) ) //case of int variable
          {
            AssignValueToNodes<>(KratosComponents< Variable<int> >::Get(mvariable_name) , mint_value);
          }
          else if( KratosComponents< Variable<bool> >::Has( mvariable_name ) ) //case of bool variable
          {
            AssignValueToNodes<>(KratosComponents< Variable<bool> >::Get(mvariable_name), mbool_value);
          }
          else
          {
            KRATOS_ERROR << "Not able to set the variable. Attempting to set variable:" << mvariable_name << std::endl;
          }

        }
        else if( mEntity == CONDITIONS ) {

          if( KratosComponents< Variable<double> >::Has( mvariable_name ) ) //case of double variable
          {
            AssignValueToConditions<>(KratosComponents< Variable<double> >::Get(mvariable_name), mdouble_value);
          }
          else if( KratosComponents< Variable<int> >::Has( mvariable_name ) ) //case of int variable
          {
            AssignValueToConditions<>(KratosComponents< Variable<int> >::Get(mvariable_name) , mint_value);
          }
          else if( KratosComponents< Variable<bool> >::Has( mvariable_name ) ) //case of bool variable
          {
            AssignValueToConditions<>(KratosComponents< Variable<bool> >::Get(mvariable_name), mbool_value);
          }
          else
          {
            KRATOS_ERROR << "Not able to set the variable. Attempting to set variable:" << mvariable_name << std::endl;
          }

        }
        else if( mEntity == ELEMENTS ) {

          if( KratosComponents< Variable<double> >::Has( mvariable_name ) ) //case of double variable
          {
            AssignValueToElements<>(KratosComponents< Variable<double> >::Get(mvariable_name), mdouble_value);
          }
          else if( KratosComponents< Variable<int> >::Has( mvariable_name ) ) //case of int variable
          {
            AssignValueToElements<>(KratosComponents< Variable<int> >::Get(mvariable_name) , mint_value);
          }
          else if( KratosComponents< Variable<bool> >::Has( mvariable_name ) ) //case of bool variable
          {
            AssignValueToElements<>(KratosComponents< Variable<bool> >::Get(mvariable_name), mbool_value);
          }
          else
          {
            KRATOS_ERROR << "Not able to set the variable. Attempting to set variable:" << mvariable_name << std::endl;
          }

        }

        KRATOS_CATCH("");

    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
    }

    /// this function is designed for being execute once before the solution loop but after all of the
    /// solvers where built
    void ExecuteBeforeSolutionLoop() override
    {
    }


    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    void ExecuteFinalizeSolutionStep() override
    {
    }


    /// this function will be executed at every time step BEFORE  writing the output
    void ExecuteBeforeOutputStep() override
    {
    }


    /// this function will be executed at every time step AFTER writing the output
    void ExecuteAfterOutputStep() override
    {
    }


    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    void ExecuteFinalize() override
    {

        KRATOS_TRY

        if( mEntity == CONDITIONS ){

          if( KratosComponents< Variable<double> >::Has( mvariable_name ) ) //case of double variable
          {
            double double_value = 0;
            AssignValueToConditions<>(KratosComponents< Variable<double> >::Get(mvariable_name), double_value);
          }
          else if( KratosComponents< Variable<int> >::Has( mvariable_name ) ) //case of int variable
          {
            int int_value = 0;
            AssignValueToConditions<>(KratosComponents< Variable<int> >::Get(mvariable_name), int_value);
          }
          else if( KratosComponents< Variable<bool> >::Has( mvariable_name ) ) //case of bool variable
          {
            bool bool_value = !mbool_value;
            AssignValueToConditions<>(KratosComponents< Variable<bool> >::Get(mvariable_name), bool_value);
          }
          else
          {
            KRATOS_ERROR << "Not able to set the variable. Attempting to set variable:" << mvariable_name << std::endl;
          }
        }

        KRATOS_CATCH("")

    }


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
    std::string Info() const override
    {
        return "AssignScalarVariableToEntitiesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AssignScalarVariableToEntitiesProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }


    ///@}
    ///@name Friends
    ///@{
    ///@}

protected:

    ///@name Protected static Member Variables
    ///@{

    ModelPart& mrModelPart;
    std::string mvariable_name;

    EntityType mEntity;
    AssignmentType mAssignment;

    ///@}
    ///@name Protected member Variables
    ///@{
    ///@}
    ///@name Protected Operators
    ///@{

    /// Copy constructor.
    AssignScalarVariableToEntitiesProcess(AssignScalarVariableToEntitiesProcess const& rOther);

    ///@}
    ///@name Protected Operations
    ///@{

    // set assignment method

    void SetAssignmentType(std::string method, AssignmentType& rAssignment)
    {
        //compound_assignment:

        //implemented:
        //  = direct
        // += addition
        // -= subtraction
        // *= multiplication
        // /= division

        if( method == "direct" ){
          rAssignment = DIRECT;
        }
        else if(  method == "addition" ){
          rAssignment = ADDITION;
        }
        else if(  method == "subtraction" ){
          rAssignment = SUBTRACTION;
        }
        else if(  method == "multiplication" ){
          rAssignment = MULTIPLICATION;
        }
        else if(  method == "division" ){
          rAssignment = DIVISION;
        }
        else{
          KRATOS_ERROR <<" Assignment type "<< method << " is not supported " << std::endl;
        }

    }

    // nodes

    template< class TVarType, class TDataType >
    void DirectAssignValue(ModelPart::NodeType& rNode, const TVarType& rVariable, const TDataType& value)
    {
      rNode.FastGetSolutionStepValue(rVariable) = value;
    }

    template< class TVarType, class TDataType >
    void AddAssignValue(ModelPart::NodeType& rNode, const TVarType& rVariable, const TDataType& value)
    {
      rNode.FastGetSolutionStepValue(rVariable) += value;
    }

    template< class TVarType, class TDataType >
    void SubtractAssignValue(ModelPart::NodeType& rNode, const TVarType& rVariable, const TDataType& value)
    {
      rNode.FastGetSolutionStepValue(rVariable) -= value;
    }

    template< class TVarType, class TDataType >
    void MultiplyAssignValue(ModelPart::NodeType& rNode, const TVarType& rVariable, const TDataType& value)
    {
      rNode.FastGetSolutionStepValue(rVariable) *= value;
    }

    template< class TVarType, class TDataType >
    void DivideAssignValue(ModelPart::NodeType& rNode, const TVarType& rVariable, const TDataType& value)
    {
      if(value!=0)
        rNode.FastGetSolutionStepValue(rVariable) /= value;
    }

    // elements and conditions

    template< class TEntityType, class TVarType, class TDataType >
    void DirectAssignValue(TEntityType& rEntity, const TVarType& rVariable, const TDataType& value)
    {
      rEntity.SetValue(rVariable,value);
    }

    template< class TEntityType, class TVarType, class TDataType >
    void AddAssignValue(TEntityType& rEntity, const TVarType& rVariable, const TDataType& value)
    {
      TDataType AddedValue = rEntity.GetValue(rVariable)+value;
      rEntity.SetValue(rVariable,AddedValue);
    }

    template< class TEntityType, class TVarType, class TDataType >
    void SubtractAssignValue(TEntityType& rEntity, const TVarType& rVariable, const TDataType& value)
    {
      TDataType SubtractedValue = rEntity.GetValue(rVariable)-value;
      rEntity.SetValue(rVariable,SubtractedValue);
    }

    template< class TEntityType, class TVarType, class TDataType >
    void MultiplyAssignValue(TEntityType& rEntity, const TVarType& rVariable, const TDataType value)
    {
      TDataType MultipliedValue = rEntity.GetValue(rVariable)*value;
      rEntity.SetValue(rVariable,MultipliedValue);
    }

    template< class TEntityType, class TVarType, class TDataType >
    void DivideAssignValue(TEntityType& rEntity, const TVarType& rVariable, const TDataType& value)
    {
      TDataType DividedValue = rEntity.GetValue(rVariable);
      if(value!=0)
        DividedValue/=value;
      rEntity.SetValue(rVariable,DividedValue);
    }

    template< class TEntityType >
    void MultiplyAssignValue(TEntityType& rEntity, const Variable<Vector>& rVariable, const Vector& value)
    {
      Vector MultipliedValue = rEntity.GetValue(rVariable);
      for(unsigned int i=0; i<MultipliedValue.size(); ++i)
      {
        MultipliedValue[i]*=value[i];
      }
      rEntity.SetValue(rVariable,MultipliedValue);
    }

    template< class TEntityType >
    void DivideAssignValue(TEntityType& rEntity, const Variable<Vector>& rVariable, const Vector& value)
    {
      Vector DividedValue = rEntity.GetValue(rVariable);
      for(unsigned int i=0; i<DividedValue.size(); ++i)
      {
        if(value[i]!=0)
          DividedValue[i]/=value[i];
      }
      rEntity.SetValue(rVariable,DividedValue);
    }

    template< class TEntityType >
    void MultiplyAssignValue(TEntityType& rEntity, const Variable<array_1d<double,3>>& rVariable, const array_1d<double,3>& value)
    {
      Vector MultipliedValue = rEntity.GetValue(rVariable);
      for(unsigned int i=0; i<3; ++i)
      {
        MultipliedValue[i]*=value[i];
      }
      rEntity.SetValue(rVariable,MultipliedValue);
    }

    template< class TEntityType >
    void DivideAssignValue(TEntityType& rEntity, const Variable<array_1d<double,3>>& rVariable, const array_1d<double,3>& value)
    {
      Vector DividedValue = rEntity.GetValue(rVariable);
      for(unsigned int i=0; i<3; ++i)
      {
        if(value[i]!=0)
          DividedValue[i]/=value[i];
      }
      rEntity.SetValue(rVariable,DividedValue);
    }




    template< class TMethodPointerType >
    TMethodPointerType GetAssignmentMethod()
    {
        TMethodPointerType AssignmentMethod = nullptr;
        switch( mAssignment )
        {
          case DIRECT:
            AssignmentMethod = &AssignScalarVariableToEntitiesProcess::DirectAssignValue;
            break;
          case ADDITION:
            AssignmentMethod = &AssignScalarVariableToEntitiesProcess::AddAssignValue;
            break;
          case SUBTRACTION:
            AssignmentMethod = &AssignScalarVariableToEntitiesProcess::SubtractAssignValue;
            break;
          case MULTIPLICATION:
            AssignmentMethod = &AssignScalarVariableToEntitiesProcess::MultiplyAssignValue;
            break;
          case DIVISION:
            AssignmentMethod = &AssignScalarVariableToEntitiesProcess::DivideAssignValue;
            break;
          default:
            KRATOS_ERROR << "Unexpected value for Assignment method: " << mAssignment << std::endl;
        }
        return AssignmentMethod;
    }

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

    bool   mbool_value;
    int    mint_value;
    double mdouble_value;

    ///@}
    ///@name Private Operators
    ///@{

    template< class TVarType, class TDataType >
    void AssignValueToNodes(TVarType& rVariable, const TDataType value)
    {
        const int nnodes = mrModelPart.Nodes().size();

        typedef void (AssignScalarVariableToEntitiesProcess::*AssignmentMethodPointer) (ModelPart::NodeType&, const TVarType&, const TDataType&);

        AssignmentMethodPointer AssignmentMethod = this->GetAssignmentMethod<AssignmentMethodPointer>();

        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh().NodesBegin();

            #pragma omp parallel for
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                (this->*AssignmentMethod)(*it, rVariable, value);
            }
        }
    }

    template< class TVarType, class TDataType >
    void AssignValueToConditions(TVarType& rVariable, const TDataType value)
    {
        const int nconditions = mrModelPart.GetMesh().Conditions().size();

        typedef void (AssignScalarVariableToEntitiesProcess::*AssignmentMethodPointer) (ModelPart::ConditionType&, const TVarType&, const TDataType&);

        AssignmentMethodPointer AssignmentMethod = this->GetAssignmentMethod<AssignmentMethodPointer>();

        if(nconditions != 0)
        {
            ModelPart::ConditionsContainerType::iterator it_begin = mrModelPart.GetMesh().ConditionsBegin();

            #pragma omp parallel for
            for(int i = 0; i<nconditions; i++)
            {
                ModelPart::ConditionsContainerType::iterator it = it_begin + i;

                (this->*AssignmentMethod)(*it, rVariable, value);
            }
        }
    }

    template< class TVarType, class TDataType >
    void AssignValueToElements(TVarType& rVariable, const TDataType value)
    {
        const int nelements = mrModelPart.GetMesh().Elements().size();

        typedef void (AssignScalarVariableToEntitiesProcess::*AssignmentMethodPointer) (ModelPart::ElementType&, const TVarType&, const TDataType&);

        AssignmentMethodPointer AssignmentMethod = this->GetAssignmentMethod<AssignmentMethodPointer>();

        if(nelements != 0)
        {
            ModelPart::ElementsContainerType::iterator it_begin = mrModelPart.GetMesh().ElementsBegin();

            #pragma omp parallel for
            for(int i = 0; i<nelements; i++)
            {
                ModelPart::ElementsContainerType::iterator it = it_begin + i;

                (this->*AssignmentMethod)(*it, rVariable, value);
            }
        }
    }



    ///@}
    ///@name Private Operations
    ///@{
    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    AssignScalarVariableToEntitiesProcess& operator=(AssignScalarVariableToEntitiesProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class AssignScalarVariableToEntitiesProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignScalarVariableToEntitiesProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignScalarVariableToEntitiesProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ASSIGN_SCALAR_VARIABLES_TO_ENTITIES_PROCESS_H_INCLUDED  defined