//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined(KRATOS_QHULL_TRIANGULAR_MESH_2D_MESHER_H_INCLUDED )
#define  KRATOS_QHULL_TRIANGULAR_MESH_2D_MESHER_H_INCLUDED


#if !defined(KRATOS_QHULL_EXTERNAL_H_INCLUDED)
#define  KRATOS_QHULL_EXTERNAL_H_INCLUDED
#include "libqhull/libqhull.h"
#include "libqhull/qset.h"
#endif

// System includes

// Project includes
#include "geometries/triangle_2d_3.h"
#include "custom_meshers/mesher.hpp"

#if __cplusplus
extern "C" {
  int isatty(int);
}

#elif _MSC_VER
#include <io.h>
#define isatty _isatty
#else
int isatty(int);
#endif


///VARIABLES used:
//Data:
//StepData:
//Flags:    (checked)
//          (set)
//          (modified)
//          (reset)


namespace Kratos
{

    // extern "C" {
    //   void runQhull3D(const std::vector<vec3> &points, const char* args);
    //   void runQhull(const PointCoordinates &points, const char *qhullCommand2);
    // }


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
class KRATOS_API(DELAUNAY_MESHING_APPLICATION) QhullTriangularMesh2DMesher
  : public Mesher
{
protected:

    enum class QhullError {INPUT_MEMORY_ERROR=1, INTERNAL_ERROR=2, INVALID_GEOMETRY_ERROR=3};

public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of TriGenCDT
    KRATOS_CLASS_POINTER_DEFINITION( QhullTriangularMesh2DMesher );

    typedef MesherUtilities::MeshingInfoParameters              InfoParametersType;
    typedef MesherUtilities::MeshingParameters               MeshingParametersType;
    typedef MesherUtilities::RefiningParameters               RefineParametersType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    QhullTriangularMesh2DMesher(): Mesher() {

      QHULL_LIB_CHECK /* Check for compatible library */

      int argc;
      char *argv[] = {" arguments "};
      char hidden_options[] = " d n ";
      int curlong, totlong; /* used !qh_NOmem */
      int exitcode, numpoints, dim;
      coordT *points;
      boolT ismalloc;

      qh_init_A(stdin, stdout, stderr, argc, argv);

      exitcode= setjmp(qh errexit); /* simple statement for CRAY J916 */


      if (!exitcode) {
        qh NOerrexit = False;
        qh_option("delaunay  Qbbound-last", NULL, NULL);
        qh DELAUNAY= True;     /* 'd'   */
        qh SCALElast= True;    /* 'Qbb' */
        qh KEEPcoplanar= True; /* 'Qc', to keep coplanars in 'p' */
        qh_checkflags(qh qhull_command, hidden_options);
        qh_initflags(qh qhull_command);
        points= qh_readpoints(&numpoints, &dim, &ismalloc);
        if (dim >= 5) {
          qh_option("Qxact_merge", NULL, NULL);
          qh MERGEexact= True; /* 'Qx' always */
        }
        qh_init_B(points, numpoints, dim, ismalloc);
        qh_qhull();
        qh_check_output();
        qh_produce_output();
        if (qh VERIFYoutput && !qh FORCEoutput && !qh STOPpoint && !qh STOPcone)
          qh_check_points();
        exitcode= qh_ERRnone;
      }
    } //

    /// Copy constructor.
    QhullTriangularMesh2DMesher(QhullTriangularMesh2DMesher const& rOther): Mesher(rOther) {}

    /// Destructor.
    virtual ~QhullTriangularMesh2DMesher() {}

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

    /// Turn back information as a string.
    std::string Info() const override
    {
	return "";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override{}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override{}


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

    /// Assignment operator.
   QhullTriangularMesh2DMesher& operator=(QhullTriangularMesh2DMesher const& rOther);

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
    ///@name Unaccessible methods
    ///@{

    ///@}

}; // Class QhullTriangularMesh2DMesher

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
    inline std::istream& operator >> (std::istream& rIStream,
				      QhullTriangularMesh2DMesher& rThis);

/// output stream function
    inline std::ostream& operator << (std::ostream& rOStream,
				      const QhullTriangularMesh2DMesher& rThis)
    {
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
    }
///@}


}  // namespace Kratos.

#endif // KRATOS_QHULL_TRIANGULAR_MESH_2D_MESHER_H_INCLUDED  defined
