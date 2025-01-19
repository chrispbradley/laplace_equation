/*
 * 
 * This is an example program to solve Laplace's equation using OpenCMISS calls from C.
 * by Chris Bradley
 *
 */
#include <stdlib.h>
#include <stdio.h>

#include "opencmiss/opencmiss.h"

#define STRING_SIZE 255

#define HEIGHT 1.0
#define WIDTH 2.0
#define LENGTH 3.0

#define CONTEXT_USER_NUMBER 1
#define COORDINATE_SYSTEM_USER_NUMBER 2
#define REGION_USER_NUMBER 3
#define BASIS_USER_NUMBER 4
#define GENERATED_MESH_USER_NUMBER 5
#define MESH_USER_NUMBER 6
#define DECOMPOSITION_USER_NUMBER 7
#define DECOMPOSER_USER_NUMBER 8
#define GEOMETRIC_FIELD_USER_NUMBER 9
#define EQUATIONS_SET_USER_NUMBER 10
#define EQUATIONS_SET_FIELD_USER_NUMBER 11
#define DEPENDENT_FIELD_USER_NUMBER 12
#define PROBLEM_USER_NUMBER 13

#define MAX_COORDINATES 3

int main(int argc, char *argv[])
{
  OC_BasisType basis = (OC_BasisType)NULL;
  OC_BoundaryConditionsType boundaryConditions=(OC_BoundaryConditionsType)NULL;
  OC_ComputationEnvironmentType computationEnvironment=(OC_ComputationEnvironmentType)NULL;
  OC_ContextType context=(OC_ContextType)NULL;
  OC_CoordinateSystemType coordinateSystem=(OC_CoordinateSystemType)NULL;
  OC_DecompositionType decomposition=(OC_DecompositionType)NULL;
  OC_DecomposerType decomposer=(OC_DecomposerType)NULL;
  OC_EquationsType equations=(OC_EquationsType)NULL;
  OC_EquationsSetType equationsSet=(OC_EquationsSetType)NULL;
  OC_FieldsType fields=(OC_FieldsType)NULL;
  OC_FieldType geometricField=(OC_FieldType)NULL,dependentField=(OC_FieldType)NULL,equationsSetField=(OC_FieldType)NULL;
  OC_GeneratedMeshType generatedMesh=(OC_GeneratedMeshType)NULL;
  OC_MeshType mesh=(OC_MeshType)NULL;
  OC_ProblemType problem=(OC_ProblemType)NULL;
  OC_RegionType region=(OC_RegionType)NULL,worldRegion=(OC_RegionType)NULL;
  OC_SolverType solver=(OC_SolverType)NULL;
  OC_SolverEquationsType solverEquations=(OC_SolverEquationsType)NULL;
  OC_WorkGroupType worldWorkGroup=(OC_WorkGroupType)NULL;

  int numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements,interpolationType;
  
  int numberOfComputationalNodes,computationalNodeNumber;
  int decompositionIndex,equationsSetIndex;
  int firstNodeNumber,lastNodeNumber;
  int firstNodeDomain,lastNodeDomain;

  int basisInterpolation[MAX_COORDINATES];
  int numberOfDimensions;
  int numberOfGauss[MAX_COORDINATES];
  int numberOfGaussXi;
  int numberXiElements[MAX_COORDINATES];
  int controlLoopIdentifier[1];
  double meshExtent[MAX_COORDINATES];

  int equationsSetSpecification[3];
  int problemSpecification[3];

  char filename[STRING_SIZE],format[STRING_SIZE],regionName[STRING_SIZE];

  int err;

  controlLoopIdentifier[0]=OC_CONTROL_LOOP_NODE;

  if(argc >= 4)
    {
      numberOfGlobalXElements = atoi(argv[1]);
      numberOfGlobalYElements = atoi(argv[2]);
      numberOfGlobalZElements = atoi(argv[3]);
      interpolationType = atoi(argv[4]);
    }
  else
    {
      numberOfGlobalXElements = 1;
      numberOfGlobalYElements = 3;
      numberOfGlobalZElements = 1;
      interpolationType = OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION;
   }

  err = OC_Initialise();
  OPENCMISS_CHECK_ERROR(err,"Initialising OpenCMISS");
  err = OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR);
  sprintf(filename,"%s_%dx%dx%d_%d","Laplace",numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements,interpolationType);
  err = OC_OutputSetOn(STRING_SIZE,filename);

  err = OC_Context_Initialise(&context);
  OPENCMISS_CHECK_ERROR(err,"Initialising context");
  err = OC_Context_Create(CONTEXT_USER_NUMBER,context);
  OPENCMISS_CHECK_ERROR(err,"Creating context");

  err = OC_Region_Initialise(&worldRegion);
  err = OC_Context_WorldRegionGet(context,worldRegion);
  

  /* Get the computational nodes information */
  err = OC_ComputationEnvironment_Initialise(&computationEnvironment);
  err = OC_Context_ComputationEnvironmentGet(context,computationEnvironment);

  err = OC_WorkGroup_Initialise(&worldWorkGroup);
  err = OC_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup);
  err = OC_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,&numberOfComputationalNodes);
  err = OC_WorkGroup_GroupNodeNumberGet(worldWorkGroup,&computationalNodeNumber);

  /* Start the creation of a new RC coordinate system */
  err = OC_CoordinateSystem_Initialise(&coordinateSystem);
  err = OC_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem);
  if(numberOfGlobalZElements == 0)
    {
      /* Set the coordinate system to be 2D */
      numberOfDimensions = 2;
    }
  else
    {
      /* Set the coordinate system to be 3D */
      numberOfDimensions = 3;
    }
  err = OC_CoordinateSystem_DimensionSet(coordinateSystem,numberOfDimensions);
  /* Finish the creation of the coordinate system */
  err = OC_CoordinateSystem_CreateFinish(coordinateSystem);

  /* Start the creation of the region */
  sprintf(regionName,"%s","LaplaceEquation");
  err = OC_Region_Initialise(&region);
  err = OC_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region);
  /* Set the regions coordinate system to the 2D RC coordinate system that we have created */
  err = OC_Region_CoordinateSystemSet(region,coordinateSystem);
  /* Set the regions name */
  err = OC_Region_LabelSet(region,STRING_SIZE,regionName);
  /* Finish the creation of the region */
  err = OC_Region_CreateFinish(region);

  /* Start the creation of a basis (default is trilinear lagrange) */
  err = OC_Basis_Initialise(&basis);
  err = OC_Basis_CreateStart(BASIS_USER_NUMBER,context,basis);
  switch(interpolationType)
    {
    case 1:
    case 2:
    case 3:
    case 4:
      err = OC_Basis_TypeSet(basis,OC_BASIS_LAGRANGE_HERMITE_TP_TYPE);
      break;
    case 7:
    case 8:
    case 9:
      err = OC_Basis_TypeSet(basis,OC_BASIS_SIMPLEX_TYPE);
      break;
    default:
      OPENCMISS_CHECK_ERROR(OC_FORCE_ERROR,"Invalid interpolation type");
    }
  switch(interpolationType)
    {
    case 1:
      numberOfGaussXi = 2;
      break;
    case 2:
      numberOfGaussXi = 3;
      break;
    case 3:
    case 4:
      numberOfGaussXi = 4;
      break;
    default:
      numberOfGaussXi = 0;
    }
  basisInterpolation[0] = interpolationType;
  basisInterpolation[1] = interpolationType;
  numberOfGauss[0] = numberOfGaussXi;
  numberOfGauss[1] = numberOfGaussXi;
  if(numberOfGlobalZElements != 0)
    {
      basisInterpolation[2] = interpolationType;
      numberOfGauss[2] = numberOfGaussXi;
    }
  err = OC_Basis_NumberOfXiSet(basis,numberOfDimensions);
  err = OC_Basis_InterpolationXiSet(basis,numberOfDimensions,basisInterpolation);
  err = OC_Basis_QuadratureNumberOfGaussXiSet(basis,numberOfDimensions,numberOfGauss);
  /* Finish the creation of the basis */
  err = OC_Basis_CreateFinish(basis);

  /* Start the creation of a generated mesh in the region */
  err = OC_GeneratedMesh_Initialise(&generatedMesh);
  err = OC_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh);
  /* Set up a regular x*y*z mesh */
  err = OC_GeneratedMesh_TypeSet(generatedMesh,OC_GENERATED_MESH_REGULAR_MESH_TYPE);
  /* Set the default basis */
  err = OC_GeneratedMesh_BasisSet(generatedMesh,1,&basis);
  /* Define the mesh on the region */
  meshExtent[0] = WIDTH;
  meshExtent[1] = HEIGHT;
  numberXiElements[0] = numberOfGlobalXElements;
  numberXiElements[1] = numberOfGlobalYElements;
  if(numberOfGlobalZElements != 0)
    {
      meshExtent[2] = LENGTH;
      numberXiElements[2] = numberOfGlobalZElements;
    }
  err = OC_GeneratedMesh_ExtentSet(generatedMesh,MAX_COORDINATES,meshExtent);
  err = OC_GeneratedMesh_NumberOfElementsSet(generatedMesh,MAX_COORDINATES,numberXiElements);
  /* Finish the creation of a generated mesh in the region */
  err = OC_Mesh_Initialise(&mesh);
  /* Finish the creation of a generated mesh in the region */
  err = OC_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh);

  /* Create a decomposition */
  err = OC_Decomposition_Initialise(&decomposition);
  err = OC_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition);
  /* Finish the decomposition */
  err = OC_Decomposition_CreateFinish(decomposition);

  /* Decompose */
  err = OC_Decomposer_Initialise(&decomposer);
  err = OC_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer);
  /* Add in the decomposition */
  err = OC_Decomposer_DecompositionAdd(decomposer,decomposition,&decompositionIndex);
  /* Finish the decomposer */
  err = OC_Decomposer_CreateFinish(decomposer);

  /* Start to create a default (geometric) field on the region */
  err = OC_Field_Initialise(&geometricField);
  err = OC_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region,geometricField);
  /* Set the decomposition to use */
  err = OC_Field_DecompositionSet(geometricField,decomposition);
  /* Set the domain to be used by the field components. */
  err = OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,1,1);
  err = OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,2,1);
  if(numberOfGlobalZElements != 0)
    {
      err = OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,3,1);
    }
  /* Finish creating the field */
  err = OC_Field_CreateFinish(geometricField);

  /* Update the geometric field parameters */
  err = OC_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField);

  /* Create the equations_set */
  err = OC_EquationsSet_Initialise(&equationsSet);
  err = OC_Field_Initialise(&equationsSetField);
  equationsSetSpecification[0] = OC_EQUATIONS_SET_CLASSICAL_FIELD_CLASS;
  equationsSetSpecification[1] = OC_EQUATIONS_SET_LAPLACE_EQUATION_TYPE;
  equationsSetSpecification[2] = OC_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE;
  err = OC_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,geometricField, \
    3,equationsSetSpecification,EQUATIONS_SET_FIELD_USER_NUMBER, \
    equationsSetField,equationsSet);
  /* Finish creating the equations set */
  err = OC_EquationsSet_CreateFinish(equationsSet);

  /* Create the equations set dependent field variables */
  err = OC_Field_Initialise(&dependentField);
  err = OC_EquationsSet_DependentCreateStart(equationsSet,DEPENDENT_FIELD_USER_NUMBER,dependentField);
  /* Finish the equations set dependent field variables */
  err = OC_EquationsSet_DependentCreateFinish(equationsSet);

  /* Create the equations set equations */
  err = OC_Equations_Initialise(&equations);
  err = OC_EquationsSet_EquationsCreateStart(equationsSet,equations);
  /* Set the equations matrices sparsity type */
  err = OC_Equations_SparsityTypeSet(equations,OC_EQUATIONS_SPARSE_MATRICES);
  /* Set the equations set output */
  /* err = OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_NO_OUTPUT); */
  err = OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_TIMING_OUTPUT);
  /* err = OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_MATRIX_OUTPUT); */
  /* err = OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_ELEMENT_MATRIX_OUTPUT); */
  /* Finish the equations set equations */
  err = OC_EquationsSet_EquationsCreateFinish(equationsSet);

  /* Start the creation of a problem, setting the problem to be a standard Laplace problem. */
  err = OC_Problem_Initialise(&problem);
  problemSpecification[0] = OC_PROBLEM_CLASSICAL_FIELD_CLASS;
  problemSpecification[1] = OC_PROBLEM_LAPLACE_EQUATION_TYPE;
  problemSpecification[2] = OC_PROBLEM_STANDARD_LAPLACE_SUBTYPE;
  err = OC_Problem_CreateStart(PROBLEM_USER_NUMBER,context,3,problemSpecification,problem);
  /* Finish the creation of a problem. */
  err = OC_Problem_CreateFinish(problem);

  /* Start the creation of the problem control loop */
  err = OC_Problem_ControlLoopCreateStart(problem);
  /* Finish creating the problem control loop */
  err = OC_Problem_ControlLoopCreateFinish(problem);

  /* Start the creation of the problem solvers */
  err = OC_Solver_Initialise(&solver);
  err = OC_Problem_SolversCreateStart(problem);
  err = OC_Problem_SolverGet(problem,1,controlLoopIdentifier,1,solver);
  /* err = OC_Solver_OutputTypeSet(solver,OC_SOLVER_NO_OUTPUT); */
  /* err = OC_Solver_OutputTypeSet(solver,OC_SOLVER_PROGRESS_OUTPUT); */
  /* err = OC_Solver_OutputTypeSet(solver,OC_SOLVER_TIMING_OUTPUT); */
  /* err = OC_Solver_OutputTypeSet(solver,OC_SOLVER_SOLVER_OUTPUT); */
  err = OC_Solver_OutputTypeSet(solver,OC_SOLVER_MATRIX_OUTPUT);
  err = OC_Solver_LinearTypeSet(solver,OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE);
  err = OC_Solver_LibraryTypeSet(solver,OC_SOLVER_MUMPS_LIBRARY);
  /* Finish the creation of the problem solver */
  err = OC_Problem_SolversCreateFinish(problem);

  /* Start the creation of the problem solver equations */
  solver=(OC_SolverType)NULL;
  err = OC_Solver_Initialise(&solver);
  err = OC_SolverEquations_Initialise(&solverEquations);
  err = OC_Problem_SolverEquationsCreateStart(problem);
  /* Get the solve equations */
  err = OC_Problem_SolverGet(problem,1,controlLoopIdentifier,1,solver);
  err = OC_Solver_SolverEquationsGet(solver,solverEquations);
  /* Set the solver equations sparsity */
  err = OC_SolverEquations_SparsityTypeSet(solverEquations,OC_SOLVER_SPARSE_MATRICES);
  /* err = OC_SolverEquations_SparsityTypeSet(solverEquations,OC_SOLVER_FULL_MATRICES);  */
  /* Add in the equations set */
  err = OC_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,&equationsSetIndex);
  /* Finish the creation of the problem solver equations */
  err = OC_Problem_SolverEquationsCreateFinish(problem);

  /* Start the creation of the equations set boundary conditions */
  err = OC_BoundaryConditions_Initialise(&boundaryConditions);
  err = OC_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions);
  /* Set the first node to 0.0 and the last node to 1.0 */
  firstNodeNumber = 1;
  if(numberOfGlobalZElements == 0)
    {
      lastNodeNumber = (numberOfGlobalXElements+1)*(numberOfGlobalYElements+1);
    }
  else
    {
      lastNodeNumber = (numberOfGlobalXElements+1)*(numberOfGlobalYElements+1)*(numberOfGlobalZElements+1);
    }
  err = OC_Decomposition_NodeDomainGet(decomposition,1,firstNodeNumber,&firstNodeDomain);
  err = OC_Decomposition_NodeDomainGet(decomposition,1,lastNodeNumber,&lastNodeDomain);
  if(firstNodeDomain == computationalNodeNumber)
    {
      err = OC_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,firstNodeNumber,1, \
        OC_BOUNDARY_CONDITION_FIXED,0.0);
    }
  if(lastNodeDomain == computationalNodeNumber)
    {
      err = OC_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,lastNodeNumber,1, \
        OC_BOUNDARY_CONDITION_FIXED,1.0);
    }
  /* Finish the creation of the equations set boundary conditions */
  err = OC_SolverEquations_BoundaryConditionsCreateFinish(solverEquations);

  /* Solve the problem */
  err = OC_Problem_Solve(problem);

  /* Output results */
  sprintf(filename,"%s","LaplaceEquation");
  sprintf(format,"%s","FORTRAN");
  err = OC_Fields_Initialise(&fields);
  err = OC_Fields_CreateRegion(region,fields);
  err = OC_Fields_NodesExport(fields,STRING_SIZE,filename,STRING_SIZE,format);
  err = OC_Fields_ElementsExport(fields,STRING_SIZE,filename,STRING_SIZE,format);
  err = OC_Fields_Finalise(&fields);

  /* Destroy the context */
  err = OC_Context_Destroy(context);
  /* Finalise OpenCMISS */
  err = OC_Finalise();

  return err;
}
