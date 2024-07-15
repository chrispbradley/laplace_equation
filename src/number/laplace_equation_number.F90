PROGRAM LaplaceEquation

  USE OpenCMISS
  
  IMPLICIT NONE

  !-----------------------------------------------------------------------------------------------------------
  ! PROGRAM VARIABLES AND TYPES
  !-----------------------------------------------------------------------------------------------------------

  !Test program parameters
  REAL(OC_RP), PARAMETER :: HEIGHT=1.0_OC_RP
  REAL(OC_RP), PARAMETER :: WIDTH=2.0_OC_RP
  REAL(OC_RP), PARAMETER :: LENGTH=3.0_OC_RP
 
  INTEGER(OC_Intg), PARAMETER :: CONTEXT_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=2
  INTEGER(OC_Intg), PARAMETER :: REGION_USER_NUMBER=3
  INTEGER(OC_Intg), PARAMETER :: BASIS_USER_NUMBER=4
  INTEGER(OC_Intg), PARAMETER :: GENERATED_MESH_USER_NUMBER=5
  INTEGER(OC_Intg), PARAMETER :: MESH_USER_NUMBER=6
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSITION_USER_NUMBER=7
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSER_USER_NUMBER=8
  INTEGER(OC_Intg), PARAMETER :: GEOMETRIC_FIELD_USER_NUMBER=9
  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER=10
  INTEGER(OC_Intg), PARAMETER :: DEPENDENT_FIELD_USER_NUMBER=11
  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=12
  INTEGER(OC_Intg), PARAMETER :: PROBLEM_USER_NUMBER=13

  !Program types

  !Program variables
  INTEGER(OC_Intg) :: numberOfArguments,argumentLength,status
  INTEGER(OC_Intg) :: numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements, &
    & interpolationType,numberOfGaussXi
  CHARACTER(LEN=255) :: commandArgument,filename

  !CMISS variables

  !Generic CMISS variables
  INTEGER(OC_Intg) :: worldRegionUserNumber
  INTEGER(OC_Intg) :: numberOfComputationalNodes,computationalNodeNumber
  INTEGER(OC_Intg) :: decompositionIndex,equationsSetIndex
  INTEGER(OC_Intg) :: firstNodeNumber,lastNodeNumber
  INTEGER(OC_Intg) :: firstNodeDomain,lastNodeDomain
  INTEGER(OC_Intg) :: err

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM CONTROL PANEL
  !-----------------------------------------------------------------------------------------------------------

  numberOfArguments = COMMAND_ARGUMENT_COUNT()
  IF(numberOfArguments >= 4) THEN
    !If we have enough arguments then use the first four for setting up the problem. The subsequent arguments may be used to
    !pass flags to, say, PETSc.
    CALL GET_COMMAND_ARGUMENT(1,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 1.")
    READ(commandArgument(1:argumentLength),*) numberOfGlobalXElements
    IF(numberOfGlobalXElements<=0) CALL HandleError("Invalid number of X elements.")
    CALL GET_COMMAND_ARGUMENT(2,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 2.")
    READ(commandArgument(1:argumentLength),*) numberOfGlobalYElements
    IF(numberOfGlobalYElements<=0) CALL HandleError("Invalid number of Y elements.")
    CALL GET_COMMAND_ARGUMENT(3,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 3.")
    READ(commandArgument(1:argumentLength),*) numberOfGlobalZElements
    IF(numberOfGlobalZElements<0) CALL HandleError("Invalid number of Z elements.")
    CALL GET_COMMAND_ARGUMENT(4,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 4.")
    READ(commandArgument(1:argumentLength),*) interpolationType
    IF(interpolationType<=0) CALL HandleError("Invalid Interpolation specification.")
  ELSE
    !If there are not enough arguments default the problem specification
    numberOfGlobalXElements=1
    numberOfGlobalYElements=3
    numberOfGlobalZElements=1
    interpolationType=OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION
!    interpolationType=OC_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
!    interpolationType=OC_BASIS_CUBIC_LAGRANGE_INTERPOLATION
  ENDIF

  !Intialise OpenCMISS
  CALL OC_Initialise(err)
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,err)
  !CALL OC_DiagnosticsSetOn(OC_IN_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["Laplace_FiniteElementCalculate"],err)
  WRITE(filename,'(A,"_",I0,"x",I0,"x",I0,"_",I0)') "Laplace",numberOfGlobalXElements,numberOfGlobalYElements, &
    & numberOfGlobalZElements,interpolationType
  CALL OC_OutputSetOn(filename,err)
  !Create a context
  CALL OC_Context_Create(CONTEXT_USER_NUMBER,err)
  CALL OC_Context_WorldRegionGet(CONTEXT_USER_NUMBER,worldRegionUserNumber,err)
  CALL OC_Context_RandomSeedsSet(CONTEXT_USER_NUMBER,9999,err)

  !Get the computational nodes information
  CALL OC_WorkGroup_NumberOfGroupNodesGet(CONTEXT_USER_NUMBER,0,numberOfComputationalNodes,err)
  CALL OC_WorkGroup_GroupNodeNumberGet(CONTEXT_USER_NUMBER,0,computationalNodeNumber,err)

  !-----------------------------------------------------------------------------------------------------------
  ! COORDINATE SYSTEM
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a new RC coordinate system
  CALL OC_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,CONTEXT_USER_NUMBER,err)
  IF(numberOfGlobalZElements==0) THEN
    !Set the coordinate system to be 2D
    CALL OC_CoordinateSystem_DimensionSet(CONTEXT_USER_NUMBER,COORDINATE_SYSTEM_USER_NUMBER,2,err)
  ELSE
    !Set the coordinate system to be 3D
    CALL OC_CoordinateSystem_DimensionSet(CONTEXT_USER_NUMBER,COORDINATE_SYSTEM_USER_NUMBER,3,err)
  ENDIF
  !Finish the creation of the coordinate system
  CALL OC_CoordinateSystem_CreateFinish(CONTEXT_USER_NUMBER,COORDINATE_SYSTEM_USER_NUMBER,err)

  !-----------------------------------------------------------------------------------------------------------
  ! REGION
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the region
  CALL OC_Region_CreateStart(REGION_USER_NUMBER,CONTEXT_USER_NUMBER,worldRegionUserNumber,err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL OC_Region_CoordinateSystemSet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,COORDINATE_SYSTEM_USER_NUMBER,err)
  !Set the region label
  CALL OC_Region_LabelSet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,"LaplaceEquation",err)
  !Finish the creation of the region
  CALL OC_Region_CreateFinish(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,err)

  !-----------------------------------------------------------------------------------------------------------
  ! BASIS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a basis (default is trilinear lagrange)
  CALL OC_Basis_CreateStart(BASIS_USER_NUMBER,CONTEXT_USER_NUMBER,err)
  SELECT CASE(interpolationType)
  CASE(1,2,3,4)
    CALL OC_Basis_TypeSet(CONTEXT_USER_NUMBER,BASIS_USER_NUMBER,OC_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
  CASE(7,8,9)
    CALL OC_Basis_TypeSet(CONTEXT_USER_NUMBER,BASIS_USER_NUMBER,OC_BASIS_SIMPLEX_TYPE,err)
  CASE DEFAULT
    CALL HandleError("Invalid interpolation type.")
  END SELECT
  SELECT CASE(interpolationType)
  CASE(1)
    numberOfGaussXi=2
  CASE(2)
    numberOfGaussXi=3
  CASE(3,4)
    numberOfGaussXi=4
  CASE DEFAULT
    numberOfGaussXi=0 !Don't set number of Gauss points for tri/tet
  END SELECT
  IF(numberOfGlobalZElements==0) THEN
    !Set the basis to be a bi-interpolation basis
    CALL OC_Basis_NumberOfXiSet(CONTEXT_USER_NUMBER,BASIS_USER_NUMBER,2,err)
    CALL OC_Basis_InterpolationXiSet(CONTEXT_USER_NUMBER,BASIS_USER_NUMBER,[interpolationType,interpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL OC_Basis_QuadratureNumberOfGaussXiSet(CONTEXT_USER_NUMBER,BASIS_USER_NUMBER, &
        & [numberOfGaussXi,numberOfGaussXi],err)
    ENDIF
  ELSE
    !Set the basis to be a tri-interpolation basis
    CALL OC_Basis_NumberOfXiSet(CONTEXT_USER_NUMBER,BASIS_USER_NUMBER,3,err)
    CALL OC_Basis_InterpolationXiSet(CONTEXT_USER_NUMBER,BASIS_USER_NUMBER, &
      & [interpolationType,interpolationType,interpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL OC_Basis_QuadratureNumberOfGaussXiSet(CONTEXT_USER_NUMBER,BASIS_USER_NUMBER, &
        & [numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],err)
    ENDIF
  ENDIF
  !Finish the creation of the basis
  CALL OC_Basis_CreateFinish(CONTEXT_USER_NUMBER,BASIS_USER_NUMBER,err)

  !-----------------------------------------------------------------------------------------------------------
  ! MESH
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a generated mesh in the region
  CALL OC_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,CONTEXT_USER_NUMBER,REGION_USER_NUMBER,err)
  !Set up a regular x*y*z mesh
  CALL OC_GeneratedMesh_TypeSet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,GENERATED_MESH_USER_NUMBER, &
    & OC_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  CALL OC_GeneratedMesh_BasisSet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,GENERATED_MESH_USER_NUMBER,BASIS_USER_NUMBER,err)
  !Define the mesh on the region
  IF(numberOfGlobalZElements==0) THEN
    CALL OC_GeneratedMesh_ExtentSet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,GENERATED_MESH_USER_NUMBER,[WIDTH,HEIGHT],err)
    CALL OC_GeneratedMesh_NumberOfElementsSet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,GENERATED_MESH_USER_NUMBER, &
      & [numberOfGlobalXElements,numberOfGlobalYElements],err)
  ELSE
    CALL OC_GeneratedMesh_ExtentSet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,GENERATED_MESH_USER_NUMBER,[WIDTH,HEIGHT,LENGTH],err)
    CALL OC_GeneratedMesh_NumberOfElementsSet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,GENERATED_MESH_USER_NUMBER, &
      & [numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements],err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL OC_GeneratedMesh_CreateFinish(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,GENERATED_MESH_USER_NUMBER,MESH_USER_NUMBER,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DECOMPOSiTION
  !-----------------------------------------------------------------------------------------------------------

  !Create a decomposition
  CALL OC_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,CONTEXT_USER_NUMBER,REGION_USER_NUMBER, &
    & MESH_USER_NUMBER,err)
! Set the decomposition to be a general decomposition with the specified number of domains
  CALL OC_Decomposition_TypeSet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,MESH_USER_NUMBER,DECOMPOSITION_USER_NUMBER, &
    & OC_DECOMPOSITION_CALCULATED_TYPE,err)
  !Finish the decomposition
  CALL OC_Decomposition_CreateFinish(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,MESH_USER_NUMBER,DECOMPOSITION_USER_NUMBER,err)


  !-----------------------------------------------------------------------------------------------------------
  ! DECOMPOSER
  !-----------------------------------------------------------------------------------------------------------
  
  !Decompose
  CALL OC_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,CONTEXT_USER_NUMBER,REGION_USER_NUMBER,0,err)
  !Add in the decomposition
  CALL OC_Decomposer_DecompositionAdd(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,DECOMPOSER_USER_NUMBER,MESH_USER_NUMBER, &
    & DECOMPOSITION_USER_NUMBER,decompositionIndex,err)
  !Finish the decomposer
  CALL OC_Decomposer_CreateFinish(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,DECOMPOSER_USER_NUMBER,err)
  
  !Destory the mesh now that we have decomposed it
! CALL OC_Mesh_Destroy(mesh,err)

  !-----------------------------------------------------------------------------------------------------------
  ! GEOMETRIC FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Start to create a default (geometric) field on the region
  CALL OC_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,CONTEXT_USER_NUMBER,REGION_USER_NUMBER,err)
  !Set the decomposition to use
  CALL OC_Field_DecompositionSet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,GEOMETRIC_FIELD_USER_NUMBER, &
    & MESH_USER_NUMBER,DECOMPOSITION_USER_NUMBER,err)
  !Set the domain to be used by the field components.
  CALL OC_Field_ComponentMeshComponentSet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,GEOMETRIC_FIELD_USER_NUMBER, &
    & OC_FIELD_U_VARIABLE_TYPE,1,1,err)
  CALL OC_Field_ComponentMeshComponentSet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,GEOMETRIC_FIELD_USER_NUMBER, &
    & OC_FIELD_U_VARIABLE_TYPE,2,1,err)
  IF(numberOfGlobalZElements/=0) THEN
    CALL OC_Field_ComponentMeshComponentSet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,GEOMETRIC_FIELD_USER_NUMBER, &
      & OC_FIELD_U_VARIABLE_TYPE,3,1,err)
  ENDIF
  !Finish creating the field
  CALL OC_Field_CreateFinish(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,GEOMETRIC_FIELD_USER_NUMBER,err)

  !Update the geometric field parameters
  CALL OC_GeneratedMesh_GeometricParametersCalculate(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,GENERATED_MESH_USER_NUMBER, &
    & GEOMETRIC_FIELD_USER_NUMBER,err)

  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS SETS
  !-----------------------------------------------------------------------------------------------------------

  !Create the Standard Laplace equations set
  CALL OC_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,CONTEXT_USER_NUMBER,REGION_USER_NUMBER,GEOMETRIC_FIELD_USER_NUMBER, &
    & [OC_EQUATIONS_SET_CLASSICAL_FIELD_CLASS,OC_EQUATIONS_SET_LAPLACE_EQUATION_TYPE, &
    & OC_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE],EQUATIONS_SET_FIELD_USER_NUMBER,err)
  !Finish creating the equations set
  CALL OC_EquationsSet_CreateFinish(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DEPENDENT FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set dependent field variables
  CALL OC_EquationsSet_DependentCreateStart(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER, &
    & DEPENDENT_FIELD_USER_NUMBER,err)
  !Set the DOFs to be contiguous across components
  CALL OC_Field_DOFOrderTypeSet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,DEPENDENT_FIELD_USER_NUMBER, &
    & OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_SEPARATED_COMPONENT_DOF_ORDER,err)
  CALL OC_Field_DOFOrderTypeSet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,DEPENDENT_FIELD_USER_NUMBER, &
    & OC_FIELD_DELUDELN_VARIABLE_TYPE,OC_FIELD_SEPARATED_COMPONENT_DOF_ORDER,err)
  !Finish the equations set dependent field variables
  CALL OC_EquationsSet_DependentCreateFinish(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER,err)

  !Initialise the field with an initial guess
  CALL OC_Field_ComponentValuesInitialise(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,DEPENDENT_FIELD_USER_NUMBER, &
    & OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE,1,0.5_OC_RP,err)

  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set equations
  CALL OC_EquationsSet_EquationsCreateStart(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER,err)
  !Set the equations matrices sparsity type
  CALL OC_Equations_SparsityTypeSet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER, &
    & OC_EQUATIONS_SPARSE_MATRICES,err)
  !CALL OC_Equations_SparsityTypeSet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER, &
  !  & OC_EQUATIONS_FULL_MATRICES,err)
  !Set the equations set output
  CALL OC_Equations_OutputTypeSet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER, &
    & OC_EQUATIONS_NO_OUTPUT,err)
  !CALL OC_Equations_OutputTypeSet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER, &
  !  & OC_EQUATIONS_TIMING_OUTPUT,err)
  !CALL OC_Equations_OutputTypeSet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER, &
  !  & OC_EQUATIONS_MATRIX_OUTPUT,err)
  !CALL OC_Equations_OutputTypeSet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER, &
  !  & OC_EQUATIONS_ELEMENT_MATRIX_OUTPUT,err)
  !Finish the equations set equations
  CALL OC_EquationsSet_EquationsCreateFinish(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER,err)

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a problem.
  CALL OC_Problem_CreateStart(PROBLEM_USER_NUMBER,CONTEXT_USER_NUMBER,[OC_PROBLEM_CLASSICAL_FIELD_CLASS, &
    & OC_PROBLEM_LAPLACE_EQUATION_TYPE,OC_PROBLEM_STANDARD_LAPLACE_SUBTYPE],err)
  !Finish the creation of a problem.
  CALL OC_Problem_CreateFinish(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,err)

  !-----------------------------------------------------------------------------------------------------------
  ! CONTROL LOOPS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the problem control loop
  CALL OC_Problem_ControlLoopCreateStart(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,err)
  !Finish creating the problem control loop
  CALL OC_Problem_ControlLoopCreateFinish(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVER
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the problem solvers
  CALL OC_Problem_SolversCreateStart(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,err)
  !CALL OC_Solver_OutputTypeSet(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,OC_CONTROL_LOOP_NODE,1, &
  !  & OC_SOLVER_NO_OUTPUT,err)
  !CALL OC_Solver_OutputTypeSet(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,OC_CONTROL_LOOP_NODE,1, &
  !  & OC_SOLVER_PROGRESS_OUTPUT,err)
  !CALL OC_Solver_OutputTypeSet(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,OC_CONTROL_LOOP_NODE,1, &
  !  & OC_SOLVER_TIMING_OUTPUT,err)
  !CALL OC_Solver_OutputTypeSet(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,OC_CONTROL_LOOP_NODE,1, &
  !  & OC_SOLVER_SOLVER_OUTPUT,err)
  CALL OC_Solver_OutputTypeSet(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,OC_CONTROL_LOOP_NODE,1, &
    & OC_SOLVER_MATRIX_OUTPUT,err)
  
  CALL OC_Solver_LinearTypeSet(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,OC_CONTROL_LOOP_NODE,1, &
    & OC_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,err)
  CALL OC_Solver_LinearIterativeAbsoluteToleranceSet(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,OC_CONTROL_LOOP_NODE,1, &
    & 1.0E-12_OC_RP,err)
  CALL OC_Solver_LinearIterativeRelativeToleranceSet(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,OC_CONTROL_LOOP_NODE,1, &
    & 1.0E-12_OC_RP,err)
  
  !CALL OC_Solver_LinearTypeSet(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,OC_CONTROL_LOOP_NODE,1, &
  !  & OC_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
  !CALL OC_Solver_LibraryTypeSet(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,OC_CONTROL_LOOP_NODE,1, &
  !  & OC_SOLVER_MUMPS_LIBRARY,err)
  !CALL OC_Solver_LibraryTypeSet(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,OC_CONTROL_LOOP_NODE,1, &
  !  & OC_SOLVER_LAPACK_LIBRARY,err)
  !CALL OC_Solver_LibraryTypeSet(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,OC_CONTROL_LOOP_NODE,1, &
  !  & OC_SOLVER_SUPERLU_LIBRARY,err)
  !CALL OC_Solver_LibraryTypeSet(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,OC_CONTROL_LOOP_NODE,1, &
  !  & OC_SOLVER_PASTIX_LIBRARY,err)
  !Finish the creation of the problem solver
  CALL OC_Problem_SolversCreateFinish(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVER EQUATIONS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the problem solver equations
  CALL OC_Problem_SolverEquationsCreateStart(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,err)
  !Set the solver equations sparsity
  CALL OC_SolverEquations_SparsityTypeSet(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,OC_CONTROL_LOOP_NODE,1, &
    & OC_SOLVER_SPARSE_MATRICES,err)
  !CALL OC_SolverEquations_SparsityTypeSet(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,OC_CONTROL_LOOP_NODE,1, &
  !  & OC_SOLVER_FULL_MATRICES,err)
  !Add in the equations set
  CALL OC_SolverEquations_EquationsSetAdd(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,OC_CONTROL_LOOP_NODE,1, &
    & REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER,equationsSetIndex,err)
  !Finish the creation of the problem solver equations
  CALL OC_Problem_SolverEquationsCreateFinish(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,err)

  !-----------------------------------------------------------------------------------------------------------
  ! BOUNDARY CONDITIONS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the equations set boundary conditions
  CALL OC_SolverEquations_BoundaryConditionsCreateStart(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,OC_CONTROL_LOOP_NODE,1,err)
  !Set the first node to 0.0 and the last node to 1.0
  firstNodeNumber=1
  CALL OC_Nodes_NumberOfNodesGet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,lastNodeNumber,err)
  CALL OC_Decomposition_NodeDomainGet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,MESH_USER_NUMBER,DECOMPOSITION_USER_NUMBER, &
    & firstNodeNumber,1,firstNodeDomain,err)
  CALL OC_Decomposition_NodeDomainGet(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,MESH_USER_NUMBER,DECOMPOSITION_USER_NUMBER, &
    & lastNodeNumber,1,lastNodeDomain,err)
  IF(firstNodeDomain==computationalNodeNumber) THEN
    CALL OC_BoundaryConditions_SetNode(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,OC_CONTROL_LOOP_NODE,1, &
      & REGION_USER_NUMBER,DEPENDENT_FIELD_USER_NUMBER,OC_FIELD_U_VARIABLE_TYPE,1,1,firstNodeNumber,1, &
      & OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
  ENDIF
  IF(lastNodeDomain==computationalNodeNumber) THEN
    CALL OC_BoundaryConditions_SetNode(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,OC_CONTROL_LOOP_NODE,1, &
      & REGION_USER_NUMBER,DEPENDENT_FIELD_USER_NUMBER,OC_FIELD_U_VARIABLE_TYPE,1,1,lastNodeNumber,1, &
      & OC_BOUNDARY_CONDITION_FIXED,1.0_OC_RP,err)
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL OC_SolverEquations_BoundaryConditionsCreateFinish(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,OC_CONTROL_LOOP_NODE,1,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVE
  !-----------------------------------------------------------------------------------------------------------

  !Solve the problem
  CALL OC_Problem_Solve(CONTEXT_USER_NUMBER,PROBLEM_USER_NUMBER,err)

  !-----------------------------------------------------------------------------------------------------------
  ! OUTPUT
  !-----------------------------------------------------------------------------------------------------------

  !Export results
  CALL OC_Fields_NodesExport(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,"LaplaceEquation","FORTRAN",err)
  CALL OC_Fields_ElementsExport(CONTEXT_USER_NUMBER,REGION_USER_NUMBER,"LaplaceEquation","FORTRAN",err)

  !Destroy the context
  CALL OC_Context_Destroy(CONTEXT_USER_NUMBER,err)
  !Finialise OpenCMISS
  CALL OC_Finalise(err)
  
  WRITE(*,'("Program successfully completed.")')
  STOP

CONTAINS

  SUBROUTINE HandleError(errorString)
    CHARACTER(LEN=*), INTENT(IN) :: errorString
    WRITE(*,'(">>ERROR: ",A)') errorString(1:LEN_TRIM(errorString))
    STOP
  END SUBROUTINE HandleError

END PROGRAM LaplaceEquation
