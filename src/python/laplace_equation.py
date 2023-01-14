#!/usr/bin/env python

import sys

# Intialise OpenCMISS-Iron
from opencmiss.iron import iron

#-----------------------------------------------------------------------------------------------------------
# SET PROBLEM PARAMETERS
#-----------------------------------------------------------------------------------------------------------

HEIGHT = 1.0
WIDTH = 1.0
LENGTH = 1.0

(CONTEXT_USER_NUMBER,
 COORDINATE_SYSTEM_USER_NUMBER,
 REGION_USER_NUMBER,
 BASIS_USER_NUMBER,
 GENERATED_MESH_USER_NUMBER,
 MESH_USER_NUMBER,
 DECOMPOSITION_USER_NUMBER,
 DECOMPOSER_USER_NUMBER,
 GEOMETRIC_FIELD_USER_NUMBER,
 EQUATIONS_SET_FIELD_USER_NUMBER,
 DEPENDENT_FIELD_USER_NUMBER,
 EQUATIONS_SET_USER_NUMBER,
 PROBLEM_USER_NUMBER) = range(1,14)

NUMBER_OF_GAUSS_XI = 2

numberOfGlobalXElements = 1
numberOfGlobalYElements = 3
numberOfGlobalZElements = 1

# Override with command line arguments if need be
if len(sys.argv) > 1:
    if len(sys.argv) > 4:
        sys.exit('ERROR: too many arguments- currently only accepting up to 3 options: numberOfGlobalXElements numberOfGlobalYElements numberOfGlobalZElements')
        numberOfGlobalXElementsnumberOfGlobalXElements = int(sys.argv[1])
    if len(sys.argv) > 2:
        numberOfGlobalYElements = int(sys.argv[2])
    if len(sys.argv) > 3:
        numberOfGlobalZElements = int(sys.argv[3])

if (numberOfGlobalZElements >= 0):
    if (numberOfGlobalYElements >= 0):
        if (numberOfGlobalXElements >= 0):
            if (numberOfGlobalZElements == 0):
                if(numberOfGlobalYElements == 0):
                    numberOfDimensions = 1
                else:
                    numberOfDimensions = 2
            else:
                numberOfDimensions = 3
        else:
            sys.exit('ERROR: number of global X elements must be greater than 0.')
    else:
        sys.exit('ERROR: number of global Y elements must be greater than 0.')
else:
    sys.exit('ERROR: number of global Z elements must be greater than 0.')

#-----------------------------------------------------------------------------------------------------------
# DIAGNOSTICS AND COMPUTATIONAL NODE INFORMATION
#-----------------------------------------------------------------------------------------------------------

# Create a context for the example
context = iron.Context()
context.Create(CONTEXT_USER_NUMBER)

# Get the world region
worldRegion = iron.Region()
context.WorldRegionGet(worldRegion)

#iron.DiagnosticsSetOn(iron.DiagnosticTypes.IN,[1,2,3,4,5],"Diagnostics",["Laplace_FiniteElementCalculate"])

# Get the computational nodes information
computationEnvironment = iron.ComputationEnvironment()
context.ComputationEnvironmentGet(computationEnvironment)

worldWorkGroup = iron.WorkGroup()
computationEnvironment.WorldWorkGroupGet(worldWorkGroup)
numberOfComputationalNodes = worldWorkGroup.NumberOfGroupNodesGet()
computationalNodeNumber = worldWorkGroup.GroupNodeNumberGet()

#-----------------------------------------------------------------------------------------------------------
#COORDINATE SYSTEM
#-----------------------------------------------------------------------------------------------------------

coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context)
coordinateSystem.DimensionSet(numberOfDimensions)
coordinateSystem.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
#REGION
#-----------------------------------------------------------------------------------------------------------
region = iron.Region()
region.CreateStart(REGION_USER_NUMBER,worldRegion)
region.LabelSet("LaplaceEquation")
region.CoordinateSystemSet(coordinateSystem)
region.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
#BASIS
#-----------------------------------------------------------------------------------------------------------

basis = iron.Basis()
basis.CreateStart(BASIS_USER_NUMBER,context)
basis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
basis.NumberOfXiSet(numberOfDimensions)
basis.InterpolationXiSet([iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfDimensions)
basis.QuadratureNumberOfGaussXiSet([NUMBER_OF_GAUSS_XI]*numberOfDimensions)
basis.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
#MESH
#-----------------------------------------------------------------------------------------------------------
generatedMesh = iron.GeneratedMesh()
generatedMesh.CreateStart(GENERATED_MESH_USER_NUMBER,region)
generatedMesh.TypeSet(iron.GeneratedMeshTypes.REGULAR)
generatedMesh.BasisSet([basis])
if (numberOfDimensions == 1):
    generatedMesh.ExtentSet([WIDTH])
    generatedMesh.NumberOfElementsSet([numberOfGlobalXElements])
elif (numberOfDimensions == 2):
    generatedMesh.ExtentSet([WIDTH,HEIGHT])
    generatedMesh.NumberOfElementsSet([numberOfGlobalXElements,numberOfGlobalYElements])
elif (numberOfDimensions == 3):
    generatedMesh.ExtentSet([WIDTH,HEIGHT,LENGTH])
    generatedMesh.NumberOfElementsSet([numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements])
else:
    sys.exit('ERROR: invalid number of dimensions.')
mesh = iron.Mesh()
generatedMesh.CreateFinish(MESH_USER_NUMBER,mesh)

#-----------------------------------------------------------------------------------------------------------
#MESH DECOMPOSITION
#-----------------------------------------------------------------------------------------------------------

decomposition = iron.Decomposition()
decomposition.CreateStart(DECOMPOSITION_USER_NUMBER,mesh)
decomposition.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
#DECOMPOSER
#-----------------------------------------------------------------------------------------------------------

decomposer = iron.Decomposer()
decomposer.CreateStart(DECOMPOSER_USER_NUMBER,worldRegion,worldWorkGroup)
decompositionIndex = decomposer.DecompositionAdd(decomposition)
decomposer.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
#GEOMETRIC FIELD
#-----------------------------------------------------------------------------------------------------------

geometricField = iron.Field()
geometricField.CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region)
geometricField.DecompositionSet(decomposition)
for dimensionIdx in range(1,numberOfDimensions+1):
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,dimensionIdx,1)
geometricField.CreateFinish()

# Set geometry from the generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

#-----------------------------------------------------------------------------------------------------------
#EQUATION SETS
#-----------------------------------------------------------------------------------------------------------

# Create standard Laplace equations set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
        iron.EquationsSetTypes.LAPLACE_EQUATION,
        iron.EquationsSetSubtypes.STANDARD_LAPLACE]
equationsSet.CreateStart(EQUATIONS_SET_USER_NUMBER,region,geometricField,
        equationsSetSpecification,EQUATIONS_SET_FIELD_USER_NUMBER,equationsSetField)
equationsSet.CreateFinish()

#-----------------------------------------------------------------------------------------------------------
#DEPENDENT FIELD
#-----------------------------------------------------------------------------------------------------------

dependentField = iron.Field()
equationsSet.DependentCreateStart(DEPENDENT_FIELD_USER_NUMBER,dependentField)
dependentField.DOFOrderTypeSet(iron.FieldVariableTypes.U,iron.FieldDOFOrderTypes.SEPARATED)
dependentField.DOFOrderTypeSet(iron.FieldVariableTypes.DELUDELN,iron.FieldDOFOrderTypes.SEPARATED)
equationsSet.DependentCreateFinish()

# Initialise dependent field
dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.5)

#-----------------------------------------------------------------------------------------------------------
# EQUATIONS
#-----------------------------------------------------------------------------------------------------------

equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
#equations.OutputTypeSet(iron.EquationsOutputTypes.NONE)
#equations.OutputTypeSet(iron.EquationsOutputTypes.MATRIX)
equations.OutputTypeSet(iron.EquationsOutputTypes.ELEMENT_MATRIX)
equationsSet.EquationsCreateFinish()

#-----------------------------------------------------------------------------------------------------------
#PROBLEM
#-----------------------------------------------------------------------------------------------------------

problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.CLASSICAL_FIELD,
        iron.ProblemTypes.LAPLACE_EQUATION,
        iron.ProblemSubtypes.STANDARD_LAPLACE]
problem.CreateStart(PROBLEM_USER_NUMBER,context,problemSpecification)
problem.CreateFinish()

# Create control loops
problem.ControlLoopCreateStart()
problem.ControlLoopCreateFinish()

#-----------------------------------------------------------------------------------------------------------
#SOLVER
#-----------------------------------------------------------------------------------------------------------

# Create problem solver
solver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
#solver.OutputTypeSet(iron.SolverOutputTypes.SOLVER)
solver.OutputTypeSet(iron.SolverOutputTypes.MATRIX)
solver.LinearTypeSet(iron.LinearSolverTypes.ITERATIVE)
solver.LinearIterativeAbsoluteToleranceSet(1.0E-12)
solver.LinearIterativeRelativeToleranceSet(1.0E-12)
problem.SolversCreateFinish()

#-----------------------------------------------------------------------------------------------------------
#SOLVER EQUATIONS
#-----------------------------------------------------------------------------------------------------------

# Create solver equations and add equations set to solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.SparsityTypeSet(iron.SolverEquationsSparsityTypes.SPARSE)
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()


#-----------------------------------------------------------------------------------------------------------
#BOUNDARY CONDITIONS
#-----------------------------------------------------------------------------------------------------------

# Create boundary conditions and set first and last nodes to 0.0 and 1.0
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
firstNodeNumber=1
nodes = iron.Nodes()
region.NodesGet(nodes)
lastNodeNumber = nodes.numberOfNodes
firstNodeDomain = decomposition.NodeDomainGet(firstNodeNumber,1)
lastNodeDomain = decomposition.NodeDomainGet(lastNodeNumber,1)
if firstNodeDomain == computationalNodeNumber:
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,firstNodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
if lastNodeDomain == computationalNodeNumber:
    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,lastNodeNumber,1,iron.BoundaryConditionsTypes.FIXED,1.0)
solverEquations.BoundaryConditionsCreateFinish()

#-----------------------------------------------------------------------------------------------------------
#SOLVE
#-----------------------------------------------------------------------------------------------------------

problem.Solve()

#-----------------------------------------------------------------------------------------------------------
#OUTPUT
#-----------------------------------------------------------------------------------------------------------

# Export results
baseName = "LaplaceEquation"
dataFormat = "PLAIN_TEXT"

fml = iron.FieldMLIO()
fml.OutputCreate(mesh, "", baseName, dataFormat)
fml.OutputAddFieldNoType(baseName+".geometric", dataFormat, geometricField,iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
fml.OutputAddFieldNoType(baseName+".phi", dataFormat, dependentField,iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
fml.OutputWrite("LaplaceEquation.xml")
fml.Finalise()

fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("LaplaceEquation","FORTRAN")
fields.ElementsExport("LaplaceEquation","FORTRAN")
fields.Finalise()

# Destroy the context
context.Destroy()
# Finalise OpenCMISS
iron.Finalise()
