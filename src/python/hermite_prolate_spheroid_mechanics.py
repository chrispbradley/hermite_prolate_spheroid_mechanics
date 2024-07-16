#!/usr/bin/env python
import os
from numpy import pi
from opencmiss.opencmiss import OpenCMISS_Python as oc

import prolate_spheroid_geometry

# Prolate spheroid geometry parameters:
cutoffAngle = 120.0 * pi / 180.0
focus = 37.5  # mm
endocardiumLambda = 0.38
epicardiumLambda = 0.69

# Fibre angles in radians:
epicardiumFibreAngle = 45.0 * pi / 180.0
endocardiumFibreAngle = -90.0 * pi / 180.0
# Assume constant sheet angle
sheetAngle = 90.0 * pi / 180.0

# Simulation parameters:
cavityPressure = 0.8
numIncrements = 1

# Geometric and hydrostatic pressure interpolation: 
interpolations = ['cubic_hermite', 'linear']
# Other interpolations should also work, eg:
#interpolations = ['quadratic', 'linear']
nodalPressure = True
geometricMeshComponent = 1
pressureMeshComponent = 2
hasDerivatives = interpolations[0] == 'cubic_hermite'
# Circumferential, longitudinal, transmural number of elements:
numberGlobalElements = [4, 2, 1]

# Set up a ProlateSpheroid object, which calculates the mesh geometry:
geometry = prolate_spheroid_geometry.ProlateSpheroid(focus, endocardiumLambda, epicardiumLambda, cutoffAngle,
        numberGlobalElements, endocardiumFibreAngle, epicardiumFibreAngle, sheetAngle,
        interpolations)

# Constitutive relation setup
# Guccione constitutive relation:
constitutiveRelation = oc.EquationsSetSubtypes.TRANSVERSE_ISOTROPIC_EXPONENTIAL
constitutiveParameters = [0.88, 0.0, 18.5, 3.58, 3.26]
initialHydrostaticPressure = 0.0

# User numbers for identifying OpenCMISS objects
contextUserNumber = 1
coordinateSystemUserNumber = 1
regionUserNumber = 1
meshUserNumber = 1
decompositionUserNumber = 1
decomposerUserNumber = 1
(geometricFieldUserNumber,
    fibreFieldUserNumber,
    materialFieldUserNumber,
    dependentFieldUserNumber,
    deformedFieldUserNumber,
    equationsSetFieldUserNumber) = range(1, 7)
equationsSetUserNumber = 1
problemUserNumber = 1

#quit()

context = oc.Context()
context.Create(contextUserNumber)

worldRegion = oc.Region()
context.WorldRegionGet(worldRegion)

# Get the number of computational nodes and this computational node number
# for when running in parallel with MPI
computationEnvironment = oc.ComputationEnvironment()
context.ComputationEnvironmentGet(computationEnvironment)

worldWorkGroup = oc.WorkGroup()
computationEnvironment.WorldWorkGroupGet(worldWorkGroup)
numberOfComputationalNodes = worldWorkGroup.NumberOfGroupNodesGet()
computationalNodeNumber = worldWorkGroup.GroupNodeNumberGet()

# Create a 3D rectangular cartesian coordinate system
coordinateSystem = oc.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber,context)
coordinateSystem.DimensionSet(3)
coordinateSystem.CreateFinish()

# Create a region within the world region and
# assign the coordinate system to the region
region = oc.Region()
region.CreateStart(regionUserNumber, worldRegion)
region.LabelSet("ProlateSpheroid")
region.CoordinateSystemSet(coordinateSystem)
region.CreateFinish()

# Create mesh from the prolate spheroid geometry
mesh = geometry.generateMesh(context, region)

# Create a decomposition for the mesh
# This breaks the mesh up into multiple decomposition
# domains when running in parallel
decomposition = oc.Decomposition()
decomposition.CreateStart(decompositionUserNumber, mesh)
# Have to enable face calculation for the decomposition
# in order to be able to use pressure boundary conditions
decomposition.CalculateFacesSet(True)
decomposition.CreateFinish()

# Decompose 
decomposer = oc.Decomposer()
decomposer.CreateStart(decomposerUserNumber,worldRegion,worldWorkGroup)
decompositionIndex = decomposer.DecompositionAdd(decomposition)
decomposer.CreateFinish()

# Create a field for the geometry
geometricField = oc.Field()
geometricField.CreateStart(geometricFieldUserNumber, region)
geometricField.DecompositionSet(decomposition)
geometricField.TypeSet(oc.FieldTypes.GEOMETRIC)
geometricField.VariableLabelSet(oc.FieldVariableTypes.U, "Geometry")
# Set the x, y and z components to use the first mesh component:
for component in range(1, 4):
    geometricField.ComponentMeshComponentSet(
            oc.FieldVariableTypes.U, component, geometricMeshComponent)
geometricField.ScalingTypeSet(oc.FieldScalingTypes.UNIT)
geometricField.CreateFinish()

# Use the prolate spheroid geometry to set the geometric field parameters
geometry.setGeometry(computationEnvironment,geometricField)

# Create a fibre field and attach it to the geometric field
# This has three components describing fibre orientations as
# rotations in radians about the z, y' and x'' base vectors.
fibreField = oc.Field()
fibreField.CreateStart(fibreFieldUserNumber, region)
fibreField.TypeSet(oc.FieldTypes.FIBRE)
fibreField.DecompositionSet(decomposition)
fibreField.GeometricFieldSet(geometricField)
fibreField.VariableLabelSet(oc.FieldVariableTypes.U, "Fibre")
fibreField.ScalingTypeSet(oc.FieldScalingTypes.UNIT)
# Use linear mesh component for fibre field
for component in range(1, 4):
    fibreField.ComponentMeshComponentSet(
                oc.FieldVariableTypes.U, component, pressureMeshComponent)
fibreField.CreateFinish()

# Use the prolate spheroid geometry to set up the fibre field values
geometry.setFibres(computationEnvironment,fibreField)

# Create the equations set and equations set field
# This defines the type of equations to solve
equationsSetField = oc.Field()
equationsSet = oc.EquationsSet()
equationsSetSpecification = [oc.EquationsSetClasses.ELASTICITY, 
                             oc.EquationsSetTypes.FINITE_ELASTICITY,
                             constitutiveRelation]
equationsSet.CreateStart(equationsSetUserNumber, region, fibreField,
    equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

# Create the material field, used for setting constitutive parameters
# This has an interpolation type of constant by default
materialField = oc.Field()
equationsSet.MaterialsCreateStart(materialFieldUserNumber, materialField)
materialField.VariableLabelSet(oc.FieldVariableTypes.U, "Material")
equationsSet.MaterialsCreateFinish()

# Set constant material parameters:
for (component, value) in enumerate(constitutiveParameters, 1):
    materialField.ComponentValuesInitialise(
            oc.FieldVariableTypes.U, oc.FieldParameterSetTypes.VALUES,
            component, value)

# Create the dependent field for storing the solution
# This has one U variable for displacement values and hydrostatic pressure,
# and one T variable for the nodal forces
# Sensible defaults are defined by the equations set, but we can
# override some settings here, eg. to change the interpolation type
# for the hydrostatic pressure to node based
dependentField = oc.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
dependentField.VariableLabelSet(oc.FieldVariableTypes.U, "Dependent")
if nodalPressure:
    dependentField.ComponentInterpolationSet(
            oc.FieldVariableTypes.U, 4,
            oc.FieldInterpolationTypes.NODE_BASED)
    dependentField.ComponentInterpolationSet(
            oc.FieldVariableTypes.T, 4,
            oc.FieldInterpolationTypes.NODE_BASED)
    dependentField.ComponentMeshComponentSet(
            oc.FieldVariableTypes.U, 4, pressureMeshComponent)
    dependentField.ComponentMeshComponentSet(
            oc.FieldVariableTypes.T, 4, pressureMeshComponent)
else:
    dependentField.ComponentInterpolationSet(
            oc.FieldVariableTypes.U, 4,
            oc.FieldInterpolationTypes.ELEMENT_BASED)
    dependentField.ComponentInterpolationSet(
            oc.FieldVariableTypes.T, 4,
            oc.FieldInterpolationTypes.ELEMENT_BASED)
dependentField.ScalingTypeSet(oc.FieldScalingTypes.UNIT)
equationsSet.DependentCreateFinish()

# Initialise dependent field position values from the undeformed geometry
for component in [1, 2, 3]:
    geometricField.ParametersToFieldParametersComponentCopy(
        oc.FieldVariableTypes.U,
        oc.FieldParameterSetTypes.VALUES, component,
        dependentField, oc.FieldVariableTypes.U,
        oc.FieldParameterSetTypes.VALUES, component)
    dependentField.ComponentValuesInitialise(
        oc.FieldVariableTypes.T,
        oc.FieldParameterSetTypes.VALUES,
        component, 0.0)

# Set the initial hydrostatic pressure
dependentField.ComponentValuesInitialise(
    oc.FieldVariableTypes.U,
    oc.FieldParameterSetTypes.VALUES,
    4, initialHydrostaticPressure)
dependentField.ComponentValuesInitialise(
    oc.FieldVariableTypes.T,
    oc.FieldParameterSetTypes.VALUES,
    4, 0.0)

# Create a deformed geometry field, as cmgui doesn't like displaying
# deformed fibres from the dependent field because it isn't a geometric field.
deformedField = oc.Field()
deformedField.CreateStart(deformedFieldUserNumber, region)
deformedField.DecompositionSet(decomposition)
deformedField.TypeSet(oc.FieldTypes.GEOMETRIC)
deformedField.VariableLabelSet(oc.FieldVariableTypes.U, "DeformedGeometry")
for component in [1, 2, 3]:
    deformedField.ComponentMeshComponentSet(
            oc.FieldVariableTypes.U, component,
            geometricMeshComponent)
deformedField.ScalingTypeSet(oc.FieldScalingTypes.UNIT)
deformedField.CreateFinish()

# Create equations
equations = oc.Equations()
equationsSet.EquationsCreateStart(equations)
equations.SparsityTypeSet(oc.EquationsSparsityTypes.SPARSE)
equations.OutputTypeSet(oc.EquationsOutputTypes.NONE)
equationsSet.EquationsCreateFinish()

# Define the problem
# The problem defines how the equations are solved by
# setting up control loops and solvers
problem = oc.Problem()
problemSpecification = [oc.ProblemClasses.ELASTICITY,
        oc.ProblemTypes.FINITE_ELASTICITY,
        oc.ProblemSubtypes.STATIC_FINITE_ELASTICITY]
problem.CreateStart(problemUserNumber,context,problemSpecification)
problem.CreateFinish()

# Create the problem control loops
# For static finite elasticity, there is just a
# single load increment control loop, and we set
# the number of load increments by setting the number of
# iterations of this control loop
problem.ControlLoopCreateStart()
controlLoop = oc.ControlLoop()
problem.ControlLoopGet([oc.ControlLoopIdentifiers.NODE], controlLoop)
controlLoop.MaximumIterationsSet(numIncrements)
controlLoop.OutputTypeSet(oc.ControlLoopOutputTypes.PROGRESS)
problem.ControlLoopCreateFinish()

# Create the problem solver
solver = oc.Solver()
problem.SolversCreateStart()
problem.SolverGet([oc.ControlLoopIdentifiers.NODE], 1, solver)
solver.OutputTypeSet(oc.SolverOutputTypes.MONITOR)
solver.NewtonJacobianCalculationTypeSet(oc.JacobianCalculationTypes.FD)
solver.NewtonRelativeToleranceSet(1.0e-14)
solver.NewtonAbsoluteToleranceSet(1.0e-14)
solver.NewtonSolutionToleranceSet(1.0e-14)
# Adjust settings for the line search solver
linesearchSolver = oc.Solver()
solver.NewtonLinearSolverGet(linesearchSolver)
linesearchSolver.LinearTypeSet(oc.LinearSolverTypes.DIRECT)
problem.SolversCreateFinish()

# Create solver equations for the problem and add
# the single equations set to solver equations.
# For coupled problems there may be multiple equations sets
# solved within one solver equations
solver = oc.Solver()
solverEquations = oc.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([oc.ControlLoopIdentifiers.NODE], 1, solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.SparsityTypeSet(oc.SolverEquationsSparsityTypes.SPARSE)
# Only a single equations set to add, the finite elasticity equations:
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Prescribe boundary conditions on the solver equations
boundaryConditions = oc.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

def getDomainNodes(computationEnvironment, geometry, decomposition, component):
    component_name = interpolations[component - 1]
    computationalNodeNumber = computationEnvironment.WorldNodeNumberGet()
    nodes = geometry.componentNodes(component_name)
    meshComponent = geometry.meshComponent(component_name)
    return set(node for node in nodes
        if decomposition.NodeDomainGet(node, meshComponent) == computationalNodeNumber)
geometricDomainNodes = getDomainNodes(computationEnvironment, geometry, decomposition, geometricMeshComponent)

# Fix epicardium nodes at the base:
baseNodes = set(geometry.nodeGroup('base'))
externalNodes = set(geometry.nodeGroup('external'))
fixedNodes = baseNodes.intersection(externalNodes)
for node in fixedNodes.intersection(geometricDomainNodes):
    for component in (1, 2, 3):
        derivative, version = 1, 1
        # AddNode is used, as this fixes the nodal value at an
        # increment (in this case 0) from the current value
        boundaryConditions.AddNode(dependentField, oc.FieldVariableTypes.U,
                version, derivative, node, component,
                oc.BoundaryConditionsTypes.FIXED, 0.0)
        if hasDerivatives:
            derivative = oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1
            boundaryConditions.AddNode(dependentField, oc.FieldVariableTypes.U,
                    version, derivative, node, component,
                    oc.BoundaryConditionsTypes.FIXED, 0.0)

# Constrain degrees of freedom at the apex to collapse faces:
constrainedNodeSets = geometry.constrainedNodes()
for nodes in constrainedNodeSets:
    in_domain = [n in geometricDomainNodes for n in nodes]
    if all(in_domain):
        version = 1
        for component in (1, 2, 3):
            # Map nodal values + xi_3 derivative to same value,
            if hasDerivatives:
                mappedDerivatives = [
                    oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,
                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S3]
            else:
                mappedDerivatives = [
                    oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV]
            for derivative in mappedDerivatives:
                boundaryConditions.ConstrainNodeDofsEqual(dependentField, oc.FieldVariableTypes.U,
                        version, derivative, component, nodes, 1.0)
            # Fix xi_1 derivative to be zero
            if hasDerivatives:
                for node in nodes:
                    for derivative in [
                            oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,
                            oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,
                            oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S3,
                            oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2_S3,
                            ]:
                        boundaryConditions.SetNode(dependentField, oc.FieldVariableTypes.U,
                                version, derivative, node, component,
                                oc.BoundaryConditionsTypes.FIXED, 0.0)
    elif any(in_domain):
        raise RuntimeError("Mesh decomposition has DOFs "
                "that must be constrained to be equal in separate domains")

# Apply incremented cavity pressure on the endocardial surface:
internalNodes = set(geometry.nodeGroup('internal'))
for node in internalNodes.intersection(geometricDomainNodes):
    derivative, version = 1, 1
    # xi_3 is the transmural direction
    xiDirection = 3
    # For pressure/force boundary conditions, the T field variable is
    # constrained rather than the U field variable
    boundaryConditions.SetNode(dependentField, oc.FieldVariableTypes.T,
            version, derivative, node, xiDirection,
            oc.BoundaryConditionsTypes.PRESSURE_INCREMENTED, -cavityPressure)

solverEquations.BoundaryConditionsCreateFinish()

# Solve the problem
problem.Solve()

# Copy deformed geometry into deformed field
for component in [1, 2, 3]:
    dependentField.ParametersToFieldParametersComponentCopy(
        oc.FieldVariableTypes.U,
        oc.FieldParameterSetTypes.VALUES, component,
        deformedField, oc.FieldVariableTypes.U,
        oc.FieldParameterSetTypes.VALUES, component)

if not os.path.exists("./results"):
    os.makedirs("./results")

# Export results to exnode/exelem files
fields = oc.Fields()
fields.CreateRegion(region)
fields.NodesExport("./results/prolate_spheroid", "FORTRAN")
fields.ElementsExport("./results/prolate_spheroid", "FORTRAN")
fields.Finalise()
