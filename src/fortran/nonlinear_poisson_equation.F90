PROGRAM NonlinearPoissonEquation

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  !-----------------------------------------------------------------------------------------------------------
  ! PROGRAM VARIABLES AND TYPES
  !-----------------------------------------------------------------------------------------------------------

  !Program parameters
  REAL(CMISSRP), PARAMETER :: HEIGHT=0.5_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=0.5_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=1.0_CMISSRP

  REAL(CMISSRP), PARAMETER :: PI=3.141592653589793_CMISSRP
  
  INTEGER(CMISSIntg), PARAMETER :: CONTEXT_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: REGION_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: BASIS_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: GENERATED_MESH_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: MESH_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: DECOMPOSITION_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: DECOMPOSER_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: GEOMETRIC_FIELD_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER=2
  INTEGER(CMISSIntg), PARAMETER :: DEPENDENT_FIELD_USER_NUMBER=3
  INTEGER(CMISSIntg), PARAMETER :: MATERIALS_FIELD_USER_NUMBER=4
  INTEGER(CMISSIntg), PARAMETER :: ANALYTIC_FIELD_USER_NUMBER=5
  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: PROBLEM_USER_NUMBER=1

  !Program variables
  INTEGER(CMISSIntg) :: argumentLength,componentIdx,computationalNodeNumber,decompositionIndex,equationsSetIndex,err, &
    & gaussOrder,interpolationType,nodeDomain,nodeNumber,numberOfArguments,numberOfComponents,numberOfComputationalNodes, &
    & numberOfDimensions,numberOfGaussXi,numberOfGlobalXElements,numberOfGlobalXNodes,numberOfGlobalYElements, &
    & numberOfGlobalYNodes,numberOfGlobalZElements,numberOfGlobalZNodes,numberOfNodesXi,parameterIdx,status,yNodeIdx,zNodeIdx
  CHARACTER(LEN=255) :: commandArgument,filename
  LOGICAL :: exportField

  !OpenCMISS variables
  TYPE(cmfe_BasisType) :: basis
  TYPE(cmfe_BoundaryConditionsType) :: boundaryConditions
  TYPE(cmfe_ComputationEnvironmentType) :: computationEnvironment
  TYPE(cmfe_ContextType) :: context
  TYPE(cmfe_CoordinateSystemType) :: coordinateSystem
  TYPE(cmfe_DecompositionType) :: decomposition
  TYPE(cmfe_DecomposerType) :: decomposer
  TYPE(cmfe_EquationsType) :: equations
  TYPE(cmfe_EquationsSetType) :: equationsSet
  TYPE(cmfe_FieldType) :: analyticField,dependentField,equationsSetField,geometricField,materialsField
  TYPE(cmfe_FieldsType) :: fields
  TYPE(cmfe_GeneratedMeshType) :: generatedMesh
  TYPE(cmfe_MeshType) :: mesh
  TYPE(cmfe_ProblemType) :: problem
  TYPE(cmfe_RegionType) :: region,worldRegion
  TYPE(cmfe_SolverType) :: linearSolver,solver
  TYPE(cmfe_SolverEquationsType) :: solverEquations
  TYPE(cmfe_WorkGroupType) :: worldWorkGroup

  !-----------------------------------------------------------------------------------------------------------
  ! GET PROGRAM ARGUMENTS
  !-----------------------------------------------------------------------------------------------------------
  
  !Get input arguments
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
    IF(numberOfGlobalYElements<0) CALL HandleError("Invalid number of Y elements.")
    CALL GET_COMMAND_ARGUMENT(3,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 3.")
    READ(commandArgument(1:argumentLength),*) numberOfGlobalZElements
    IF(numberOfGlobalZElements<0) CALL HandleError("Invalid number of Z elements.")
    CALL GET_COMMAND_ARGUMENT(4,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 4.")
    READ(commandArgument(1:argumentLength),*) interpolationType
    IF(interpolationType<=0) CALL HandleError("Invalid Interpolation specification.")
    IF(numberOfGlobalZElements>0) THEN
      numberOfDimensions=3
    ELSEIF(numberOfGlobalYElements>0) THEN
      numberOfDimensions=2
    ELSE
      numberOfDimensions=1
    ENDIF
  ELSE
    !If there are not enough arguments default the problem specification
    numberOfDimensions=2
    numberOfGlobalXElements=5
    numberOfGlobalYElements=5
    numberOfGlobalZElements=0
    interpolationType=1
  ENDIF

  !-----------------------------------------------------------------------------------------------------------
  ! INITITALISE
  !-----------------------------------------------------------------------------------------------------------  

  !Intialise OpenCMISS
  CALL cmfe_Initialise(err)
  !Trap all errors
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,err)
  !Output to a file
  WRITE(filename,'(A,"_",I0,"x",I0,"x",I0,"_",I0)') "NonlinearPoisson",numberOfGlobalXElements,numberOfGlobalYElements, &
    & numberOfGlobalZElements,interpolationType
  CALL cmfe_OutputSetOn(filename,err)
  !Create a context
  CALL cmfe_Context_Initialise(context,err)
  CALL cmfe_Context_Create(CONTEXT_USER_NUMBER,context,err)
  !Get the world region
  CALL cmfe_Region_Initialise(worldRegion,err)
  CALL cmfe_Context_WorldRegionGet(context,worldRegion,err)

  !Get the computational nodes information
  CALL cmfe_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL cmfe_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL cmfe_WorkGroup_Initialise(worldWorkGroup,err)
  CALL cmfe_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL cmfe_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
  CALL cmfe_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)

  !-----------------------------------------------------------------------------------------------------------
  ! COORDINATE SYSTEM
  !-----------------------------------------------------------------------------------------------------------  

  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL cmfe_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
  !Set the coordinate system number of dimensions
  CALL cmfe_CoordinateSystem_DimensionSet(coordinateSystem,numberOfDimensions,err)
  !Set the origin
  CALL cmfe_CoordinateSystem_OriginSet(coordinateSystem,[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP],err)
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(coordinateSystem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! REGION
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the region
  CALL cmfe_Region_Initialise(region,err)
  CALL cmfe_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
  CALL cmfe_Region_LabelSet(region,"NonlinearPoisson",err)
  !Set the regions coordinate system to the RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(region,coordinateSystem,err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(region,err)

  !-----------------------------------------------------------------------------------------------------------
  ! BASIS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a basis (default is trilinear lagrange)
  CALL cmfe_Basis_Initialise(basis,err)
  CALL cmfe_Basis_CreateStart(BASIS_USER_NUMBER,context,basis,err)
  CALL cmfe_Basis_NumberOfXiSet(basis,numberOfDimensions,err)
  SELECT CASE(interpolationType)
  CASE(CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
    & CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
    & CMFE_BASIS_CUBIC_LAGRANGE_INTERPOLATION, &
    & CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION)
    CALL cmfe_Basis_TypeSet(basis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
  CASE(CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
    & CMFE_BASIS_QUADRATIC_SIMPLEX_INTERPOLATION, &
    & CMFE_BASIS_CUBIC_SIMPLEX_INTERPOLATION)
    CALL cmfe_Basis_TypeSet(basis,CMFE_BASIS_SIMPLEX_TYPE,err)
  CASE DEFAULT
    CALL HandleError("Invalid interpolation type.")
  END SELECT
  SELECT CASE(interpolationType)
  CASE(CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION)
    numberOfGaussXi=2
    numberOfNodesXi=2
    gaussOrder=0
  CASE(CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION)
    numberOfGaussXi=3
    numberOfNodesXi=3
    gaussOrder=0
  CASE(CMFE_BASIS_CUBIC_LAGRANGE_INTERPOLATION)
    numberOfGaussXi=4
    numberOfNodesXi=4
    gaussOrder=0
  CASE(CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION)
    numberOfGaussXi=4
    numberOfNodesXi=2
    gaussOrder=0
  CASE(CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION)
    numberOfGaussXi=0
    numberOfNodesXi=2    
    gaussOrder=3
  CASE(CMFE_BASIS_QUADRATIC_SIMPLEX_INTERPOLATION)
    numberOfGaussXi=0
    numberOfNodesXi=3
    gaussOrder=4
  CASE(CMFE_BASIS_CUBIC_SIMPLEX_INTERPOLATION)
    numberOfGaussXi=0
    numberOfNodesXi=4
    gaussOrder=5
  CASE DEFAULT
    CALL HandleError("Invalid interpolation type.")
  END SELECT
  IF(numberOfDimensions==1) THEN
    CALL cmfe_Basis_InterpolationXiSet(basis,[interpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi],err)
    ELSE
      CALL cmfe_Basis_QuadratureOrderSet(basis,gaussOrder,err)
    ENDIF
  ELSE IF(numberOfDimensions==2) THEN
    CALL cmfe_Basis_InterpolationXiSet(basis,[interpolationType,interpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi],err)
    ELSE
      CALL cmfe_Basis_QuadratureOrderSet(basis,gaussOrder,err)
    ENDIF
  ELSE
    CALL cmfe_Basis_InterpolationXiSet(basis,[interpolationType,interpolationType,interpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],err)
    ELSE
      CALL cmfe_Basis_QuadratureOrderSet(basis,gaussOrder,err)      
    ENDIF
  ENDIF
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(basis,err)

  !Compute the number of nodes in each global direction based on the interpolation type used.
  numberOfGlobalXNodes=(numberOfNodesXi-1)*numberOfGlobalXElements+1
  numberOfGlobalYNodes=(numberOfNodesXi-1)*numberOfGlobalYElements+1
  numberOfGlobalZNodes=(numberOfNodesXi-1)*numberOfGlobalZElements+1

  !-----------------------------------------------------------------------------------------------------------
  ! MESH
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(generatedMesh,err)
  CALL cmfe_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(generatedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(generatedMesh,basis,err)
  !Define the mesh on the region
  IF(numberOfDimensions==1) THEN
    CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH],err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements],err)
  ELSEIF(numberOfDimensions==2) THEN
    CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT],err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements],err)
  ELSE
    CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT,LENGTH],err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements, &
      & numberOfGlobalZElements],err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(mesh,err)
  CALL cmfe_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DECOMPOSITION
  !-----------------------------------------------------------------------------------------------------------

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(decomposition,err)
  CALL cmfe_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(decomposition,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DECOMPOSER
  !-----------------------------------------------------------------------------------------------------------
  
  CALL cmfe_Decomposer_Initialise(decomposer,err)
  CALL cmfe_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL cmfe_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL cmfe_Decomposer_CreateFinish(decomposer,err)

  !-----------------------------------------------------------------------------------------------------------
  ! GEOMETRIC FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(geometricField,err)
  CALL cmfe_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region,geometricField,err)
  !Set the decomposition to use
  CALL cmfe_Field_DecompositionSet(geometricField,decomposition,err)
  !Set the domain to be used by the field components.
  DO componentIdx=1,numberOfDimensions
    CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,componentIdx,1,err)
  ENDDO
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(geometricField,err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)

  !-----------------------------------------------------------------------------------------------------------
  ! MATERIAL FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create a material field for the scalars a and b and the Voigt components of the rank-2 conductivity tensor.
  CALL cmfe_Field_Initialise(materialsField,err)
  CALL cmfe_Field_CreateStart(MATERIALS_FIELD_USER_NUMBER,region,materialsField,err)
  !Set the decomposition to use
  CALL cmfe_Field_DecompositionSet(materialsField,decomposition,err)
  !Set the type
  CALL cmfe_Field_TypeSet(materialsField,CMFE_FIELD_MATERIAL_TYPE,err)
  CALL cmfe_Field_NumberOfVariablesSet(materialsField,1,err)
  CALL cmfe_Field_VariableTypesSet(materialsField,[CMFE_FIELD_U_VARIABLE_TYPE],err)
  numberOfComponents=2+numberOfDimensions+MAX(numberOfDimensions-1,0)+MAX(numberOfDimensions-2,0)
  CALL cmfe_Field_NumberOfComponentsSet(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,numberOfComponents,err)
  DO componentIdx=1,numberOfComponents
    CALL cmfe_Field_ComponentMeshComponentSet(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,componentIdx,1,err)
    CALL cmfe_Field_ComponentInterpolationSet(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,componentIdx, &
      & CMFE_FIELD_CONSTANT_INTERPOLATION,err)
  ENDDO !componentIdx
  !Set associated geometric field
  CALL cmfe_Field_GeometricFieldSet(materialsField,geometricField,err)
  !Set the label of the field
  CALL cmfe_Field_VariableLabelSet(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,"Parameters",err)
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(materialsField,err)

  !Update material field components
  !a parameter:
  CALL cmfe_Field_ComponentValuesInitialise(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,-1.0_CMISSRP,err)
  !b parameter:
  CALL cmfe_Field_ComponentValuesInitialise(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,0.5_CMISSRP,err)
  !Rank-2 conductivity tensor diagonal components in Voigt form:
  DO componentIdx=1,numberOfDimensions
    CALL cmfe_Field_ComponentValuesInitialise(materialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 2+componentIdx,1.0_CMISSRP,err)
  ENDDO

  !-----------------------------------------------------------------------------------------------------------
  ! DEPENDENT FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create dependent fields
  CALL cmfe_Field_Initialise(dependentField,err)
  CALL cmfe_Field_CreateStart(DEPENDENT_FIELD_USER_NUMBER,region,dependentField,err)
  !Set the decomposition to use
  CALL cmfe_Field_DecompositionSet(dependentField,decomposition,err)
  !Set the type
  CALL cmfe_Field_TypeSet(dependentField,CMFE_FIELD_GENERAL_TYPE,err)
  CALL cmfe_Field_DependentTypeSet(dependentField,CMFE_FIELD_DEPENDENT_TYPE,err)
  !Set associated geometric field
  CALL cmfe_Field_GeometricFieldSet(dependentField,geometricField,err)
  !Two dependent variables: primary variable 'U' and secondary variable 'DelUDelN'
  CALL cmfe_Field_NumberOfVariablesSet(dependentField,2,err)
  CALL cmfe_Field_DimensionSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,err)
  CALL cmfe_Field_DimensionSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,err)
  CALL cmfe_Field_DOFOrderTypeSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SEPARATED_COMPONENT_DOF_ORDER,err)
  CALL cmfe_Field_DOFOrderTypeSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,CMFE_FIELD_SEPARATED_COMPONENT_DOF_ORDER,err)
  !Set appropriate labels
  CALL cmfe_Field_VariableLabelSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,"U",err)
  CALL cmfe_Field_VariableLabelSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,"DelUDelN",err)
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(dependentField,err)

  !Update dependent fields (initial guess for nonlinear solver)
  CALL cmfe_Field_ComponentValuesInitialise(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,-5.0_CMISSRP,err)

  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS SETS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations_sets
  CALL cmfe_EquationsSet_Initialise(equationsSet,err)
  CALL cmfe_Field_Initialise(equationsSetField,err)
  !Poisson equation with exponential (nonlinear) source term
  CALL cmfe_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,geometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_POISSON_EQUATION_TYPE,CMFE_EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE], &
    & EQUATIONS_SET_FIELD_USER_NUMBER,equationsSetField,equationsSet,err)
  !Finish creating eqiations sets
  CALL cmfe_EquationsSet_CreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DEPENDENT FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Attach dependent field already created to equation sets
  CALL cmfe_EquationsSet_DependentCreateStart(equationsSet,DEPENDENT_FIELD_USER_NUMBER,dependentField,err)
  CALL cmfe_EquationsSet_DependentCreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! MATERIAL FIELD
  !-----------------------------------------------------------------------------------------------------------
  
  !Attach material field already created to equation sets
  CALL cmfe_EquationsSet_MaterialsCreateStart(equationsSet,MATERIALS_FIELD_USER_NUMBER,materialsField,err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! ANALYTIC FIELD
  !-----------------------------------------------------------------------------------------------------------
  
  !Create the equations set analytic field variables
  IF(numberOfDimensions==2) THEN
    CALL cmfe_Field_Initialise(analyticField,err)
    CALL cmfe_EquationsSet_AnalyticCreateStart(equationsSet,CMFE_EQUATIONS_SET_EXPONENTIAL_POISSON_EQUATION_TWO_DIM_1, &
      & ANALYTIC_FIELD_USER_NUMBER,analyticField,Err)
    !Finish the equations set analytic field variables
    CALL cmfe_EquationsSet_AnalyticCreateFinish(equationsSet,err)
    !Set the analytic solution constants
    !A
    CALL cmfe_Field_ComponentValuesInitialise(analyticField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 1,PI/(100.0_CMISSRP*LENGTH),err)
    !B
    CALL cmfe_Field_ComponentValuesInitialise(analyticField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & 2,PI/(100.0_CMISSRP*HEIGHT),err)
  ENDIF
  
  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(equations,err)
  CALL cmfe_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
  !Set the equations matrices sparsity type
  !CALL cmfe_Equations_SparsityTypeSet(equations,CMFE_EQUATIONS_FULL_MATRICES,err)
  CALL cmfe_Equations_SparsityTypeSet(equations,CMFE_EQUATIONS_SPARSE_MATRICES,err)
  !Set the equations set output
  !CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_NO_OUTPUT,err)
  CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM
  !-----------------------------------------------------------------------------------------------------------  

  !Start the creation of a problem.
  CALL cmfe_Problem_Initialise(problem,err)
  CALL cmfe_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[CMFE_PROBLEM_CLASSICAL_FIELD_CLASS, &
    & CMFE_PROBLEM_POISSON_EQUATION_TYPE,CMFE_PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE],problem,err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(problem,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! CONTROL LOOP
  !-----------------------------------------------------------------------------------------------------------  

  !Start the creation of the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(problem,err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVER
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the problem solvers
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_Solver_Initialise(linearSolver,err)
  CALL cmfe_Problem_SolversCreateStart(problem,err)
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
  !Set the solver output
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_NO_OUTPUT,err)
  CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_PROGRESS_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_TIMING_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_SOLVER_OUTPUT,err)
  CALL cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_MATRIX_OUTPUT,err)
  !Set the Jacobian type to be either calculated analytically or via finite differences
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(solver,CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,err)
  !CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(solver,CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,err)
  !Set the nonlinear Newton solver tolerances etc.
  CALL cmfe_Solver_NewtonAbsoluteToleranceSet(solver,1.0E-8_CMISSRP,err)
  CALL cmfe_Solver_NewtonRelativeToleranceSet(solver,1.0E-8_CMISSRP,err)
  CALL cmfe_Solver_NewtonMaximumIterationsSet(solver,100000,err)
  !Get the associated linear solver
  CALL cmfe_Solver_NewtonLinearSolverGet(solver,linearSolver,err)
  !Set the linear iterative solver tolerances etc.
  CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(linearSolver,1.0E-8_CMISSRP,err)
  CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(linearSolver,1.0E-8_CMISSRP,err)
  CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(linearSolver,10000,err)
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVER EQUATIONS
  !-----------------------------------------------------------------------------------------------------------  

  !Start the creation of the problem solver equations
  CALL cmfe_Solver_Initialise(solver,err)
  CALL cmfe_SolverEquations_Initialise(solverEquations,err)
  CALL cmfe_Problem_SolverEquationsCreateStart(problem,err)
  !Get the solve equations
  CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,err)
  CALL cmfe_Solver_SolverEquationsGet(solver,solverEquations,err)
  !Set the solver equations sparsity
  !CALL cmfe_SolverEquations_SparsityTypeSet(solverEquations,CMFE_SOLVER_SPARSE_MATRICES,err)
  CALL cmfe_SolverEquations_SparsityTypeSet(solverEquations,CMFE_SOLVER_FULL_MATRICES,err)
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! BOUNDARY CONDITIONS
  !-----------------------------------------------------------------------------------------------------------

  !Set up the boundary conditions
  CALL cmfe_BoundaryConditions_Initialise(boundaryConditions,err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)
  IF(numberOfDimensions==2) THEN
    !Use analytic boundary conditions
    CALL cmfe_SolverEquations_BoundaryConditionsAnalytic(solverEquations,err)
  ELSE
    !Set the fixed boundary conditions on opposide sides
    DO zNodeIdx=1,numberOfGlobalZNodes
      DO yNodeIdx=1,numberOfGlobalYNodes
        !Left side
        nodeNumber=1+(yNodeIdx-1)*numberOfGlobalXNodes+(zNodeIdx-1)*numberOfGlobalXNodes*numberOfGlobalYNodes
        !Check what computational node the node is on
        CALL cmfe_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
        IF(nodeDomain==computationalNodeNumber) THEN
          !If the node is on my computational node then set left hand side BC to 0.0
          CALL cmfe_BoundaryConditions_SetNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,1, &
            & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
        ENDIF
        !Right side
        nodeNumber=numberOfGlobalXNodes+(yNodeIdx-1)*numberOfGlobalXNodes+(zNodeIdx-1)*numberOfGlobalXNodes*numberOfGlobalYNodes
        !Check what computational node the node is on
        CALL cmfe_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
        IF(nodeDomain==computationalNodeNumber) THEN
          !If the node is on my computational node then set left hand side BC to 0.0
          CALL cmfe_BoundaryConditions_SetNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,1, &
            & CMFE_BOUNDARY_CONDITION_FIXED,1.0_CMISSRP,err)
        ENDIF
      ENDDO !yNodeIdx
    ENDDO !zNodeIdx
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVE
  !-----------------------------------------------------------------------------------------------------------

  !Export results
  exportField=.TRUE.
  IF(exportField) THEN
    CALL cmfe_Fields_Initialise(fields,err)
    CALL cmfe_Fields_Create(region,fields,err)
    CALL cmfe_Fields_NodesExport(fields,"NonlinearPoisson","FORTRAN",err)
    CALL cmfe_Fields_ElementsExport(fields,"NonlinearPoisson","FORTRAN",err)
    CALL cmfe_Fields_Finalise(fields,err)
  ENDIF
  
  !Solve the problem
  CALL cmfe_Problem_Solve(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! OUTPUT
  !-----------------------------------------------------------------------------------------------------------

  IF(numberOfDimensions==2) THEN
    !Output Analytic analysis
    CALL cmfe_AnalyticAnalysis_Output(dependentField,"NonlinearPoissonAnalytic",err)
  ENDIF
  
  !Export results
  exportField=.TRUE.
  IF(exportField) THEN
    CALL cmfe_Fields_Initialise(fields,err)
    CALL cmfe_Fields_Create(region,fields,err)
    CALL cmfe_Fields_NodesExport(fields,"NonlinearPoisson","FORTRAN",err)
    CALL cmfe_Fields_ElementsExport(fields,"NonlinearPoisson","FORTRAN",err)
    CALL cmfe_Fields_Finalise(fields,err)
  ENDIF

  !-----------------------------------------------------------------------------------------------------------
  ! FINALISE
  !-----------------------------------------------------------------------------------------------------------  

  !Destroy the context
  CALL cmfe_Context_Destroy(context,err)
  !Finialise OpenCMISS
  CALL cmfe_Finalise(err)

  WRITE(*,'("Program successfully completed.")')
  STOP

CONTAINS

  SUBROUTINE HandleError(errorString)

    CHARACTER(LEN=*), INTENT(IN) :: errorString

    WRITE(*,'(">>ERROR: ",A)') errorString(1:LEN_TRIM(errorString))
    STOP

  END SUBROUTINE HandleError

END PROGRAM NonlinearPoissonEquation
