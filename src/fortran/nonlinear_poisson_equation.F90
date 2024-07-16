PROGRAM NonlinearPoissonEquation

  USE OpenCMISS

  IMPLICIT NONE

  !-----------------------------------------------------------------------------------------------------------
  ! PROGRAM VARIABLES AND TYPES
  !-----------------------------------------------------------------------------------------------------------

  !Program parameters
  REAL(OC_RP), PARAMETER :: HEIGHT=0.5_OC_RP
  REAL(OC_RP), PARAMETER :: WIDTH=0.5_OC_RP
  REAL(OC_RP), PARAMETER :: LENGTH=1.0_OC_RP

  REAL(OC_RP), PARAMETER :: PI=3.141592653589793_OC_RP
  
  INTEGER(OC_Intg), PARAMETER :: CONTEXT_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: REGION_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: BASIS_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: GENERATED_MESH_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: MESH_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSITION_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: DECOMPOSER_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: GEOMETRIC_FIELD_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER=2
  INTEGER(OC_Intg), PARAMETER :: DEPENDENT_FIELD_USER_NUMBER=3
  INTEGER(OC_Intg), PARAMETER :: MATERIALS_FIELD_USER_NUMBER=4
  INTEGER(OC_Intg), PARAMETER :: ANALYTIC_FIELD_USER_NUMBER=5
  INTEGER(OC_Intg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=1
  INTEGER(OC_Intg), PARAMETER :: PROBLEM_USER_NUMBER=1

  !Program variables
  INTEGER(OC_Intg) :: argumentLength,componentIdx,computationalNodeNumber,decompositionIndex,equationsSetIndex,err, &
    & gaussOrder,interpolationType,nodeDomain,nodeNumber,numberOfArguments,numberOfComponents,numberOfComputationalNodes, &
    & numberOfDimensions,numberOfGaussXi,numberOfGlobalXElements,numberOfGlobalXNodes,numberOfGlobalYElements, &
    & numberOfGlobalYNodes,numberOfGlobalZElements,numberOfGlobalZNodes,numberOfNodesXi,parameterIdx,status,yNodeIdx,zNodeIdx
  CHARACTER(LEN=255) :: commandArgument,filename
  LOGICAL :: exportField

  !OpenCMISS variables
  TYPE(OC_BasisType) :: basis
  TYPE(OC_BoundaryConditionsType) :: boundaryConditions
  TYPE(OC_ComputationEnvironmentType) :: computationEnvironment
  TYPE(OC_ContextType) :: context
  TYPE(OC_CoordinateSystemType) :: coordinateSystem
  TYPE(OC_DecompositionType) :: decomposition
  TYPE(OC_DecomposerType) :: decomposer
  TYPE(OC_EquationsType) :: equations
  TYPE(OC_EquationsSetType) :: equationsSet
  TYPE(OC_FieldType) :: analyticField,dependentField,equationsSetField,geometricField,materialsField
  TYPE(OC_FieldsType) :: fields
  TYPE(OC_GeneratedMeshType) :: generatedMesh
  TYPE(OC_MeshType) :: mesh
  TYPE(OC_ProblemType) :: problem
  TYPE(OC_RegionType) :: region,worldRegion
  TYPE(OC_SolverType) :: linearSolver,solver
  TYPE(OC_SolverEquationsType) :: solverEquations
  TYPE(OC_WorkGroupType) :: worldWorkGroup

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
  CALL OC_Initialise(err)
  !Trap all errors
  CALL OC_ErrorHandlingModeSet(OC_ERRORS_TRAP_ERROR,err)
  !Output to a file
  WRITE(filename,'(A,"_",I0,"x",I0,"x",I0,"_",I0)') "NonlinearPoisson",numberOfGlobalXElements,numberOfGlobalYElements, &
    & numberOfGlobalZElements,interpolationType
  CALL OC_OutputSetOn(filename,err)
  !Create a context
  CALL OC_Context_Initialise(context,err)
  CALL OC_Context_Create(CONTEXT_USER_NUMBER,context,err)
  !Get the world region
  CALL OC_Region_Initialise(worldRegion,err)
  CALL OC_Context_WorldRegionGet(context,worldRegion,err)

  !Get the computational nodes information
  CALL OC_ComputationEnvironment_Initialise(computationEnvironment,err)
  CALL OC_Context_ComputationEnvironmentGet(context,computationEnvironment,err)
  
  CALL OC_WorkGroup_Initialise(worldWorkGroup,err)
  CALL OC_ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err)
  CALL OC_WorkGroup_NumberOfGroupNodesGet(worldWorkGroup,numberOfComputationalNodes,err)
  CALL OC_WorkGroup_GroupNodeNumberGet(worldWorkGroup,computationalNodeNumber,err)

  !-----------------------------------------------------------------------------------------------------------
  ! COORDINATE SYSTEM
  !-----------------------------------------------------------------------------------------------------------  

  !Start the creation of a new RC coordinate system
  CALL OC_CoordinateSystem_Initialise(coordinateSystem,err)
  CALL OC_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,context,coordinateSystem,err)
  !Set the coordinate system number of dimensions
  CALL OC_CoordinateSystem_DimensionSet(coordinateSystem,numberOfDimensions,err)
  !Set the origin
  CALL OC_CoordinateSystem_OriginSet(coordinateSystem,[0.0_OC_RP,0.0_OC_RP,0.0_OC_RP],err)
  !Finish the creation of the coordinate system
  CALL OC_CoordinateSystem_CreateFinish(coordinateSystem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! REGION
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the region
  CALL OC_Region_Initialise(region,err)
  CALL OC_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region,err)
  CALL OC_Region_LabelSet(region,"NonlinearPoisson",err)
  !Set the regions coordinate system to the RC coordinate system that we have created
  CALL OC_Region_CoordinateSystemSet(region,coordinateSystem,err)
  !Finish the creation of the region
  CALL OC_Region_CreateFinish(region,err)

  !-----------------------------------------------------------------------------------------------------------
  ! BASIS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a basis (default is trilinear lagrange)
  CALL OC_Basis_Initialise(basis,err)
  CALL OC_Basis_CreateStart(BASIS_USER_NUMBER,context,basis,err)
  CALL OC_Basis_NumberOfXiSet(basis,numberOfDimensions,err)
  SELECT CASE(interpolationType)
  CASE(OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
    & OC_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
    & OC_BASIS_CUBIC_LAGRANGE_INTERPOLATION, &
    & OC_BASIS_CUBIC_HERMITE_INTERPOLATION)
    CALL OC_Basis_TypeSet(basis,OC_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
  CASE(OC_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
    & OC_BASIS_QUADRATIC_SIMPLEX_INTERPOLATION, &
    & OC_BASIS_CUBIC_SIMPLEX_INTERPOLATION)
    CALL OC_Basis_TypeSet(basis,OC_BASIS_SIMPLEX_TYPE,err)
  CASE DEFAULT
    CALL HandleError("Invalid interpolation type.")
  END SELECT
  SELECT CASE(interpolationType)
  CASE(OC_BASIS_LINEAR_LAGRANGE_INTERPOLATION)
    numberOfGaussXi=2
    numberOfNodesXi=2
    gaussOrder=0
  CASE(OC_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION)
    numberOfGaussXi=3
    numberOfNodesXi=3
    gaussOrder=0
  CASE(OC_BASIS_CUBIC_LAGRANGE_INTERPOLATION)
    numberOfGaussXi=4
    numberOfNodesXi=4
    gaussOrder=0
  CASE(OC_BASIS_CUBIC_HERMITE_INTERPOLATION)
    numberOfGaussXi=4
    numberOfNodesXi=2
    gaussOrder=0
  CASE(OC_BASIS_LINEAR_SIMPLEX_INTERPOLATION)
    numberOfGaussXi=0
    numberOfNodesXi=2    
    gaussOrder=3
  CASE(OC_BASIS_QUADRATIC_SIMPLEX_INTERPOLATION)
    numberOfGaussXi=0
    numberOfNodesXi=3
    gaussOrder=4
  CASE(OC_BASIS_CUBIC_SIMPLEX_INTERPOLATION)
    numberOfGaussXi=0
    numberOfNodesXi=4
    gaussOrder=5
  CASE DEFAULT
    CALL HandleError("Invalid interpolation type.")
  END SELECT
  IF(numberOfDimensions==1) THEN
    CALL OC_Basis_InterpolationXiSet(basis,[interpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL OC_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi],err)
    ELSE
      CALL OC_Basis_QuadratureOrderSet(basis,gaussOrder,err)
    ENDIF
  ELSE IF(numberOfDimensions==2) THEN
    CALL OC_Basis_InterpolationXiSet(basis,[interpolationType,interpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL OC_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi],err)
    ELSE
      CALL OC_Basis_QuadratureOrderSet(basis,gaussOrder,err)
    ENDIF
  ELSE
    CALL OC_Basis_InterpolationXiSet(basis,[interpolationType,interpolationType,interpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL OC_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],err)
    ELSE
      CALL OC_Basis_QuadratureOrderSet(basis,gaussOrder,err)      
    ENDIF
  ENDIF
  !Finish the creation of the basis
  CALL OC_Basis_CreateFinish(basis,err)

  !Compute the number of nodes in each global direction based on the interpolation type used.
  numberOfGlobalXNodes=(numberOfNodesXi-1)*numberOfGlobalXElements+1
  numberOfGlobalYNodes=(numberOfNodesXi-1)*numberOfGlobalYElements+1
  numberOfGlobalZNodes=(numberOfNodesXi-1)*numberOfGlobalZElements+1

  !-----------------------------------------------------------------------------------------------------------
  ! MESH
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a generated mesh in the region
  CALL OC_GeneratedMesh_Initialise(generatedMesh,err)
  CALL OC_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh,err)
  !Set up a regular x*y*z mesh
  CALL OC_GeneratedMesh_TypeSet(generatedMesh,OC_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  CALL OC_GeneratedMesh_BasisSet(generatedMesh,basis,err)
  !Define the mesh on the region
  IF(numberOfDimensions==1) THEN
    CALL OC_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH],err)
    CALL OC_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements],err)
  ELSEIF(numberOfDimensions==2) THEN
    CALL OC_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT],err)
    CALL OC_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements],err)
  ELSE
    CALL OC_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT,LENGTH],err)
    CALL OC_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements, &
      & numberOfGlobalZElements],err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL OC_Mesh_Initialise(mesh,err)
  CALL OC_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DECOMPOSITION
  !-----------------------------------------------------------------------------------------------------------

  !Create a decomposition
  CALL OC_Decomposition_Initialise(decomposition,err)
  CALL OC_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition,err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL OC_Decomposition_TypeSet(decomposition,OC_DECOMPOSITION_CALCULATED_TYPE,err)
  !Finish the decomposition
  CALL OC_Decomposition_CreateFinish(decomposition,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DECOMPOSER
  !-----------------------------------------------------------------------------------------------------------
  
  CALL OC_Decomposer_Initialise(decomposer,err)
  CALL OC_Decomposer_CreateStart(DECOMPOSER_USER_NUMBER,region,worldWorkGroup,decomposer,err)
  !Add in the decomposition
  CALL OC_Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err)
  !Finish the decomposer
  CALL OC_Decomposer_CreateFinish(decomposer,err)

  !-----------------------------------------------------------------------------------------------------------
  ! GEOMETRIC FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Start to create a default (geometric) field on the region
  CALL OC_Field_Initialise(geometricField,err)
  CALL OC_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region,geometricField,err)
  !Set the decomposition to use
  CALL OC_Field_DecompositionSet(geometricField,decomposition,err)
  !Set the domain to be used by the field components.
  DO componentIdx=1,numberOfDimensions
    CALL OC_Field_ComponentMeshComponentSet(geometricField,OC_FIELD_U_VARIABLE_TYPE,componentIdx,1,err)
  ENDDO
  !Finish creating the field
  CALL OC_Field_CreateFinish(geometricField,err)

  !Update the geometric field parameters
  CALL OC_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,err)

  !-----------------------------------------------------------------------------------------------------------
  ! MATERIAL FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create a material field for the scalars a and b and the Voigt components of the rank-2 conductivity tensor.
  CALL OC_Field_Initialise(materialsField,err)
  CALL OC_Field_CreateStart(MATERIALS_FIELD_USER_NUMBER,region,materialsField,err)
  !Set the decomposition to use
  CALL OC_Field_DecompositionSet(materialsField,decomposition,err)
  !Set the type
  CALL OC_Field_TypeSet(materialsField,OC_FIELD_MATERIAL_TYPE,err)
  CALL OC_Field_NumberOfVariablesSet(materialsField,1,err)
  CALL OC_Field_VariableTypesSet(materialsField,[OC_FIELD_U_VARIABLE_TYPE],err)
  numberOfComponents=2+numberOfDimensions+MAX(numberOfDimensions-1,0)+MAX(numberOfDimensions-2,0)
  CALL OC_Field_NumberOfComponentsSet(materialsField,OC_FIELD_U_VARIABLE_TYPE,numberOfComponents,err)
  DO componentIdx=1,numberOfComponents
    CALL OC_Field_ComponentMeshComponentSet(materialsField,OC_FIELD_U_VARIABLE_TYPE,componentIdx,1,err)
    CALL OC_Field_ComponentInterpolationSet(materialsField,OC_FIELD_U_VARIABLE_TYPE,componentIdx, &
      & OC_FIELD_CONSTANT_INTERPOLATION,err)
  ENDDO !componentIdx
  !Set associated geometric field
  CALL OC_Field_GeometricFieldSet(materialsField,geometricField,err)
  !Set the label of the field
  CALL OC_Field_VariableLabelSet(materialsField,OC_FIELD_U_VARIABLE_TYPE,"Parameters",err)
  !Finish creating the field
  CALL OC_Field_CreateFinish(materialsField,err)

  !Update material field components
  !a parameter:
  CALL OC_Field_ComponentValuesInitialise(materialsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 1,-1.0_OC_RP,err)
  !b parameter:
  CALL OC_Field_ComponentValuesInitialise(materialsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 2,0.5_OC_RP,err)
  !Rank-2 conductivity tensor diagonal components in Voigt form:
  DO componentIdx=1,numberOfDimensions
    CALL OC_Field_ComponentValuesInitialise(materialsField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
      & 2+componentIdx,1.0_OC_RP,err)
  ENDDO

  !-----------------------------------------------------------------------------------------------------------
  ! DEPENDENT FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create dependent fields
  CALL OC_Field_Initialise(dependentField,err)
  CALL OC_Field_CreateStart(DEPENDENT_FIELD_USER_NUMBER,region,dependentField,err)
  !Set the decomposition to use
  CALL OC_Field_DecompositionSet(dependentField,decomposition,err)
  !Set the type
  CALL OC_Field_TypeSet(dependentField,OC_FIELD_GENERAL_TYPE,err)
  CALL OC_Field_DependentTypeSet(dependentField,OC_FIELD_DEPENDENT_TYPE,err)
  !Set associated geometric field
  CALL OC_Field_GeometricFieldSet(dependentField,geometricField,err)
  !Two dependent variables: primary variable 'U' and secondary variable 'DelUDelN'
  CALL OC_Field_NumberOfVariablesSet(dependentField,2,err)
  CALL OC_Field_DimensionSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_SCALAR_DIMENSION_TYPE,err)
  CALL OC_Field_DimensionSet(dependentField,OC_FIELD_DELUDELN_VARIABLE_TYPE,OC_FIELD_SCALAR_DIMENSION_TYPE,err)
  CALL OC_Field_DOFOrderTypeSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_SEPARATED_COMPONENT_DOF_ORDER,err)
  CALL OC_Field_DOFOrderTypeSet(dependentField,OC_FIELD_DELUDELN_VARIABLE_TYPE,OC_FIELD_SEPARATED_COMPONENT_DOF_ORDER,err)
  !Set appropriate labels
  CALL OC_Field_VariableLabelSet(dependentField,OC_FIELD_U_VARIABLE_TYPE,"U",err)
  CALL OC_Field_VariableLabelSet(dependentField,OC_FIELD_DELUDELN_VARIABLE_TYPE,"DelUDelN",err)
  !Finish creating the field
  CALL OC_Field_CreateFinish(dependentField,err)

  !Update dependent fields (initial guess for nonlinear solver)
  CALL OC_Field_ComponentValuesInitialise(dependentField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
    & 1,-5.0_OC_RP,err)

  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS SETS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations_sets
  CALL OC_EquationsSet_Initialise(equationsSet,err)
  CALL OC_Field_Initialise(equationsSetField,err)
  !Poisson equation with exponential (nonlinear) source term
  CALL OC_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,geometricField,[OC_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & OC_EQUATIONS_SET_POISSON_EQUATION_TYPE,OC_EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE], &
    & EQUATIONS_SET_FIELD_USER_NUMBER,equationsSetField,equationsSet,err)
  !Finish creating eqiations sets
  CALL OC_EquationsSet_CreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DEPENDENT FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Attach dependent field already created to equation sets
  CALL OC_EquationsSet_DependentCreateStart(equationsSet,DEPENDENT_FIELD_USER_NUMBER,dependentField,err)
  CALL OC_EquationsSet_DependentCreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! MATERIAL FIELD
  !-----------------------------------------------------------------------------------------------------------
  
  !Attach material field already created to equation sets
  CALL OC_EquationsSet_MaterialsCreateStart(equationsSet,MATERIALS_FIELD_USER_NUMBER,materialsField,err)
  CALL OC_EquationsSet_MaterialsCreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! ANALYTIC FIELD
  !-----------------------------------------------------------------------------------------------------------
  
  !Create the equations set analytic field variables
  IF(numberOfDimensions==2) THEN
    CALL OC_Field_Initialise(analyticField,err)
    CALL OC_EquationsSet_AnalyticCreateStart(equationsSet,OC_EQUATIONS_SET_EXPONENTIAL_POISSON_EQUATION_TWO_DIM_1, &
      & ANALYTIC_FIELD_USER_NUMBER,analyticField,Err)
    !Finish the equations set analytic field variables
    CALL OC_EquationsSet_AnalyticCreateFinish(equationsSet,err)
    !Set the analytic solution constants
    !A
    CALL OC_Field_ComponentValuesInitialise(analyticField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
      & 1,PI/(100.0_OC_RP*LENGTH),err)
    !B
    CALL OC_Field_ComponentValuesInitialise(analyticField,OC_FIELD_U_VARIABLE_TYPE,OC_FIELD_VALUES_SET_TYPE, &
      & 2,PI/(100.0_OC_RP*HEIGHT),err)
  ENDIF
  
  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set equations
  CALL OC_Equations_Initialise(equations,err)
  CALL OC_EquationsSet_EquationsCreateStart(equationsSet,equations,err)
  !Set the equations matrices sparsity type
  !CALL OC_Equations_SparsityTypeSet(equations,OC_EQUATIONS_FULL_MATRICES,err)
  CALL OC_Equations_SparsityTypeSet(equations,OC_EQUATIONS_SPARSE_MATRICES,err)
  !Set the equations set output
  !CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_NO_OUTPUT,err)
  CALL OC_Equations_OutputTypeSet(equations,OC_EQUATIONS_ELEMENT_MATRIX_OUTPUT,err)
  !Finish the equations set equations
  CALL OC_EquationsSet_EquationsCreateFinish(equationsSet,err)

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM
  !-----------------------------------------------------------------------------------------------------------  

  !Start the creation of a problem.
  CALL OC_Problem_Initialise(problem,err)
  CALL OC_Problem_CreateStart(PROBLEM_USER_NUMBER,context,[OC_PROBLEM_CLASSICAL_FIELD_CLASS, &
    & OC_PROBLEM_POISSON_EQUATION_TYPE,OC_PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE],problem,err)
  !Finish the creation of a problem.
  CALL OC_Problem_CreateFinish(problem,err)
  
  !-----------------------------------------------------------------------------------------------------------
  ! CONTROL LOOP
  !-----------------------------------------------------------------------------------------------------------  

  !Start the creation of the problem control loop
  CALL OC_Problem_ControlLoopCreateStart(problem,err)
  !Finish creating the problem control loop
  CALL OC_Problem_ControlLoopCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVER
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the problem solvers
  CALL OC_Solver_Initialise(solver,err)
  CALL OC_Solver_Initialise(linearSolver,err)
  CALL OC_Problem_SolversCreateStart(problem,err)
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,solver,err)
  !Set the solver output
  !CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_NO_OUTPUT,err)
  CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_PROGRESS_OUTPUT,err)
  !CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_TIMING_OUTPUT,err)
  !CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_SOLVER_OUTPUT,err)
  CALL OC_Solver_OutputTypeSet(solver,OC_SOLVER_MATRIX_OUTPUT,err)
  !Set the Jacobian type to be either calculated analytically or via finite differences
  CALL OC_Solver_NewtonJacobianCalculationTypeSet(solver,OC_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,err)
  !CALL OC_Solver_NewtonJacobianCalculationTypeSet(solver,OC_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,err)
  !Set the nonlinear Newton solver tolerances etc.
  CALL OC_Solver_NewtonAbsoluteToleranceSet(solver,1.0E-8_OC_RP,err)
  CALL OC_Solver_NewtonRelativeToleranceSet(solver,1.0E-8_OC_RP,err)
  CALL OC_Solver_NewtonMaximumIterationsSet(solver,100000,err)
  !Get the associated linear solver
  CALL OC_Solver_NewtonLinearSolverGet(solver,linearSolver,err)
  !Set the linear iterative solver tolerances etc.
  CALL OC_Solver_LinearIterativeRelativeToleranceSet(linearSolver,1.0E-8_OC_RP,err)
  CALL OC_Solver_LinearIterativeAbsoluteToleranceSet(linearSolver,1.0E-8_OC_RP,err)
  CALL OC_Solver_LinearIterativeMaximumIterationsSet(linearSolver,10000,err)
  !Finish the creation of the problem solver
  CALL OC_Problem_SolversCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVER EQUATIONS
  !-----------------------------------------------------------------------------------------------------------  

  !Start the creation of the problem solver equations
  CALL OC_Solver_Initialise(solver,err)
  CALL OC_SolverEquations_Initialise(solverEquations,err)
  CALL OC_Problem_SolverEquationsCreateStart(problem,err)
  !Get the solve equations
  CALL OC_Problem_SolverGet(problem,OC_CONTROL_LOOP_NODE,1,solver,err)
  CALL OC_Solver_SolverEquationsGet(solver,solverEquations,err)
  !Set the solver equations sparsity
  !CALL OC_SolverEquations_SparsityTypeSet(solverEquations,OC_SOLVER_SPARSE_MATRICES,err)
  CALL OC_SolverEquations_SparsityTypeSet(solverEquations,OC_SOLVER_FULL_MATRICES,err)
  !Add in the equations set
  CALL OC_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err)
  !Finish the creation of the problem solver equations
  CALL OC_Problem_SolverEquationsCreateFinish(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! BOUNDARY CONDITIONS
  !-----------------------------------------------------------------------------------------------------------

  !Set up the boundary conditions
  CALL OC_BoundaryConditions_Initialise(boundaryConditions,err)
  CALL OC_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err)
  IF(numberOfDimensions==2) THEN
    !Use analytic boundary conditions
    CALL OC_SolverEquations_BoundaryConditionsAnalytic(solverEquations,err)
  ELSE
    !Set the fixed boundary conditions on opposide sides
    DO zNodeIdx=1,numberOfGlobalZNodes
      DO yNodeIdx=1,numberOfGlobalYNodes
        !Left side
        nodeNumber=1+(yNodeIdx-1)*numberOfGlobalXNodes+(zNodeIdx-1)*numberOfGlobalXNodes*numberOfGlobalYNodes
        !Check what computational node the node is on
        CALL OC_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
        IF(nodeDomain==computationalNodeNumber) THEN
          !If the node is on my computational node then set left hand side BC to 0.0
          CALL OC_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,1, &
            & OC_BOUNDARY_CONDITION_FIXED,0.0_OC_RP,err)
        ENDIF
        !Right side
        nodeNumber=numberOfGlobalXNodes+(yNodeIdx-1)*numberOfGlobalXNodes+(zNodeIdx-1)*numberOfGlobalXNodes*numberOfGlobalYNodes
        !Check what computational node the node is on
        CALL OC_Decomposition_NodeDomainGet(decomposition,nodeNumber,1,nodeDomain,err)
        IF(nodeDomain==computationalNodeNumber) THEN
          !If the node is on my computational node then set left hand side BC to 0.0
          CALL OC_BoundaryConditions_SetNode(boundaryConditions,dependentField,OC_FIELD_U_VARIABLE_TYPE,1,1,nodeNumber,1, &
            & OC_BOUNDARY_CONDITION_FIXED,1.0_OC_RP,err)
        ENDIF
      ENDDO !yNodeIdx
    ENDDO !zNodeIdx
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL OC_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVE
  !-----------------------------------------------------------------------------------------------------------

  !Export results
  exportField=.TRUE.
  IF(exportField) THEN
    CALL OC_Fields_Initialise(fields,err)
    CALL OC_Fields_Create(region,fields,err)
    CALL OC_Fields_NodesExport(fields,"NonlinearPoisson","FORTRAN",err)
    CALL OC_Fields_ElementsExport(fields,"NonlinearPoisson","FORTRAN",err)
    CALL OC_Fields_Finalise(fields,err)
  ENDIF
  
  !Solve the problem
  CALL OC_Problem_Solve(problem,err)

  !-----------------------------------------------------------------------------------------------------------
  ! OUTPUT
  !-----------------------------------------------------------------------------------------------------------

  IF(numberOfDimensions==2) THEN
    !Output Analytic analysis
    CALL OC_AnalyticAnalysis_Output(dependentField,"NonlinearPoissonAnalytic",err)
  ENDIF
  
  !Export results
  exportField=.TRUE.
  IF(exportField) THEN
    CALL OC_Fields_Initialise(fields,err)
    CALL OC_Fields_Create(region,fields,err)
    CALL OC_Fields_NodesExport(fields,"NonlinearPoisson","FORTRAN",err)
    CALL OC_Fields_ElementsExport(fields,"NonlinearPoisson","FORTRAN",err)
    CALL OC_Fields_Finalise(fields,err)
  ENDIF

  !-----------------------------------------------------------------------------------------------------------
  ! FINALISE
  !-----------------------------------------------------------------------------------------------------------  

  !Destroy the context
  CALL OC_Context_Destroy(context,err)
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

END PROGRAM NonlinearPoissonEquation
