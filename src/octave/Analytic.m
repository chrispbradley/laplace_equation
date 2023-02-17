# Constant declarations

global LINEAR_LAGRANGE_INTERPOLATION = 1;
global QUADRATIC_LAGRANGE_INTERPOLATION = 2;
global CUBIC_LAGRANGE_INTERPOLATION = 3;

global OFF_BOUNDARY_NODE = 0;
global ON_BOUNDARY_NODE = 1;

global DOF_FIXED = 0;
global DOF_FREE = 1;

global STANDARD_LAPLACE_1D_ANALYTIC_TYPE_1 = 1;
global STANDARD_LAPLACE_2D_ANALYTIC_TYPE_1 = 2;
global STANDARD_LAPLACE_2D_ANALYTIC_TYPE_2 = 3;
global STANDARD_LAPLACE_3D_ANALYTIC_TYPE_1 = 4;

# Example parameters and options

global diagnosticsLevel = 1;

global analyticType = STANDARD_LAPLACE_2D_ANALYTIC_TYPE_2;

#global L = 1.0;
global L = 2.0;
global H = 1.0;
global W = 1.0;

global baseNumXElements = 3;
global baseNumYElements = 2;
global baseNumZElements = 0;

global refinementLevel = 5;

interpolationType = LINEAR_LAGRANGE_INTERPOLATION;

# Should not need to change below here

global numDimensions = gt( baseNumXElements, 0 ) + gt( baseNumYElements, 0 ) + gt( baseNumZElements, 0 );

global numXElements = baseNumXElements * refinementLevel
global numYElements = baseNumYElements * refinementLevel
global numZElements = baseNumZElements * refinementLevel
global numElements = numXElements * eq( numDimensions, 1 ) + numXElements * numYElements * eq( numDimensions, 2 ) + numXElements * numYElements * numZElements * eq( numDimensions, 3 )
global numNodesPerXi = (interpolationType + 1)
global numElementDOFS = numNodesPerXi ** numDimensions
global numXNodes = (numNodesPerXi-1)*numXElements + 1
global numYNodes = (numNodesPerXi-1)*numYElements + 1
global numZNodes = (numNodesPerXi-1)*numZElements + 1
global numNodes = numXNodes*numYNodes*numZNodes
global numDOFS = numNodes
global numFixedDOFS = 2 * eq( numDimensions, 1 ) + (2 * numXNodes + 2 * (numYNodes-2)) * eq( numDimensions, 2 ) + (2 * numXNodes * numYNodes + 4 * numYNodes * (numZNodes - 2)) * eq( numDimensions, 3 );
global numFreeDOFS = numDOFS - numFixedDOFS;

function [ nodePositions, boundaryNodes ] = computeGeometry( interpolationType )

  global DIAGNOSTICS_OFF;
  global diagnosticsLevel;
  global ON_BOUNDARY_NODE;
  global OFF_BOUNDARY_NODE;
  global LINEAR_LAGRANGE_INTERPOLATION;
  global QUADRATIC_LAGRANGE_INTERPOLATION;
  global CUBIC_LAGRANGE_INTERPOLATION;
  global L;
  global H;
  global W;
  global numDimensions;
  global numXElements;
  global numYElements;
  global numZElements;

  nodeIdx = 0;
  switch( numDimensions )
    case 1
      for i = 1:numXElements+1
	nodeIdx = nodeIdx + 1;
	boundaryNodes(nodeIdx,1) = OFF_BOUNDARY_NODE;
	nodePositions(nodeIdx,1) = (i-1)*L/numXElements;
	if( or( i==1, i==numXElements+1 ) )
          boundaryNodes(nodeIdx,1) = ON_BOUNDARY_NODE;
	endif
      endfor
    case 2
      for j = 1:numYElements+1
	for i = 1:numXElements+1
	  nodeIdx = nodeIdx + 1;
	  boundaryNodes(nodeIdx,1) = OFF_BOUNDARY_NODE;
	  nodePositions(nodeIdx,1) = (i-1)*L/numXElements;
	  nodePositions(nodeIdx,2) = (j-1)*H/numYElements;
	  if( or( i==1, i==numXElements+1, j==1, j==numYElements+1 ) )
            boundaryNodes(nodeIdx,1) = ON_BOUNDARY_NODE;
	  endif
	endfor
      endfor
    case 3
      for k = 1:numZElements+1
	for j = 1:numYElements+1
	  for i = 1:numXElements+1
	    nodeIdx = nodeIdx + 1;
	    boundaryNodes(nodeIdx,1) = OFF_BOUNDARY_NODE;
	    nodePositions(nodeIdx,1) = (i-1)*L/numXElements;
	    nodePositions(nodeIdx,2) = (j-1)*H/numYElements;
	    nodePositions(nodeIdx,3) = (k-1)*W/numZElements;
	    if( or( i==1, i==numXElements+1, j==1, j==numYElements+1, k==1, k==numZElements+1 ) )
              boundaryNodes(nodeIdx,1) = ON_BOUNDARY_NODE;
	    endif
	  endfor
	endfor
      endfor
    otherwise
      printf("ERROR: The number of dimensions of %d is invalid.\n",numDimensions);
  endswitch

  if( diagnosticsLevel > DIAGNOSTICS_OFF )
    nodePositions
    boundaryNodes
  endif

endfunction

function [ elementDOFNumbers ] = computeElementDOFs( interpolationType, elementIdx )

  global DIAGNOSTICS_OFF;
  global diagnosticsLevel;
  global LINEAR_LAGRANGE_INTERPOLATION;
  global QUADRATIC_LAGRANGE_INTERPOLATION;
  global CUBIC_LAGRANGE_INTERPOLATION;
  global numDimensions;
  global numXNodes;
  global numYNodes;
  global numZNodes;

  switch( interpolationType )
    case LINEAR_LAGRANGE_INTERPOLATION
      switch( numDimensions )
	case 1
	  elementDOFNumbers(1) = elementIdx(1);
	  elementDOFNumbers(2) = 1 + elementIdx(1);
	case 2
	  elementDOFNumbers(1) = elementIdx(1) + (elementIdx(2) - 1) * numXNodes;
	  elementDOFNumbers(2) = 1 + elementIdx(1) + (elementIdx(2) - 1) * numXNodes;
	  elementDOFNumbers(3) = elementIdx(1) + elementIdx(2) * numXNodes;
	  elementDOFNumbers(4) = 1 + elementIdx(1) + elementIdx(2) * numXNodes;
	case 3
	  elementDOFNumbers(1) = elementIdx(1) + (elementIdx(2) - 1) * numXNodes + (elementIdx(3) - 1) * (numXNodes * numYNodes);
	  elementDOFNumbers(2) = 1 + elementIdx(1) + (elementIdx(2) - 1) * numXNodes + (elementIdx(3) - 1) * (numXNodes * numYNodes);
	  elementDOFNumbers(3) = elementIdx(1) + elementIdx(2) * numXNodes + (elementIdx(3) - 1) * (numXNodes * numYNodes);
	  elementDOFNumbers(4) = 1 + elementIdx(1) + elementIdx(2) * numXNodes + (elementIdx(3) - 1) * (numXNodes * numYNodes);
	  elementDOFNumbers(5) = elementIdx(1) + (elementIdx(2) - 1) * numXNodes + elementIdx(3) * (numXNodes * numYNodes);
	  elementDOFNumbers(6) = 1 + elementIdx(1) + (elementIdx(2) - 1) * numXNodes + elementIdx(3) * (numXNodes * numYNodes);
	  elementDOFNumbers(7) = elementIdx(1) + elementIdx(2) * numXNodes + elementIdx(3) * (numXNodes * numYNodes);
	  elementDOFNumbers(8) = 1 + elementIdx(1) + elementIdx(2) * numXNodes + elementIdx(3) * (numXNodes * numYNodes);
	otherwise
	  printf("ERROR: The number of dimensions of %d is invalid.\n",numDimensions);
      endswitch
    case QUADRATIC_LAGRANGE_INTERPOLATION
      printf("ERROR: Quadratic Lagrange interpolation not implemented.\n");
      quit;
    case CUBIC_LAGRANGE_INTERPOLATION
      printf("ERROR: Cubic Lagrange interpolation not implemented.\n");
      quit;
    otherwise
      printf("ERROR: The interpolation type of %d is invalid.\n",interpolationType);
      quit;
  endswitch	  

  if( diagnosticsLevel > DIAGNOSTICS_OFF )
    elementDOFNumbers
  endif
  
endfunction

function [ Ke, fe ] = computeElementStiffnessMatrix( interpolationType )

  global DIAGNOSTICS_OFF;
  global diagnosticsLevel;
  global LINEAR_LAGRANGE_INTERPOLATION;
  global QUADRATIC_LAGRANGE_INTERPOLATION;
  global CUBIC_LAGRANGE_INTERPOLATION;
  global L;
  global H;
  global W;
  global numDimensions;
  global numXElements;
  global numYElements;
  global numZElements;

  switch( interpolationType )
    case LINEAR_LAGRANGE_INTERPOLATION
      switch( numDimensions )
	case 1
	  printf("ERROR: Linear Lagrange interpolation not implemented.\n");
	case 2      
	  Ke(1,1) = ( 2.0*H*H*numXElements*numXElements + 2.0*L*L*numYElements*numYElements)/(6.0*L*H*numXElements*numYElements);
	  Ke(1,2) = (-2.0*H*H*numXElements*numXElements + 1.0*L*L*numYElements*numYElements)/(6.0*L*H*numXElements*numYElements);
	  Ke(1,3) = ( 1.0*H*H*numXElements*numXElements - 2.0*L*L*numYElements*numYElements)/(6.0*L*H*numXElements*numYElements);
	  Ke(1,4) = (-1.0*H*H*numXElements*numXElements - 1.0*L*L*numYElements*numYElements)/(6.0*L*H*numXElements*numYElements);
	  Ke(2,1) = Ke(1,2);
	  Ke(2,2) = ( 2.0*H*H*numXElements*numXElements + 2.0*L*L*numYElements*numYElements)/(6.0*L*H*numXElements*numYElements);
	  Ke(2,3) = (-1.0*H*H*numXElements*numXElements - 1.0*L*L*numYElements*numYElements)/(6.0*L*H*numXElements*numYElements);
	  Ke(2,4) = ( 1.0*H*H*numXElements*numXElements - 2.0*L*L*numYElements*numYElements)/(6.0*L*H*numXElements*numYElements);
	  Ke(3,1) = Ke(1,3);
	  Ke(3,2) = Ke(2,3);
	  Ke(3,3) = ( 2.0*H*H*numXElements*numXElements + 2.0*L*L*numYElements*numYElements)/(6.0*L*H*numXElements*numYElements);
	  Ke(3,4) = (-2.0*H*H*numXElements*numXElements + 1.0*L*L*numYElements*numYElements)/(6.0*L*H*numXElements*numYElements);
	  Ke(4,1) = Ke(1,4);
	  Ke(4,2) = Ke(2,4);
	  Ke(4,3) = Ke(3,4);
	  Ke(4,4) = ( 2.0*H*H*numXElements*numXElements + 2.0*L*L*numYElements*numYElements)/(6.0*L*H*numXElements*numYElements);
	  fe(1,1) = 0.0;
	  fe(2,1) = 0.0;
	  fe(3,1) = 0.0;
	  fe(4,1) = 0.0;
	case 3
	  printf("ERROR: Tri-Linear Lagrange interpolation not implemented.\n");
	otherwise
	  printf("ERROR: The number of dimensions of %d is invalid.\n",numDimensions);
      endswitch
    case QUADRATIC_LAGRANGE_INTERPOLATION
      switch( numDimensions )
	case 1
	  printf("ERROR: Quadratic Lagrange interpolation not implemented.\n");
	case 2      
	  printf("ERROR: Bi-Quadratic Lagrange interpolation not implemented.\n");
	case 3
	  printf("ERROR: Tri-Quadric Lagrange interpolation not implemented.\n");
	otherwise
	  printf("ERROR: The number of dimensions of %d is invalid.\n",numDimensions);
      endswitch
      quit;
    case CUBIC_LAGRANGE_INTERPOLATION
      switch( numDimensions )
	case 1
	  printf("ERROR: Cubic Lagrange interpolation not implemented.\n");
	case 2      
	  printf("ERROR: Bi-Cubic Lagrange interpolation not implemented.\n");
	case 3
	  printf("ERROR: Tri-Cubic Lagrange interpolation not implemented.\n");
	otherwise
	  printf("ERROR: The number of dimensions of %d is invalid.\n",numDimensions);
      endswitch
      quit;
    otherwise
      printf("ERROR: The interpolation type of %d is invalid.\n",interpolationType);
      quit;
  endswitch	  

  if( diagnosticsLevel > DIAGNOSTICS_OFF )
    Ke
    fe
  endif
  
endfunction

function [ K, f ] = computeLinearMatrices( interpolationType )

  global DIAGNOSTICS_OFF;
  global diagnosticsLevel;
  global numDimensions;
  global numXElements;
  global numYElements;
  global numZElements;
  global numElementDOFS;
  global numDOFS;

  K = zeros(numDOFS,numDOFS);
  f = zeros(numDOFS,1);

  Ke = zeros(numElementDOFS,numElementDOFS);
  fe = zeros(numElementDOFS,1);
  elementDOFNumbers = zeros(numElementDOFS,1);
  
  [ Ke, fe ] = computeElementStiffnessMatrix( interpolationType )

  elementIdx = 0;
  for zElementIdx = 1:max(numZElements,1)
    for yElementIdx = 1:max(numYElements,1)
      for xElementIdx = 1:max(numXElements,1)
	elementIdx=elementIdx + 1;
	[ elementDOFNumbers ] = computeElementDOFs( interpolationType, [ xElementIdx, yElementIdx, zElementIdx ] );
	for rowDOFIdx = 1:numElementDOFS
	  rowNumber = elementDOFNumbers( rowDOFIdx );
	  for columnDOFIdx = 1:numElementDOFS
	    columnNumber = elementDOFNumbers( columnDOFIdx );
	    K(rowNumber,columnNumber) = K(rowNumber,columnNumber) + Ke(rowDOFIdx,columnDOFIdx);	
	  endfor
	  f(rowNumber,1) = f(rowNumber,1) + fe(rowDOFIdx,1);
	endfor
      endfor
    endfor
  endfor
  
  if( diagnosticsLevel > DIAGNOSTICS_OFF )
    K
    f
  endif
  
endfunction

function [ analyticU, gradAnalyticU, hessAnalyticU ] = analytic( analyticType, nodePositions )

  global DIAGNOSTICS_OFF;
  global diagnosticsLevel;
  global GENERALISED_LAPLACE_ANALYTIC_TYPE;
  global STANDARD_LAPLACE_1D_ANALYTIC_TYPE_1;
  global STANDARD_LAPLACE_2D_ANALYTIC_TYPE_1;
  global STANDARD_LAPLACE_2D_ANALYTIC_TYPE_2;
  global STANDARD_LAPLACE_3D_ANALYTIC_TYPE_1;
  global numNodes;

  switch( analyticType )
    case STANDARD_LAPLACE_1D_ANALYTIC_TYPE_1
      print("ERROR: Standard Laplace 1D analytic type 1 is not implemented.\n");
      quit;
    case STANDARD_LAPLACE_2D_ANALYTIC_TYPE_1
      for nodeIdx = 1:numNodes
	analyticU(nodeIdx,1) = nodePositions(nodeIdx,1)*nodePositions(nodeIdx,1)+2*nodePositions(nodeIdx,1)*nodePositions(nodeIdx,2)-nodePositions(nodeIdx,2)*nodePositions(nodeIdx,2);
	gradAnalyticU(nodeIdx,1) = 2.0*nodePositions(nodeIdx,1)+2.0*nodePositions(nodeIdx,2);
	gradAnalyticU(nodeIdx,1) = 2.0*nodePositions(nodeIdx,1)-2.0*nodePositions(nodeIdx,2);
	hessAnalyticU(nodeIdx,1,1) = 0.0;
	hessAnalyticU(nodeIdx,1,2) = 0.0;
	hessAnalyticU(nodeIdx,2,1) = 0.0;
	hessAnalyticU(nodeIdx,2,2) = 0.0;
      endfor
    case STANDARD_LAPLACE_2D_ANALYTIC_TYPE_2
      for nodeIdx = 1:numNodes
	analyticU(nodeIdx,1) = cos(nodePositions(nodeIdx,1))*cosh(nodePositions(nodeIdx,2));
	gradAnalyticU(nodeIdx,1) = -sin(nodePositions(nodeIdx,1))*cosh(nodePositions(nodeIdx,2));
	gradAnalyticU(nodeIdx,1) = cos(nodePositions(nodeIdx,1))*sinh(nodePositions(nodeIdx,2));
	hessAnalyticU(nodeIdx,1,1) = -cos(nodePositions(nodeIdx,1))*cosh(nodePositions(nodeIdx,2));
	hessAnalyticU(nodeIdx,1,2) = -sin(nodePositions(nodeIdx,1))*sinh(nodePositions(nodeIdx,2));
	hessAnalyticU(nodeIdx,2,1) = -sin(nodePositions(nodeIdx,1))*sinh(nodePositions(nodeIdx,2));
	hessAnalyticU(nodeIdx,2,2) = cos(nodePositions(nodeIdx,1))*cosh(nodePositions(nodeIdx,2));
      endfor
    case STANDARD_LAPLACE_3D_ANALYTIC_TYPE_1
      print("ERROR: Standard Laplace 3D analytic type 1 is not implemented.\n");
      quit;
    otherwise
      printf("ERROR: The analytic type of %d is invalid.\n",analyticType);
      quit;
  endswitch
  
  if( diagnosticsLevel > DIAGNOSTICS_OFF )
    analyticU
    #gradAnalyticU
    #hessAnalyticU
  endif
  
endfunction

function [ u, dofMap, freeDOFS, fixedDOFS ] = setBoundaryConditions( boundaryNodes, analyticU )

  global DIAGNOSTICS_OFF;
  global diagnosticsLevel;
  global ON_BOUNDARY_NODE;
  global DOF_FIXED;
  global DOF_FREE;
  global numNodes;

  dofIdx = 0;
  freeDOFIdx = 0;
  fixedDOFIdx = 0;
  for nodeIdx = 1:numNodes
    dofIdx = dofIdx + 1;
    if( boundaryNodes( nodeIdx ) == ON_BOUNDARY_NODE )
      fixedDOFIdx = fixedDOFIdx + 1;
      dofMap( dofIdx, 1 ) = DOF_FIXED;
      dofMap( dofIdx, 2 ) = fixedDOFIdx;
      fixedDOFS( fixedDOFIdx ) = dofIdx;
      u( nodeIdx, 1 ) = analyticU( nodeIdx, 1 );
    else
      freeDOFIdx = freeDOFIdx + 1;
      dofMap( dofIdx, 1 ) = DOF_FREE;
      dofMap( dofIdx, 2 ) = freeDOFIdx;
      freeDOFS( freeDOFIdx ) = dofIdx;
      u( nodeIdx, 1 ) = 0.0;
    endif
  endfor
  
  if( diagnosticsLevel > DIAGNOSTICS_OFF )
    u
    dofMap
    freeDOFS
    fixedDOFS
  endif
  
endfunction

function [ A, b ] = reduceGlobalSystem( dofMap, freeDOFS, fixedDOFS, K, u, f )
  
  global DIAGNOSTICS_OFF;
  global diagnosticsLevel;
  global DOF_FIXED;
  global DOF_FREE;
  global numDOFS;
  global numFreeDOFS;
  global numFixedDOFS;
  
  for rowDOFIdx = 1:numFreeDOFS
    globalRowNumber = freeDOFS( rowDOFIdx );
    # Form reduced A matrix
    for columnDOFIdx = 1:numFreeDOFS
      globalColumnNumber = freeDOFS( columnDOFIdx );
      A( rowDOFIdx, columnDOFIdx ) = K( globalRowNumber, globalColumnNumber );      
    endfor
    # Form RHS vector
    b( rowDOFIdx, 1 ) = f( globalRowNumber, 1 );
    for columnDOFIdx = 1:numFixedDOFS
      globalColumnNumber = fixedDOFS( columnDOFIdx );
      b( rowDOFIdx, 1 ) = b( rowDOFIdx, 1 ) - K( globalRowNumber, globalColumnNumber )*u( globalColumnNumber, 1 );
    endfor
  endfor

  if( diagnosticsLevel > DIAGNOSTICS_OFF )
    A
    b
    u
  endif
  
endfunction

function [ x ] = solveSystem( A, b );

  global DIAGNOSTICS_OFF;
  global diagnosticsLevel;

  x = linsolve( A, b );

  if( diagnosticsLevel > DIAGNOSTICS_OFF )
    x
  endif

endfunction

function [ u, f ] = updateDOFValues( dofMap, freeDOFS, fixedDOFS, x, K, u, f )
  
  global DIAGNOSTICS_OFF;
  global diagnosticsLevel;
  global DOF_FIXED;
  global DOF_FREE;
  global numDOFS;
  global numFreeDOFS;
  global numFixedDOFS;
  global tempRHS=zeros(numDOFS,1)

  # Update u with solved values
  for dofIdx = 1:numFreeDOFS
    globalDOFIdx = freeDOFS( dofIdx );
    u( globalDOFIdx, 1 ) = x( dofIdx, 1 );
    f( globalDOFIdx, 1 ) = f( globalDOFIdx, 1 );
  endfor
  # Backsubstitute
  tempRHS = K*u;
  for dofIdx = 1:numFixedDOFS
    globalDOFIdx = fixedDOFS( dofIdx );
    u( globalDOFIdx, 1 ) = u( globalDOFIdx );
    f( globalDOFIdx, 1 ) = -1.0*tempRHS( globalDOFIdx, 1 );
  endfor
  
  if( diagnosticsLevel > DIAGNOSTICS_OFF )
    u
    f
  endif
  
endfunction

function [ err, percentErr, rmsErr ] = computeErrors( u, analyticU )

  global DIAGNOSTICS_OFF;
  global diagnosticsLevel;
  global numDOFS;
  
  sum = 0.0;
  printf("  row    value analytic      err    %%err \n");
  for dofIdx = 1:numDOFS
    err(dofIdx,1) = (u(dofIdx,1)-analyticU(dofIdx,1));
    if( abs(analyticU(dofIdx,1)) > 0.00000001 )
      percentErr(dofIdx,1) = 100.0*(u(dofIdx,1)-analyticU(dofIdx,1))/analyticU(dofIdx,1);
    else
      percentErr(dofIdx,1) = 0.0;
    endif
    sum = sum + err(dofIdx,1)*err(dofIdx,1);
    printf("%5d %8.5f %8.5f %8.5f %7.2f\n",dofIdx,u(dofIdx,1),analyticU(dofIdx,1),err(dofIdx,1),percentErr(dofIdx,1));
  endfor
  rmsErr = sqrt(sum/numDOFS)
  
  if( diagnosticsLevel > DIAGNOSTICS_OFF )
  endif

endfunction

# Start of main code

global nodePositions = zeros(numNodes,numDimensions);
global boundaryNodes = zeros(numNodes,1);
global u = zeros(numNodes,1);
global analyticU = zeros(numNodes,1);
global gradAnalyticU = zeros(numNodes,numDimensions);
global hessAnalyticU = zeros(numNodes,numDimensions,numDimensions);
global K = zeros(numDOFS,numDOFS);
global f = zeros(numDOFS,1);
global dofMap = zeros(numDOFS,2);
global freeDOFS = zeros(numFreeDOFS,1);
global fixedDOFS = zeros(numFixedDOFS,1);
global A = zeros(numFreeDOFS,numFreeDOFS);
global x = zeros(numFreeDOFS,1);
global b = zeros(numFreeDOFS,1);
global uError = zeros(numNodes,1);
global uPerError = zeros(numNodes,1);
global uRMSError = 0.0;


[ nodePositions, boundaryNodes ] = computeGeometry( interpolationType );
[ analyticU, gradAnalyticU, hessAnalyticU ] = analytic( analyticType, nodePositions );
[ K, f ] = computeLinearMatrices( interpolationType );
[ u, dofMap, freeDOFS, fixedDOFS ] = setBoundaryConditions( boundaryNodes, analyticU );
[ A, b ] = reduceGlobalSystem( dofMap, freeDOFS, fixedDOFS, K, u, f );
[ x ] = solveSystem( A, b );
[ u, f ] = updateDOFValues( dofMap, freeDOFS, fixedDOFS, x, K, u, f );
[ err, percentErr, rmsErr ] = computeErrors( u, analyticU );

