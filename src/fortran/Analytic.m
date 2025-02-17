global L = 1.0;
global H = 1.0;

global sigma = [ 1.0, 0.0; 0.0, 1.0 ];
global a = -1.0;
global b = 0.5;

global A;
global B;
global C = 0.0;
global alpha = -a;
global beta = b;

global numXElements = 2;
global numYElements = 2;

global numElements = numXElements * numYElements; 
global numNodes = (numXElements + 1) * (numYElements + 1);


function x = computeGeometry( )

  global L;
  global H;
  global numXElements;
  global numYElements;
  global numNodes;
  global x = zeros(numNodes,1);

  nodeIdx = 0;
  for j = 1:numYElements+1
    for i = 1:numXElements+1
      nodeIdx = nodeIdx + 1;
      x(nodeIdx,1) = (i-1)*L/numXElements;
    endfor
    x(nodeIdx,2) = (j-1)*H/numYElements;
  endfor

endfunction	    

function [ K, f ] = computeLinearMatrices( )

endfunction

function [ analyticU ] = analytic( x )

  global A;
  global B;
  global C;
  global alpha;
  global beta;
  global numNodes;

  for i = 1:numNodes
    analyticU(i) = (1.0/beta)*ln(2.0*(A*A+B*B)/(alpha*beta*cos(A*x(i,1)+B*x(i,2)+C)*cos(A*x(i,1)+B*x(i,2)+C))
  endfor
  
endfunction

function [ err, percentErr, rmsErr ] = error( u, analyticU )

  [ numRows, numCols ] = size( u );
  err = zeros(numRows,1);
  percentErr = zeros(numRows,1);
  sum = 0.0;
  printf("  row    value analytic      err    %%err \n");
  for i = 1:numRows
    err(i) = (u(i)-analyticU(i));
    if( abs(analyticU(i)) > 0.00000001 )
      percentErr(i) = 100.0*(u(i)-analyticU(i))/analyticU(i);
    else
      percentErr(i) = 0.0;
    endif
    sum = sum + err(i)*err(i);
    printf("%5d %8.5f %8.5f %8.5f %7.2f\n",i,u(i),analyticU(i),err(i),percentErr(i));
  endfor
  rmsErr = sqrt(sum/numRows)

endfunction


function [ reducedR, reducedJ ] = reducedFunction ( reducedAlpha )

   global alpha;
   global numNodes;

   myAlpha = zeros(numNodes,1); 
  
   myAlpha(1) = alpha(1);
   myAlpha(2:numNodes-1) = reducedAlpha(1:numNodes-2);
   myAlpha(numNodes) = alpha(numNodes)

   r = Residual( myAlpha )

   reducedR = r(2:numNodes-1)

   normReducedR = norm(reducedR)

   if (nargout == 2)
   
      J = Jacobian( myAlpha )

      reducedJ = J(2:numNodes-1,2:numNodes-1)

   endif

endfunction


global x = zeros(numNodes,2);
global u = zeros(numNodes,1);
global analyticU = zeros(numNodes,1);
global K = zeros(numNodes,numNodes);
global resid = zeros(numNodes,1);

x = computeGeometry( );

[ reducedAlpha, fval, info ] = fsolve( @reducedFunction, startAlpha );
#[ reducedAlpha, fval, info ] = fsolve( @reducedFunction, startAlpha, optimset ("jacobian", "on") )
