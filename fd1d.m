function [u,h,xx,dirval1,dirval2,qflux1,qflux2] = ...
    fd1d (nx,a,b,val1,val2,flag1,flag2,flagf) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FD1D 
%% fd1d(nx,a,b,val1,val2,flag1,flag2,flagf) 
%% 1d FD cell-centered solution to Poisson equation on (a,b)
%% user must code permfun, rhsfun, exfun (at bottom of code)
%% val1, val2 are values of bdary data on left and right endpoints
%% flag1/2 == 0: val1/2 is a Dirichlet value
%% flag1/2 == -1: val1/2 is a Neumann value
%% flagf == 0: Set rhs = 0 (optional argument)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs:
%% u is the numerical solution at xx, the cell centers.
%% h is the maximum step size between nodes
%% dirval1/2 are the dirichlet data
%% qflux1/2 are the boundary fluxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% examples
%% [u,h,x,d1,d2,f1,f2]=fd1d(1000,0,1,0,0,0,0);
%% [u,h,x,d1,d2,f1,f2]=fd1d(30,0,1,0,-pi,0,-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin <= 7 
    flagf = 1;
end;
dx = (b-a)/nx*ones(nx,1); x0 = a; 

%%%%% create data structures
x = dx; x(1)=x0; for i=2:nx x(i)=x(i-1)+dx(i);end;

%% permeability coefficient
perm = ones(nx,1);  %%or use hard-coded permfun below

%%%%% compute transmissibilities
tx = zeros(nx+1,1); 
for i=2:nx 
       tx(i)=2/(dx(i-1)/perm(i-1)+dx(i)/perm(i));
end;
 
if flag1 == 0     %% Dirichlet contributions to the transmissibilities
  i = 1;         tx(i)=2/(dx(i)/perm(i));
end;
if flag2 == 0
  i = nx + 1;    tx(i)=2/(dx(i-1)/perm(i-1));
end

stiff = sparse (nx,nx);
for i=2:nx                                       
  gl = i-1;                                             
  gr = i;
  stiff(gl,gl) = stiff(gl,gl) + tx(i);                                    
  stiff(gl,gr) = stiff(gl,gr) - tx(i);                                                              
  stiff(gr,gl) = stiff(gr,gl) - tx(i);                                    
  stiff(gr,gr) = stiff(gr,gr) + tx(i);     
end;

%% compute rhs

if flagf == 0
    q = zeros(nx,1);
else
    q = dx.* rhsfun (x+dx/2);
end;

%% contributions of bdary conditions to the matrix and rhs 
if flag1 == 0
    i = 1; gr = 1;       
    stiff(gr,gr) = stiff(gr,gr) + tx(i);   
    q(1) = q(1) + tx(1) * val1;
else
    q(1) = q(1) + val1;
end;

if flag2 == 0
    i = nx + 1; gl = nx; 
    stiff(gl,gl) = stiff(gl,gl) + tx(i);  
    q(nx) = q(nx) + tx(nx+1) * val2;
else    
    q(nx) = q(nx) + val2;
end

nsol = stiff \ q;

xplot = x + dx/2;

plot(xplot, exfun(xplot),xplot,nsol,'r*');

%%%%%% post-processing: extract boundary fluxes and Dirichlet data 

if flag1 == 0
    dirval1 = val1;
    qflux1 = - tx(1)*(nsol(1)-val1);
else 
    i = 1;  tx(i)=2/(dx(i)/perm(i));
    dirval1 = nsol(1) - tx(2)/tx(1)*(nsol(2)-nsol(1));
    qflux1 = val1;
end;

if flag2 == 0
    dirval2 = val2;
    qflux2 = - tx(nx+1)*(nsol(nx)-val2);
else
    i = nx + 1;    tx(i)=2/(dx(i-1)/perm(i-1));
    dirval2 = nsol(nx) - tx(nx)/tx(nx+1)*(nsol(nx-1)-nsol(nx));
    qflux2 = val2;
end;

xx = xplot;
h = max(dx);
u = nsol;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comment/Uncomment related examples below, or code/create your own. 
%% Ensure input parameters are consistent with choices below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = permfun(x)
v = 1;

function v = rhsfun(x)
v = pi * pi * sin(pi*x);
%v = pi * pi * cos(pi*x);
%v = 1;

function v = exfun(x)
v = sin(pi*x);     %% with dirichlet 0,0 first order, with flux condition pi,-pi quadratic
%v = cos(pi*x);    %% with Neumann conditions
%v = 1/2*x.*(1-x); %% with dirichlet: 0,0 first order, with flux superconvergent
