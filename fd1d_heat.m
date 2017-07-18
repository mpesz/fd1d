function [xplot] = fd1d_heat (nx,dt,bdaryflag,outflag) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fd1d_heat (nx,dt,bdaryflag,outflag)
%% 1D FD cell-centered solution to heat equation on (xbeg,xend)
%% user must code the functions permfun, porfun, rhsfunt, 
%%    initfun, exfun, dexfun INSIDE the code (below)
%% nx is the number of sub-intervals of (0,1) for uniform grid
%% dt is the time step
%% dt == 0: steady state solution only
%% bdaryflag == 0: Dirichlet (nonhomogeneous) BC on both ends
%% bdaryflag != 0: Neumann (nonhomogeneous) BC on both ends
%% note: user must code exfun/dexfun to deliver the BC values 
%% outflag == 0: only plot solution in time 
%% outflag != 0: consider exact solution, compute error etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs:
%% xplot : the cell-centers
%% nsols : the numerical solution
%%   the code will also plot the numerical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xbeg = 0; xend = 1;               %% or change to another interval
dx = (xend-xbeg)/nx * ones(nx,1); %% or code the nonuniform grid 


t1 = 0; t2 = 1;     %% or use a different time interval
if dt == 0 
    nt = 1; transient = false;  dt = 1;
else 
     transient = true; nt = (t2-t1) / dt; 
end;

%%%%% create data structures
x = linspace(xbeg,xend-dx(nx),nx)'; %% x is the left end of each subinterval
xplot = x + dx/2;  %% xplot is the center of each subinterval

%% permeability coefficient
perm = 0*xplot; perm = permfun(xplot); %#ok<*NASGU>

%% porosity coefficient
if transient, por = 0*x; por= porfun(xplot);end

%%%%% compute transmissibilities
tx = zeros(nx+1,1); 
for j=2:nx 
       tx(j)=2/(dx(j-1)/perm(j-1)+dx(j)/perm(j));
end
 
if bdaryflag == 0     %% Dirichlet contributions to the transmissibilities
  j = 1;         tx(j)=1/(dx(j)/perm(j));
  j = nx + 1;    tx(j)=1/(dx(j-1)/perm(j-1));
end
tx = tx*dt;
    
stiff = sparse (nx,nx);
for j=2:nx                                       
  gl = j-1;                                             
  gr = j;
  stiff(gl,gl) = stiff(gl,gl) + tx(j);                                    
  stiff(gl,gr) = stiff(gl,gr) - tx(j);                                                              
  stiff(gr,gl) = stiff(gr,gl) - tx(j);                                    
  stiff(gr,gr) = stiff(gr,gr) + tx(j);     
end;
    
%% Dirichlet contributions to matrix
if bdaryflag == 0
    j = 1; gr = 1;       stiff(gr,gr) = stiff(gr,gr) + tx(j);     
    j = nx + 1; gl = nx; stiff(gl,gl) = stiff(gl,gl) + tx(j);  
end

%
errvec = zeros(nt,1);

if transient
    for j=1:nx
        stiff(j,j) = stiff(j,j) + por(j)*dx(j); %#ok<SPRIX>
    end
    nsol = initfun ( xplot);    
end

%nsols = [];
for n = 1 : nt %  time step loop
    t = t1 + n*dt;
    
    if ~transient 
        q = dx.* rhsfunt (x+dx/2,0);
    else
        q = dx.* (dt*rhsfunt ( x+dx/2, t) + por.*nsol);
    end;
    
    if bdaryflag == 0
        val1 = exfun(xbeg,t-dt/2);
        val2 = exfun(xend,t-dt/2);
    else
        val1 = - dexfun(xbeg,t-dt/2);
        val2 =  dexfun(xend,t-dt/2);
    end;
    
    %% contributions of bdary conditions to the rhs 
    if bdaryflag == 0
        q(1) = q(1) + tx(1) * val1;
        q(nx) = q(nx) + tx(nx+1) * val2;
    else
        q(1) = q(1) + dt*val1;
        q(nx) = q(nx) + dt*val2;
    end

%%%%% find numerical solution
    nsol = stiff \ q;
 
    %% nsols = [nsols, nsol];
    if outflag == 0
        plot(xplot,nsol,'r*-'),
        axis([0 1 -1 1]);
        pause(.2);
    else
        if transient, sol = exfun(xplot,t); else sol = exfun(xplot,0); end;
        
        plot(xplot,sol,'b',xplot,nsol,'r*-'),
        axis ([0 1 0 0.5]);
        pause(.1);
        errvec (n)= norm(sol-nsol,inf);
    end
end
if outflag ~= 0 
    errornorm = max(errvec);
    fprintf('dx = %g dt = %g error = %g \n',max(dx),dt,errornorm);
end

end  %% function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comment/Uncomment related examples below, or code/create your own. 
%% Ensure input parameters are consistent with choices below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = permfun(x)
%%
v = zeros(size(x,1),1) + 1;
%% uncomment to have a low permeability zone 
% v ( find (abs(x-0.3)<0.2)) = 1e-1;
end

function v = porfun(x)
%% porosity 
v = 1 + 0*x;
%% realistic porosities are < 1.0
end

function v = rhsfunt(x,t)
v = 2*(t) + x.*(1-x);
%v = zeros(size(x));
end

function v = initfun(x)
v = exfun(x,0);
%v = zeros(size(x));
end

function v = exfun(x,t)
v = (t)*x.*(1-x); %% with dirichlet: 0,0 first order, with flux -1,-1 quadratic (superconvergent)
%v = zeros(size(x));
%v = x;
end

function v = dexfun (x,t)
v = (t)*(1-2*x);
%v = zeros(size(x));
%v = 1;
end
