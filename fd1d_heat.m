function [xplot,nsols] = fd1d_heat (nx,dt,bdaryflag,outflag) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fd1d (nx,dt,bdaryflag,outflag)
%% 1D FD cell-centered solution to heat equation on (0,1)
%% user must code the functions permfun, rhsfun, exfun,dexfun,initfun INSIDE the code
%% bdaryflag == 0: Dirichlet (nonhomogeneous) values on both ends
%% bdaryflag !=0: Neumann (nonhomogeneous) values on both ends
%% note: user must code exfun. dexfun to deliver these values
%% dt ==0 : steady state solution only 
%% outflag =0: only plot solution in time 
%% outflag > 0: consider exact solution, compute error etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = 1/nx * ones(nx,1); x0 = 0; 

t1 = 0; t2 = 1;
if dt == 0 
    nt = 1; transient = false;  
else 
     transient = true; nt = (t2-t1) / dt; 
end;

%%%%% create data structures
x = dx; x(1)=x0; for i=2:nx x(i)=x(i-1)+dx(i);end;
xplot = x + dx/2;

%% permeability coefficient
perm = x; perm(:,1) = permfun(x(:,1)+ dx(:,1)/2);
%% porosity coefficient
if transient, por = x; por(:,1) = porfun(x(:,1)+ dx(:,1)/2);end;

%%%%% compute transmissibilities
tx = zeros(nx+1,1); 
for i=2:nx 
       tx(i)=2/(dx(i-1)/perm(i-1)+dx(i)/perm(i));
end;
 
if bdaryflag == 0     %% Dirichlet contributions to the transmissibilities
  i = 1;         tx(i)=1/(dx(i)/perm(i));
  i = nx + 1;    tx(i)=1/(dx(i-1)/perm(i-1));
end

tx = tx*dt;
    
stiff = sparse (nx,nx);

for i=2:nx                                       
  gl = i-1;                                             
  gr = i;
  stiff(gl,gl) = stiff(gl,gl) + tx(i);                                    
  stiff(gl,gr) = stiff(gl,gr) - tx(i);                                                              
  stiff(gr,gl) = stiff(gr,gl) - tx(i);                                    
  stiff(gr,gr) = stiff(gr,gr) + tx(i);     
end;
    
%% Dirichlet contributions to matrix
if bdaryflag == 0
    i = 1; gr = 1;       stiff(gr,gr) = stiff(gr,gr) + tx(i);     
    i = nx + 1; gl = nx; stiff(gl,gl) = stiff(gl,gl) + tx(i);  
end
errvec = [];

if transient
    for i=1:nx
        stiff(i,i) = stiff(i,i) + por(i)*dx(i);
    end;
    nsol = initfun ( x + dx/2 );    
end;

%nsols = zeros(nt,nx);
nsols = [];
for n = 1 : nt %  time step loop
    t = t1 + n*dt;
    
    if ~transient 
        q = dx.* rhsfunt (x+dx/2,0);
    else
        q = dx.* (dt*rhsfunt ( x+dx/2, t) + por.*nsol);
    end;
    
    if bdaryflag == 0
        val1 = exfun(0,t-dt/2);
        val2 = exfun(1,t-dt/2);
    else
        val1 = - dexfun(0,t-dt/2);
        val2 =  dexfun(1,t-dt/2);
    end;
    
    %% contributions of bdary conditions to the rhs 
    if bdaryflag == 0
        q(1) = q(1) + tx(1) * val1;
        q(nx) = q(nx) + tx(nx+1) * val2;
    else
        q(1) = q(1) + dt*val1;
        q(nx) = q(nx) + dt*val2;
    end

% condest(stiff)
%%%%% find numerical solution
    nsol = stiff \ q;
 
    nsols = [nsols, nsol];
    if outflag == 0
        fprintf('\ntime step %d at t = %g',n,t1+n*dt);
        plot(xplot,nsol,'r*-'),
        hold on;
        %plot(xplot,xplot-nsol,'b-'),
        hold off;
        axis([0 1 -1 1]);
        pause(.2);
    else
        if transient, sol = exfun(xplot,t); else sol = exfun(xplot,0); end;
        
        plot(xplot,sol,'b',xplot,nsol,'r*-'),
        axis ([0 1 0 0.5]);
        pause(.1);
        %plot(xplot,sol,'b*',xplot,xplot,'r');axis([0,1,-1,1]);
        errvec = [errvec; norm(sol-nsol,inf)];
    end
% return;
% pause
end
errornorm = max(errvec);
fprintf('dx = %g dt = %g error = %g ',max(dx),dt,errornorm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function v = permfun(x)
v = zeros(size(x,1),1);
for j=1:size(x,1)
    if abs(x(j)-.3) < 0.2 v(j) = 1e-1;
    else v(j) = 1e0;
    end;
end;

function v = porfun(x)
v = 1;

function v = rhsfunt(x,t)
%v = 2*(t) + x.*(1-x);
v = zeros(size(x));

function v = initfun(x)
%v = exfun(x,0);
%v = ones(size(x));
v = zeros(size(x));

function v = exfun(x,t)
%v = (t)*x.*(1-x); %% with dirichlet: 0,0 first order, with flux -1,-1 quadratic (superconvergent)
%v = zeros(size(x));
v = x;

function v = dexfun (x,t)
v = (t)*(1-2*x);
