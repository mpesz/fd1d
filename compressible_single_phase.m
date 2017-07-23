function compressible_single_phase (nx,dfac,dt0,t1,t2,...
    bdary1,bdary2,val1,val2,...
    implicit_explicit,...
    pmin,pmax, ...
    ifpause)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1D FD cell-centered solution to compressible flow equation 
%% PRESSURE UNKNOWN
%% <nx, dfac, dt0, t1, t2>: physical domain
%%    nx: Number of sub-intervals of (a,b) (hard-coded(0,0.6))
%%    dfac: sin(angle), angle=0 horizontal, angle= pi/2 vertical
%%    dt0: time step over the time interval (t1,t2).
%% dfac == 0: depth = 0,    dfac == 1: depth = x;
%% bdary1/2 == 0: Dirichlet (nonhomogeneous) values val1/2
%% bdary1/2 == 1: Neumann flux values val1/2
%% implicit_explicit == 'explicit': Solve sequentially
%% <pmin, pmax, ifpause>: parameters for plotting
%% Note: This is hard-coded with SI units, K with [m^2], 
%%		pressure with [Pa].  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs:
%%   the code will plot the numerical solution over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Examples:
%%  No-flow at the bottom, vertical
%%   compressible_single_phase(10,1,0.1,0,1,0,1,2000,0,0,0,5e4,1)
%%  Horizontal case, Dirichlet conditions
%%   compressible_single_phase(10,0,0.1,0,1,0,0,0,1000,0,0,1e3,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo_input = 1;
if echo_input
    fprintf('***********************\n');
    fprintf ('nx=%d\n',nx);
    fprintf('dfacl=%g\n',dfac);
    fprintf('dt0=%g T in [%g,%g]\n',dt0,t1,t2);
    fprintf('initialization option: impl_expl=%\n',implicit_explicit);
    fprintf('bdary conditions: type %d %d values %g %g\n',bdary1,bdary2,val1,val2);
    fprintf('display options: plot p in  [%g,%g] \n', pmin,pmax);
    fprintf('pause parameter=%d\n',  ifpause) ;
    fprintf('***********************\n');
end
%%%%%%%%%%%%%%%%%%%%%% rock/sand data

%%% constants to be used for permeability and porosity
permeability = 1e-10;  %%[m^2]
porosity     = 0.4;

%%%%%%%%%%%%%%%%%%%%%% tolerance of Newton solver 
tol = 1e-4; atol = 1e-8; maxiter = 20;
if strcmp(implicit_explicit,'explicit'), maxiter = 1;end %% only one Newton iteration is taken: sequential algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% discretization parameters
%% depth of the reservoir
dx = (0.6-0)/nx * ones(nx,1); %% uniform grid (code does not depend on it)
x0 = 0; x = zeros(size(dx)); 
x(1) = x0; for j = 2:nx, x(j)=x(j-1)+dx(j);end
xplot = x + dx/2;

%%%%% dfac is the sin(angle), where angle is=0 for horizontal, angle =pi/2
%%%%% for vertical
depth = dfac*xplot;

%% constant porosity coefficient
por = ones(size(x))*porosity; 

%%%%%% hydraulic conductivity coefficient: use Darcy permeability divided by visc
visc_w    = 1e-3;
Ks        = permeability/visc_w;
perm      = ones(size(x))*Ks;

%%%%%%%%%%%%%%%%% gravity premultiplied by density
pref = 0; rhowref = 1e3; compr_w = 1e-5;
density_w = @(p)(rhowref*exp(compr_w*(p-pref)));
ddensity_w = @(p)(compr_w*rhowref*exp(compr_w*(p-pref)));

density_r = @(p)(rhorref*exp(1e-6*(p-pref)));
grav = 9.8066; 

%%%%% compute transmissibilities: just geometry
txx = zeros(nx+1,1); 
for j=2:nx 
       txx(j)=2.0/(dx(j-1)/perm(j-1)+dx(j)/perm(j));
end
 
if bdary1 == 0     %% Dirichlet contributions to the transmissibilities
        j = 1;         txx(j) = 2/(dx(j)/perm(j));
        depthl = 3/2*depth(j)-1/2*depth(2);
end
if bdary2 == 0     %% Dirichlet contributions to the transmissibilities
        j = nx + 1;    txx(j) = 2/(dx(j-1)/perm(j-1));
        depthr = 3/2*depth(nx)-1/2*depth(nx-1);
end

%%%%%%%%%%%%% initialize pressures and saturations
p = zeros(size(x));

%% initialize the pressure to be approximately hydrostatic
p = pref + rhowref*grav*depth;

if ifpause
    plot(xplot,p); axis([0 0.6 pmin pmax]);
    title('Initial condition: pressure');
    pause;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
fprintf('Boundary conditions: pressure: %g %g\n',val1,val2);
%%
porv = por.*dx;

%%%%%%%%%%%%%%%%%%% 
dt = dt0;
t = t1 ;
up = 0;
down = 0;
n = 0;
control = 0;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TIME STEP LOOP
while  t< t2 %% && err_level == 0 %  time step loop
    n = n + 1;
    olddeltat = dt;
    
    %%% primitive time step control: select dt based on number of iters in previous time step;   
    if control == 1
        if n>1 && iter >= maxiter/2+1, dt = dt/2;
            fprintf('n=%d iter =%d dt=%g <- %g\n',n,iter, dt*2,dt);down = down +1;end;   
        if n>1 && iter <= 3, dt = dt*2;
            fprintf('n=%d iter= %d dt=%g -> %g\n',n,iter,dt/2,dt);up = up+1;end;   
    end
    %%
    if t + dt > t2, dt = t2 - t; end;
    %%
    t = t + dt;
    tx = txx*dt;    
    
    %%%%%% old time step values: they can be used for extrapolation
    %%%%%% example how boundary conditions can change in time 
    if t> 10000
        pval1 = 0.05;
    end
    
    %%%%%% p 
    pold = p; rhoold = density_w(pold);
    
    %%%%%% solve for new time step value of p using Newton iteration
    pguess = p; 
    %%% start Newton loop
    iter = 0; 
    while iter < maxiter 
        iter = iter + 1;
        %%%%%% given pguess compute properties (p is primary unknown)  
        rhoguess = density_w(pguess);
        drhodp = ddensity_w(pguess);
            
        jac = sparse (nx,nx);
        r   = zeros(size(pguess));

        %%% symmetric part of jacobian and residual
        for j = 1:nx
            r(j) = porv(j)*(rhoguess(j)-rhoold(j));
            jac (j,j) = porv(j)*drhodp(j); 
        end
               
        %%%% compute non-symmetric part of residual and jacobian: 
        %%%% loop over the cell edges
        
        for j=2:nx 
            delpp = pguess(j)-pguess(j-1);
            rhoedge  = (rhoguess(j-1)+rhoguess(j))/2;      
            deld = - grav *(depth(j)-depth(j-1));
            delp =  delpp + rhoedge*deld;
            
            
            r(j-1) = r(j-1) - tx(j)*rhoedge*delp;
            r(j)   = r(j)   + tx(j)*rhoedge*delp;
            %%             
            jac(j,j) = jac(j,j) + tx(j) * (drhodp(j)/2*delp + ...
                rhoedge * (1 + drhodp(j)/2*deld) );
            jac(j,j-1) = jac(j,j-1) + tx(j) * (drhodp(j-1)/2*delp + ...
                rhoedge * (-1 + drhodp(j)/2*deld) );
            %%
            jac(j-1,j) = jac(j-1,j) - tx(j) * (drhodp(j)/2*delp + ...
                rhoedge * (1 + drhodp(j)/2*deld) );
            jac(j-1,j-1) = jac(j-1,j-1) - tx(j) * (drhodp(j-1)/2*delp + ...
                rhoedge * (-1 + drhodp(j)/2*deld) );
                    
        end
        
         if bdary1 == 0 %% Dirichlet condition for pressure
            pressl = val1; 
            rhoedge = density_w(pressl); %% no averaging for simplicity
            
            delp_l = pguess(1)-pressl-grav*rhoedge*(depth(1)-depthl);
            
            r(1) = r(1) + tx(1)*rhoedge*delp_l;     
            jac(1,1) = jac(1,1) +  tx(1) * rhoedge;
            
            flux_l = tx(1)*rhoedge*delp_l/dt;
         end
        
        if bdary2 == 0 
            pressr = val2; 
            rhoedge = density_w(pressr);
            delp_r = pressr - pguess(nx)-grav*rhoedge*(depthr-depth(nx));
            
            r(nx) = r(nx) - tx(nx+1) * rhoedge * delp_r;
            %%
            jac(nx,nx) = jac(nx,nx) -  tx(nx+1) * rhoedge * (-1);
             
            flux_r = - tx(nx+1)*rhoedge*delp_r / dt;
        end
        
        if bdary1 == 1 %% Neumann flux 
            r(1) = r(1) + val1*dt;
            flux_l = val1;
        end
        if bdary2 == 1 %% Neumann flux 
            r(nx) = r(nx) + val2*dt;
            flux_r = val2;
        end
  
        
        %%% decide if we need to solve or not 
        rn = norm(r,inf);
        if iter == 1, rn0 = rn;end;    
        if rn < tol*rn0 || rn < atol, break; end;

        %%%%% solve for pressure and update
        dvar =  jac \ r;
    
        pguess = pguess - dvar;
    
    end %% Newtonian loop
    %% report on Newtonian iteration
    if ifpause
        fprintf('Step=%d time=%g nits=%d res = %g mass=%g maxp=%g minp=%g flux=%g %g\n',...
            n,t,iter,rn,sum(rhoguess.*porv),max(pguess),min(pguess),flux_l,flux_r);
    end
    %% save new time step value
    p = pguess;
    
    if ifpause == 1 
        titp = sprintf('Pressure   at t=%g step = %d',t,n);
      
        plot(xplot,p);axis([0 0.6 pmin pmax]);
        title(titp);
        pause(0.05),
    end    
    
end   %% end time loop
fprintf('Number of up=%d down=%d original dt=%g final=%g\n',up,down,dt0,dt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comment/Uncomment examples below, or code/create your own, if 
%% non-constant permeability and/or porosity is desired. 
%% Ensure input parameters are consistent with choices below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%function v = permfun(x) %Permeability
%v = ones(size(x));
%end

%function v = porfun(x) %Porosity
%v = ones(size(x));
%end



