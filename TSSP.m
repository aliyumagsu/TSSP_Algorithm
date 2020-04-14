function [iter,fval,norm_Fx0] = TSSP(f,x0)
% This MATLAB code implement the TSSP algorithm presented in: 
% A.M. Awwal, L. Wang, P. Kumam, and H. Mohammad
% A two-step spectral gradient projection method for system of nonlinear 
%monotone equations and image deblurring problems
% Symmetry MDPI 2020 (to appear)
% Written by: A.M. Awwal and H. Mohammad (March, 2020)
% Input:  x0 is the initial point and f is the function handle
% Output: iter = number of iterations, fval = number of function evaluations
% norm_Fx0 = norm of the residual function at the solution
% Note: This code does not come with guarantee of any sort.
% Version 1.0 : Last modified by A.M. Awwal (April, 2020)
tic;
%Step 0 Initialization
iter=0; fval=0; bck=0;%inalpha=1; 
Fx0=feval(f,x0); % evaluating F(x0);
fval=fval+1;
%Compute ||F(x0)||
norm_Fx0= norm(Fx0);
%Set d(x0)=-F(x0) for k=0
BB1=1;
dx0=-BB1*Fx0;
while(iter<=maxiter && norm_Fx0>tol)
    alpha=1/(iter+1)^2;
    %alpha=1/(iter+1);
     r=0.01; 
 % Update the iterate using the first step
    w0=x0+alpha*dx0;
% Compute F(w0)
    Fw0=feval(f,w0);
    %Compute sk and yk w.r.t. w0 and x0
    sw0=w0-x0;
    yw0=Fw0-Fx0+r*sw0;
    %Compute BB1 and search direction w.r.t. w0 and x0
    BB2 = (sw0'*yw0)/(yw0'*yw0);
    dw0=-BB2*Fx0;
    % Prepare the line search
    sig=1e-2;  beeta=0.5; t=1;
    while (fval<maxfval)
        if -(feval(f,x0+t*dw0))'*dw0 < sig*t*(dw0'*dw0)*(norm(feval(f,x0+t*dw0)))^(0.5)
            t=beeta*t;
        else
            break
        end
        bck=bck+1;        
    end
    beeta=t;
    % Update the iterate using the second step
    zk=x0+beeta*dw0;
    %Compute ||F(zk)||
    Fzk=feval(f,zk);
    norm_Fzk= norm(Fzk);
    fval=fval+1;
    if (feval(proj,(zk),l,u)==zk & norm_Fzk<tol)
        x0=zk;
        Fx0=Fzk;
        norm_Fx0= norm(Fx0);
        disp('zk is in the convex set and its the solution at iteration number')
        disp([num2str(iter)])
        break
    else
        zetak=Fzk'*(x0-zk)/(Fzk'*Fzk); % computing zetak
        P=feval(proj,(x0-zetak*Fzk),l,u); % projection on convex set
        xk=P;
        Fxk=feval(f,xk);
        %Compute sk and yk w.r.t. xk and x0
        sxk=xk-x0;
        yxk=Fxk-Fx0+r*sxk;
        %Compute BB1 and search direction w.r.t. xk and w0
        BB1 = (sxk'*sxk)/(sxk'*yxk);
        dxk = -BB1*Fxk;
    end
    dx0=dxk;
    x0=xk;
    Fx0=Fxk;
    norm_Fx0= norm(Fx0);
    iter=iter+1; 
end
toc