%--- University of Washington, Department of Aeronautics & Astronautics ---
%---------- Advanced Dynamics, Validation & Control Research Lab ----------
%
% Runs COSY continuation on equation with varying parameter designations to
% produce multiple data sets for multi-parameter hybrid stability automata
%
% Currently works for 1D and 2D systems
%
% Author: Peter Uth
% Created: March 2017
%
% Inputs:
% f     - original function, e.g. f = @(x,u,p) p(1)+p(2)*x+p(3)*x^2
% PD    - parameter designation matrix
% x0,p0 - initial solution guesses for states and parameters
% plim  - active parameter limits
%
% Outputs:
% F     - system equations, column vector with first entry as original eq
%         e.g. F = {@(x,u,p)p(1)+p(2)*x^2;@(x,u,p1)p1+x^2;@(x,u,p2)p2*x^2}
% B, BP - solution branch and branch point data sets from COSY
%--------------------------------------------------------------------------
function [F,B,BP] = dataMP(f,PD,x0,p0,plim)

F(1,:) = {f};
B = cell(size(PD,1),1);
BP = cell(size(PD,1),1);
for i = 1:size(PD,1)
    ftemp = func2str(f);
    
    % Determine sub-equation
    for j = 1:size(PD,2)
        if isnan(PD(i,j))
            ftemp = strrep(ftemp,sprintf('p(%.0f)',j),sprintf('p%0.f',j));
            ftemp = strrep(ftemp,'@(x,u,p)',sprintf('@(x,u,p%.0f)',j));
            pvar = sprintf('p%0.f',j);
        else
            ftemp = strrep(ftemp,sprintf('p(%.0f)',j),sprintf('%.0f',PD(i,j)));
        end
    end
    ftemp = str2func(ftemp);
     
    % Run COSY on sub-equation
    c = cosy(ftemp,size(x0,2),0,1);
    c.cm.max_steps = 15000;
    c.p_min = plim(i,1);
    c.p_max = plim(i,2);
    [b,bp] = c.trace_solution_branches(x0(i,:)',[],p0(i,1));   
        
    % Simplify function for display if symbolic toolbox is available
    try
        if size(x0,2) == 1
            syms x u 
            ftemp = sprintf('@(x,u,%s)%s',pvar,char(vpa(ftemp(x,u,sym(pvar)),4)));
            ftemp = str2func(ftemp);
        elseif size(x0,2) == 2
            syms x1 x2 u       
            ftemp = strrep(strrep(char(ftemp),'x(1)','x1'),'x(2)','x2');
            ftemp = strrep(ftemp,'@(x,u,','@(x1,x2,u,');
            ftemp = str2func(ftemp);
            ft = ftemp(x1,x2,u,sym(pvar));
            ftemp = sprintf('@(x,u,%s) [%s ; %s]',pvar,char(ft(1)),char(ft(2)));
            ftemp = strrep(strrep(char(ftemp),'x1','x(1)'),'x2','x(2)');
            ftemp = str2func(ftemp);
        end
    catch 
        disp('Symbolic Toolbox not available - function not simplified');
    end
    
    F(i+1,:) = {ftemp};
    B(i,:) = {b};
    BP(i,:) = {bp};
end
