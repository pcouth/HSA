%--- University of Washington, Department of Aeronautics & Astronautics ---
%---------- Advanced Dynamics, Validation & Control Research Lab ----------
%
% Hybrid Stability Automaton (HSA) generator. Extracts stability mode data
% from COSY numberical continuation tool output, then creates and executes
% a GraphViz file to produce a PNG of the HSA. Also creates mode data and 
% an executable function m-file that provides stability information.
%
% Multi-dimensional version (currently only handles 2D)
%
% Author: Peter Uth
% Created: April 2017
%
% Inputs:
% label - Name as a string, used for labeling output files
% F     - System equation(s) in form @(x,u,p)p+x^2
%         If multiple, column vector with first entry as original equation
%         e.g. F = {@(x,u,p)p(1)+p(2)*x^2;@(x,u,p1)p1+x^2;@(x,u,p2)p2*x^2}
% B, BP - Solution branch and branch point data from COSY
% P     - User defined p value, leave blank [] for general case
% opt   - Optional inputs as struct with the following fields:
%       PD      - parameter designation matrix
%       dot2tex - 'dot2tex' as string to enable LaTeX-style output format
%       texlab  - replace x and p with others symbols while using dot2tex
%       ztol    - zero tolerance, default 0.00001
%
% Output:
% node  - Resulting mode data as struct
%
% Requires folder containing GraphViz "dot.exe" in environmental variables
%--------------------------------------------------------------------------
function [node] = hsaMD(label,F,B,BP,P,opt)

% Set optional inputs
if exist('opt','var')

    % Set parameter designation matrix if PD is opt field or input is a double
    if isa(opt,'struct') && isfield(opt,'PD')
        PD = opt.PD;
    elseif isa(opt,'double')
        PD = opt;
    end

    % Set latex replace if texlab is opt field or input is a cell
    if (isa(opt,'struct') && isfield(opt,'texlab')) || isa(opt,'cell')
        if isa(opt,'struct') && isfield(opt,'texlab')
            texlab = opt.texlab;
        elseif isa(opt,'cell')
            texlab = opt;
        end
    end
end

% Set zero tolerance if ztol is opt field or use default
if exist('opt','var') && isa(opt,'struct') && isfield(opt,'ztol')
    ztol = opt.ztol;
else
    ztol = 0.00001;
end

% Identify function(s), create multiple equation formats for later
if size(F,1) > 1
    fun.orig = func2str(F{1});
    F = F(2:end);
else
    fun.orig = func2str(F);
    F = {F};
    B = {B};
    BP = {BP};
end
fun.origgv = strrep(fun.orig,'@(x,u,p)','x<SUP>''</SUP> = ');
fun.origgvtex = strrep(fun.orig,'@(x,u,p)','\begin{bmatrix}\dot{x}_1\\\\\dot{x}_2\end{bmatrix}$=$');
fun.origgvtex = strrep(fun.origgvtex,'*','');
if size(F,1) > 1
    for i = 1:size(PD,2)
        fun.origgvtex = strrep(fun.origgvtex,sprintf('p(%.0f)',i),sprintf('p_%.0f',i));
    end
end
fun.origgvtex = strrep(fun.origgvtex,'x(1)','x_1');
fun.origgvtex = strrep(fun.origgvtex,'x(2)','x_2');
fun.origgvtex = strrep(fun.origgvtex,'[','\begin{bmatrix}');
fun.origgvtex = strrep(fun.origgvtex,']','\end{bmatrix}');
fun.origgvtex = strrep(fun.origgvtex,';','\\\');

% Set standard widths and heights for GraphViz nodes
node_width = 2.5;
node_height = 1.5;
node_width2 = 2;
node_height2 = 1.25;

% Generate DOT file with header, fid: GraphViz, fid2: dot2tex, fid3: texlab
dotfile = [label '.dot'];
dotfile_tex = [label 'Tex' '.dot'];
dotfile_texR = [label 'TexR' '.dot'];
fid = fopen(dotfile,'w');
fid2 = fopen(dotfile_tex,'w');
if exist('texlab','var')
	fid3 = fopen(dotfile_texR,'w');
    fclose(fid3);
end
fprintf(fid,'digraph {\n');
fprintf(fid2,'digraph {\n');
fprintf(fid,'compound=true;\n');
fprintf(fid2,'compound=true;\n');
fprintf(fid,'graph[style="rounded,filled" fillcolor=wheat]\n');
fprintf(fid2,'graph[style="rounded" fillcolor=wheat]\n');
fprintf(fid,'nodesep=1\n');
fprintf(fid2,'nodesep=1\n');
fprintf(fid,'label=<%s <br/>> \n',fun.origgv);
fprintf(fid2,'label="%s" \n',fun.origgvtex);
fclose(fid);
fclose(fid2);

% Preallocate output node size, run main for-loop
node(size(F,1)) = struct();
for I = 1:size(F,1)

    % Clear temporary variables, activate relevant function and data
    clearvars -except F B BP P label ztol I dotfile fun PD PD2 node dotfile_tex...
        node_width node_height node_width2 node_height2 opt dotfile_texR texlab
    node(I).func = F{I};
    node(I).name = label;
    b = B{I};
    bp = BP{I};

    % Identify active parameter and create label for inactive parameters
    if exist('PD','var')
        active_parameter = find(isnan(PD(I,:)));
        inactive_parameters = char();
        inactive_parameters_tex = char();
        for i = 1:size(PD,2)
            if i == active_parameter
                continue
            else
                if isempty(inactive_parameters)
                    inactive_parameters = strcat(inactive_parameters,sprintf('p%0.f=%.2g',i,PD(I,i)));
                    inactive_parameters_tex = strcat(inactive_parameters_tex,sprintf('p_%0.f=%.2g',i,PD(I,i)));
                else
                    inactive_parameters = strcat(inactive_parameters,{', '},sprintf('p%0.f=%.2g',i,PD(I,i)));
                    inactive_parameters_tex = strcat(inactive_parameters_tex,sprintf(',\\;p_%0.f=%.2g',i,PD(I,i)));
                end
            end
        end
        inactive_parameters = char(inactive_parameters);
        node(I).PD = inactive_parameters;
    else
        active_parameter = [];
    end

    % Create top-modes with function labels
    if size(F,1) > 1
        fid = fopen(dotfile,'a');
        fid2 = fopen(dotfile_tex,'a');
        fprintf(fid,'subgraph cluster%.0f {\n',I);
        fprintf(fid2,'subgraph cluster%.0f {\n',I);
        fun.gv = strrep(func2str(node(I).func),sprintf('@(x,u,p%.0f)',active_parameter),'x<SUP>''</SUP> = ');
        fprintf(fid,'label=<%s<br/>%s<br/>> \n',fun.gv,inactive_parameters);
        fun.gvtex = strrep(func2str(node(I).func),sprintf('@(x,u,p%.0f)',active_parameter),'\begin{bmatrix}\dot{x}_1\\\dot{x}_2\end{bmatrix}=');
        fun.gvtex = strrep(fun.gvtex,'x(1)','x_1');
        fun.gvtex = strrep(fun.gvtex,'x(2)','x_2');
        fun.gvtex = strrep(fun.gvtex,'[','\begin{bmatrix}');
        fun.gvtex = strrep(fun.gvtex,']','\end{bmatrix}');
        fun.gvtex = strrep(fun.gvtex,';','\\');
        fun.gvtex = strrep(fun.gvtex,sprintf('p%.0f',active_parameter),sprintf('p_%.0f',active_parameter));
        fun.gvtex = strrep(fun.gvtex,'*','');
        fprintf(fid2,'texlbl="$\\begin{matrix}%s\\\\%s\\end{matrix}$" \n',fun.gvtex,inactive_parameters_tex);
        fclose(fid);
        fclose(fid2);
    else
        fun.gv = strrep(fun.orig,'@(x,u,p)','x<SUP>''</SUP> = ');
    end

    % Separate stable/unstable data within same solution branches
    event = [];
    split = 0;
    for i = 1:numel(b)
        b(1,i).x(abs(b(1,i).x) < ztol) = 0;
        b(1,i).p(abs(b(1,i).p) < ztol) = 0;
        b(1,i).eigval(abs(b(1,i).eigval) < ztol) = 0;
        branch(i+split).x = b(1,i).x'; %#ok<*AGROW>
        branch(i+split).p = b(1,i).p';
        branch(i+split).eigval = b(1,i).eigval';
        if max(sign(real(branch(i+split).eigval(2:end-1,1)))) * min(sign(real(branch(i+split).eigval(2:end-1,1)))) < 0 || ...
                max(sign(real(branch(i+split).eigval(2:end-1,2)))) * min(sign(real(branch(i+split).eigval(2:end-1,2)))) < 0
            LOC = [];
            for k = 2:size(branch(i+split).x,1) - 2
                if sign(real(branch(i+split).eigval(k,1))) * sign(real(branch(i+split).eigval(k+1,1))) < 0 || ...
                        sign(real(branch(i+split).eigval(k,2))) * sign(real(branch(i+split).eigval(k+1,2))) < 0
                    LOC = [LOC;k];
                    event = [event;branch(i+split).p(k,:) branch(i+split).x(k,:)];
                end
            end
            btemp = branch(i+split);
            for k = 1:size(event,1)
                if k == 1
                    branch(i+split).x = btemp.x(1:LOC(k,:),:);
                    branch(i+split).p = btemp.p(1:LOC(k,:),:);
                    branch(i+split).eigval = btemp.eigval(1:LOC(k,:),:);
                else
                    branch(i+split).x = btemp.x(LOC(k-1):LOC(k,:),:);
                    branch(i+split).p = btemp.p(LOC(k-1):LOC(k,:),:);
                    branch(i+split).eigval = btemp.eigval(LOC(k-1,:):LOC(k,:),:);
                end
                split = split + 1;
                if k == size(event,1)
                    branch(i+split).x = btemp.x(LOC(k,:):end,:);
                    branch(i+split).p = btemp.p(LOC(k,:):end,:);
                    branch(i+split).eigval = btemp.eigval(LOC(k,:):end,:);
                else
                    branch(i+split).x = btemp.x(LOC(k,:):LOC(k+1),:);
                    branch(i+split).p = btemp.p(LOC(k,:):LOC(k+1),:);
                    branch(i+split).eigval = btemp.eigval(LOC(k,:):LOC(k+1),:);
                end
            end
        end
    end
    num_branch = numel(branch);

    % Add COSY bp branch point output data to event data
    for i = 1:numel(bp)
        bp(i).x(abs(bp(i).x) < ztol) = 0;
        bp(i).p(abs(bp(i).p) < ztol) = 0;
        event(i+split,:) = [bp(i).p bp(i).x'];
    end

    % Total number of events (COSY branch points + number of branch splits)
    num_event = numel(bp) + split;

    % Create p-domain and p-label matrices
    if num_event == 0
        pdom = [-Inf Inf];
        labp = {sprintf('-Inf < p%.0f < Inf',active_parameter)};
        node(I).p(1).pdom = [-Inf Inf];
    else
        [~,order] = sort(event(:,1));
        event = event(order,:);
        pdom(1,:) = [-Inf event(1,1)];
        pdom(num_event+1,:) = [event(end,1) Inf];
        node(I).p(1).pdom = [-Inf event(1,1)];
        node(I).p(num_event+1).pdom = [event(end,1) Inf];
        labp{1} = sprintf('p%.0f < %.4g',active_parameter,event(1,1));
        labp{num_event+1} = sprintf('p%.0f > %.4g',active_parameter,event(end,1));
        if num_event > 1
            for i = 2:num_event
                pdom(i,:) = [event(i-1,1) event(i,1)];
                node(I).p(i).pdom = [event(i-1,1) event(i,1)];
                labp{i} = sprintf('%.4g < p%.0f < %.4g',event(i-1,1),active_parameter,event(i,1));
            end
        end
    end

    % Separate solution branch data into p-domains
    for i = 1:num_event + 1
        n = 0;
        for j = 1:num_branch
            xtemp = branch(j).x(pdom(i,1) < branch(j).p & branch(j).p < pdom(i,2),:);
            ptemp = branch(j).p(pdom(i,1) < branch(j).p & branch(j).p < pdom(i,2),:);
            eigtemp = branch(j).eigval(pdom(i,1) < branch(j).p & branch(j).p < pdom(i,2),:);
            if ~isempty(ptemp)
                n = n + 1;
                pdat(i,n).x = xtemp;
                pdat(i,n).p = ptemp;
                pdat(i,n).eigval = eigtemp;
                node(I).p(i).x(:,n) = {xtemp};
                node(I).p(i).p(:,n) = {ptemp};
                node(I).p(i).eigval(:,n) = {eigtemp};
            end
        end
        num_branch_p(i,1) = n;
    end

    % Determine branch x-ordering for each p-domain
    for i = 1:size(num_branch_p,1)
        stabp(i).p = [];
        for j = 1:num_branch_p(i,1)
            stabp(i).p(j,:) = [median(pdat(i,j).p) median(pdat(i,j).x) median(pdat(i,j).eigval)];
        end
    end

    % Determine branch stability
    ks = 1;
    ku = 1;
    for i = 1:num_event + 1
        if isempty(stabp(i).p)
            node(I).p(i).stab = [];
        else
            for j = 1:size(stabp(i).p,1)
                if stabp(i).p(j,4) < 0 && stabp(i).p(j,5) < 0
                    node(I).p(i).stab{1,j} = 'stable';
                    stabp(i).dom(j,:) = sprintf('Ds%.0f%.0f',I,ks);
                    ks = ks + 1;
                else
                    node(I).p(i).stab{1,j} = 'unstable';
                    stabp(i).dom(j,:) = sprintf('Du%.0f%.0f',I,ku);
                    ku = ku + 1;
                end
            end
        end
    end

    %% General Case -------------------------------------------------------
    if isempty(P)

        % Determine solution equations
        test = strrep(strrep(strrep(func2str(node(I).func),'x(1)','x1'),'x(2)','x2'),'@(x,u,','@(x1,x2,u,');
        soltemp = solve(eval(test));
        soltemp = [soltemp.x1 soltemp.x2];
        sol = [];
        for i = 1:size(soltemp,1)
        	sol = [sol;soltemp(i,:)];
        end

        % Use average branch values to match equations to nodes
        for i = 1:num_event + 1
            if size(stabp(i).p,1) == 0
                stabp(i).eq = [];
            end
            for j = 1:size(stabp(i).p,1)
                eval(sprintf('p%.0f=stabp(i).p(j,1);',active_parameter))
                soltemp = eval(sol);
                [~,loc] = min(abs(stabp(i).p(j,2)-soltemp(:,1)));
                stabp(i).eq(j,:) = sol(loc,:);
            end
            node(I).p(i).eq = stabp(i).eq;
        end

        % Add to DOT file
        fid = fopen(dotfile,'a');
        fid2 = fopen(dotfile_tex,'a');
        for i = 1:num_event + 1
            fprintf(fid,'subgraph cluster%.0f%.0f {\n fillcolor=white; \n',I,i);
            fprintf(fid2,'subgraph cluster%.0f%.0f {\n fillcolor=white; \n',I,i);
            fprintf(fid,'label=<%s > \n',strrep(strrep(labp{i},'<','&lt;'),'>','&gt;'));
            if size(F,1) > 1
                fprintf(fid2,'texlbl="$%s$"\n',strrep(strrep(labp{i},'p','p_'),'Inf','\infty'));
            else
                fprintf(fid2,'texlbl="$%s$"\n',strrep(labp{i},'Inf','\infty'));
            end
            fprintf(fid,'node [style=filled] \n');
            fprintf(fid2,'node [style=filled] \n');

            % Create DOT file nodes
            if isempty(stabp(i).p)
                fprintf(fid,'"eq%.0f%.0f%.0f" [fillcolor=grey label="no equilibrium" ]; \n',I,i,1);
                fprintf(fid2,'"eq%.0f%.0f%.0f" [style=dashed label="\\textbf{no equilibrium}" width=%.2f height=%.2f]; \n',I,i,1,node_width,node_height);
            else
                for j = 1:size(stabp(i).p,1)
                    if stabp(i).p(j,4) < 0 && stabp(i).p(j,5) < 0
                        fprintf(fid,'"eq%.0f%.0f%.0f" [fillcolor=green label="stable \n x1* = %s \n x2* = %s \n x(Ds)" ]; \n',I,i,j,char(stabp(i).eq(j,1)),char(stabp(i).eq(j,2)));
                        if size(F,1) > 1
                            fprintf(fid2,'"eq%.0f%.0f%.0f" [fillcolor=white color=green4 texlbl="$\\begin{matrix}\\textbf{stable}\\\\x_1^{*_%.0f}=%s\\\\x_2^{*_%.0f}=%s\\\\x\\in %s\\end{matrix}$" width=%.2f height=%.2f];\n',I,i,j,j,strrep(strrep(latex(stabp(i).eq(j,1)),'p','p_'),'\mathrm',''),j,strrep(strrep(latex(stabp(i).eq(j,2)),'p','p_'),'\mathrm',''),[strrep(stabp(i).dom(j,:),'D','D_{') '}'],node_width+round(0.02*size(char(stabp(i).eq(j,1)),2),1),node_height);
                        else
                            fprintf(fid2,'"eq%.0f%.0f%.0f" [fillcolor=white color=green4 texlbl="$\\begin{matrix}\\textbf{stable}\\\\x_1^{*_%.0f}=%s\\\\x_2^{*_%.0f}=%s\\\\x\\in %s\\end{matrix}$" width=%.2f height=%.2f];\n',I,i,j,j,strrep(latex(stabp(i).eq(j,1)),'\mathrm',''),j,strrep(latex(stabp(i).eq(j,2)),'\mathrm',''),[strrep(stabp(i).dom(j,:),'D','D_{') '}'],node_width+round(0.02*size(char(stabp(i).eq(j,1)),2),1),node_height);
                        end
                    else
                        fprintf(fid,'subgraph cluster%0.f%0.f%.0f {\n label="Unknown" \n graph[style="dashed,rounded"] \n',I,i,1);
                        fprintf(fid,'"eq%.0f%.0f%.0f" [fillcolor=red label="unstable \n x1* = %s \n x2* = %s" ]; \n',I,i,j,char(stabp(i).eq(j,1)),char(stabp(i).eq(j,2)));
                        fprintf(fid,'}\n');
                        if size(F,1) > 1
                            fprintf(fid2,'subgraph cluster%0.f%0.f%.0f {\n texlbl="$\\begin{matrix}\\textbf{Unknown}\\\\x\\in D_u \\end{matrix}$" \n graph[style="dashed,rounded"] \n',I,i,1);
                            fprintf(fid2,'"eq%.0f%.0f%.0f" [fillcolor=white color=red texlbl="$\\begin{matrix}\\textbf{unstable}\\\\x_1^{*_%.0f}=%s\\\\x_2^{*_%.0f}=%s\\end{matrix}$" width=%.2f height=%.2f];\n',I,i,j,j,strrep(strrep(latex(stabp(i).eq(j,1)),'p','p_'),'\mathrm',''),j,strrep(strrep(latex(stabp(i).eq(j,2)),'p','p_'),'\mathrm',''),node_width+round(0.02*size(char(stabp(i).eq(j,1)),2),1),node_height);
                            fprintf(fid2,'}\n');
                        else
                            fprintf(fid2,'subgraph cluster%0.f%0.f%.0f {\n texlbl="$\\begin{matrix}\\textbf{Unknown}\\\\x\\in D_u \\end{matrix}$" \n graph[style="dashed,rounded"] \n',I,i,1);
                            fprintf(fid2,'"eq%.0f%.0f%.0f" [fillcolor=white color=red texlbl="$\\begin{matrix}\\textbf{unstable}\\\\x_1^{*_%.0f}=%s\\\\x_2^{*_%.0f}=%s\\end{matrix}$" width=%.2f height=%.2f];\n',I,i,j,j,strrep(latex(stabp(i).eq(j,1)),'\mathrm',''),j,strrep(latex(stabp(i).eq(j,2)),'\mathrm',''),node_width+round(0.02*size(char(stabp(i).eq(j,1)),2),1),node_height);
                            fprintf(fid2,'}\n');
                        end
                    end
                end
            end
            fprintf(fid,'}\n');
            fprintf(fid2,'}\n');

            % Create DOT file x-mode connection arrows
            if size(stabp(i).p,1) > 1
                for j = 1:size(stabp(i).p,1)
                    for h = 1:size(stabp(i).p,1)
                        if j == h
                            continue
                        else
                            if strcmp(node(I).p(i).stab(h),'unstable') && strcmp(node(I).p(i).stab(j),'unstable')
                                fprintf(fid,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[color=white] \n',I,i,j,I,i,h);
                                fprintf(fid2,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[color=white] \n',I,i,j,I,i,h);
                            elseif strcmp(node(I).p(i).stab(h),'unstable')
                                fprintf(fid,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[lhead=cluster%0.f%.0f%.0f,label="x(%s)"]; \n',I,i,j,I,i,h,I,i,1,stabp(i).dom(h,:));
                                fprintf(fid2,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[lhead=cluster%0.f%.0f%.0f,label="x\\in %s"]; \n',I,i,j,I,i,h,I,i,1,'D_u');
                            elseif strcmp(node(I).p(i).stab(j),'unstable')
                                fprintf(fid,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[ltail=cluster%0.f%.0f%.0f,label="x(%s)"]; \n',I,i,j,I,i,h,I,i,1,stabp(i).dom(h,:));
                                fprintf(fid2,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[ltail=cluster%0.f%.0f%.0f,label="x\\in %s"]; \n',I,i,j,I,i,h,I,i,1,[strrep(stabp(i).dom(h,:),'D','D_{') '}']);
                            else
                                fprintf(fid,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[label="x(%s)"] \n',I,i,j,I,i,h,stabp(i).dom(h,:));
                                fprintf(fid2,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[label="x\\in %s"] \n',I,i,j,I,i,h,[strrep(stabp(i).dom(h,:),'D','D_{') '}']);
                            end
                        end
                    end
                end
            end
        end
        fclose(fid);
        fclose(fid2);

        % Create DOT file p-mode connection arrows
        fid = fopen(dotfile,'a');
        fid2 = fopen(dotfile_tex,'a');
        if num_event > 0
            for i = 1:num_event
                fprintf(fid,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[ltail=cluster%.0f%.0f,lhead=cluster%.0f%.0f,label="%s"]; \n',I,i,1,I,i+1,1,I,i,I,i+1,labp{i+1});
                if size(F,1) > 1
                    fprintf(fid2,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[ltail=cluster%.0f%.0f,lhead=cluster%.0f%.0f,label="%s"]; \n',I,i,1,I,i+1,1,I,i,I,i+1,strrep(labp{i+1},'p','p_'));
                    fprintf(fid2,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[ltail=cluster%.0f%.0f,lhead=cluster%.0f%.0f,label="%s"]; \n',I,i+1,1,I,i,1,I,i+1,I,i,strrep(labp{i},'p','p_'));
                else
                    fprintf(fid2,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[ltail=cluster%.0f%.0f,lhead=cluster%.0f%.0f,label="%s"]; \n',I,i,1,I,i+1,1,I,i,I,i+1,labp{i+1});
                    fprintf(fid2,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[ltail=cluster%.0f%.0f,lhead=cluster%.0f%.0f,label="%s"]; \n',I,i+1,1,I,i,1,I,i+1,I,i,labp{i});
                end
                fprintf(fid,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[ltail=cluster%.0f%.0f,lhead=cluster%.0f%.0f,label="%s"]; \n',I,i+1,1,I,i,1,I,i+1,I,i,labp{i});
            end
        end
        fprintf(fid,'}\n');
        fprintf(fid2,'}\n');
        fclose(fid);
        fclose(fid2);

    %% Fixed-p Case -------------------------------------------------------
    else

        % Set fixed-p value, find relevant p-domain and number of branches
        p = P;
        fun.fp = strrep(fun.gv,sprintf('p%.0f',active_parameter),sprintf('%.4g',p));
        pdom_idx = find(p>pdom(:,1) & p<pdom(:,2));
        num_branch_fp = num_branch_p(pdom_idx);

        % Interpolate to obtain x-values for fixed-p value
        xval = zeros(num_branch_fp,2);
        pval = zeros(num_branch_fp,1);
        eigval = zeros(num_branch_fp,2);
        for i = 1:num_branch_fp
            branch2 = [pdat(pdom_idx,i).p pdat(pdom_idx,i).x pdat(pdom_idx,i).eigval];
            [~,order] = sort(branch2(:,1));
            branch2 = branch2(order,:);
            pdiff = (branch2(:,1)-p);
            [ptemp,ploc] = min(abs(pdiff));
            if ptemp < ztol
                xval(i,1:2) = branch2(ploc,2:3);
                pval(i,1) = branch2(ploc,1);
                eigval(i,1:2) = branch2(ploc,4:5);
            else
                if (ploc == 1 && p < branch2(ploc,1)) || (ploc == size(branch2,1) && p > branch2(ploc,1))
                    xval(i,1:2) = NaN;
                    pval(i,1) = NaN;
                    eigval(i,1:2) = NaN;
                else
                    if p > branch2(ploc,1)
                        ploc2 = ploc + 1;
                    elseif p < branch2(ploc,1)
                        ploc2 = ploc - 1;
                    end
                    pval(i,1) = p;
                    xval(i,1) = interp1([branch2(ploc,1) branch2(ploc2,1)],[branch2(ploc,2) branch2(ploc2,2)],p);
                    xval(i,2) = interp1([branch2(ploc,1) branch2(ploc2,1)],[branch2(ploc,3) branch2(ploc2,3)],p);
                    eigval(i,1) = interp1([branch2(ploc,1) branch2(ploc2,1)],[branch2(ploc,4) branch2(ploc2,4)],p);
                    eigval(i,2) = interp1([branch2(ploc,1) branch2(ploc2,1)],[branch2(ploc,5) branch2(ploc2,5)],p);
                end
            end
        end

        % Create equilibrium data for fixed-p value, remove unused branches
        equil = [pval xval eigval];
        equil = equil(~any(isnan(equil),2),:);
        numeq = size(equil,1);
        node(I).p(pdom_idx).equil = equil;
        node(I).p(pdom_idx).pd = active_parameter;
        node(I).fp_idx = pdom_idx;
        node(I).fp = p;

        % Add to DOT file
        fid = fopen(dotfile,'a');
        fid2 = fopen(dotfile_tex,'a');
        fprintf(fid,'subgraph cluster%0.f%0.f {\n fillcolor=white; \n',I,I);
        fprintf(fid2,'subgraph cluster%0.f%0.f {\n fillcolor=white; \n',I,I);
        fprintf(fid,'label=<%s <br/>p%.0f = %.4g > \n;',strrep(strrep(labp{pdom_idx},'<','&lt;'),'>','&gt;'),active_parameter,p);
        if size(F,1) > 1
            fprintf(fid2,'texlbl="$\\begin{matrix}%s\\\\p_%.0f=%.4g\\end{matrix}$"\n',strrep(strrep(labp{pdom_idx},'p','p_'),'Inf','\infty'),active_parameter,p);
        else
            fprintf(fid2,'texlbl="$\\begin{matrix}%s\\\\p=%.4g\\end{matrix}$"\n',strrep(labp{pdom_idx},'Inf','\infty'),p);
        end
        fprintf(fid,'node [style=filled] \n');
        fprintf(fid2,'node [style=filled] \n');

        % Create DOT file nodes
        if numeq == 0
            fprintf(fid,'"eq%.0f1%.0f" [fillcolor=grey label="no equilbrium" ]; \n',I,1);
            fprintf(fid2,'"eq%.0f1%.0f" [style=dashed label="\\textbf{no equilibrium}" width=%.2f height=%.2f]; \n',I,1,node_width,node_height);
        else
            for i = 1:numeq
                if equil(i,4) > 0 || equil(i,5) > 0
                    fprintf(fid,'subgraph cluster%0.f%0.f {\n label="Unknown" \n graph[style="dashed,rounded"] \n',I,11);
                    fprintf(fid,'"eq%0.f1%0.f" [fillcolor=red label="unstable \n x1* = %.4g \n x2* = %.4g \n"]\n',I,i,equil(i,2),equil(i,3));
                    fprintf(fid,'}\n');
                    fprintf(fid2,'subgraph cluster%0.f%0.f {\n texlbl="$\\begin{matrix}\\textbf{Unknown}\\\\x\\in D_u \\end{matrix}$" \n graph[style="dashed,rounded"] \n',I,11);
                    fprintf(fid2,'"eq%0.f1%0.f" [fillcolor=white color=red texlbl="$\\begin{matrix}\\textbf{unstable}\\\\x_1^{*_%.0f}=%.4g\\\\x_2^{*_%.0f}=%.4g\\end{matrix}$" width=%.2f height=%.2f];\n',I,i,i,equil(i,2),i,equil(i,3),node_width2,node_height2);
                    fprintf(fid2,'}\n');
                elseif equil(i,4) <= 0 && equil(i,5) <= 0
                    fprintf(fid,'"eq%0.f1%0.f" [fillcolor=green label="stable \n x1* = %.4g \n x2* = %.4g \n x (Ds) \n"]\n',I,i,equil(i,2),equil(i,3));
                    fprintf(fid2,'"eq%0.f1%0.f" [fillcolor=white color=green4 texlbl="$\\begin{matrix}\\textbf{stable}\\\\x_1^{*_%.0f}=%.4g\\\\x_2^{*_%.0f}=%.4g\\\\x\\in %s\\end{matrix}$" width=%.2f height=%.2f];\n',I,i,i,equil(i,2),i,equil(i,3),[strrep(stabp(pdom_idx).dom(i,:),'D','D_{') '}'],node_width,node_height);
                end
            end
        end
        fprintf(fid,'}\n');
        fprintf(fid2,'}\n');

        % Create DOT file x-mode connection arrows
        if numeq > 1
            for j = 1:numeq
                for h = 1:numeq
                    if j == h
                        continue
                    else
                        if strcmp(node(I).p(pdom_idx).stab(h),'unstable') && strcmp(node(I).p(pdom_idx).stab(j),'unstable')
                            fprintf(fid,'"eq%.0f1%.0f" -> "eq%.0f1%.0f"[color=white] \n',I,j,I,h);
                            fprintf(fid2,'"eq%.0f1%.0f" -> "eq%.0f1%.0f"[color=white] \n',I,j,I,h);
                        elseif strcmp(node(I).p(pdom_idx).stab(h),'unstable')
                            fprintf(fid,'"eq%.0f1%.0f" -> "eq%.0f1%.0f"[lhead=cluster%0.f%.0f,label="x(%s)"]; \n',I,j,I,h,I,11,stabp(pdom_idx).dom(h,:));
                            fprintf(fid2,'"eq%.0f1%.0f" -> "eq%.0f1%.0f"[lhead=cluster%0.f%.0f,label="x\\in %s"]; \n',I,j,I,h,I,11,'D_u');
                        elseif strcmp(node(I).p(pdom_idx).stab(j),'unstable')
                            fprintf(fid,'"eq%.0f1%.0f" -> "eq%.0f1%.0f"[ltail=cluster%0.f%.0f,label="x(%s)"]; \n',I,j,I,h,I,11,stabp(pdom_idx).dom(h,:));
                            fprintf(fid2,'"eq%.0f1%.0f" -> "eq%.0f1%.0f"[ltail=cluster%0.f%.0f,label="x\\in %s"]; \n',I,j,I,h,I,11,[strrep(stabp(pdom_idx).dom(h,:),'D','D_{') '}']);
                        else
                            fprintf(fid,'"eq%.0f1%.0f" -> "eq%.0f1%.0f"[label="x(%s)"] \n',I,j,I,h,stabp(pdom_idx).dom(h,:));
                            fprintf(fid2,'"eq%.0f1%.0f" -> "eq%.0f1%.0f"[label="x\\in %s"] \n',I,j,I,h,[strrep(stabp(pdom_idx).dom(h,:),'D','D_{') '}']);
                        end
                    end
                end
            end
        end
        fprintf(fid,'}\n');
        fprintf(fid2,'}\n');
        fclose(fid);
        fclose(fid2);
    end
end

% Create DOT file top-mode connection arrows (all top-modes attached)
fid = fopen(dotfile,'a');
fid2 = fopen(dotfile_tex,'a');
if size(F,1) > 1
    for i = 1:size(F,1)
        for k = 1:size(F,1)
            if k == i
                continue
            else
                fprintf(fid,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[ltail=cluster%0.f,lhead=cluster%0.f,label="%s"]; \n',i,1,1,k,1,1,i,k,node(k).PD);
                fprintf(fid2,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[ltail=cluster%0.f,lhead=cluster%0.f,label="%s"]; \n',i,1,1,k,1,1,i,k,strrep(strrep(node(k).PD,'p','\\;p_'),'=','$=$'));
            end
        end
    end
    fprintf(fid,'}\n');
    fprintf(fid2,'}\n');
end
fclose(fid);
fclose(fid2);

%% Run GraphViz on DOT file, open resulting image file --------------------
imageFile = [label '.png'];
system(sprintf('dot.exe -Tpng -Gsize="10,10" "%s" -o"%s"',dotfile,imageFile));
system(imageFile);

%% Run GraphViz with dot2tex on DOT file, open resulting image file -------
if exist('opt','var')
    if exist('texlab','var')
        fid2 = fopen(dotfile_tex,'r');
        fid3 = fopen(dotfile_texR,'a');
        tline = fgetl(fid2);
        while ischar(tline)
            states = [];
            for j = 1:size(texlab,1)
                if ~isempty(strfind(texlab{j,1},'x'))
                    states = [states texlab{j,2} ', '];
                end
            end
            tline = strrep(tline,'ex','eee');
            tline = strrep(tline,'ix','iii');
            tline = strrep(tline,'ap','aaa');
            tline = strrep(tline,'mp','mmm'); tline = strrep(tline,'ep','EEE');
            for i = 1:size(texlab,1)
                tline = strrep(tline,texlab{i,1},texlab{i,2});
                tline = strrep(tline,[texlab{i,1}(1) '}' texlab{i,1}(2:end)],[texlab{i,2} '}' ]);
                tline = strrep(tline,'x\in',[states(1:end-2) '\in']);
            end
            tline = strrep(tline,'eee','ex');
            tline = strrep(tline,'iii','ix');
            tline = strrep(tline,'aaa','ap');
            tline = strrep(tline,'mmm','mp'); tline = strrep(tline,'EEE','ep');
            fprintf(fid3,'%s\n',tline);
            tline = fgetl(fid2);
        end
        fclose(fid2);
        fclose(fid3);
        system(sprintf('dot2tex -tmath --preproc --crop --figpreamble="\\Large" --autosize %sTexR.dot > %sTexR.tex',label,label));
        system(sprintf('pdflatex %sTexR.tex',label))
        winopen(sprintf('%sTexR.pdf',label))
    elseif (isa(opt,'struct') && isfield(opt,'dot2tex') && ~isfield(opt,'texlab')) || isa(opt,'char')
        system(sprintf('dot2tex -tmath --preproc --crop --figpreamble="\\Large" --autosize %sTex.dot > %sTex.tex',label,label));
        system(sprintf('pdflatex %sTex.tex',label))
        winopen(sprintf('%sTex.pdf',label))
    end
end

%% Create executable function file ----------------------------------------
funFile = [label '.m'];
fid4 = fopen(funFile,'w');
fprintf(fid4,'function dxdt = %s(x,u,p)\n\n',label);

% Load basin information, if available
if ~isempty(P)
    fprintf(fid4,'try load %s_basins.mat\nend\n\n',label);
end

if size(F,1) == 1
    fprintf(fid4,'dxdt = %s;\n',strrep(func2str(node(I).func),'@(x,u,p)',''));
    fprintf(fid4,'end\n');
else
    for i = 1:size(PD,1)

        % Print parameter conditions
        pvar = find(isnan(PD(i,:)));
        cond = strrep(node(i).PD,' p',' && p');
        cond = strrep(cond,'=','==');
        cond = strrep(cond,',','');
        for j = 1:size(PD,2)
            cond = strrep(cond,sprintf('p%.0f',j),sprintf('p(%.0f)',j));
        end

        % Print subset dynamics
        fprintf(fid4,'if %s\n',cond);
        fprintf(fid4,'\tdxdt = %s;\n',strrep(strrep(func2str(node(i).func),sprintf('@(x,u,p%.0f)',pvar),''),sprintf('p%.0f',pvar),sprintf('p(%.0f)',pvar)));
        
        % Print stability basin requirements
        if ~isempty(P)
            fprintf(fid4,'\tif %s==%g && exist(''basin'',''var'')\n',sprintf('p(%.0f)',pvar),P);
            fprintf(fid4,'\t\tif isfield(basin(%.0f),''Ds'')\n',i);
            fprintf(fid4,'\t\t\tif inpolygon(x(1),x(2),basin(%.0f).xb(:,1),basin(%.0f).xb(:,2))\n',i,i);
            fprintf(fid4,'\t\t\t\tfprintf(''Stable: converges to (%%.2f, %%.2f)\\n'',basin(%.0f).eq(1),basin(%.0f).eq(2))\n',i,i);
            fprintf(fid4,'\t\t\telse\n\t\t\t\tfprintf(''Unknown behavior\\n'')\n\t\t\tend\n');
            fprintf(fid4,'\t\telse\n\t\t\tfprintf(''No stable equilibria'')\n');
            fprintf(fid4,'\t\tend\n\tend\n');
        end
        fprintf(fid4,'end\n');
    end
    fprintf(fid4,'if ~exist(''dxdt'',''var'')\n\tdisp(''No data for this parameter designation'')\nend');
    fprintf(fid4,'\nend\n');
end
fclose(fid4);
