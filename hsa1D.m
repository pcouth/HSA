%--- University of Washington, Department of Aeronautics & Astronautics ---
%---------- Advanced Dynamics, Validation & Control Research Lab ----------
%
% Hybrid Stability Automaton (HSA) generator. Extracts stability mode data
% from COSY numberical continuation tool output, then creates and executes
% a GraphViz file to produce a PNG of the HSA. Also creates mode data and 
% an executable function m-file that provides stability information.
%
% One-dimensional version
%
% Author: Peter Uth
% Created: February 2017
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
function [node] = hsa1D(label,F,B,BP,P,opt)

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
fun.origgvtex = strrep(fun.orig,'@(x,u,p)','\dot{x}$=$');
fun.origgvtex = strrep(fun.origgvtex,'*','');
if size(F,1) > 1
    for i = 1:size(PD,2)
        fun.origgvtex = strrep(fun.origgvtex,sprintf('p(%.0f)',i),sprintf('p_%.0f',i));
    end
end

% Set standard width and height for GraphViz nodes
node_width = 2.75;
node_height = 1.5;

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
        node_width node_height opt dotfile_texR texlab
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
        fun.gvtex = strrep(func2str(node(I).func),sprintf('@(x,u,p%.0f)',active_parameter),'\dot{x} = ');
        fun.gvtex = strrep(fun.gvtex,sprintf('p%.0f',active_parameter),sprintf('p_%.0f',active_parameter));
        fun.gvtex = strrep(fun.gvtex,'*','');
        fprintf(fid,'label=<%s<br/>%s<br/>> \n',fun.gv,inactive_parameters);
        fprintf(fid2,'texlbl="$\\begin{matrix}%s\\\\%s\\end{matrix}$"\n',fun.gvtex,inactive_parameters_tex);
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
        if max(sign(branch(i+split).eigval)) * min(sign(branch(i+split).eigval)) < 0
            for k = 1:numel(branch(i+split).x) - 1
                if sign(branch(i+split).eigval(k)) * sign(branch(i+split).eigval(k+1)) < 0
                    [~,loc] = min([abs(b(i+split).eigval(k)) abs(b(i+split).eigval(k+1))]);
                    event = [event;b(i+split).p(k+loc-1) b(i+split).x(k+loc-1)];
                end
            end
            evtemp(1,:) = [-Inf -Inf];
            evtemp(2:size(event,1)+1,:) = event;
            evtemp(end+1,:) = [Inf Inf];
            for k = 1:size(event,1)
                xtemp = branch(i+split).x;
                ptemp = branch(i+split).p;
                eigtemp = branch(i+split).eigval;
                branch(i+split).x = xtemp(evtemp(k,2) < xtemp & xtemp < evtemp(k+1,2));
                branch(i+split).p = ptemp(evtemp(k,2) < xtemp & xtemp < evtemp(k+1,2));
                branch(i+split).eigval = eigtemp(evtemp(k,2) < xtemp & xtemp < evtemp(k+1,2));
                split = split + 1;
                branch(i+split).x = xtemp(evtemp(k+1,2) < xtemp & xtemp < evtemp(k+2,2));
                branch(i+split).p = ptemp(evtemp(k+1,2) < xtemp & xtemp < evtemp(k+2,2));
                branch(i+split).eigval = eigtemp(evtemp(k+1,2) < xtemp & xtemp < evtemp(k+2,2));
            end
        end
    end
    num_branch = numel(branch);

    % Add COSY bp branch point output data to event data
    for i = 1:numel(bp)
        bp(i).x(abs(bp(i).x) < ztol) = 0;
        bp(i).p(abs(bp(i).p) < ztol) = 0;
        event(i+split,:) = [bp(i).p bp(i).x];
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
            xtemp = branch(j).x(pdom(i,1) < branch(j).p & branch(j).p < pdom(i,2));
            ptemp = branch(j).p(pdom(i,1) < branch(j).p & branch(j).p < pdom(i,2));
            eigtemp = branch(j).eigval(pdom(i,1) < branch(j).p & branch(j).p < pdom(i,2));
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
            stabp(i).p(j,:) = [median(pdat(i,j).x) median(pdat(i,j).p) median(pdat(i,j).eigval)];
        end
        if ~isempty(stabp(i).p)
            [~,order] = sort(stabp(i).p(:,1));
            stabp(i).p = stabp(i).p(order,:);
        end
    end

    % Determine x-domains for each p-domain
    for i = 1:num_event + 1
        if size(stabp(i).p,1) == 1
            stabp(i).dom{1,:} = {'-Inf' 'Inf'};
        elseif size(stabp(i).p,1) == 0
            stabp(i).dom{1,:} = [];
        else
            for j = 1:size(stabp(i).p,1)
                if stabp(i).p(j,3) > 0
                    if j == 1
                        stabp(i).dom{j,:} = {'-Inf' sprintf('xeq%.0f',j)};
                    elseif j == size(stabp(i).p,1)
                        stabp(i).dom{j,:} = {sprintf('xeq%.0f',j) 'Inf'};
                    else
                        stabp(i).dom{j,:} = {sprintf('xeq%.0f',j) sprintf('xeq%.0f',j)};
                    end
                elseif stabp(i).p(j,3) < 0
                    if j == 1
                        stabp(i).dom{j,:} = {'-Inf' sprintf('xeq%.0f',j+1)};
                    elseif j == size(stabp(i).p,1)
                        stabp(i).dom{j,:} = {sprintf('xeq%.0f',j-1) 'Inf'};
                    else
                        stabp(i).dom{j,:} = {sprintf('xeq%.0f',j-1) sprintf('xeq%.0f',j+1)};
                    end
                end
            end
        end
        node(I).p(i).xdom = stabp(i).dom';
    end
    
    % Determine solution equations
    soltemp = solve(node(I).func);
    sol = [];
    for i = 1:numel(soltemp)
        if isempty(strfind(strrep(char(soltemp(i,1)),'pi',''),'i'))
            sol = [sol;soltemp(i,1)];
        end
    end

    % Use average branch values to match equations to nodes
    for i = 1:num_event + 1
        if size(stabp(i).p,1) == 0
            stabp(i).eq = [];
        end
        for j = 1:size(stabp(i).p,1)
            eval(sprintf('p%.0f=stabp(i).p(j,2);',active_parameter))
            soltemp = eval(sol);
            [~,loc] = min(abs(stabp(i).p(j,1)-soltemp));
            stabp(i).eq(1,j) = sol(loc);
        end
        node(I).p(i).eq = stabp(i).eq;

        % Determine branch stability
        if isempty(stabp(i).p)
            node(I).p(i).stab = [];
        else
            for j = 1:size(stabp(i).p,1)
                if stabp(i).p(j,3) < 0
                    node(I).p(i).stab{1,j} = 'stable';
                elseif stabp(i).p(j,3) > 0
                    node(I).p(i).stab{1,j} = 'unstable';
                end
            end
        end
    end

    %% General Case -------------------------------------------------------
    if isempty(P)
        
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
                    if stabp(i).p(j,3) < 0
                        fprintf(fid,'"eq%.0f%.0f%.0f" [fillcolor=green label="stable \n x* = %s \n x(%s,%s)" ]; \n',I,i,j,char(stabp(i).eq(j)),stabp(i).dom{j,1}{1},stabp(i).dom{j,1}{2});
                        if size(F,1) > 1
                            fprintf(fid2,'"eq%.0f%.0f%.0f" [fillcolor=white color=green4 texlbl="$\\begin{matrix}\\textbf{stable}\\\\x^{*_%.0f}=%s\\\\x\\in(%s,%s)\\end{matrix}$" width=%.2f height=%.2f];\n',I,i,j,j,strrep(strrep(latex(stabp(i).eq(j)),'p','p_'),'\mathrm',''),strrep([strrep(stabp(i).dom{j,1}{1},'xeq','x^{*_') '}'],'Inf}','\infty'),strrep([strrep(stabp(i).dom{j,1}{2},'xeq','x^{*_') '}'],'Inf}','\infty'),node_width,node_height);
                        else
                            fprintf(fid2,'"eq%.0f%.0f%.0f" [fillcolor=white color=green4 texlbl="$\\begin{matrix}\\textbf{stable}\\\\x^{*_%.0f}=%s\\\\x\\in(%s,%s)\\end{matrix}$" width=%.2f height=%.2f];\n',I,i,j,j,strrep(latex(stabp(i).eq(j)),'\mathrm',''),strrep([strrep(stabp(i).dom{j,1}{1},'xeq','x^{*_') '}'],'Inf}','\infty'),strrep([strrep(stabp(i).dom{j,1}{2},'xeq','x^{*_') '}'],'Inf}','\infty'),node_width,node_height);
                        end
                    elseif stabp(i).p(j,3) > 0
                        fprintf(fid,'"eq%.0f%.0f%.0f" [fillcolor=red label="unstable \n x* = %s \n x(%s,%s)" ]; \n',I,i,j,char(stabp(i).eq(j)),stabp(i).dom{j,1}{1},stabp(i).dom{j,1}{2});
                        if size(F,1) > 1
                            fprintf(fid2,'"eq%.0f%.0f%.0f" [fillcolor=white color=red texlbl="$\\begin{matrix}\\textbf{unstable}\\\\x^{*_%.0f}=%s\\\\x\\in(%s,%s)\\end{matrix}$" width=%.2f height=%.2f];\n',I,i,j,j,strrep(strrep(latex(stabp(i).eq(j)),'p','p_'),'\mathrm',''),strrep([strrep(stabp(i).dom{j,1}{1},'xeq','x^{*_') '}'],'Inf}','\infty'),strrep([strrep(stabp(i).dom{j,1}{2},'xeq','x^{*_') '}'],'Inf}','\infty'),node_width,node_height);
                        else
                            fprintf(fid2,'"eq%.0f%.0f%.0f" [fillcolor=white color=red texlbl="$\\begin{matrix}\\textbf{unstable}\\\\x^{*_%.0f}=%s\\\\x\\in(%s,%s)\\end{matrix}$" width=%.2f height=%.2f];\n',I,i,j,j,strrep(latex(stabp(i).eq(j)),'\mathrm',''),strrep([strrep(stabp(i).dom{j,1}{1},'xeq','x^{*_') '}'],'Inf}','\infty'),strrep([strrep(stabp(i).dom{j,1}{2},'xeq','x^{*_') '}'],'Inf}','\infty'),node_width,node_height);
                        end
                    end
                end
            end
            fprintf(fid,'}\n');
            fprintf(fid2,'}\n');

            % Create DOT file x-mode connection arrows
            if size(stabp(i).p,1) > 1
                for j = 1:size(stabp(i).p,1)
                    if strcmp(stabp(i).dom{j,1}{1},'-Inf')
                        fprintf(fid,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[label="x<%s"] \n',I,i,j+1,I,i,j,stabp(i).dom{j,1}{2});
                        fprintf(fid2,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[label="x<%s"] \n',I,i,j+1,I,i,j,strrep([stabp(i).dom{j,1}{2} '}'],'xeq','x^{*_'));
                    elseif strcmp(stabp(i).dom{j,1}{2},'Inf')
                        fprintf(fid,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[label="x>%s"] \n',I,i,j-1,I,i,j,stabp(i).dom{j,1}{1});
                        fprintf(fid2,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[label="x>%s"] \n',I,i,j-1,I,i,j,strrep([stabp(i).dom{j,1}{1} '}'],'xeq','x^{*_'));
                    elseif strcmp(stabp(i).dom{j,1}{1},stabp(i).dom{j,1}{2})
                        fprintf(fid,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[label="x=%s"] \n',I,i,j+1,I,i,j,stabp(i).dom{j,1}{1});
                        fprintf(fid2,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[label="x=%s"] \n',I,i,j+1,I,i,j,strrep([stabp(i).dom{j,1}{1} '}'],'xeq','x^{*_'));
                        fprintf(fid,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[label="x=%s"] \n',I,i,j-1,I,i,j,stabp(i).dom{j,1}{1});
                        fprintf(fid2,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[label="x=%s"] \n',I,i,j-1,I,i,j,strrep([stabp(i).dom{j,1}{1} '}'],'xeq','x^{*_'));
                    else
                        fprintf(fid,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[label="%s<x<%s"] \n',I,i,j-1,I,i,j,stabp(i).dom{j,1}{1},stabp(i).dom{j,1}{2});
                        fprintf(fid2,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[label="%s<x<%s"] \n',I,i,j-1,I,i,j,stabp(i).dom{j,1}{1},strrep([stabp(i).dom{j,1}{2} '}'],'xeq','x^{*_'));
                        fprintf(fid,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[label="%s<x<%s"] \n',I,i,j+1,I,i,j,stabp(i).dom{j,1}{1},stabp(i).dom{j,1}{2});
                        fprintf(fid2,'"eq%.0f%.0f%.0f" -> "eq%.0f%.0f%.0f"[label="%s<x<%s"] \n',I,i,j+1,I,i,j,stabp(i).dom{j,1}{1},strrep([stabp(i).dom{j,1}{2} '}'],'xeq','x^{*_'));
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
        pdom_idx = find(p > pdom(:,1) & p < pdom(:,2));
        num_branch_fp = num_branch_p(pdom_idx);

        % Interpolate to obtain x-values for fixed-p value
        xval = zeros(num_branch_fp,1);
        pval = zeros(num_branch_fp,1);
        eigval = zeros(num_branch_fp,1);
        for i = 1:num_branch_fp
            branch2 = [pdat(pdom_idx,i).x pdat(pdom_idx,i).p pdat(pdom_idx,i).eigval];
            [~,order] = sort(branch2(:,2));
            branch2 = branch2(order,:);
            pdiff = (branch2(:,2)-p);
            [ptemp,ploc] = min(abs(pdiff));
            if ptemp == 0
                xval(i,1) = branch2(ploc,1);
                pval(i,1) = branch2(ploc,2);
                eigval(i,1) = branch2(ploc,3);
            else
                if (ploc == 1 && p < branch2(ploc,2)) || (ploc == size(branch2,1) && p > branch2(ploc,2))  
                    xval(i,1) = NaN;
                    pval(i,1) = NaN;
                    eigval(i,1) = NaN;
                else
                    if p > branch2(ploc,2)
                        ploc2 = ploc + 1;
                    elseif p < branch2(ploc,2)
                        ploc2 = ploc - 1;
                    end
                    pval(i,1) = p;
                    xval(i,1) = interp1([branch2(ploc,2) branch2(ploc2,2)],[branch2(ploc,1) branch2(ploc2,1)],p);
                    eigval(i,1) = interp1([branch2(ploc,2) branch2(ploc2,2)],[branch2(ploc,3) branch2(ploc2,3)],p);
                end
            end
        end

        % Create equilibrium data for fixed-p value, remove unused branches
        equil = [xval pval eigval];
        equil = equil(~any(isnan(equil),2),:);
        [~,order] = sort(equil(:,1));
        equil = equil(order,:);
        numeq = size(equil,1);

        % Determine x-domain for each equilibrium
        domain = zeros(numeq,2);
        if numeq == 1
            domain = [-Inf Inf];
        else
            for i = 1:numeq
                if equil(i,3) > 0
                    if i == 1
                        domain(i,:) = [-Inf equil(i,1)];
                    elseif i == numeq
                        domain(i,:) = [equil(i,1) Inf];
                    else
                        domain(i,:) = [equil(i,1) equil(i,1)];
                    end
                elseif equil(i,3) < 0
                    if i == 1
                        domain(i,:) = [-Inf equil(i+1,1)];
                    elseif i == numeq
                        domain(i,:) = [equil(i-1,1) Inf];
                    else
                        domain(i,:) = [equil(i-1,1) equil(i+1,1)];
                    end
                end
            end
        end

        % Add to DOT file
        fid = fopen(dotfile,'a');
        fid2 = fopen(dotfile_tex,'a');
        fprintf(fid,'subgraph cluster%0.f%0.f {\n fillcolor=white; \n',I,i);
        fprintf(fid2,'subgraph cluster%0.f%0.f {\n fillcolor=white; \n',I,i);
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
                if equil(i,3) > 0
                    fprintf(fid,'"eq%0.f1%0.f" [fillcolor=red label="unstable \n x* = %.4g \n x (%.4g,%.4g) \n"]\n',I,i,equil(i,1),domain(i,1),domain(i,2));
                    fprintf(fid2,'"eq%0.f1%0.f" [fillcolor=white color=red texlbl="$\\begin{matrix}\\textbf{unstable}\\\\x^{*_%.0f}=%.4g\\\\x\\in(%s,%s)\\end{matrix}$" width=%.2f height=%.2f];\n',I,i,i,equil(i,1),strrep(num2str(domain(i,1)),'Inf','\infty'),strrep(num2str(domain(i,2)),'Inf','\infty'),node_width,node_height);
                elseif equil(i,3) < 0
                    fprintf(fid,'"eq%0.f1%0.f" [fillcolor=green label="stable \n x* = %.4g \n x (%.4g,%.4g) \n"]\n',I,i,equil(i,1),domain(i,1),domain(i,2));
                    fprintf(fid2,'"eq%0.f1%0.f" [fillcolor=white color=green4 texlbl="$\\begin{matrix}\\textbf{stable}\\\\x^{*_%.0f}=%.4g\\\\x\\in(%s,%s)\\end{matrix}$" width=%.2f height=%.2f];\n',I,i,i,equil(i,1),strrep(num2str(domain(i,1)),'Inf','\infty'),strrep(num2str(domain(i,2)),'Inf','\infty'),node_width,node_height);
                else 
                    fprintf(fid,'"eq%0.f1%0.f" [label="marginally stable"]\n',I,i);
                    fprintf(fid2,'"eq%0.f1%0.f" [label="\\textbf{marginally stable}" width=%.2f height=%.2f]; \n',I,i,node_width,node_height);
                end
            end
        end
        fprintf(fid,'}\n');
        fprintf(fid2,'}\n');

        % Create DOT file x-mode connection arrows
        if numeq > 1
            for i = 1:numeq
                if domain(i,1) == -Inf
                    fprintf(fid,'"eq%0.f1%0.f" -> "eq%0.f1%0.f"[label="x<%.4g"] \n',I,i+1,I,i,domain(i,2));
                    fprintf(fid2,'"eq%0.f1%0.f" -> "eq%0.f1%0.f"[label="x<%.4g"] \n',I,i+1,I,i,domain(i,2));
                elseif domain(i,2) == Inf
                    fprintf(fid,'"eq%0.f1%0.f" -> "eq%0.f1%0.f"[label="x>%.4g"] \n',I,i-1,I,i,domain(i,1));
                    fprintf(fid2,'"eq%0.f1%0.f" -> "eq%0.f1%0.f"[label="x>%.4g"] \n',I,i-1,I,i,domain(i,1));
                elseif domain(i,1) == domain(i,2)
                    fprintf(fid,'"eq%0.f1%0.f" -> "eq%0.f1%0.f"[label="x=%.4g"] \n',I,i+1,I,i,domain(i,1));
                    fprintf(fid2,'"eq%0.f1%0.f" -> "eq%0.f1%0.f"[label="x=%.4g"] \n',I,i+1,I,i,domain(i,1));
                    fprintf(fid,'"eq%0.f1%0.f" -> "eq%0.f1%0.f"[label="x=%.4g"] \n',I,i-1,I,i,domain(i,1));
                    fprintf(fid2,'"eq%0.f1%0.f" -> "eq%0.f1%0.f"[label="x=%.4g"] \n',I,i-1,I,i,domain(i,1));
                else
                    fprintf(fid,'"eq%0.f1%0.f" -> "eq%0.f1%0.f"[label="%.4g<x<%.4g"] \n',I,i-1,I,i,domain(i,1),domain(i,2));
                    fprintf(fid2,'"eq%0.f1%0.f" -> "eq%0.f1%0.f"[label="%.4g<x<%.4g"] \n',I,i-1,I,i,domain(i,1),domain(i,2));
                    fprintf(fid,'"eq%0.f1%0.f" -> "eq%0.f1%0.f"[label="%.4g<x<%.4g"] \n',I,i+1,I,i,domain(i,1),domain(i,2));
                    fprintf(fid2,'"eq%0.f1%0.f" -> "eq%0.f1%0.f"[label="%.4g<x<%.4g"] \n',I,i+1,I,i,domain(i,1),domain(i,2));
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
            tline = strrep(tline,'ex','eee');
            tline = strrep(tline,'ix','iii');
            tline = strrep(tline,'ap','aaa');
            tline = strrep(tline,'mp','mmm'); tline = strrep(tline,'ep','EEE');
            for i = 1:size(texlab,1)
                tline = strrep(tline,texlab{i,1},texlab{i,2});
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

if size(F,1) == 1
    fprintf(fid4,'dxdt = %s;\n',strrep(func2str(node(I).func),'@(x,u,p)',''));
    for j = 1:size(node.p,2)
        if j == 1
            fprintf(fid4,'if %s < p && p < %s\n',num2str(node.p(j).pdom(1)),num2str(node.p(j).pdom(2)));
        else
            fprintf(fid4,'elseif %s < p && p < %s\n',num2str(node.p(j).pdom(1)),num2str(node.p(j).pdom(2)));
        end
        for k = 1:size(node.p(j).eq,2)
            fprintf(fid4,'\txeq%.0f = %s;\n',k,char(node.p(j).eq(k)));
        end
        for k = 1:size(node.p(j).xdom,2)
            if ~isempty(node.p(j).xdom{k})
                fprintf(fid4,'\tif %s < x && x < %s\n',node.p(j).xdom{k}{1},node.p(j).xdom{k}{2});
                fprintf(fid4,'\t\txf = %s;\n',char(node.p(j).eq(k)));
                if strcmp('stable',node.p(j).stab{k})
                    fprintf(fid4,'\t\tfprintf(''Stable: converges to %%s\\n'',num2str(xf))\n');
                elseif strcmp('unstable',node.p(j).stab{k})
                    fprintf(fid4,'\t\tfprintf(''Unstable: diverges from %%s\\n'',num2str(xf))\n');
                end
                fprintf(fid4,'\tend\n');
            else
                fprintf(fid4,'\txf = [];\n\tdisp(''no equilibria'')\n');
            end
        end
        fprintf(fid4,'\tif ~exist(''xf'',''var'')\n\t\tdisp(''On unstable equilibrium'')\n\tend\n');
    end
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
        
        % Print stability information
        for j = 1:size(node(i).p,2)
            if j == 1
                fprintf(fid4,'\tif %s < p(%.0f) && p(%.0f) < %s\n',num2str(node(i).p(j).pdom(1)),pvar,pvar,num2str(node(i).p(j).pdom(2)));
            else
                fprintf(fid4,'\telseif %s < p(%.0f) && p(%.0f) < %s\n',num2str(node(i).p(j).pdom(1)),pvar,pvar,num2str(node(i).p(j).pdom(2)));
            end
            for k = 1:size(node(i).p(j).eq,2)
                fprintf(fid4,'\t\txeq%.0f = %s;\n',k,strrep(char(node(i).p(j).eq(k)),sprintf('p%.0f',pvar),sprintf('p(%.0f)',pvar)));
            end
            for k = 1:size(node(i).p(j).xdom,2)
                if ~isempty(node(i).p(j).xdom{k})
                    fprintf(fid4,'\t\tif %s < x && x < %s\n',node(i).p(j).xdom{k}{1},node(i).p(j).xdom{k}{2});
                    fprintf(fid4,'\t\t\txf = %s;\n',strrep(char(node(i).p(j).eq(k)),sprintf('p%.0f',pvar),sprintf('p(%.0f)',pvar)));
                    if strcmp('stable',node(i).p(j).stab{k})
                        fprintf(fid4,'\t\t\tfprintf(''Stable: converges to %%s\\n'',num2str(xf))\n');
                    elseif strcmp('unstable',node(i).p(j).stab{k})
                        fprintf(fid4,'\t\t\tfprintf(''Unstable: diverges from %%s\\n'',num2str(xf))\n');
                    end
                    fprintf(fid4,'\t\tend\n');
                else
                    fprintf(fid4,'\t\txf = [];\n\t\tdisp(''no equilibria'')\n');
                end
            end
            fprintf(fid4,'\t\tif ~exist(''xf'',''var'')\n\t\t\tdisp(''On unstable equilibrium'')\n\t\tend\n');
        end
        fprintf(fid4,'\telse\n\t\tdisp(''On bifurcation point'')\n\tend\n');
        fprintf(fid4,'end\n\n');
    end
    fprintf(fid4,'if ~exist(''dxdt'',''var'')\n\tdisp(''No data for this parameter designation'')\nend');
    fprintf(fid4,'\nend\n');
end
fclose(fid4);
