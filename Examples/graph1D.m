%--- University of Washington, Department of Aeronautics & Astronautics ---
%---------- Advanced Dynamics, Validation & Control Research Lab ----------
% Hybrid graph representation generator. Extracts stability node data from
% COSY numerical continuation tool output, then creates and executes a
% GraphViz file to produce a graph representation PNG. [created Feb 2017]

% Inputs 
% func  - system equation in the form @(x,u,p)p+x^2
% b, bp -  solution branch and branch point data from COSY
% p     - user defined p value, leave blank [] for general case
% label - string used to name output PNG file
% ztol  - abs() below this value will be =0, optional

% Note: add folder containing GraphViz "dot.exe" to environmental variables

function graph1D(func,b,bp,p,label,ztol)

fun.orig = func2str(func);                                  
fun.gv = strrep(fun.orig,'@(x,u,p)','x<SUP>''</SUP> = '); 	

%% extract solution branches and event points info from COSY outputs b, bp
if ~exist('ztol','var')
    ztol = 0.00001;     % default zero-tolerance
end

% separate unstable/stable within same branches ---------------------------
event = [];         
j = 0;                                              
for i = 1:numel(b)    
    b(1,i).x(abs(b(1,i).x)<ztol) = 0;               
    b(1,i).p(abs(b(1,i).p)<ztol) = 0;               
    b(1,i).eigval(abs(b(1,i).eigval)<ztol) = 0;     
    branch(i+j).x = b(1,i).x';
    branch(i+j).p = b(1,i).p';
    branch(i+j).eigval = b(1,i).eigval';
    if max(sign(branch(i+j).eigval))*min(sign(branch(i+j).eigval)) < 0   
        for k = 1:numel(branch(i+j).x)-1
            if sign(branch(i+j).eigval(k))*sign(branch(i+j).eigval(k+1)) < 0 
                [~,loc] = min([abs(b(i+j).eigval(k)) abs(b(i+j).eigval(k+1))]);
                event = [event;b(i+j).p(k+loc-1) b(i+j).x(k+loc-1)];
            end
        end
        evtemp(1,:) = [-inf -inf];
        evtemp(2:size(event,1)+1,:) = event;
        evtemp(end+1,:) = [inf inf];
        for k = 1:size(event,1)
            xtemp = branch(i+j).x; ptemp = branch(i+j).p; eigtemp = branch(i+j).eigval;	
            branch(i+j).x = xtemp(evtemp(k,2)<xtemp & xtemp<evtemp(k+1,2));
            branch(i+j).p = ptemp(evtemp(k,2)<xtemp & xtemp<evtemp(k+1,2));
            branch(i+j).eigval = eigtemp(evtemp(k,2)<xtemp & xtemp<evtemp(k+1,2));
            j = j+1;                              	
            branch(i+j).x = xtemp(evtemp(k+1,2)<xtemp & xtemp<evtemp(k+2,2));
            branch(i+j).p = ptemp(evtemp(k+1,2)<xtemp & xtemp<evtemp(k+2,2));
            branch(i+j).eigval = eigtemp(evtemp(k+1,2)<xtemp & xtemp<evtemp(k+2,2));
        end
    end
end
num_branch = numel(branch); 

% add COSY bp branch point output data to event data ----------------------
for i = 1:numel(bp)
    bp(i).x(abs(bp(i).x)<ztol) = 0;             
    bp(i).p(abs(bp(i).p)<ztol) = 0;           	
    event(i+j,:) = [bp(i).p bp(i).x];
end
num_event = numel(bp)+j;        

% create p-domain and p-label matrices ------------------------------------
if num_event == 0
    pdom = [-inf inf];
    labp = {'-Inf < p < Inf'};
else
    [~,order] = sort(event(:,1));       
    event = event(order,:);	
    pdom(1,:) = [-inf event(1,1)];
    pdom(num_event+1,:) = [event(end,1) inf];
    labp{1} = sprintf('p < %.4g',event(1,1));
    labp{num_event+1} = sprintf('p > %.4g',event(end,1));
    if num_event > 1
        for i = 2:num_event
            pdom(i,:) = [event(i-1,1) event(i,1)];
            labp{i} = sprintf('%.4g < p < %.4g',event(i-1,1),event(i,1));
        end
    end
end

% separate solution branch data into p-domains ----------------------------
for i = 1:num_event+1
    n = 0;
    for j = 1:num_branch
        xtemp = branch(j).x(pdom(i,1)<branch(j).p & branch(j).p<pdom(i,2));
        ptemp = branch(j).p(pdom(i,1)<branch(j).p & branch(j).p<pdom(i,2));
        eigtemp = branch(j).eigval(pdom(i,1)<branch(j).p & branch(j).p<pdom(i,2));
        if ~isempty(ptemp)
            n = n+1;
            pdat(i,n).x = xtemp;
            pdat(i,n).p = ptemp;
            pdat(i,n).eigval = eigtemp;
        end
    end
    pbnum(i,1) = n;     
    pfun{i,1} = strrep(func2str(func),'@(x,u,p)','');
end

% determine branch stability ordering for each p-domain -------------------
for i = 1:size(pbnum,1)
    stabp(i).p = [];
    for j = 1:pbnum(i,1)
        stabp(i).p(j,:) = [median(pdat(i,j).x) median(pdat(i,j).p) median(pdat(i,j).eigval)];
    end
    if ~isempty(stabp(i).p)
        [~,order] = sort(stabp(i).p(:,1));      
        stabp(i).p = stabp(i).p(order,:);       
    end    
end

% determine x domains for each p-domain -----------------------------------
for i = 1:num_event+1
    if size(stabp(i).p,1) == 1                  
    	stabp(i).dom{1,:} = {'-inf' 'inf'};
    else
        for j = 1:size(stabp(i).p,1)   
            if stabp(i).p(j,3) > 0              
                if j == 1
                    stabp(i).dom{j,:} = {'-inf' sprintf('xeq%.0f',j)};
                elseif j == size(stabp(i).p,1)
                    stabp(i).dom{j,:} = {sprintf('xeq%.0f',j) 'inf'};
                else
                    stabp(i).dom{j,:} = {sprintf('xeq%.0f',j) sprintf('xeq%.0f',j)};
                end     
            elseif stabp(i).p(j,3) < 0         
                if j == 1
                    stabp(i).dom{j,:} = {'-inf' sprintf('xeq%.0f',j+1)};
                elseif j == size(stabp(i).p,1)
                    stabp(i).dom{j,:} = {sprintf('xeq%.0f',j-1) 'inf'};
                else
                    stabp(i).dom{j,:} = {sprintf('xeq%.0f',j-1) sprintf('xeq%.0f',j+1)};
                end
            end
        end
    end
end

%% General-p case
if isempty(p)
    soltemp = solve(func);      
    sol = [];
    for i = 1:numel(soltemp)
        if isempty(strfind(char(soltemp(i,1)),'i'))
            sol = [sol;soltemp(i,1)];
        end
    end
            
    % match equations to nodes --------------------------------------------
    for i = 1:num_event+1
        for j = 1:size(stabp(i).p,1)   
            p = stabp(i).p(j,2);
            soltemp = eval(sol);
            [~,loc] = min(abs(stabp(i).p(j,1)-soltemp));
            stabp(i).eq(j,1) = sol(loc);
        end
    end

    % Create DOT file  ----------------------------------------------------
    dotFile = [label '.dot'];       
    fid = fopen(dotFile,'w');                   
    fprintf(fid,'digraph {\n');
    fprintf(fid,'compound=true;\n');
    fprintf(fid,'graph[style="rounded"]\n');
    fprintf(fid,'nodesep=0.75\n');
    fprintf(fid,'label=<%s <br/>> \n',fun.gv);
         
    for i = 1:num_event+1
        fprintf(fid,'subgraph cluster%.0f {\n',i);
        fprintf(fid,'label=<%s <br/>%s> \n',strrep(strrep(labp{i},'<','&lt;'),'>','&gt;'),fun.gv);
        fprintf(fid,'node [style=filled] \n');

    % DOT file plot nodes  ------------------------------------------------
        if isempty(stabp(i).p)
            fprintf(fid,'"eq%.0f%.0f" [fillcolor=grey label="no equilbrium" ]; \n',i,1);
        else
            for j = 1:size(stabp(i).p,1)   
                if stabp(i).p(j,3) < 0
                    fprintf(fid,'"eq%.0f%.0f" [fillcolor=green label="stable \n x* = %s \n x(%s,%s)" ]; \n',i,j,char(stabp(i).eq(j)),stabp(i).dom{j,1}{1},stabp(i).dom{j,1}{2});
                elseif stabp(i).p(j,3) > 0
                    fprintf(fid,'"eq%.0f%.0f" [fillcolor=red label="unstable \n x* = %s \n x(%s,%s)" ]; \n',i,j,char(stabp(i).eq(j)),stabp(i).dom{j,1}{1},stabp(i).dom{j,1}{2});
                end
            end
        end
        fprintf(fid,'}\n');
    end
  
    % Plot x-node transitions ---------------------------------------------
    for i = 1:num_event+1
        if size(stabp(i).p,1) > 1
            for j = 1:size(stabp(i).p,1)   
                if strcmp(stabp(i).dom{j,1}{1},'-inf')
                    fprintf(fid,'"eq%.0f%.0f" -> "eq%.0f%.0f"[label="x<%s"] \n',i,j+1,i,j,stabp(i).dom{j,1}{2});
                elseif strcmp(stabp(i).dom{j,1}{2},'inf')
                    fprintf(fid,'"eq%.0f%.0f" -> "eq%.0f%.0f"[label="x>%s"] \n',i,j-1,i,j,stabp(i).dom{j,1}{1});
                elseif strcmp(stabp(i).dom{j,1}{1},stabp(i).dom{j,1}{2})
                    fprintf(fid,'"eq%.0f%.0f" -> "eq%.0f%.0f"[label="x=%s"] \n',i,j+1,i,j,stabp(i).dom{j,1}{1});
                    fprintf(fid,'"eq%.0f%.0f" -> "eq%.0f%.0f"[label="x=%s"] \n',i,j-1,i,j,stabp(i).dom{j,1}{1});
                else
                    fprintf(fid,'"eq%.0f%.0f" -> "eq%.0f%.0f"[label="%s<x<%s"] \n',i,j-1,i,j,stabp(i).dom{j,1}{1},stabp(i).dom{j,1}{2});
                    fprintf(fid,'"eq%.0f%.0f" -> "eq%.0f%.0f"[label="%s<x<%s"] \n',i,j+1,i,j,stabp(i).dom{j,1}{1},stabp(i).dom{j,1}{2});
                end
            end
        end       
    end
    
    % Plot p-node transitions ---------------------------------------------
    if num_event > 0
        for i = 1:num_event
            fprintf(fid,'"eq%.0f%.0f" -> "eq%.0f%.0f"[ltail=cluster%.0f,lhead=cluster%.0f,label="%s"]; \n',i,1,i+1,1,i,i+1,labp{i+1});
            fprintf(fid,'"eq%.0f%.0f" -> "eq%.0f%.0f"[ltail=cluster%.0f,lhead=cluster%.0f,label="%s"]; \n',i+1,1,i,1,i+1,i,labp{i});
        end
    end
    fprintf(fid,'}\n');
    fclose(fid);
                                   
else
%% Fixed-p case
    fun.fp = strrep(fun.gv,'p',sprintf('%.4g',p));
    pdom_idx = find(p>pdom(:,1) & p<pdom(:,2));         
    numf_branch = pbnum(pdom_idx);                     
    xval = zeros(numf_branch,1); pval = zeros(numf_branch,1); eigval = zeros(numf_branch,1);  
       
    % Interpolates to obtain x-values for specified-p  --------------------
    for i = 1:numf_branch     
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
                    ploc2 = ploc+1;
                elseif p < branch2(ploc,2)
                    ploc2 = ploc-1;
                end
                pval(i,1) = p;
                xval(i,1) = interp1([branch2(ploc,2) branch2(ploc2,2)],[branch2(ploc,1) branch2(ploc2,1)],p);
                eigval(i,1) = interp1([branch2(ploc,2) branch2(ploc2,2)],[branch2(ploc,3) branch2(ploc2,3)],p);        
            end
        end
    end   
    equil = [xval pval eigval];                 
    equil = equil(~any(isnan(equil),2),:);      
    [~,order] = sort(equil(:,1));               
    equil = equil(order,:);                     
    numeq = size(equil,1);

    % Determine domain for each equilibrium branch  -----------------------
    domain = zeros(numeq,2);    
    if numeq == 1                               
        domain = [-inf inf];
    else
        for i = 1:numeq
            if equil(i,3) > 0                   
                if i == 1
                    domain(i,:) = [-inf equil(i,1)];
                elseif i == numeq
                    domain(i,:) = [equil(i,1) inf];
                else
                    domain(i,:) = [equil(i,1) equil(i,1)];
                end     
            elseif equil(i,3) < 0               
                if i == 1
                    domain(i,:) = [-inf equil(i+1,1)];
                elseif i == numeq
                    domain(i,:) = [equil(i-1,1) inf];
                else
                    domain(i,:) = [equil(i-1,1) equil(i+1,1)];
                end
            end
        end
    end

    % Create DOT file  ----------------------------------------------------
    dotFile = [label '.dot'];       
    fid = fopen(dotFile,'w');                   
    fprintf(fid,'digraph {\n');
    fprintf(fid,'compound=true;\n');
    fprintf(fid,'graph[style="rounded"]\n');
    fprintf(fid,'nodesep=0.75\n');
    fprintf(fid,'label=<%s <br/>p = %.4g> \n',fun.gv,p);
    
    fprintf(fid,'subgraph cluster1 {\n');
    fprintf(fid,'label=<%s <br/>%s> \n;',strrep(strrep(labp{pdom_idx},'<','&lt;'),'>','&gt;'),fun.fp);
    fprintf(fid,'node [style=filled]');

    % DOT file plot nodes  ------------------------------------------------
    if numeq == 0   
        fprintf(fid,'"no equilibria" [label="no equilibria"]\n');
    else
        for i = 1:numeq     
            if equil(i,3) > 0 
                fprintf(fid,'"eq%s" [fillcolor=red label="unstable \n x* = %.4g \n x (%.4g,%.4g) \n"]\n',i,equil(i,1),domain(i,1),domain(i,2));
            elseif equil(i,3) < 0
                fprintf(fid,'"eq%s" [fillcolor=green label="stable \n x* = %.4g \n x (%.4g,%.4g) \n"]\n',i,equil(i,1),domain(i,1),domain(i,2));
            else 
                fprintf(fid,'"eq%s" [label="marginally stable"]\n',i);
            end
        end
    end
    fprintf(fid,'}\n');
    
    % Plot x-node transitions ---------------------------------------------
    if numeq > 1
        for i = 1:numeq
            if domain(i,1) == -inf
                fprintf(fid,'"eq%s" -> "eq%s"[label="x<%.4g"] \n',i+1,i,domain(i,2));
            elseif domain(i,2) == inf
                fprintf(fid,'"eq%s" -> "eq%s"[label="x>%.4g"] \n',i-1,i,domain(i,1));
            elseif domain(i,1) == domain(i,2)
                fprintf(fid,'"eq%s" -> "eq%s"[label="x=%.4g"] \n',i+1,i,domain(i,1));
                fprintf(fid,'"eq%s" -> "eq%s"[label="x=%.4g"] \n',i-1,i,domain(i,1));
            else
                fprintf(fid,'"eq%s" -> "eq%s"[label="%.4g<x<%.4g"] \n',i-1,i,domain(i,1),domain(i,2));
                fprintf(fid,'"eq%s" -> "eq%s"[label="%.4g<x<%.4g"] \n',i+1,i,domain(i,1),domain(i,2));
            end
        end
    end       
    fprintf(fid,'}\n');
    fclose(fid);                
    
end

% Run GraphViz on DOT file, open resulting image file ---------------------
imageFile = [label '.png'];                         
system(sprintf('dot.exe -Tpng -Gsize="10,10" "%s" -o"%s"',dotFile,imageFile)); 
system(imageFile);     

%% Create executable function file
% funFile = [label '.m'];       
% fid = fopen(funFile,'w');      
% 
% fprintf(fid,'function dxdt = %s(x,u,p)\n\n',label);
% 
% if num_event == 0
%     fprintf(fid,'\ndxdt = %s;\n',strrep(func2str(func),'@(x,u,p)',''));
% else
%     for i = 1:num_event+1
%         if i == 1
%             fprintf(fid,'if %s\n\tdxdt = %s;\n',labp{i},pfun{1});
%         elseif i == num_event+1
%             fprintf(fid,'elseif %s\n\tdxdt = %s;\nelse\n\tdxdt = %s;\n\tdisp(''on bifurcation point!'')\nend\n',labp{i},pfun{1},strrep(func2str(func),'@(x,u,p)',''));
%         else
%             fprintf(fid,'elseif %s && %s\n\tdxdt = %s;\n',labp{i}(1:strfind(labp{2},'p')),labp{i}(strfind(labp{2},'p'):end),pfun{1});
%         end
%     end
% end
% fclose(fid);


