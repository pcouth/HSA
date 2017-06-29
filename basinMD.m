%--- University of Washington, Department of Aeronautics & Astronautics ---
%---------- Advanced Dynamics, Validation & Control Research Lab ----------
%
% Basin of attraction identification tool for use with multi-dimensional,
% fixed-p hybrid stability automata (HSA). Takes HSA node data and 
% approximates basins of attraction for stable equilibria. Produces basin 
% data and visualizations.
%
% Currently works for 2D systems
%
% Basin approximation type: 
%       "Simulation"  - repeated simulation, fixed-step sampling
%       "Monte Carlo" - repeated simulation, random sampling
% 
% Author: Peter Uth
% Created: May 2017
%--------------------------------------------------------------------------
function basin = basinMD(node,type,xlimit,segs,tspan)

% Options
ztol = 0.4;         % zero tolerance for value at final simulation time
time_plots = 0;     % 1 / 0 - plot / do not plot time simulations

% Preallocate basin size, run main for-loop
basin(1,size(node,2)) = struct();
for I = 1:size(node,2)
    
    % Clear temporary variables, identify stable and unstable modes
    clearvars -except node type xlimit segs tspan I ztol basin time_plots
    stable_idx = (strcmp(node(I).p(node(I).fp_idx).stab,'stable'));
    unstable_idx = (strcmp(node(I).p(node(I).fp_idx).stab,'unstable'));
    [~,loc] = find(unstable_idx == 1);
    unstable = node(I).p(node(I).fp_idx).equil(loc,2:3);
       
    % Run for-loop for each stable equilibrium
    basin(I).x = [];
    r = 0;
    for i = 1:numel(stable_idx)
        if stable_idx(i) == 1
            r = r + 1;
            stable = node(I).p(node(I).fp_idx).equil(i,2:3);
            
            %% "Simulation" method - repeated simulation, fixed-step sampling
            if strcmp(type,'Simulation')
                for j = 1:segs + 1
                    size1 = (xlimit(1,2) - xlimit(1,1)) / segs;
                    x01 = xlimit(1,1) + (j - 1) * size1;
                    for h = 1:segs + 1
                        size2 = (xlimit(2,2) - xlimit(2,1)) / segs;
                        x02 = xlimit(2,1) + (h - 1) * size2;
                        func = strrep(func2str(node(I).func),sprintf('@(x,u,p%.0f)',node(I).p(node(I).fp_idx).pd),'@(t,x)');
                        func = strrep(func,sprintf('p%.0f',node(I).p(node(I).fp_idx).pd),sprintf('(%.4g)',node(I).fp));
                        func = str2func(func);
                        [t,x] = ode45(func,tspan,[x01 x02]);

                        figure(I)
                        title({sprintf('Basin for $D_{s%.0f%.0f}$',I,r) ; ''},'interpreter','latex','fontsize',14)
                        basin(I).Ds = sprintf('Ds%.0f%.0f',I,r);
                        hold on
                        leg2(1) = scatter(stable(1),stable(2),'k','filled');
                        leg2(2) = scatter(unstable(1),unstable(2),'k');
                        if abs(stable(1) - x(end,1)) < ztol && abs(stable(2) - x(end,2)) < ztol
                            leg2(3) = plot(x01,x02,'.g','markers',20);
                            basin(I).x = [basin(I).x;x01 x02];
%                             leg2(3) = plot(x(:,1),x(:,2),'g','linewidth',1); % alternative line style
                        else
                            leg2(4) = plot(x01,x02,'.r','markers',20);
%                             leg2(4) = plot(x(:,1),x(:,2),'r','linewidth',1); % alternative line style
                        end                        
                        
                        if time_plots == 1
                            figure(I+2*size(node,2))
                            title({sprintf('Basin for $D_{s%.0f%.0f}$',I,r) ; ''},'interpreter','latex','fontsize',14)
                            hold on
                            leg1(1) = plot3([stable(1) stable(1)],[tspan(1) tspan(2)],[stable(2) stable(2)],'--k','linewidth',1.5);
                            leg1(2) = plot3([unstable(1) unstable(1)],[tspan(1) tspan(2)],[unstable(2) unstable(2)],'--r','linewidth',1.5);
                            if abs(stable(1)-x(end,1)) < ztol && abs(stable(2)-x(end,2)) < ztol
                                leg1(3) = plot3(x(:,1),t,x(:,2),'g');
                                leg1(4) = plot3(x01,0,x02,'.g','markers',15);
                            else
                                leg1(5) = plot3(x(:,1),t,x(:,2),'r');
                                leg1(6) = plot3(x01,0,x02,'.r','markers',15);
                            end
                        end                        
                    end
                end
            end
            
            %% "Monte Carlo" method - repeated simulation, random sampling
            if strcmp(type,'MonteCarlo')
                samples(:,1) = (xlimit(1,2) - xlimit(1,1)) .* rand(segs,1) + xlimit(1,1);
                samples(:,2) = (xlimit(2,2) - xlimit(2,1)) .* rand(segs,1) + xlimit(2,1);
                for j = 1:segs
                    x01 = samples(j,1);
                    x02 = samples(j,2);
                    func = strrep(func2str(node(I).func),sprintf('@(x,u,p%.0f)',node(I).p(node(I).fp_idx).pd),'@(t,x)');
                    func = strrep(func,sprintf('p%.0f',node(I).p(node(I).fp_idx).pd),sprintf('(%.4g)',node(I).fp));
                    func = str2func(func);
                    [t,x] = ode45(func,tspan,[x01 x02]);
                    
                    figure(I)
                    title({sprintf('Basin for $D_{s%.0f%.0f}$',I,r) ; ''},'interpreter','latex','fontsize',14)
                    basin(I).Ds = sprintf('Ds%.0f%.0f',I,r);
                    hold on
                    leg2(1) = scatter(stable(1),stable(2),'k','filled');
                    leg2(2) = scatter(unstable(1),unstable(2),'k');
                    if abs(stable(1) - x(end,1)) < ztol && abs(stable(2) - x(end,2)) < ztol
                        leg2(3) = plot(x01,x02,'.g','markers',20);
%                         basin(:,:,I) = [basin;x01 x02];
                        basin(I).x = [basin(I).x;x01 x02];
                    else
                        leg2(4) = plot(x01,x02,'.r','markers',20);
                    end
                    
                    if time_plots == 1
                        figure(I+2*size(node,2))
                        title({sprintf('Basin for $D_{s%.0f%.0f}$',I,r) ; ''},'interpreter','latex','fontsize',14)
                        hold on
                        leg1(1) = plot3([stable(1) stable(1)],[tspan(1) tspan(2)],[stable(2) stable(2)],'--k','linewidth',1.5);
                        leg1(2) = plot3([unstable(1) unstable(1)],[tspan(1) tspan(2)],[unstable(2) unstable(2)],'--r','linewidth',1.5);
                        if abs(stable(1)-x(end,1)) < ztol && abs(stable(2)-x(end,2)) < ztol
                            leg1(3) = plot3(x(:,1),t,x(:,2),'g');
                            leg1(4) = plot3(x01,0,x02,'.g','markers',15);
                        else
                            leg1(5) = plot3(x(:,1),t,x(:,2),'r');
                            leg1(6) = plot3(x01,0,x02,'.r','markers',15);
                        end                
                    end                   
                end
            end
            
        else
            continue
        end
    end
    
    % Create basin boundary from results
    D = boundary(basin(I).x(:,1),basin(I).x(:,2));
    basin(I).xb = [basin(I).x(D,1) basin(I).x(D,2)];
    basin(I).xlimit = xlimit;
    basin(I).eq = stable;
    
    % Plot sample results
    figure(I)
    xlim(xlimit(1,:))
    ylim(xlimit(2,:))
    xlabel('$x_1$','interpreter','latex','fontsize',18)
    ylabel('$x_2$','interpreter','latex','fontsize',18)
    grid on
    box on
    legend2 = legend(leg2,'Stable Equilibrium','Unstable Equilibrium','Stable IC','Unstable IC','orientation','horizontal','interpreter','latex');
    set(legend2,'box','off')
    
    % Plot filled basin of attraction
    figure(I+size(node,2))
    fill(basin(I).xb(:,1),basin(I).xb(:,2),'g')
    hold on
    plot(basin(I).xb(:,1),basin(I).xb(:,2),'k')
    scatter(stable(1),stable(2),'k','filled')
    title(sprintf('$D_{s%.0f%.0f}$ Basin of Attraction',I,r),'interpreter','latex','fontsize',14)
    xlim(xlimit(1,:))
    ylim(xlimit(2,:))
    xlabel('$x_1$','interpreter','latex','fontsize',18)
    ylabel('$x_2$','interpreter','latex','fontsize',18)
    grid on
    box on
    
    % Plot time trajectories
    if time_plots == 1
        figure(I+2*size(node,2))
        hold on
        xlim(xlimit(1,:))
        ylim(tspan)
        zlim(xlimit(2,:))
        xlabel('$x_1$','interpreter','latex','fontsize',18)
        ylabel('time','interpreter','latex','fontsize',18)
        zlabel('$x_2$','interpreter','latex','fontsize',18)
        grid on
        box on
        legend1 = legend(leg1,'Stable Equilibrium','Unstable Equilibrium','Stable Trajectory','Stable IC','Unstable Trajectory','Unstable IC','interpreter','latex','orientation','horizontal');
        set(legend1,'box','off')
    end    
end

% Save basin data for external use
save(sprintf('%s_basins.mat',node(1).name),'basin','basin')
