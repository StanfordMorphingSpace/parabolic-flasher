% plots in 3d the origami represented by nodes, edges

function fig = plot3dNodesEdges(nodes, edges, angles, varargin)
    if size(varargin) == 0
        fig     = figure('Color', [1 1 1]);
    elseif size(varargin) == 1
        fig     = figure(varargin{1});
    end
    
    nEdges = size(edges,1);
        
    hold on
    for i = 1:nEdges        
        if isnan(angles(i))
            color = [0 0 0];        
        elseif angles(i) > 0 % valley
            color = [0 0 1]; % + angles(i)/pi()*[0 1 1];
        else                 % mountain
            color = [1 0 0]; % - angles(i)/pi()*[1 1 0];
        end
        
        %color = 'k';
        
        p1 = nodes(:,edges(i,1));
        p2 = nodes(:,edges(i,2));
        lh = plot3([p1(1) p2(1)],[p1(2) p2(2)], [p1(3) p2(3)], 'Color', color, 'LineWidth', 2);
        if ~isnan(angles(i))
            lh.Color = [lh.Color 1-abs(angles(i)/pi)]; % set alpha value for origami simulator
        end
    end
    
    % plot nodes
   % plot3(nodes(1,:), nodes(2,:), nodes(3,:), 'LineStyle', 'none', 'Marker', '.', 'MarkerEdgeColor', 'k', 'MarkerSize', 3)

    % % plot node index numbers
    % for i=1:size(nodes,2) 
    %     text(nodes(1,i)+0.01/1000,nodes(2,i)-0.01/1000,nodes(3,i),num2str(i),'Fontsize',8);
    % end
       
    hold off
    
    view(38,30)
    axis equal
    axis vis3d
    axis off
    ax = gca;               % get the current axis
    ax.Clipping = 'off';    % turn clipping off
    
end