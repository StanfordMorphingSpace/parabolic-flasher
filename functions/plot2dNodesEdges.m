% plots in 2d the origami represented by nodes, edges

function fig = plot2dNodesEdges(nodes, edges, angles, varargin)
    if size(varargin) == 0
        fig     = figure('Color', [1 1 1]);
    elseif size(varargin) == 1
        fig     = figure(varargin{1});
    end
    
    nEdges  = size(edges, 1);    
    hold on
    
    for i = 1:nEdges                
        p1 = nodes(:,edges(i,1));
        p2 = nodes(:,edges(i,2));
        
        
        if isnan(angles(i))
            color = [0 0 0];        
        elseif angles(i) > 0
            color = [0 0 1]; % + angles(i)/pi()*[0 1 1];
        else
            color = [1 0 0]; % - angles(i)/pi()*[1 1 0];
        end
        
        %color = 'k';
        
        % line([p1(1) p2(1)],[p1(2) p2(2)], 'Color', color, 'LineWidth', 1)
        lh = plot([p1(1) p2(1)],[p1(2) p2(2)], 'Color', color, 'LineWidth', 1.5);
        if ~isnan(angles(i))
            lh.Color = [lh.Color 1-abs(angles(i)/pi)];
        end
    end
    
    plot(nodes(1,:), nodes(2,:), 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 1, 'MarkerEdgeColor', 'k')

    % % plot node index numbers
    % for i=1:size(nodes,2) 
    %     text(nodes(1,i)+0.01,nodes(2,i)-0.01,num2str(i),'Fontsize',8);
    % end
    
    hold off    
    axis equal
    axis tight
    axis off
    
    ax = gca;
    ax.Clipping = 'off';
end