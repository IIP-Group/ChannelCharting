% =========================================================================
% -- Functions to Plot Scenario and Channel Charts
% -------------------------------------------------------------------------
% -- (c) 2016-2021 Christoph Studer, Emre Gonultas, and Said Medjkouh
% -- e-mail: studer@ethz.ch, eg566@cornell.edu, sm2685@cornell.edu
% =========================================================================

function plotMapping(mappedX,par,plotChoice)

switch (plotChoice)
    
    case ('channel') % plot channel scenario
        
        h = figure();
        scatter3(par.bs.x,par.bs.y,par.bs.z,50,'+')
        hold on
        plot3([0,0],[0,0],[0,par.bs.z(1)],'k--')
        scatter3(par.ue.x,par.ue.y,par.ue.z,50,par.colorMap,'filled')
        
        % show boundary
        if (par.boundary)
            plot3(par.ue.x(par.boundIdx),par.ue.y(par.boundIdx),par.ue.z(par.boundIdx),'k-','LineWidth',1)
        end
        
        % color shape1
        if (par.Circ >0)
            colorCirc = [repmat([0 0],par.Circ,1) linspace(0.2,1,par.Circ)'];
            scatter3(par.ue.x(1:par.Circ),par.ue.y(1:par.Circ),par.ue.z(1:par.Circ),50,colorCirc,'filled');
        end
        if (par.shape>0)
            colorShape = [linspace(0.2,1,par.shape)' repmat([0 0],par.shape,1)];
            scatter3(par.ue.x(par.coordShape),par.ue.y(par.coordShape),par.ue.z(par.coordShape),300,colorShape,'filled');
        end
        
        % set axes
        hold off
        xlabel('x [m]', 'fontsize', 14);
        ylabel('y [m]', 'fontsize', 14);
        zlabel('z [m]', 'fontsize', 14);
        view(45,45)
        axis tight
        
        % save figure (in color and with a reasonable bounding box)
        print(h,'-loose','-dpng','./output/scenario.png')
        
    case('embedding2D') % plot channel chart
        
        %normalize before plotting
        mappedX(:,1)  =  (  mappedX(:,1) - mean(mappedX(:,1)) ) / std(mappedX(:,1));
        mappedX(:,2)  =  (  mappedX(:,2) - mean(mappedX(:,2)) ) / std(mappedX(:,2));
        
        %% plot results
        figure;
        
        % plot boundary
        if(par.boundary)
            plot(mappedX(par.boundIdx,1),mappedX(par.boundIdx,2),'k-','LineWidth',0.5)
        end
        
        hold on
        scatter(mappedX(:,1),mappedX(:,2),50,par.colorMap,'filled')
        
        % show the circle
        colorCirc = [repmat([0 0],par.Circ,1) linspace(0.2,1,par.Circ)'];
        scatter(mappedX(1:par.Circ,1),mappedX(1:par.Circ,2),70,colorCirc,'filled')
        plot(mappedX(1:par.Circ,1),mappedX(1:par.Circ,2),'b-','LineWidth',0.5);
        
        % show the shape
        colorShape = [linspace(0.2,1,par.shape)' repmat([0 0],par.shape,1)];
        scatter(mappedX(par.coordShape,1),mappedX(par.coordShape,2),200,colorShape,'filled')
        
        plot(mappedX(par.coordShape,1),mappedX(par.coordShape,2),'r-')
        
        hold off
        grid on
        xlabel('coordinate 1')
        ylabel('coordiante 2')
        str = sprintf('%s with U=%d and B=%d',par.nameParam, par.U,par.B);
        title(strrep(str,'_',' '))
        
        % save channel chart as png
        nameFig = sprintf('CC_%s_U%d_B%d',par.nameParam, par.U,par.B);
        print(['./output/',nameFig,'.png'],'-dpng')
        
end