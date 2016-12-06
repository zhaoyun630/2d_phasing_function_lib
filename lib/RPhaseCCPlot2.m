function RPhaseCCPlot2(Rms_cc_RFactor,runID)

%*************************************************************************
% Used to plot figure of merits such as rms, R factor, phase CC, etc
% Please add argument on Rms_cc_RFactor size. Or set it as a class with 
% three attributes.
%*************************************************************************


% rms_cc = importdata('metrics_shrinkwrap_hole_v4_1.mat');
[nIter,nVar] = size(Rms_cc_RFactor);
x = linspace(1,nIter,nIter);

figure(1)


switch nVar
    case 1
        plot(x,Rms_cc_RFactor(:,1),'b-*')
        x1=xlabel('interation');
        set(x1,'FontSize', 15);
        set(x1,'FontWeight','bold');
        y1 = ylabel('R factor');
        set(y1,'FontSize',15);
        set(y1,'FontWeight','bold')
        fleg = legend('rms','CC','R factor');
        set(fleg,'FontSize',15);
        set(fleg,'FontWeight','bold');
        figName = strcat('rms_runID_',int2str(runID),'.tif');
        saveas(gcf,figName);
    case 2
        plot(x,Rms_cc_RFactor(:,1),'b-*',x,Rms_cc_RFactor(:,2),'r-x')
        x1=xlabel('interation');
        set(x1,'FontSize', 15);
        set(x1,'FontWeight','bold');
        y1 = ylabel('R factor and rms');
        set(y1,'FontSize',15);
        set(y1,'FontWeight','bold')
        fleg = legend('rms','CC','R factor');
        set(fleg,'FontSize',15);
        set(fleg,'FontWeight','bold');
        figName = strcat('rms_cc_runID_',int2str(runID),'.tif');
        saveas(gcf,figName);
    case 3
%         plot(x,Rms_cc_RFactor(:,1),'b-*',x,Rms_cc_RFactor(:,2),'r-x',x,Rms_cc_RFactor(:,3),'g-x')
        plot(x,Rms_cc_RFactor(:,1),'b-*',x,Rms_cc_RFactor(:,3),'g-x')
        x1=xlabel('interation');
        set(x1,'FontSize', 15);
        set(x1,'FontWeight','bold');
%         y1 = ylabel('rms, CC and R factor');
        y1 = ylabel('rms and R factor');
        set(y1,'FontSize',15);
        set(y1,'FontWeight','bold')
%         fleg = legend('rms','CC','R factor');
        fleg = legend('rms','R factor');
        set(fleg,'FontSize',15);
        set(fleg,'FontWeight','bold');
        figName = strcat('rms_cc_rfactor_runID_',int2str(runID),'.tif');
        saveas(gcf,figName);
    otherwise
        fprintf('Wrong option\n');
end
       