function RPhaseCCPlot3(Rms_cc_RFactor,runID)

%*************************************************************************
% Used to plot figure of merits such as rms, R factor, phase CC, etc
% Please add argument on Rms_cc_RFactor size. Or set it as a class with 
% three attributes.
%*************************************************************************

%% all in one figure

% rms_cc = importdata('metrics_shrinkwrap_hole_v4_1.mat');
% [nIter,~] = size(Rms_cc_RFactor);
% x = linspace(1,nIter,nIter);
% 
% 
% plot(x,Rms_cc_RFactor(:,1),'b-*',x,Rms_cc_RFactor(:,2),'r-x',x,Rms_cc_RFactor(:,3),'g-x',...
%     x,Rms_cc_RFactor(:,4),'m-*',x,Rms_cc_RFactor(:,5),'k-x');
% %         plot(x,Rms_cc_RFactor(:,1),'b-*',x,Rms_cc_RFactor(:,3),'g-x')
% x1=xlabel('interation');
% set(x1,'FontSize', 15);
% set(x1,'FontWeight','bold');
% y1 = ylabel('rms, CC and R factor');
% %         y1 = ylabel('rms and R factor');
% set(y1,'FontSize',15);
% set(y1,'FontWeight','bold')
% fleg = legend('rms','CC','R factor','phase err ave','phase err weighted');
% %         fleg = legend('rms','R factor');
% set(fleg,'FontSize',15);
% set(fleg,'FontWeight','bold');
% figName = strcat('rms_cc_rfactor_runID_',int2str(runID),'.tif');
% saveas(gcf,figName);


%% plot two figures using subplot

[nIter,~] = size(Rms_cc_RFactor);
x = linspace(1,nIter,nIter);

subplot(2,1,1)
plot(x,Rms_cc_RFactor(:,1),'b',x,Rms_cc_RFactor(:,2),'r',x,Rms_cc_RFactor(:,3),'g','LineWidth',2)
    
%         plot(x,Rms_cc_RFactor(:,1),'b-*',x,Rms_cc_RFactor(:,3),'g-x')
% x1=xlabel('iteration');
% set(x1,'FontSize', 15);
% set(x1,'FontWeight','bold');
y1 = ylabel('rms, CC and R factor');
%         y1 = ylabel('rms and R factor');
set(y1,'FontSize',15);
set(y1,'FontWeight','bold')
fleg = legend('rms','CC','R factor');
%         fleg = legend('rms','R factor');
set(fleg,'FontSize',15);
set(fleg,'FontWeight','bold');

subplot(2,1,2)
plot(x,Rms_cc_RFactor(:,4)*180/pi,'m',x,Rms_cc_RFactor(:,5)*180/pi,'b','LineWidth',2);
xlabel('iteration','FontSize', 15,'FontWeight','bold');
ylabel('phase error','FontSize', 15,'FontWeight','bold');
fleg = legend('average phase error','weighted phase error');
%         fleg = legend('rms','R factor');
set(fleg,'FontSize',15);
set(fleg,'FontWeight','bold');


figName = strcat('rms_cc_rfactor_runID_',int2str(runID),'.tif');
saveas(gcf,figName);