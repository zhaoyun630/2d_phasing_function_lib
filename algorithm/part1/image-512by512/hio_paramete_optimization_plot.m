% This piece of code is specifically used for plots of
% hio_parameter_optimization.

% First created by Yun Zhao on Aug 28, 2014
% Last editted by Yun Zhao on Aug 28, 2014

prompt = 'please give the output file from hio_parameter_optimization: \n The file is probably phase_error.mat \n';
file_name = input(prompt,'s');
phase_summary = importdata(file_name);
[length,width] = size(phase_summary);

prompt = 'please give the matrix with true phase: \n The name is diff_phase_true.mat \n';
file_phase = input(prompt,'s');
true_phase = importdata(file_phase);
[a,b,c] = size(true_phase);
phi_100  = true_phase((a+3)/2,(b+1)/2,(c+1)/2)*ones(length,1);
phi_010  = true_phase((a+1)/2,(b+3)/2,(c+1)/2)*ones(length,1);
phi_m100 = true_phase((a-1)/2,(b+1)/2,(c+1)/2)*ones(length,1);
phi_0m10 = true_phase((a+1)/2,(b-1)/2,(c+1)/2)*ones(length,1);
%plot(phase_summary(:,1),phase_summary(:,2),phase_summary(:,1),phase_summary(:,3),phase_summary(:,1),phase_summary(:,4),phase_summary(:,1),phase_summary(:,5),phase_summary(:,1),phase_summary(:,6),phase_summary(:,1),phase_summary(:,6));
plot(phase_summary(:,1),phase_summary(:,2),'m',phase_summary(:,1),phase_summary(:,3),'r');
hold on
plot(phase_summary(:,1),phase_summary(:,4),'g',phase_summary(:,1),phase_summary(:,5),'b');
hold on
plot(phase_summary(:,1),phase_summary(:,6),'c');
hold on
plot(phase_summary(:,1),phase_summary(:,7),'Color',[0.1,0.4,0.6]);
hold on
plot(phase_summary(:,1),phi_100,'m-o',phase_summary(:,1),phi_010,'r-o');
hold on
plot(phase_summary(:,1),phi_m100,'g-o',phase_summary(:,1),phi_0m10,'b-o');

x1=xlabel('feedback parameter    \beta','FontSize', 15,'FontWeight','bold');
% x1=xlabel('feedback parameter \beta');
% set(x1,'FontSize', 15,'FontWeight','bold');
y1 = ylabel('phase_{hio}, phase_{true}, rms and cc','FontSize', 15,'FontWeight','bold');

t1 = title('HIO parameter optimization','FontSize', 15,'FontWeight','bold'); 

fleg = legend;
set(fleg,'Location','NorthEastOutside');
saveas(gcf,'hio_parameter_optimization.tif');
