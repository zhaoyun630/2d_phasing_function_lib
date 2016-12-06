% show results from each iteration from hio.
% YZ. Mar 18, 2014.
% Last editted by YZ May 7, 2014


% for w=1:2
%  for j = 1:40                      
%    b1=['w_' int2str(w) 'j_' int2str(j)];  
%    load(b1);
%    imshow(b1_1,[]);
%    pause(0.3);
%  end
% end
rms_cc_temp = importdata('rms_cc_0.mat');
rms_cc_res = rms_cc_temp(160,:);

for i = 1:26
    n_info = ['rms_cc_' int2str(i) '.mat'];
    rms_cc_temp = importdata(n_info);
    rms_cc_res = [rms_cc_res; rms_cc_temp(160,:)];
end

save('rms_cc_res.mat','rms_cc_res');

figure(2)
rms_cc_res = importdata('rms_cc_res.mat');
L = size(rms_cc_res);
x = linspace(0,L(1)-1,L(1))/80;
plot(x,rms_cc_res(:,1),'b',x,rms_cc_res(:,2),'r')
t1 = title('rms and CC value at convergence'); 
x1=xlabel('k(1/A)');
set(x1,'FontSize', 15);
set(x1,'FontWeight','bold');
y1 = ylabel('CC and rms');
set(y1,'FontSize',15);
set(y1,'FontWeight','bold')
fleg = legend('rms','CC');
set(fleg,'FontSize',15);
set(fleg,'FontWeight','bold');
saveas(gcf,'rms_cc_res.png');
