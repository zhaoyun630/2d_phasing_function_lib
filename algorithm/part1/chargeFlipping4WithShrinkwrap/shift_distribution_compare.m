
x_I_abs1=x_I_abs +density_shift;
figure(2);
plt2=plot(x_I_abs,n_count,'b',x_I_abs1,n_count,'r');
xlabel('charge density','FontSize', 15,'FontWeight','bold');
ylabel('voxel counts','FontSize', 15,'FontWeight','bold');
title('charge density distribution within unit cell','FontSize', 15,'FontWeight','bold');
lg1_text1 = 'shifted density distribution';
lg1_text2 = 'original density distribution';
leg = legend(lg1_text1,lg1_text2);
set(leg,'FontSize',15,'FontWeight','bold')