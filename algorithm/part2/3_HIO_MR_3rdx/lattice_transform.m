% lattice transform plot

n1 = 5;
n2 = 10;
n3 = 100;
N = 601;
x = linspace(-3,3,N);

y1 = abs((sin(pi*n1*x)./sin(pi*x)));
y2 = abs((sin(pi*n2*x)./sin(pi*x)));
y3 = abs((sin(pi*n3*x)./sin(pi*x)));
for i=1:N
    if rem(abs(x(i)),1)== 0
        y1(i)=n1;
        y2(i)=n2;
        y3(i)=n3;
    end
end

% plot(x,y1,'r',x,y2,'b',x,y3,'y');
subplot(3,1,1);
plot(x,y1,'b');
title('N=5');
subplot(3,1,2);
plot(x,y2,'b');
title('N=10');
ylabel('|sin(N\pix)|/|sin(\pix)|');
subplot(3,1,3);
plot(x,y3,'b');
title('N=100');
xlabel('x');
