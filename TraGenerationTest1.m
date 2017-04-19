%test program for "minimun snap trajectory generation and control for
%quadrotors"
%for three points
%����汾û���н��߽��ٶȵ�����
clear;
%;����
fixPoints = [0,0; 2,3; 3,4];
%��ʱ��
alpha = 60;

[m, n] = size(fixPoints);
totalLen = zeros(m,1);
for i = 2:m
    totalLen(i-1) = sqrt((fixPoints(i,1)-fixPoints(i-1,1))^2 + (fixPoints(i,2)-fixPoints(i-1,2))^2);
    totalLen(m) = totalLen(m) + totalLen(i-1); 
end
traTime = zeros(m-1,1);
for i = 1:(size(fixPoints(:,1))-1)
    traTime(i) = alpha / totalLen(m) * totalLen(i);
end

%0~1��·��
H = zeros(8,8);
for i=5:8
   for j=5:8
      H(i,j) = (i-1)*(i-2)*(i-3)*(i-4)*(j-1)*(j-2)*(j-3)*(j-4)/(i+j-7);
   end
end
%A = [];
%b = [];
%x����
Aeq = [1, 1, 1, 1, 1, 1, 1, 1;      %t1ʱ��λ������1
       %0, 1, 2, 3, 4, 5, 6, 7;      %t1ʱ���ٶ����ƣ�1�ε�
       %0, 0, 2, 6, 12, 20, 30, 42;  %t1ʱ�̼��ٶ����ƣ�2�ε�
       %0, 0, 0, 6, 24, 60, 120, 210;%t1ʱ�̼Ӽ��ٶ����ƣ�3�ε�
       1, 0, 0, 0, 0, 0, 0, 0;      %t0ʱ��λ������0
       0, 1, 0, 0, 0, 0, 0, 0;      %t0ʱ���ٶ�����0��1�ε�
       0, 0, 2, 0, 0, 0, 0, 0;      %t0ʱ�̼��ٶ�����0��2�ε�
       0, 0, 0, 6, 0, 0, 0, 0     %t0ʱ�̼Ӽ��ٶ�����0��3�ε�
      ];
 
if ((fixPoints(2,1)-fixPoints(1,1)) * (fixPoints(3,1)-fixPoints(2,1)) > 0)
    Vtempx = (fixPoints(3,1)-fixPoints(1,1))/alpha;
else
    Vtempx = 0;
end
      
Beq = [1;
       %Vtempx/(fixPoints(2,1)-fixPoints(1,1))*traTime(1);       %11/9����ע������·��λ�����ٶȼ��ٶȵȵ�ת��
       %0;
       %0;
       0;
       0;
       0;
       0];
sigma_ones(:,1) = quadprog(H, [], [], [], Aeq, Beq);

syms t;
temp = 0;
for i=1:8
    temp = temp + sigma_ones(i,1) * (t/traTime(1))^(i-1);
end
x(1,1) = temp * (fixPoints(2,1)-fixPoints(1,1)) + fixPoints(1,1);

%y����
Aeq = [1, 1, 1, 1, 1, 1, 1, 1;      %t1ʱ��λ������1
       %0, 1, 2, 3, 4, 5, 6, 7;      %t1ʱ���ٶ����ƣ�1�ε�
       %0, 0, 2, 6, 12, 20, 30, 42;  %t1ʱ�̼��ٶ����ƣ�2�ε�
       %0, 0, 0, 6, 24, 60, 120, 210;%t1ʱ�̼Ӽ��ٶ����ƣ�3�ε�
       1, 0, 0, 0, 0, 0, 0, 0;      %t0ʱ��λ������0
       0, 1, 0, 0, 0, 0, 0, 0;      %t0ʱ���ٶ�����0��1�ε�
       0, 0, 2, 0, 0, 0, 0, 0;      %t0ʱ�̼��ٶ�����0��2�ε�
       0, 0, 0, 6, 0, 0, 0, 0];     %t0ʱ�̼Ӽ��ٶ�����0��3�ε�

if ((fixPoints(2,2)-fixPoints(1,2)) * (fixPoints(3,2)-fixPoints(2,2)) > 0)
    Vtempy = (fixPoints(3,2)-fixPoints(1,2))/alpha;
else
    Vtempy = 0;
end
       
Beq = [1;
       %Vtempy/(fixPoints(2,2)-fixPoints(1,2))*traTime(1);%5/6;       %����ע������·��λ�����ٶȼ��ٶȵȵ�ת��
       %0;
       %0;
       0;
       0;
       0;
       0];
sigma_ones(:,2) = quadprog(H, [], [], [], Aeq, Beq);

syms t;
temp = 0;
for i=1:8
    temp = temp + sigma_ones(i,2) * (t/traTime(1))^(i-1);
end
y(1,1) = temp * (fixPoints(2,2)-fixPoints(1,2)) + fixPoints(1,2);

%1~2��·��
%x
Aeq = [1, 1, 1, 1, 1, 1, 1, 1;      %t2ʱ��λ������1
       0, 1, 2, 3, 4, 5, 6, 7;      %t2ʱ���ٶ����ƣ�1�ε�
       0, 0, 2, 6, 12, 20, 30, 42;  %t2ʱ�̼��ٶ����ƣ�2�ε�
       0, 0, 0, 6, 24, 60, 120, 210;%t2ʱ�̼Ӽ��ٶ����ƣ�3�ε�
       1, 0, 0, 0, 0, 0, 0, 0;      %t1ʱ��λ������
       0, 1, 0, 0, 0, 0, 0, 0;      %t1ʱ���ٶ����ƣ�1�ε�
       0, 0, 2, 0, 0, 0, 0, 0;      %t1ʱ�̼��ٶ����ƣ�2�ε�
       0, 0, 0, 6, 0, 0, 0, 0      %t1ʱ�̼Ӽ��ٶ����ƣ�3�ε�
       ];     
   
%��������˱߽�����
Vtempx = Aeq(2,:)*sigma_ones(:,1)*abs(fixPoints(2,1)-fixPoints(1,1))/traTime(1)*traTime(2)/abs(fixPoints(3,1)-fixPoints(2,1));
dVtempx = Aeq(3,:)*sigma_ones(:,1)*abs(fixPoints(2,1)-fixPoints(1,1))/(traTime(1)^2)*(traTime(2)^2)/abs(fixPoints(3,1)-fixPoints(2,1));
ddVtempx = Aeq(4,:)*sigma_ones(:,1)*abs(fixPoints(2,1)-fixPoints(1,1))/(traTime(1)^3)*(traTime(2)^3)/abs(fixPoints(3,1)-fixPoints(2,1));
if((fixPoints(3,1)-fixPoints(2,1))*(fixPoints(2,1)-fixPoints(1,1))<0)

    Vtempx = -Vtempx;

    dVtempx = -dVtempx;

    ddVtempx = -ddVtempx;
end
%if (Vtempx < 0)end if (dVtempx < 0)endif (ddVtempx < 0)end

Beq = [1;
       0;
       0;
       0;
       0;
       Vtempx;%11/12
       dVtempx;
       ddVtempx
       ];
sigma_ones(:,3) = quadprog(H, [], [], [], Aeq, Beq);

syms t;
temp = 0;
for i=1:8
    temp = temp + sigma_ones(i,3) * ((t-traTime(1))/traTime(2))^(i-1);
end
x(2,1) = temp * (fixPoints(3,1)-fixPoints(2,1)) + fixPoints(2,1);

%y
Aeq = [1, 1, 1, 1, 1, 1, 1, 1;      %t2ʱ��λ������1
       0, 1, 2, 3, 4, 5, 6, 7;      %t2ʱ���ٶ����ƣ�1�ε�
       0, 0, 2, 6, 12, 20, 30, 42;  %t2ʱ�̼��ٶ����ƣ�2�ε�
       0, 0, 0, 6, 24, 60, 120, 210;%t2ʱ�̼Ӽ��ٶ����ƣ�3�ε�
       1, 0, 0, 0, 0, 0, 0, 0;      %t1ʱ��λ������
       0, 1, 0, 0, 0, 0, 0, 0;      %t1ʱ���ٶ����ƣ�1�ε�
       0, 0, 2, 0, 0, 0, 0, 0;      %t1ʱ�̼��ٶ����ƣ�2�ε�
       0, 0, 0, 6, 0, 0, 0, 0      %t1ʱ�̼Ӽ��ٶ����ƣ�3�ε�
       ];     

Vtempy = Aeq(2,:)*sigma_ones(:,2)*abs(fixPoints(2,2)-fixPoints(1,2))/traTime(1)*traTime(2)/abs(fixPoints(3,2)-fixPoints(2,2));
dVtempy = Aeq(3,:)*sigma_ones(:,2)*abs(fixPoints(2,2)-fixPoints(1,2))/(traTime(1)^2)*(traTime(2)^2)/abs(fixPoints(3,2)-fixPoints(2,2));
ddVtempy = Aeq(4,:)*sigma_ones(:,2)*abs(fixPoints(2,2)-fixPoints(1,2))/(traTime(1)^3)*(traTime(2)^3)/abs(fixPoints(3,2)-fixPoints(2,2));
if ((fixPoints(3,2)-fixPoints(2,2))*(fixPoints(2,2)-fixPoints(1,2))<0)
    Vtempy = -Vtempy;

    dVtempy = -dVtempy;

    ddVtempy = -ddVtempy;
end   
Beq = [1;
       0;
       0;
       0;
       0;
       Vtempy;%10/9;
       dVtempy;
       ddVtempy
       ];
sigma_ones(:,4) = quadprog(H, [], [], [], Aeq, Beq);

syms t;
temp = 0;
for i=1:8
    temp = temp + sigma_ones(i,4) * ((t-traTime(1))/traTime(2))^(i-1);
end
y(2,1) = temp * (fixPoints(3,2)-fixPoints(2,2)) + fixPoints(2,2);

%��ʾ·�����ٶȵ�����
%{%}
%�ٶ�1
subplot(1,2,1)
fplot(diff(x(1,1)), [0,traTime(1)]);hold on;grid on;
fplot(diff(x(2,1)), [traTime(1),alpha]);
title('x�ٶ�');
subplot(1,2,2)%figure(2);
fplot(diff(y(1,1)), [0,traTime(1)]);hold on;grid on;
fplot(diff(y(2,1)), [traTime(1),alpha]);
title('y�ٶ�');

%���ٶ�2
figure(2);
subplot(1,2,1)
fplot(diff(diff(x(1,1))), [0,traTime(1)]);hold on;grid on;
fplot(diff(diff(x(2,1))), [traTime(1),alpha]);
title('x���ٶ�');
subplot(1,2,2)%figure(4);
fplot(diff(diff(y(1,1))), [0,traTime(1)]);hold on;grid on;
fplot(diff(diff(y(2,1))), [traTime(1),alpha]);
title('y���ٶ�');
%�Ӽ��ٶ�3
figure(3);
subplot(1,2,1)
fplot(diff(diff(diff(x(1,1)))), [0,traTime(1)]);hold on;grid on;
fplot(diff(diff(diff(x(2,1)))), [traTime(1),alpha]);
title('x�Ӽ��ٶ�');
subplot(1,2,2)%figure(4);
fplot(diff(diff(diff(y(1,1)))), [0,traTime(1)]);hold on;grid on;
fplot(diff(diff(diff(y(2,1)))), [traTime(1),alpha]);
title('y�Ӽ��ٶ�');
%�ӼӼ��ٶ�4

%�켣
figure(9);
t=0:0.1:traTime(1);
plot(eval(x(1,1)), eval(y(1,1)));hold on;grid on;
t=traTime(1):0.1:alpha;
plot(eval(x(2,1)), eval(y(2,1)));
title('�켣');
