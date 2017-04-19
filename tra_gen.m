function [ traLines ] = tra_gen( fixPoints, alpha )
%input: fixPoints - the points that along the road; alpha - total time
%output:the road function
%this version limits the v, dv, ddv, dv and ddv are set to 0 during the
%road points

%m - the number of points
[m, n] = size(fixPoints);
if(n <= 1) 
    return;
end
%n col are trajectory, 2n col are v, 3n col are dv, 4n col are ddv
traLines = zeros(m-1, n*4);
traLines = sym(traLines);
%totalLen - 1 col lenth of each road and the last one is the total lenth of
%all roads
totalLen = zeros(m,2);
for i = 2:m
    totalLen(i-1) = sqrt((fixPoints(i,1)-fixPoints(i-1,1))^2 + (fixPoints(i,2)-fixPoints(i-1,2))^2);
    totalLen(m) = totalLen(m) + totalLen(i-1); 
end
%traTime - time for each road
traTime = zeros(m-1,2);
for i = 1:(m-1)
    traTime(i, 1) = alpha / totalLen(m) * totalLen(i);
    traTime(i+1, 2) = traTime(i, 2) + traTime(i, 1);
end

%H - hession matric
H = zeros(8, 8);
for i=5:8
   for j=5:8
      H(i,j) = (i-1)*(i-2)*(i-3)*(i-4)*(j-1)*(j-2)*(j-3)*(j-4)/(i+j-7);
   end
end

%Aeq - for quadprog function, Aeq*sigma=beq
Aeq = [1, 1, 1, 1, 1, 1, 1, 1;      %t1时刻位置限制
       0, 1, 2, 3, 4, 5, 6, 7;      %t1时刻速度限制，1次导
       0, 0, 2, 6, 12, 20, 30, 42;  %t1时刻加速度限制，2次导
       0, 0, 0, 6, 24, 60, 120, 210;%t1时刻加加速度限制，3次导
       1, 0, 0, 0, 0, 0, 0, 0;      %t0时刻位置限制
       0, 1, 0, 0, 0, 0, 0, 0;      %t0时刻速度限制，1次导
       0, 0, 2, 0, 0, 0, 0, 0;      %t0时刻加速度限制，2次导
       0, 0, 0, 6, 0, 0, 0, 0     %t0时刻加加速度限制，3次导
      ];
  

sigma_ones = zeros(8, (m-1)*n);
%generate trajectory for x,y,z or x,y
for j = 1:n
    Vin = 0;
    Vout = 0;
    %for 1 to m-1 roads
    for i = 1:(m-1)
        %cal Vin
        if(i==1)
           Vin = 0;
        else if(Vout ~= 0)
                Vin = Aeq(2,:)*sigma_ones(:,j+(i-2)*n)*abs(fixPoints(i,j)-fixPoints(i-1,j))/traTime(i-1)*traTime(i)/abs(fixPoints(i+1,j)-fixPoints(i,j));
            else
                Vin = 0;
            end
        end
        %cal Vout
        if(i==m-1)  %last road
            Vout = 0;
        else if ((fixPoints(i+1,j)-fixPoints(i,j)) * (fixPoints(i+2,j)-fixPoints(i+1,j)) > 0)
                Vout = (fixPoints(i+2,j)-fixPoints(i,j))/alpha/(fixPoints(i+1,j)-fixPoints(i,j))*traTime(i);
            else
                Vout = 0;
            end
        end
        
        Beq = [1;
               Vout;       %11/9这里注意两段路单位化后速度加速度等的转换
               0;
               0;
               0;
               Vin;
               0;
               0];

        sigma_ones(:,j+n*(i-1)) = quadprog(H, [], [], [], Aeq, Beq);
        syms t;
        temp = 0;
        for k=1:8
            temp = temp + sigma_ones(k,j+n*(i-1)) * ((t-traTime(i, 2))/traTime(i,1))^(k-1);
        end
        traLines(i,j) = temp * (fixPoints(i+1,j)-fixPoints(i,j)) + fixPoints(i,j);
    end
end

%plot,i次导
for i=1:4
    figure(i);    
    for j=1:n  %x,y,z
        for k=1:m-1  %k段
            subplot(1,n,j);
            traLines(k, j+n*i) = diff(traLines(k,j+n*i-n));
            fplot(traLines(k,j+n*i), [traTime(k,2),traTime(k+1,2)]);hold on;grid on;
            title(['此处为',num2str(i),'阶导']);
        end
    end
end
figure(i+1);
if(n==2)
    for k=1:m-1
        t=traTime(k,2):0.1:traTime(k+1,2);
        plot(eval(traLines(k,1)), eval(traLines(k,2)));hold on;grid on;
    end
else if(n==3)
        for k=1:m-1
            t=traTime(k,2):0.1:traTime(k+1,2);
            plot3(eval(traLines(k,1)), eval(traLines(k,2)), eval(traLines(k,3)));hold on;grid on;
        end
    end
end
end


