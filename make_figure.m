close all
%% set propety
set(0, 'defaultAxesFontSize', 20); % 軸のフォントサイズ
set(0, 'defaultTextFontSize', 20); % タイトル、注釈などのフォントサイズ
set(0, 'DefaultLineLineWidth', 2.5); % 波形の線の太さ
set(0, 'DefaultAxesLineWidth', 2);
%% figure for state
figure(1)
hold on; grid on; box on;
set(1,'Position',[10 60 1400 400])
iro = ['b','r','g','c','y','k'];
for i = 1:N
    plot(data_x(:,1), data_x(:,2*i),iro(i))
end
xlim([0,8])
ylim([-2,2])
legend('agent1','agent2','agent3','agent4','agent5','agent6')
title('x1')

figure(2)
set(2,'Position',[10 550 1400 400])
hold on; grid on; box on;
for i = 1:N
    plot(data_x(:,1), data_x(:,2*i+1))
end
xlim([0,time])
ylim([-2,1])
title('x2')

%% figure for interval
figure(3)
hold on; grid on; box on;
for i = 1:N
    plot(data_instant(:,N+1), i*data_instant(:,i),'o')
end
xlim([0,time])
ylim([0,6.5])
title('triggering instants')

%% figure for input
figure(4)
hold on; grid on; box on;
set(4,'Position',[200 300 1400 400])
for i = 1:N
    plot(data_u(:,1), data_u(:,i+1),iro(i))
end
xlim([0,time])
legend('agent1','agent2','agent3','agent4','agent5','agent6')
title('input')

%% figure for eror
figure(5)
set(5,'Position',[1920 31 1920 969]) 
for i = 1:N
    subplot(3,2,i)
    hold on; grid on; box on;
    plot(data_error(:,1), data_error(:,i+1),'r')
    plot(data_error(:,1), data_error(:,N+2),'k--')
    tit = ['agent' sprintf('%.0f', i)];
    title(tit)
    xlim([0,time])
end