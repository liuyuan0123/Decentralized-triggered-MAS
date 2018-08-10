function draw_calc_time(n)
% ���݂̌v�Z�̐i�s�󋵂�\������֐�
% ����...
% 0:����;  1:�i�s�󋵂�\��;  2:figure��close
%%
persistent h
% �A���ϐ��̒�`
time_f=evalin('base', 'time'); %% receive variable from base workspace
%% draw calculation time
if n==0
    scrsz = get(0,'ScreenSize');
    figure('MenuBar','None','Toolbar','None','Name','calculating...',...
        'Position',[scrsz(3)*74/100 scrsz(4)*88/100 scrsz(3)/4 scrsz(4)/12])
    hold on
    h=fill([0 0 0/time_f 100*0/time_f],[0 1 1 0],'green','EraseMode','xor','LineWidth',2);
    set(gca,'XTick',0:25:100); set(gca,'YTick',[0 1]);
    axis([0 100 0 1])
    drawnow

%% draw calculating...
elseif n==1
    t_f=evalin('base', 't'); %% receive variable from base workspace
%     if rem(t, 10)==0
        set(h,'XData',[0 0 100*t_f/time_f 100*t_f/time_f],'YData',[0 1 1 0])
        drawnow
%     end
%% close 
else
    close(1)
end