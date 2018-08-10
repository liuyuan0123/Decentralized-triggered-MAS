clc
clear all
close all
%% Description
% Decentralized event-triggered consensus for linear multi-agent systems under general directed graphs
% ���ۋ쓮�x�[�X�̍��Ӑ���
% ���K��l�@�����˂āD
% dynamics:
% 
% 
%% setting
dt = 0.1;
time = 10;
tspan = 8.0;

global A; global B; global K;
global alpha; global c; global c1;
global A_adj; global N;
A = [0, 1;
     -1, 0];
B = [1; 1];
n=size(A,2);
node = 6; N = node;
G = [1, 2; 2, 3; 1, 3; 6, 1; 4, 1; 1, 4; 4, 5; 5, 6; 5, 1]; % graph topology
K = [-2, -1]; % eig(A+BK)<0:Hurwitz�ɂȂ�悤�ɃQ�C����ݒ�

rand('state',0)
W = 2*rand(5,6);
Y = 2*rand(6,5);
%% initial condition
 % each initial condition
x1_0 = [0.4, 0.3]';
x2_0 = [0.5, 0.2]';
x3_0 = [0.6, 0.1]';
x4_0 = [0.7, 0.0]';
x5_0 = [0.8, -0.1]';
x6_0 = [0.4, -0.2]';
x0 = [x1_0', x2_0', x3_0', x4_0', x5_0', x6_0']'; % integration
x = x0; % 
xs = x0; % ode45�V�~�����[�V�����p�ϐ�(x_simulation)�F���ۂ̏�Ԃ̐��ڂ�\��
x_t_k = x0; % �e�G�[�W�F���g���g���C�茳�Ɏ����Ă����Ԃ̏��
%% Laplacian
[L, D, A_adj] = graph(node, G, 'directed'); % 'directed' or 'undirected'
%% original parameter
c = 1.1;
c1 = 0.6;
alpha = 0.4;
Delta = W*L*Y;
Pi = kron(eye(N-1),A)+kron(c*Delta,B*K);
%% Algorithm step1: calculation gain
% �Q�C���̌v�Z�C���̏����𖞂����悤�ɐ݌v����
% ���������̏����ɂ��Ă͂܂�P�͂�����ł����邽�ߍœK���⃍�o�X�g���Ȃǂ��l���ɉ����邱�Ƃ��\�H
% ����̃V�~�����[�V�����ł̓Q�C���͗^�����Ă��邽�߂����p����D
if 0
P = sdpvar(2,2,'symmetric');
LMI_P = P > 0; % P is positive-define
LMI_1 = A*P + P*A' - 2*B*B' < 0;
LMI = [LMI_P, LMI_1];
sol = solvesdp(LMI);
P = double(P);
% K = B'*inv(P);
end
%% Algorithm step2: setting parameter
eig_lap = real(eig(L));
eig_lap(eig_lap == min(eig_lap)) = [];
c_cand = min(eig_lap);
disp(sprintf('Should be c > %d', 1/c_cand))
[V, D] = eig(L');
% 0�����s�񐬕���T��
for j = 1:size(D);
    if find(D(j,j) <= 1e-5) == 1
        col_zero = j;
    end
end
r = V(:,col_zero)';
r = r/sum(r); % �ŗL�l0�ɂ������鍶�ŗL�x�N�g��
% delta = 
%% simulation
t = 0;              % ������
% ts: ode45�V�~�����[�V�����p�ϐ�(t_simulation)
t_k = zeros(N,1);   % �g���K�����O����
t_k_min = 0;        % �g���K�����O�Ԋu
ele_i_min = [1:N]'; % �g���K���ꂽ�G�[�W�F���g
u = zeros(N,1);     % ���͓��ꕨ

data = []; data_u = [];
for h = 1:50
    tic
    % ���͂̌v�Z
    if ~isempty(t_k_min) % �g���K�����ꍇ
        t = t + t_k_min;
        for i = 1:length(ele_i_min);
            u(ele_i_min(i),1) = 0; % ���Ă͂܂����Ƃ�������
            for j = 1:N
                u(ele_i_min(i),1) = u(ele_i_min(i),1) +...
                    c*K*A_adj(ele_i_min(i),j)*(expm(A*(t-t_k_min))...
                    *x_t_k(2*ele_i_min(i)-1:2*ele_i_min(i),1) - expm(A*(t-t_k_min))*x_t_k(2*j-1:2*j,1));
            end
        end
    else % �g���K����Ȃ������ꍇ
        t = t + tspan;
    end
    %% ode45
    [ts,xs] = ode45(@(ts,xs) mas(ts,xs,u), [0, tspan], x);
    
    xs = xs'; % �����₷���悤�Ɍ�����ς���
    e_x_t_k = kron(ones(1,length(ts)),x_t_k);
    e = zeros(n*N, length(ts));
    f = zeros(N, length(ts));
    t_k = zeros(N,1);
    % triggering check
    for i = 1:N
        for k = 1:length(ts)
%             e(2*i-1:2*i,k) = expm(A*(t(k)-t_k(i)))*e_x_t_k(2*i-1:2*i,k) - x(2*i-1:2*i,k);
            e(2*i-1:2*i,k) = expm(A*ts(k))*e_x_t_k(2*i-1:2*i,k) - xs(2*i-1:2*i,k); % tk=0�Ŗ���v�Z����Ă���
            f(i,k) = norm(e(2*i-1:2*i,k)) - c1*exp(-alpha*(ts(k)-t));
            if (f(i,k) > 0)
                t_k(i) = ts(k);
                break;
            end
        end
    end
    t_k_min = min(t_k(t_k~=0)); % �g���K���ꂽ�ŏ�����
    if ~isempty(t_k_min) % �g���K����Ȃ��ꍇ�̃G���[���
        ele_i_min = find(t_k == t_k_min); % ���ꂪ�s����G�[�W�F���g
        ele_k_min = find(ts == t_k_min); % ���̎����̗v�f�ԍ�(��X�g������k�Ŏw��)
        x = xs(:,ele_k_min); % ����̏����l�ɂȂ�
        x_t_k(2*ele_i_min-1:2*ele_i_min,1) = x(ele_i_min,1);
    else
        disp('ERROR!')
        ele_k_min = length(ts);
        x = xs(:,ele_k_min); % ����̏����l�ɂȂ�
        x_t_k(2*ele_i_min-1:2*ele_i_min,1) = x(ele_i_min,1);
%         break;
    end
    data = [data; t+ts(1:ele_k_min), xs(:,1:ele_k_min)', ts(1:ele_k_min)];
%     data_u = [data_u; u', t_k_min];
toc
end
%% figure
figure(1)
hold on; grid on; box on;
iro = ['b','r','g','c','y','k'];
for i = 1:N
    plot(data(:,1), data(:,2*i),iro(i))
end
xlim([0,8])
ylim([-2,2])

figure(2)
hold on; grid on; box on;
for i = 1:N
    plot(data(:,1), data(:,2*i+1))
end
