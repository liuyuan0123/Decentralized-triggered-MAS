clc
clear all
close all
%% Description
% Decentralized event-triggered consensus for linear multi-agent systems under general directed graphs
% 事象駆動ベースの合意制御
% 練習や考察を兼ねて．
% 
%% setting
time = 8;
tspan = 2.0;

global A; global B; global K;
global alpha; global c; global c1;
global A_adj; global N;
A = [0, 1;
     -1, 0];
B = [1; 1];
n=size(A,2);
node = 6;
N = node;
G = [1, 2; 2, 3; 1, 3; 6, 1; 4, 1; 1, 4; 4, 5; 5, 6; 5, 1]; % graph topology
K = [-2, -1]; % eig(A+BK)<0:Hurwitzになるようにゲインを設定

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
x = x0; % 離散的な状態の遷移，トリガされたときの全状態を表し，odeの初期状態になる
xs = x0; % ode45シミュレーション用変数(x_simulation)：実際の状態の推移を表す
x_t_k = x0; % 各エージェントが使う，手元に持っている状態の情報
%% Laplacian
[L, D, A_adj] = graph(node, G, 'directed'); % 'directed' or 'undirected'
%% original parameter
c = 1.1;
c1 = 0.6;
alpha = 0.4;
Delta = W*L*Y;
Pi = kron(eye(N-1),A)+kron(c*Delta,B*K);
%% Algorithm step1: calculation gain
% ゲインの計算，この条件を満たすように設計する
% ただしこの条件にあてはまるPはいくらでもあるため最適性やロバスト性などを考慮に加えることが可能？
% 今回のシミュレーションではゲインは与えられているためそれを用いる．
if 0
P = sdpvar(2,2,'symmetric');
LMI_P = P > 0; % P is positive-define
LMI_1 = A*P + P*A' - 2*B*B' < 0;
LMI = [LMI_P, LMI_1];
sol = solvesdp(LMI);
P = double(P);
K = B'*inv(P);
end
%% Algorithm step2: setting parameter
eig_lap = real(eig(L));
eig_lap(eig_lap == min(eig_lap)) = [];
c_cand = min(eig_lap);
disp(sprintf('Should be c > %d', 1/c_cand))
[V, D] = eig(L');
% 0をもつ行列成分を探す
for j = 1:size(D);
    if find(D(j,j) <= 1e-5) == 1
        col_zero = j;
    end
end
r = V(:,col_zero)';
r = r/sum(r); % 固有値0にたいする左固有ベクトル
%% simulation
t = 0;              % 実時間
% ts: ode45シミュレーション用変数(t_simulation)
t_k_i = zeros(N,1); % 各エージェントのトリガリング時刻の計算用
t_k = zeros(N,1);   % 全エージェントの最も最近行われたトリガリング時刻
t_k_min = 0;        % トリガリング間隔
ele_i_min = [1:N]'; % トリガされたエージェント
% u = zeros(N,1);     % 入力入れ物

data_x = [];
data_u = [];
data_instant = [];
data_error = [];
h = 0;
draw_calc_time(0)
while (t < time)
    h = h + 1;
% for h = 1:100
    tic
    if ~isempty(t_k_min) % トリガされる場合
        t = t + t_k_min;
    else % トリガされなかった場合
        t = t + tspan;
    end
%     t
%     t_k
    %% ode45
    [ts,xs] = ode45(@(ts,xs) mas_continuous(ts,xs,t,t_k,x_t_k), [0, tspan], x);
    xs = xs'; % 扱いやすいように向きを変える
    
    % triggering check
    X_t_k = kron(ones(1,length(ts)),x_t_k); % 差分を取るためにサイズを拡大
    e = zeros(n*N, length(ts)); % 状態の誤差の入れ物
    e_norm = zeros(N, length(ts)); % 誤差のノルム plot用
    thre = zeros(length(ts), 1);
    f = zeros(N, length(ts)); % 判定関数
%     t_k = zeros(N,1); % トリガ時刻 初期化
    t_k_i = -1*ones(N,1);
    u = zeros(length(ts),N);
    for i = 1:N
        for j = 1:N
            for k = 1:length(ts)
                u(k,i) = u(k,i) + c*K*A_adj(i,j)*(expm(A*(ts(k)+t-t_k(i)))*x_t_k(2*i-1:2*i,1)...
                    - expm(A*(ts(k)+t-t_k(j)))*x_t_k(2*j-1:2*j,1));
            end
        end
        for k = 1:length(ts)
%             e(2*i-1:2*i,k) = expm(A*(t(k)-t_k(i)))*e_x_t_k(2*i-1:2*i,k) - x(2*i-1:2*i,k);
            e(2*i-1:2*i,k) = expm(A*(ts(k)+t-t_k(i)))*X_t_k(2*i-1:2*i,k) - xs(2*i-1:2*i,k);
%             f(i,k) = norm(e(2*i-1:2*i,k)) -
%             c1*exp(-alpha*(ts(k)+t-t_k(i))); % トリガリング時刻に依存するver
            e_norm(i,k) = norm(e(2*i-1:2*i,k)); % 誤差ノルム，変更可
%             thre(k) = c1*exp(-alpha*(ts(k)+t)); % しきい値，exponentail ver.
            thre(k) = 2; % しきい値，任意ver.
            f(i,k) = e_norm(i,k) - thre(k); % 収束するver
            if (f(i,k) > 0)
                t_k_i(i) = ts(k);
                break;
            end
        end
    end
    t_k_min = min(t_k_i(t_k_i~=-1)); % トリガされた最小時間
    if ~isempty(t_k_min) % トリガされない場合のエラー回避
        ele_i_min = find(t_k_i == t_k_min); % それが行われるエージェント
        ele_k_min = find(ts == t_k_min); % その時刻の要素番号(後々使うためkで指定)
        x = xs(:,ele_k_min); % 毎回の初期値になる
        t_k_only_i = zeros(N,1); % トリガリングが複数台で行われる場合
        if length(ele_i_min) == 1
            x_t_k(2*ele_i_min-1:2*ele_i_min,1) = x(2*ele_i_min-1:2*ele_i_min,1);
            t_k_only_i(ele_i_min) = t_k_i(ele_i_min);
        else
            for m = 1:length(ele_i_min)
                x_t_k(2*ele_i_min(m)-1:2*ele_i_min(m),1) = x(2*ele_i_min-1:2*ele_i_min,1);
                t_k_only_i(ele_i_min(m)) = t_k_i(ele_i_min(m));
            end
        end
    else % tspan区間でトリガが生じなかった場合
        disp('Warning! shortage of simulation time')
        ele_i_min = [];
        ele_k_min = length(ts);
        x = xs(:, ele_k_min); % 毎回の初期値になる
        x_t_k = x_t_k_m1;
    end
%         t_k_min = ts(ele_k_min+1); % zeno behavior 回避ver
    t_k(ele_i_min) = t + t_k_min;
    % tspanのうちに状態が更新されなかった場合，引き続き前回の状態を使い続ける
    % そのために状態を保存しておく必要がある
    x_t_k_m1 = x_t_k;
    draw_calc_time(1)
    
    data_x = [data_x; t+ts(1:ele_k_min), xs(:,1:ele_k_min)', ts(1:ele_k_min)];
    data_instant(h,1:N) = -0.1;
    data_instant(h,ele_i_min) = 1;
    data_instant(h,N+1) = t;
    data_u = [data_u; t+ts(1:ele_k_min), u(1:ele_k_min,:)];
    data_error = [data_error;t+ts(1:ele_k_min), e_norm(:,1:ele_k_min)', thre(1:ele_k_min)];
% toc
end
draw_calc_time(1)
%% figure
make_figure
%% display results
%
%
%
%
%