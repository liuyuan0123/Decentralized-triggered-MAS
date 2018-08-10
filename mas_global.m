function dx = mas_global(t, x, x_t_k, interval, t_k)
global A;   global B; global K;
global alpha; global c; global c1;
global A_adj
x_t_k = x_t_k';
N=6;
% dx_i = A*x_i +B*u_i
% e_i = expm(A*(t-t_k_i))*x_i_t_k - x_i
% f = norm(e_i) - c_1*exp(-alpha*t)
% dx(1:2) = A*x(1:2) - B1*u1 + D1*omega(1);

if (interval) % 伝送が行われるとき
    l = []; % 毎回初期化からのスタート
    l = find(t_k);
%     if ~(length(l) == 1)
%         l
%     end
    u = zeros(N,1);
    for m = 1:length(l)
        for j = 1:N
            u(l(m)) = u(l(m)) + c*K*A_adj(l(m),j)*(expm(A*(t-t_k(l(m))))*x(2*l(m)-1:2*l(m)) - expm(A*(t-t_k(j)))*x(2*j-1:2*j));
        end
    end
end

% dynamics
for l = 1:N
    t_k(l) = 0;
    dx(2*l-1:2*l) = A*x(2*l-1:2*l) + B*u(l);
    e(2*l-1:2*l)= expm(A*(t-t_k(l)))*x_t_k(2*l-1:2*l) - x(2*l-1:2*l);
    f(l) = norm(e(2*l-1:2*l))-c1*exp(-alpha*(t-t_k(l)));
    if (f(l) > 0)
        t_k(l) = t;
    end
end



end