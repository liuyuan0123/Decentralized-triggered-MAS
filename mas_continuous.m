function dx = mas_continuous(t, x, t_f, t_k, x_t_k)
global A;   global B; global N;
% global alpha; global c1;
global A_adj; global K; global c;
% dx_i = A*x_i +B*u_i
% e_i = expm(A*(t-t_k_i))*x_i_t_k - x_i
% f = norm(e_i) - c_1*exp(-alpha*t)

u=zeros(N,1);
for i = 1:N
    for j = 1:N
        u(i,1) = u(i,1) + c*K*A_adj(i,j)*(expm(A*(t+t_f-t_k(i)))*x_t_k(2*i-1:2*i,1) - expm(A*(t+t_f-t_k(j)))*x_t_k(2*j-1:2*j,1));
    end
end

dx = zeros(2*N,1);
 %% dynamics
for i = 1:N
    dx(2*i-1:2*i) = A*x(2*i-1:2*i) + B*u(i);
end
end