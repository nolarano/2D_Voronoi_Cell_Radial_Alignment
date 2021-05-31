function [x,y] = flow(Cdat)
global criteria outside_cells

%Minimize energy of the voronoi tensellation 
%Using Cdat to indicate voronoi's centers positions

N = size(Cdat,1);
X = Cdat(:,1:2);
X0 = reshape(X,2*N,1)';
t = 1;

error = 10*criteria;
iter = 0;
ITER = 20;
inc = 1; %0.1

while error>criteria
    if iter>ITER 
        inc = inc/10;
        iter = 0;
    end
iter = iter +1;

Options = odeset('AbsTol',10^-6);
% [tX,X] = ode45(@solver,[t t+0.1],X0,Options);
[tX,X] = ode45(@solver,[0 inc],X0,Options);
error = max(abs(X(end,:)-X0))
X0=X(end,:);

fprintf('ode # iterations: %d \n',size(X,1));

X = X(end,:);
x = X(1:N)'; % size: N x 1
y = X(N+1:end)'; % size: N x 1

end
end