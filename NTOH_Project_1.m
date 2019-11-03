clc;
clear;
Nx = 100;
Ny = 100;
hx = 0.001;
hy = 0.001;
X = zeros(Nx,1);
Y = zeros(Ny,1);
phiold = zeros(Nx,Ny);
phinew = zeros(Nx,Ny);

for i = 1:Nx
    X(i,1) = i*hx;
end
for j = 1:Ny
    Y(j,1) = j*hy;
end

for i = 1:Nx
    for j = 1:Ny
        phiold(i,j) = sin(Nx*hx*i) + sin(Ny*hy*j);
    end
end

[XX,YY]= meshgrid(X,Y);
%surfnorm(XX,YY,grid);
mesh(X,Y,phiold');




