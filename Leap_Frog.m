close all;
clc;
clear;

%Simulation grid parameters
Nx = 200;
Ny = 10;
L = 200;
d =  10;
hx = L/Nx;
hy = d/Ny;
Nt = 20000;
ht = 0.02;
X = zeros(Nx,1);
Y = zeros(Ny,1);
for i = 1:Nx
    X(i,1) = i*hx;
end
for j = 1:Ny
    Y(j,1) = j*hy;
end
itanew = zeros(Nx,Ny);
itaold = zeros(Nx,Ny);
uold = zeros(Nx,Ny);
vold = zeros(Nx,Ny);
unew = zeros(Nx,Ny);
vnew = zeros(Nx,Ny);
uolder = zeros(Nx,Ny);
volder = zeros(Nx,Ny);
itaolder = zeros(Nx,Ny);
g= 9.81;

%LAWT Parameters
A = 0.2;
Tx = 20;

depth = 2.5;
knum = (2*pi)/98.63;
omx = 2*pi/Tx;

for k = 1:Nt
    z = A*sin(-omx*ht*k);
    uold(1,:) = (pi*2*z)/(knum*Tx)/depth;
    itaold(1,:) = z;
    uold(Nx,:) = 0;
    vold(Nx,:) = vold(Nx-1,:);
    itaold(Nx,:) = itaold(Nx-1,:);
    vold(1,:) = 0;
    vold(:,Ny) = 0;
    uold(:,Ny) = uold(:,Ny-1);
    uold(:,1) = uold(:,2);
    vold(:,1) = 0;
    itaold(:,1) = itaold(:,2);
    itaold(:,Ny) = itaold(:,Ny-1);
    if k ==1
        for i = 2:Nx-1
            for j = 2:Ny-1
                itanew(i,j) = itaold(i,j) - ht*((uold(i+1,j) - uold(i-1,j))/(2*hx)) * (depth) - ht*(vold(i,j+1) - vold(i,j-1))/(2*hy) * (depth);
                unew(i,j) = -((g * ht)/(2*hx)) * (itaold(i+1,j) - itaold(i-1,j)) + uold(i,j);
                vnew(i,j) = -((g * ht)/(2*hy)) * (itaold(i,j+1) - itaold(i,j-1)) + vold(i,j);
                %itanew(i,j) = itaold(i,j) - ht*((itaold(i+1,j) - itaold(i-1,j))/(2*hx)) * uold(i,j) - ht*((uold(i+1,j) - uold(i-1,j))/(2*hx)) * (depth) - ht*(vold(i,j+1) - vold(i,j-1))/(2*hy) * (depth) - ht*((itaold(i,j+1) - itaold(i,j-1))/(2*hy))*vold(i,j);
                %unew(i,j) = -((g * ht)/(2*hx)) * (itaold(i+1,j) - itaold(i-1,j)) - (vold(i,j)*ht/(2*hy))*(uold(i,j+1) - uold(i,j-1)) - (uold(i,j)*ht/(2*hx))*(uold(i+1,j) - uold(i-1,j)) + uold(i,j);
                %vnew(i,j) = -((g * ht)/(2*hy)) * (itaold(i,j+1) - itaold(i,j-1)) - (vold(i,j)*ht/(2*hy))*(vold(i,j+1) - vold(i,j-1)) - (uold(i,j)*ht/(2*hx))*(vold(i+1,j) - vold(i-1,j)) + vold(i,j);
            
            end
        end
    else
        for i = 2:Nx-1
            for j = 2:Ny-1
                itanew(i,j) = itaolder(i,j) - 2*ht*((uold(i+1,j) - uold(i-1,j))/(2*hx)) * (depth) - 2*ht*(vold(i,j+1) - vold(i,j-1))/(2*hy) * (depth);
                unew(i,j) = -((g * 2*ht)/(2*hx)) * (itaold(i+1,j) - itaold(i-1,j)) + uolder(i,j);
                vnew(i,j) = -((g * 2*ht)/(2*hy)) * (itaold(i,j+1) - itaold(i,j-1)) + volder(i,j);
                %itanew(i,j) = itaolder(i,j) - 2*ht*((itaold(i+1,j) - itaold(i-1,j))/(2*hx)) * uold(i,j) - 2*ht*((uold(i+1,j) - uold(i-1,j))/(2*hx)) * (depth) - 2*ht*(vold(i,j+1) - vold(i,j-1))/(2*hy) * (depth) - 2*ht*((itaold(i,j+1) - itaold(i,j-1))/(2*hy))*vold(i,j);
                %unew(i,j) = -((g * 2*ht)/(2*hx)) * (itaold(i+1,j) - itaold(i-1,j)) - (vold(i,j)*2*ht/(2*hy))*(uold(i,j+1) - uold(i,j-1)) - (uold(i,j)*2*ht/(2*hx))*(uold(i+1,j) - uold(i-1,j)) + uolder(i,j);
                %vnew(i,j) = -((g * 2*ht)/(2*hy)) * (itaold(i,j+1) - itaold(i,j-1)) - (vold(i,j)*2*ht/(2*hy))*(vold(i,j+1) - vold(i,j-1)) - (uold(i,j)*2*ht/(2*hx))*(vold(i+1,j) - vold(i-1,j)) + volder(i,j);
            
            end
        end
    end
    itaolder = itaold;
    itaold = itanew;
    uolder = uold;
    uold = unew;
    volder = vold;
    vold = vnew;
    
    upr = unew;
    vpr = vnew;
    itapr = itanew;
    
    upr(:,Ny) = upr(:,Ny-1);
    upr(1,:) = upr(2,:);
    upr(:,1) = upr(:,2);
    upr(Nx,:) = upr(Nx-1,:);
    
    vpr(Nx,:) = vpr(Nx-1,:);
    vpr(1,:) = vpr(2,:);
    vpr(:,Ny) = vpr(:,Ny-1);
    vpr(:,2) = vpr(:,1);
    
    itapr(:,1) = itapr(:,2);
    itapr(:,Ny) = itapr(:,Ny-1)
    itapr(1,:) = itapr(2,:);
    itapr(Nx,:) = itapr(Nx-1,:);
    
    mesh(X,Y,itapr');
    axis([0 200 0 10 -0.3 0.3]);
    %view(0,0);
    xlabel("Length of the flume");
    ylabel("Width of the flume");
    zlabel("Wave")
    Mesh(k) = getframe;
end