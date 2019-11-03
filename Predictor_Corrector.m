close all;
clc;
clear;
%Simulation grid parameters
Nx = 200;
Ny = 10;
L = 400;
d =  10;
hx = L/Nx;
hy = d/Ny;
Nt = 4000;
ht = 0.02;
X = zeros(Nx,1);
Y = zeros(Ny,1);
for i = 1:Nx
    X(i,1) = i*hx;
end
for j = 1:Ny
    Y(j,1) = j*hy;
end

itac = zeros(Nx,Ny);
itap = zeros(Nx,Ny);
itaold = zeros(Nx,Ny);

up = zeros(Nx,Ny);
uc = zeros(Nx,Ny);
uold = zeros(Nx,Ny);

vc = zeros(Nx,Ny);
vp = zeros(Nx,Ny);
vold = zeros(Nx,Ny);

g= 9.81;
A = 0.2;
Tx = 20;
depth = 2.5;
knum = (2*pi)/98.63;
omx = 2*pi/Tx;

for k = 1:Nt
    z = A*sin(hx*knum-omx*ht*k);
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
    for i = 2:Nx-1
        for j = 2:Ny-1
            itap(i,j) = itaold(i,j) - ht*((uold(i+1,j) - uold(i-1,j))/(2*hx)) * (depth) - ht*(vold(i,j+1) - vold(i,j-1))/(2*hy) * (depth);
            up(i,j) = -((g * ht)/(2*hx)) * (itaold(i+1,j) - itaold(i-1,j)) + uold(i,j);
            vp(i,j) = -((g * ht)/(2*hy)) * (itaold(i,j+1) - itaold(i,j-1)) + vold(i,j);
            %itap(i,j) = itaold(i,j) - ht*((itaold(i+1,j) - itaold(i-1,j))/(2*hx)) * uold(i,j) - ht*((uold(i+1,j) - uold(i-1,j))/(2*hx)) * (depth) - ht*(vold(i,j+1) - vold(i,j-1))/(2*hy) * (depth) - ht*((itaold(i,j+1) - itaold(i,j-1))/(2*hy))*vold(i,j);
            %up(i,j) = -((g * ht)/(2*hx)) * (itaold(i+1,j) - itaold(i-1,j)) - (vold(i,j)*ht/(2*hy))*(uold(i,j+1) - uold(i,j-1)) - (uold(i,j)*ht/(2*hx))*(uold(i+1,j) - uold(i-1,j)) + uold(i,j);
            %vp(i,j) = -((g * ht)/(2*hy)) * (itaold(i,j+1) - itaold(i,j-1)) - (vold(i,j)*ht/(2*hy))*(vold(i,j+1) - vold(i,j-1)) - (uold(i,j)*ht/(2*hx))*(vold(i+1,j) - vold(i-1,j)) + vold(i,j);
        end
    end
    z1 = A*sin(hx*knum-omx*(ht)*(k+1/2));
    up(1,:) = (pi*2*z1)/(knum*Tx)/depth;
    itap(1,:) = z1;
    up(Nx,:) = 0;
    vp(Nx,:) = vp(Nx-1,:);
    itap(Nx,:) = itap(Nx-1,:);
    vp(1,:) = 0;
    vp(:,Ny) = 0;
    up(:,Ny) = up(:,Ny-1);
    up(:,1) = up(:,2);
    vp(:,1) = 0;
    itap(:,1) = itap(:,2);
    itap(:,Ny) = itap(:,Ny-1);
    for i = 2:Nx-1
        for j = 2:Ny-1
            itac(i,j) = itaold(i,j) - (ht/2*((uold(i+1,j) - uold(i-1,j))/(2*hx)) * (depth+z1) - ht/2*(vold(i,j+1) - vold(i,j-1))/(2*hy) * (depth+z1)) - (ht/2*((up(i+1,j) - up(i-1,j))/(2*hx)) * (depth+z1) - ht/2*(vp(i,j+1) - vp(i,j-1))/(2*hy) * (depth+z1));
            uc(i,j) = uold(i,j) - (((g * ht)/(4*hx)) * (itaold(i+1,j) - itaold(i-1,j))) - (((g * ht)/(4*hx)) * (itap(i+1,j) - itap(i-1,j)));
            vc(i,j) = vold(i,j) - (((g * ht)/(4*hy)) * (itaold(i,j+1) - itaold(i,j-1))) - (((g * ht)/(4*hy)) * (itap(i,j+1) - itap(i,j-1)));
            %itac(i,j) = itaold(i,j) - (ht/2*((itaold(i+1,j) - itaold(i-1,j))/(2*hx)) * uold(i,j) + ht/2*((uold(i+1,j) - uold(i-1,j))/(2*hx)) * (depth + z) + ht/2*(vold(i,j+1) - vold(i,j-1))/(2*hy) * (depth + z) + ht/2*((itaold(i,j+1) - itaold(i,j-1))/(2*hy))*vold(i,j)) - (ht/2*((itap(i+1,j) - itap(i-1,j))/(2*hx)) * up(i,j) + ht/2*((up(i+1,j) - up(i-1,j))/(2*hx)) * (depth+z1) + ht/2*(vp(i,j+1) - vp(i,j-1))/(2*hy) * (depth+z1) + ht/2*((itap(i,j+1) - itap(i,j-1))/(2*hy))*vp(i,j));
            %uc(i,j) = uold(i,j) - (((g * ht)/(4*hx)) * (itaold(i+1,j) - itaold(i-1,j)) - (vold(i,j)*ht/(4*hy))*(uold(i,j+1) - uold(i,j-1)) - (uold(i,j)*ht/(4*hx))*(uold(i+1,j) - uold(i-1,j))) - (((g * ht)/(4*hx)) * (itap(i+1,j) - itap(i-1,j)) - (vp(i,j)*ht/(4*hy))*(up(i,j+1) - up(i,j-1)) - (up(i,j)*ht/(4*hx))*(up(i+1,j) - up(i-1,j)));
            %vc(i,j) = vold(i,j) - (((g * ht)/(4*hy)) * (itaold(i,j+1) - itaold(i,j-1)) - (vold(i,j)*ht/(4*hy))*(vold(i,j+1) - vold(i,j-1)) - (uold(i,j)*ht/(4*hx))*(vold(i+1,j) - vold(i-1,j))) - (((g * ht)/(4*hy)) * (itap(i,j+1) - itap(i,j-1)) - (vp(i,j)*ht/(4*hy))*(vp(i,j+1) - vp(i,j-1)) - (up(i,j)*ht/(4*hx))*(vp(i+1,j) - vp(i-1,j)));
        end
    end
    
    itaold = itac;
    uold = uc;
    vold = vc;
    

    upr = uc;
    vpr = vc;
    itapr = itac;
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
    
%{    
    itapr = itac;
    upr = uc;
    vpr = vc;
    z1 = A*sin(-omx*(ht)*(k+1));
    upr(1,:) = (pi*2*z1)/(knum*Tx)/depth;
    itapr(1,:) = z1;
    upr(Nx,:) = 0;
    vpr(Nx,:) = vpr(Nx-1,:);
    itapr(Nx,:) = itapr(Nx-1,:);
    vpr(1,:) = 0;
    vpr(:,Ny) = 0;
    upr(:,Ny) = upr(:,Ny-1);
    upr(:,1) = upr(:,2);
    vpr(:,1) = 0;
    itapr(:,1) = itapr(:,2);
    itapr(:,Ny) = itapr(:,Ny-1);
%}
    mesh(X,Y,itapr');
    axis([0 200 0 10 -0.3 0.3]);
    view(0,0);
    xlabel("Length of the flume");
    ylabel("Width of the flume");
    zlabel("Wave");
    Mesh(k) = getframe;
end
                