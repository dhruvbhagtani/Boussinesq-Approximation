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
Nt = 2000;
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
g= 9.81;

%LAWT Parameters
A = 0.2;
Tx = 20;
Ty = 2;

B = 1.5;
depth = 2.5;
knum = (2*pi)/98.63;
knum2 = (2*pi)/6.25;

omx = 2*pi/Tx;
omy = 2*pi/Ty;

disp(uold);
%At t=0, wave height:
%for i = 2:Nx-1
%        itaold(i,:) = 0.2*exp((-8*(hx*i-100)^2)/(50^2));
%end
for k = 1:Nt
    for i = 2:Nx-1
        for j = 2:Ny-1
            itanew(i,j) = itaold(i,j) - ht*((itaold(i+1,j) - itaold(i-1,j))/(2*hx)) * uold(i,j) - ht*((uold(i+1,j) - uold(i-1,j))/(2*hx)) * (depth) - ht*(vold(i,j+1) - vold(i,j-1))/(2*hy) * (depth) - ht*((itaold(i,j+1) - itaold(i,j-1))/(2*hy))*vold(i,j);
            unew(i,j) = -((g * ht)/(2*hx)) * (itaold(i+1,j) - itaold(i-1,j)) - (vold(i,j)*ht/(2*hy))*(uold(i,j+1) - uold(i,j-1)) - (uold(i,j)*ht/(2*hx))*(uold(i+1,j) - uold(i-1,j)) + uold(i,j);
            vnew(i,j) = -((g * ht)/(2*hy)) * (itaold(i,j+1) - itaold(i,j-1)) - (vold(i,j)*ht/(2*hy))*(vold(i,j+1) - vold(i,j-1)) - (uold(i,j)*ht/(2*hx))*(vold(i+1,j) - vold(i-1,j)) + vold(i,j);
            itanew(i,j) = itaold(i,j) - ht*((uold(i+1,j) - uold(i-1,j))/(2*hx)) * (depth) - ht*(vold(i,j+1) - vold(i,j-1))/(2*hy) * (depth);
            unew(i,j) = -((g * ht)/(2*hx)) * (itaold(i+1,j) - itaold(i-1,j)) + uold(i,j);
            vnew(i,j) = -((g * ht)/(2*hy)) * (itaold(i,j+1) - itaold(i,j-1)) + vold(i,j);
        end
    end
    for j = 1:Ny
       z = A*sin(knum*0 - omx*ht*k);
       %unew(1,j) = ((pi*2*A)/Tx) * (cosh(knum*(depth+z))/sinh(knum*depth)) * cos(knum*0 - omx*ht*k);
       unew(1,j) = (pi*2*A/(knum*Tx)*sin(knum*0 - omx*ht*k))/depth;
       %unew(1,j) = (g/depth)^0.5 * z;
       itanew(1,j) = z;
       %unew(1,j) = 0;
       %itanew(1,j) = itanew(2,j);
       unew(Nx,j) = 0;
       vnew(Nx,j) = vnew(Nx-1,j);
       itanew(Nx,j) = itanew(Nx-1,j);
       vnew(1,j) = 0;
    end
    for i = 1:Nx
        vnew(i,Ny) = 0;
        unew(i,Ny) = unew(i,Ny-1);
        unew(i,1) = unew(i,2);
        vnew(i,1) = 0;
        itanew(i,1) = itanew(i,2);
        itanew(i,Ny) = itanew(i,Ny-1);
    end
    itaold = itanew;
    uold = unew;
    vold = vnew;
    mesh(X,Y,itanew');
    view(0,0)
    %pause(0.1);
    xlabel("Length of the flume");
    ylabel("Width of the flume");
    zlabel("Wave")
    Mesh(k) = getframe;
end


