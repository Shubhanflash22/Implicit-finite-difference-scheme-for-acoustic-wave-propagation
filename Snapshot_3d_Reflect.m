%%Solving a PDE
clear;
%Equation
% wtt= c^2 wxx + c^2 wyy + c^2 wzz + f

%%Domain
%space
Lx=10;
Ly=10;
Lz=10;
dx=0.1;
dy=dx;
dz=dx;
nx=fix(Lx/dx);
ny=fix(Ly/dy);
nz=fix(Lz/dz);
x=linspace(0,Lx,nx);
y=linspace(0,Ly,ny);
z=linspace(0,Lz,nz);

%Time
T=10;

%%Field variables
%Variables
wn=zeros(nx,ny,nz);
wnm1=wn; % w at time n-1
wnp1=wn; % w at time n+1

%Parameters
CFL=0.5;
c=1;
dt=CFL*dx/c;

%%Initial Conditions
%% Time Stepping loop
t=0;

while(t<T)
  %Reflecting boundary conditions
  wn(:,[1 end],:)=0;
  wn([1 end],:,:)=0;
  wn(:,:,[1 end])=0;
  
  %Absorbing boundary conditions
  %wnp1(1,:,:)=wn(2,:,:) + ((CFL-1)/(CFL+1))* (wnp1(2,:,:)-wn (1,:,:)); 
  %wnp1(end,:,:)=wn(end-1,:,:)+ ((CFL-1)/(CFL+1))*(wnp1(end-1,:,:)-wn(end,:,:));
  %wnp1(:,1,:)=wn(:,2,:)+ ((CFL-1)/(CFL+1))*(wnp1(:,2,:)-wn (:,1,:));
  %wnp1(:,end,:)=wn(:, end-1,:) + ((CFL-1)/(CFL+1))*(wnp1(:, end-1,:)-wn(:, end,:));
  %wnp1(:,:,1)=wn(:,:,2)+ ((CFL-1)/(CFL+1))*(wnp1(:,:,2)-wn (:,:,1));
  %wnp1(:,:,end)=wn(:,:,end-1) + ((CFL-1)/(CFL+1))*(wnp1(:,:,end-1)-wn(:,:,end));
  
  %Solution
  t=t+dt;
  wnm1=wn; wn=wnp1; %Saving current and previous arrays
  
  %Source
  wn(50,50,50)=dt^2*200*sin(30*pi*t/20);
  
  for i=2:nx-1, for j=2:ny-1, for k=2:nz-1
          wnp1(i,j,k)= 2*wn(i,j,k)-wnm1(i,j,k)...
              +CFL^2*(wn(i+1,j,k)+wn(i,j+1,k)+wn(i,j,k+1)-6*wn(i,j,k)+wn(i-1,j,k)+wn(i,j-1,k)+wn(i,j,k-1));
      end, end, end
  
  %Check convergence
  %Visualize at selected steps
  clf;
  subplot(2,1,1);
  imagesc(x,y,wn(:,:,nz/2)'); colorbar; caxis([-0.02,0.02])
  title(sprintf('t=%.2f' , t));
  subplot(2,1,2);
  mesh(x,y,wn(:,:,nz/2)'); colorbar; caxis([-0.02 0.02])
  axis([0 Lx 0 Ly 0 Lz]);
  shg; pause(0.01);
end
