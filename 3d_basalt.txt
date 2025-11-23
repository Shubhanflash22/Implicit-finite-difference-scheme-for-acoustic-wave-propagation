% Define parameters
nx = 601;           % Number of grid points in x direction
ny = 601;           % Number of grid points in y direction
nz = 355;           % Number of grid points in z direction
nt = 1000;           % Number of time steps
dx = 10;            % Grid spacing in x direction (m)
dy = 10;            % Grid spacing in y direction (m)
dz = 10;            % Grid spacing in z direction (m)
dt = 0.001;         % Time step (s)
vp = 3500;          % P-wave velocity (m/s)
rho = 1000;         % Density (kg/m^3)
f0 = 20;            % Dominant frequency of the source (Hz)
ksrc = round(nz/2);                 % source location along OZ
jsrc = round(ny/2);                 % source location along OY
isrc = round(nx/2);                 % source location along OX
xsrc = isrc + nx*(jsrc-1) + nx*ny*(ksrc-1);   % Source index in grid
src = zeros(nt,1);
for i = 1:nt
    t = i*dt;
    src(i) = exp(-pi^2*f0^2*(t-1.5*dt)^2);
end

% Initialize matrices
u = zeros(nx,ny,nz);   % Wavefield
v = zeros(nx,ny,nz);   % Velocity

% Fill velocity matrix
v(:,:,:) = vp^2*rho;

% Compute coefficients
c1 = (dt/dx)^2;
c2 = (dt/dy)^2;
c3 = (dt/dz)^2;

% Apply boundary conditions
u(1,:,:) = 0;
u(nx,:,:) = 0;
u(:,1,:) = 0;
u(:,ny,:) = 0;
u(:,:,1) = 0;
u(:,:,nz) = 0;

% Main loop
for k = 1:nt
    % Compute time derivative
    dudt = zeros(nx,ny,nz);
    dudt(2:nx-1,2:ny-1,2:nz-1) = (c1*(v(2:nx-1,2:ny-1,2:nz-1).*(u(3:nx,2:ny-1,2:nz-1)-2*u(2:nx-1,2:ny-1,2:nz-1)+u(1:nx-2,2:ny-1,2:nz-1))) + ...
                           c2*(v(2:nx-1,2:ny-1,2:nz-1).*(u(2:nx-1,3:ny,2:nz-1)-2*u(2:nx-1,2:ny-1,2:nz-1)+u(2:nx-1,1:ny-2,2:nz-1))) + ...
                            c3*(v(2:nx-1,2:ny-1,2:nz-1).*(u(2:nx-1,2:ny-1,3:nz)-2*u(2:nx-1,2:ny-1,2:nz-1)+u(2:nx-1,2:ny-1,1:nz-2))))/2;
    
    % Add source term
    dudt(isrc,jsrc,ksrc) = dudt(isrc,jsrc,ksrc) + src(k);
     
    % Update wavefield
    u = u + dt*dudt;
    
    % Display wavefield
    imagesc([squeeze(u(round(nx/2),:,:)) squeeze(u(:,round(ny/2),:)) squeeze(u(:,:,round(nz/2)))]); axis equal tight;
    title(['Time Step: ',num2str(k)])
    colorbar;
    drawnow;
end