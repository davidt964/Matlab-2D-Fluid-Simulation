% Jonny Duong and David Tran (8/21/20)
% Program: Final Project
% {
% This script allows the study of fluid behavior in 2D using smoothed particle hydrodynamics (SPH) on a spatially-hashed 
% domain. The forces driving a single particle in your simulation are relatively simple, but the combined behavior of
% the entire particle system can model some amazingly complicated fluid behaviors. This task is accomplished by using
% an array of structures, allowing each particle to have the same set of characteristics: position, velocity, force,
% density, etc. Multiple functions are created to ease readability and expedite debugging. 
% }

% These are the clear and close commands.
clc, close all, clear all;

% ============================================================================================================
%% Part 1: Data Structures
% Set the number of particles and the bins that contain them.
N = 15^2;               % Particles in our system

% Initialize structure for particles and bins.
posVal = zeros(1,2);
velVal = zeros(1,2);
forceVal = zeros(1,2);
particles = [];
adjBins = [];

% Please note that we will use this information for part 3!!
xMax = sqrt(N);
yMax = sqrt(N);
h = 1;                  % Minimum bin dimension/smoothing radius.
Nx = floor(xMax/h);     % This will create the amount of bins in the x direction.
Ny = floor(yMax/h);     % Same as above but for the y coordinate.
dx = xMax/Nx;           % Actual bin dimensions.
dy = yMax/Ny;           % Same as above but for the y coordinate.
numBins = Nx * Ny;      % The bins that contain the particles in the system.

IDs(1:N) = struct('pos', posVal, 'vel', velVal,'force',forceVal,'den',0,'neigh',[]);
BinIDs(1:numBins) = struct('part',particles,'adj',adjBins);

% ============================================================================================================
%% Part 2: Initializing Fluid Particles
% Initialize the fluid parameters, gridSize, and smoothing radius.
p_rho0 = 1000;          % Rest density
p_rho1 = 2000;          % Rest Density of Second Fluid
p_mass = p_rho0/N;      % Mass of a single particle 
p_k = 100;              % Stiffness constant
p_mu = 0.1;             % Viscosity coefficient

% Set the time-stepping parameters.
dt = 0.02;
t_0 = 0;
t_f = 0.5;
t_steps = ceil(t_f/dt);

draw_every = 1;
v = VideoWriter('Animation.mp4');
open(v);
fig = figure(1);
% Assign each of the particles, k = 1,2, . . . N, a position on the 2D grid.
    % 15x15 box of particles
    % X from 0.25 to 3.75, Y from 0.25 to 3.75
    % 0.25 unit separation
x = 0.25;
IDcounter = 1;
y = 0.25;
for i = 1:15
    x = 0.25;
    for j = 1:15
        IDs(j + 15 * (i-1)).pos(1) = x;
        IDs(j + 15 * (i-1)).pos(2) = y;
        x = x + 0.25;
        IDcounter = IDcounter + 1;
    end
    y = y + 0.25;
end

% ============================================================================================================
%% Part 3: Initializing Bins

% Calculate the bin number for every particle in your simulation using the hash function.
for k = 1:N
  binNum = (ceil(IDs(k).pos(1)/dx) - 1)*Ny + ceil((yMax - IDs(k).pos(2))/dy);
  BinIDs(binNum).part(k) = [k];
  BinIDs(binNum).part = nonzeros(BinIDs(binNum).part);
end

% ============================================================================================================
%% LAST PART (OUT OF 8 PARTS): Visualizing Results
count = 0;
for k = 1:t_steps
   
    % Plot and update the figure data
    figure(1)
    grid on
    for i = 1:N
        xVal = IDs(i).pos(1);
        yVal = IDs(i).pos(2);
        plot(xVal, yVal,'b.');
        hold on
    end
    title('Smoothed Particle Hydrodynamic Fluid Simulation')
    xlim([0,xMax])
    ylim([0,yMax])
    BC = [xlim, ylim];
    set(gca, 'xtick', 0:dx:xMax)
    set(gca, 'ytick', 0:dy:yMax)
   
    frame = getframe(fig);
    writeVideo(v,frame);
   
    % Get the neighbor of every particle
    IDs = get_neighbors(BinIDs, IDs, numBins, Nx, Ny, h);
   
    % Calculate the density
    IDs = get_density(IDs, p_mass, h);
   
    % Calculate the force
    %if k <= N/2 (Code for mixing of two fluids with different rest densities)
        IDs = get_forces(N, IDs, h, p_mass, p_k, p_rho0, p_mu);
    %else
     %   IDs = get_forces(N, IDs, h, p_mass, p_k, p_rho1, p_mu);
    %end
   
    % Update kinematics and adjust particles due to grid boundaries
    IDs = update_kinematics(dt, IDs, BC, N);
    hold off
    count = count + 1;
end

close(v);

% ============================================================================================================
%% See Part 3 (Initializing Bins)
function adj_bins = getAdjacentBins(Nx, Ny, z)
% {
% getAdjacentBins returns a vector containing the linearly indexed ID
% numbers on a grid of size [Ny x Nx] for all bins adjacent to bin z PLUS
% the bin z itself.
% }

% Construct a vector that contains the bin IDs of every adjacent bin. 
if mod(z, Ny) == 0
    % Bottom Wall
    if z == Ny
        % Bottom-left corner
        adj_bins = [z-1, z, z+Ny-1, z+Ny];
    elseif z == Nx*Ny
        % Bottom right corner
        adj_bins = [z-Ny-1, z-Ny, z-1, z];
    else
        % Bottom wall
        adj_bins = [z-Ny-1, z-Ny, z-1, z, z+Ny-1, z+Ny];
    end
elseif mod(z,Ny) == 1
    % Top Wall
    if z == 1
        % Top-left corner
        adj_bins = [z, z+1, z+Ny, z+Ny+1];
    elseif z == Nx*Ny - Ny + 1
        % Top right corner 
        adj_bins = [z-Ny, z-Ny+1, z, z+1];
    else
        % Top wall
        adj_bins = [z-Ny, z-Ny+1, z, z+1, z+Ny, z+Ny+1];
    end
    
elseif z < Ny
    % Left Wall
    adj_bins = [z-1, z, z+1, z+Ny-1, z+Ny, z+Ny+1];
elseif z > Nx*Ny - Ny + 1
    % Right Wall
    adj_bins = [z-Ny-1, z-Ny, z-Ny+1, z-1, z, z+1];
else
    % Interior
    adj_bins = [z-Ny-1, z-Ny, z-Ny+1, z-1, z, z+1, z+Ny-1, z+Ny, z+Ny+1];
end

end

% ============================================================================================================
%% See Part 4:Identifying Neighbors

function IDs = get_neighbors(BinIDs, IDs, numBins, Nx, Ny, h)
%{
Identifies a particle's neighbors by first obtaining all of the particles
in a given bin. Then, it determines that bin's neighbors and the particles
in those neighbors. Ultimately, if the particles are within a smoothing
radius h, they are considered neighbors.
%}
for z = 1:numBins
    % Get all particles in bin z.
    binParticles = BinIDs(z).part;
    p = 0;
    % if bin z is not empty.
    if ~isempty(binParticles)
        % Get vector of bin IDs adjacent to bin z.
        binsAdj = getAdjacentBins(Nx,Ny,z);
        for w = [z, binsAdj]
            % Get all particles in bin w.
            particles = BinIDs(w).part;
            % all particles in bin z
            for k = 1:length(binParticles) 
                % (x,y) coordinate of particle k in bin z
                xk = IDs(binParticles(k,1)).pos;
                % all particles in bin w
                for j = 1:length(particles)
                    % (x,y) coordinate of particle j in bin w
                    xj = IDs(particles(j,1)).pos;
                    % distance b/t particle k in bin z and particle j in bin w
                    distance = norm(xk - xj); 
                    % if dist < h and k ~= j
                    if distance < h && binParticles(k,1) ~= particles(j,1)
                        p = p + 1;
                        %k and j are neighbors
                        IDs(binParticles(k,1)).neigh(p) = particles(j,1);
                        IDs(binParticles(k,1)).neigh = nonzeros(IDs(binParticles(k,1)).neigh);
                        IDs(binParticles(k,1)).neigh = unique(IDs(binParticles(k,1)).neigh);
                    end
                end
            end
        end
    end
end
end

% ============================================================================================================
%% See Part 5: Calculating Density
function IDs = get_density(IDs, p_mass, h)
% {
% Now that the list of neighbors is calculated for each particle, we can calculate the fluid density represented 
% at  each point. The fluid density at particle k, ?k, depends on the neighbor particles’ mass and the separation 
% distances.  Please note that the  first term on the right-hand side of Equation 2 is due to particle k’s self 
% mass and the terms in the summation are due to the k’s neighbors.
% }

constant1 = (4*p_mass)/(pi*h^2);
constant2 = (4*p_mass)/(pi*h^8);
for k = 1:length(IDs)
    currPar = IDs(k);
    currParNeighbors = currPar.neigh;
    rho_neighbors = 0;
    for j = 1:length(currParNeighbors)
        currNeighbors = IDs(currParNeighbors(j));
        xk = currPar.pos;
        xj = currNeighbors.pos;
        n = norm(xk - xj);
        rho_neighbors = rho_neighbors + (h^2 - n^2)^3;
    end
    IDs(k).den = (constant2*rho_neighbors) + constant1;
end

end

% ============================================================================================================
%% See Part 6: Calculating Forces
function IDs = get_forces(N, IDs, h, p_mass, p_k, p_rho0, p_mu)
% {
% We’ll now make use of both the stored density values and the list of neighbors to calculate the total
% force on particle k, fk. Each particle, k = 1, 2, ..., N, is affected by the sum of the fluid force
% contributions from neighbors plus any external forces like gravity, wind, etc.
% }

for k = 1:N
    neighborsOfk = IDs(k).neigh;
    sum_fkj = 0;
    f_external = [0, -9.8];
    f_fountain = [0, 1e6];
    for j = 1:length(neighborsOfk)
        qkj = norm(IDs(k).pos - IDs(neighborsOfk(j,1)).pos)/h;
        % s1, s2, ... s5 are all scalars!
        s1 = p_mass/(pi*h^4*IDs(j).den);
        s2 = (1 - qkj);
        s3 = 15*p_k*(IDs(k).den + IDs(neighborsOfk(j,1)).den - 2*p_rho0); 
        s4 = (1 - qkj)/qkj;
        x_dist = IDs(k).pos - IDs(neighborsOfk(j,1)).pos;
        s5 = 40*p_mu;
        v_dist = IDs(k).vel - IDs(neighborsOfk(j,1)).vel;
        sum_fkj = sum_fkj + (s1*s2*((s3*s4*x_dist) -(s5*v_dist)));
    end
    IDs(k).force = f_external + sum_fkj + f_fountain;
end

end

% ============================================================================================================
%% Part 7: Update Kinematics
function IDs = update_kinematics(dt, IDs, BC, N)
%{
Once the total force is known, the kinematics can be updated accordingly.
The new velocity vkp1 will be calculated using explicit Euler while the new
position xkp1 will be calculated using implicit Euler. In conjunction,
this approach is the semi-implicit Euler method. If collision at a wall
occurs, only the position and velocity components in the direction normal 
to the wall are updated.
%}
for k = 1:N
IDs(k).vel(1) = dt * (IDs(k).force(1)./IDs(k).den) + IDs(k).vel(1);
IDs(k).vel(2) = dt * (IDs(k).force(2)./IDs(k).den) + IDs(k).vel(2);
IDs(k).pos(1) = dt * IDs(k).vel(1) + IDs(k).pos(1);
IDs(k).pos(2) = dt * IDs(k).vel(2) + IDs(k).pos(2);
b = 0.9;    % Damping constant
% Check for collison with wall
if IDs(k).pos(1) <= BC(1) + 0.1
    % Collision with left wall
    IDs(k).pos(1) = 2 * (BC(1)) - IDs(k).pos(1);
    IDs(k).vel(1) = -b * IDs(k).vel(1);
elseif IDs(k).pos(1) >= BC(2) - 0.1
    % Collision with right wall
    IDs(k).pos(1) = 2 * (BC(2)) - IDs(k).pos(1);
    IDs(k).vel(1) = -b * IDs(k).vel(1);
elseif IDs(k).pos(2) <= BC(3) + 0.1
    % Collision with bottom wall
    IDs(k).pos(2) = 2 * (BC(3)) - IDs(k).pos(2);
    IDs(k).vel(2) = -b * IDs(k).vel(2);
elseif IDs(k).pos(2) >= BC(4) - 0.1
    % Collision with top wall
    IDs(k).pos(2) = 2 * (BC(4)) - IDs(k).pos(2);
    IDs(k).vel(2) = -b * IDs(k).vel(2);
end
end
end

