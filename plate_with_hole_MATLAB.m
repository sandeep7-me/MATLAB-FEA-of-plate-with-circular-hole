
%  2D Plate WITH A CIRCULAR HOLE — Triangular Mesh + Plane-Stress FEM
%  Element type: CST (Constant Strain Triangle), linear elastic, plane stress

clear; clc; close all;

% GEOMETRY & MESH PARAMETERS
Lx = 100;            % plate length (x-direction in mm)
Ly = 50;             % plate height (y-direction in mm)

hole_cx = Lx/2;     % hole center x
hole_cy = Ly/2;     % hole center y 
hole_r  = 10;      % hole radius

nPerSideLong  = 100;   % boundary points along the long (x) edges
nPerSideShort = 50;   % boundary points along the short (y) edges
nHole         = 60;   % boundary points around the hole
nInteriorX    = 85;   % interior sampling grid resolution (x)
nInteriorY    = 30;   % interior sampling grid resolution (y)
holeClearance = 0.75; % keep interior points at least this many radii away

% BOUNDARY POINTS (outer rectangle)

e1x = linspace(0,  Lx, nPerSideLong);  e1x(end) = []; e1y = zeros(size(e1x));   
e2y = linspace(0,  Ly, nPerSideShort); e2y(end) = []; e2x = Lx*ones(size(e2y)); 
e3x = linspace(Lx, 0,  nPerSideLong);  e3x(end) = []; e3y = Ly*ones(size(e3x)); 
e4y = linspace(Ly, 0,  nPerSideShort); e4y(end) = []; e4x = zeros(size(e4y)); 

bx = [e1x, e2x, e3x, e4x];
by = [e1y, e2y, e3y, e4y];

outerPts = [bx(:), by(:)];
nOuter   = size(outerPts,1);
outerEdges = [(1:nOuter)', [(2:nOuter)'; 1]];      

% BOUNDARY POINTS (hole)
theta = linspace(0, 2*pi, nHole+1); theta(end) = [];
holePts = [hole_cx + hole_r*cos(theta(:)), hole_cy + hole_r*sin(theta(:))];
nHolePts = size(holePts,1);
holeEdgesLocal = [(1:nHolePts)', [(2:nHolePts)'; 1]];

% INTERIOR POINTS
margin = min(Lx,Ly)/60;
[gx, gy] = meshgrid(linspace(margin, Lx-margin, nInteriorX), ...
                     linspace(margin, Ly-margin, nInteriorY));
interiorPts = [gx(:), gy(:)];
d = hypot(interiorPts(:,1)-hole_cx, interiorPts(:,2)-hole_cy);
interiorPts(d < holeClearance*hole_r, :) = [];

% COMBINE POINTS & CONSTRAINT EDGES
P = [outerPts; holePts; interiorPts];
holeEdges = holeEdgesLocal + nOuter;   
C = [outerEdges; holeEdges];

% Remove duplicate points
[P, ~, ic] = uniquetol(P, 1e-9, 'ByRows', true);
C = ic(C);

% CONSTRAINED DELAUNAY TRIANGULATION 
DT = delaunayTriangulation(P, C);
IO = isInterior(DT);                       % keep triangles inside plate-minus-hole
elements = DT.ConnectivityList(IO, :);
nodes    = DT.Points;
numNodesRaw = size(nodes,1);
numElements = size(elements,1);

usedNodes = unique(elements(:));
remap = zeros(numNodesRaw,1);
remap(usedNodes) = 1:numel(usedNodes);
nodes = nodes(usedNodes, :);
elements = remap(elements);
numNodes = size(nodes,1);

fprintf('Mesh generated: %d nodes, %d triangular elements (plate with hole).\n', ...
        numNodes, numElements);

% Check: total mesh area should equal Lx*Ly - pi*r^2
totalArea = 0;
for k = 1:numElements
    v = nodes(elements(k,:), :);
    totalArea = totalArea + 0.5*abs((v(2,1)-v(1,1))*(v(3,2)-v(1,2)) - ...
                                     (v(3,1)-v(1,1))*(v(2,2)-v(1,2)));
end
fprintf('Mesh area = %.4f   (expected %.4f, %.2f%% diff)\n', ...
        totalArea, Lx*Ly - pi*hole_r^2, ...
        100*abs(totalArea-(Lx*Ly-pi*hole_r^2))/(Lx*Ly-pi*hole_r^2));

% VISUALIZE MESH
figure('Name','Mesh: Plate with Hole','Color','w');
triplot(elements, nodes(:,1), nodes(:,2), 'Color',[0.2 0.4 0.8]);
hold on; axis equal tight;
th = linspace(0,2*pi,80);
plot(hole_cx+hole_r*cos(th), hole_cy+hole_r*sin(th), 'r:', 'LineWidth', 1.2);
xlabel('X'); ylabel('Y');
title(sprintf('Plate with Hole Mesh (%d elements)', numElements));

%  FINITE ELEMENT ANALYSIS: CST PLANE STRESS

% MATERIAL & LOAD PROPERTIES 
E  = 200e3;      % Youngs Modulus (inMPa)
nu = 0.3;        % Poisson's ratio
t  = 1;          % plate thickness (mm)
P_total = 500;  % total AXIAL tensile load at free end (N), positive = pulling in +x

% Plane-stress constitutive matrix
D = E/(1-nu^2) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];

% 9. GLOBAL STIFFNESS ASSEMBLY
numDof = 2*numNodes;

I_idx = zeros(36*numElements,1);
J_idx = zeros(36*numElements,1);
V_val = zeros(36*numElements,1);
ptr = 1;

elemArea = zeros(numElements,1);
elemB    = cell(numElements,1);

for k = 1:numElements
    n123 = elements(k,:);
    v = nodes(n123,:);
    x1=v(1,1); y1=v(1,2); x2=v(2,1); y2=v(2,2); x3=v(3,1); y3=v(3,2);

    Ae = 0.5*((x2*y3-x3*y2) - (x1*y3-x3*y1) + (x1*y2-x2*y1));
    if Ae <= 0
        error('Element %d has non-positive area — check node ordering.', k);
    end

    b1=y2-y3; b2=y3-y1; b3=y1-y2;
    c1=x3-x2; c2=x1-x3; c3=x2-x1;

    B = (1/(2*Ae)) * [ b1  0  b2  0  b3  0;
                         0 c1   0 c2   0 c3;
                        c1 b1  c2 b2  c3 b3];

    ke = t*Ae * (B' * D * B);

    elemArea(k) = Ae;
    elemB{k} = B;

    % Element DOF map: [u1 v1 u2 v2 u3 v3]
    eDofs = zeros(1,6);
    eDofs([1 3 5]) = 2*n123 - 1;   % ux of node 1,2,3
    eDofs([2 4 6]) = 2*n123;       % uy of node 1,2,3

    [ii, jj] = ndgrid(eDofs, eDofs);
    I_idx(ptr:ptr+35) = ii(:);
    J_idx(ptr:ptr+35) = jj(:);
    V_val(ptr:ptr+35) = ke(:);
    ptr = ptr + 36;
end

K = sparse(I_idx, J_idx, V_val, numDof, numDof);

%  BOUNDARY CONDITIONS
tol = 1e-6 * max(Lx,Ly);

fixedNodes = find(abs(nodes(:,1) - 0)  < tol);   
loadNodes  = find(abs(nodes(:,1) - Lx) < tol);  

if isempty(fixedNodes) || isempty(loadNodes)
    error('Could not find fixed/load edge nodes — check mesh boundary tolerance.');
end

fixedDofs = reshape([2*fixedNodes-1, 2*fixedNodes]', [], 1);

F = zeros(numDof,1);
% Distribute total AXIAL load equally among free-end nodes (uniform nodal load,
F(2*loadNodes-1) = P_total / numel(loadNodes);   % x-direction (axial)

freeDofs = setdiff((1:numDof)', fixedDofs);

% SOLVE
U = zeros(numDof,1);
U(freeDofs) = K(freeDofs,freeDofs) \ F(freeDofs);

Ux = U(1:2:end);
Uy = U(2:2:end);
Umag = hypot(Ux,Uy);

fprintf('\nMax displacement magnitude: %.5f mm at node %d\n', ...
        max(Umag), find(Umag==max(Umag),1));
fprintf('Free-end average axial (x) elongation: %.5f mm\n', mean(Ux(loadNodes)));

% Reference 1: simple axial bar elongation WITHOUT hole (sanity check only)
A_gross = Ly*t;
delta_bar = P_total*Lx/(A_gross*E);
fprintf('Reference (no-hole, axial bar theory) elongation: %.5f mm\n', delta_bar);
fprintf(['  (Ignores the hole and the fully-clamped end effect, but should be the\n' ...
         '   same order of magnitude as a sanity check.)\n']);

% Reference 2: classical stress-concentration factor Kt for a circular hole

d_hole  = 2*hole_r;
A_gross = Ly*t;
A_net   = (Ly - d_hole)*t;
sigma_gross = P_total / A_gross;
sigma_net   = P_total / A_net;
dW = d_hole/Ly;

fprintf('Hole diameter / plate width, d/W = %.3f\n', dW);
fprintf('Nominal stress (gross section): %.3f MPa\n', sigma_gross);
fprintf('Nominal stress (net section):   %.3f MPa\n', sigma_net);

% STRESS RECOVERY
sigma_x  = zeros(numElements,1);
sigma_y  = zeros(numElements,1);
tau_xy   = zeros(numElements,1);
vonMises = zeros(numElements,1);

for k = 1:numElements
    n123 = elements(k,:);
    eDofs = zeros(1,6);
    eDofs([1 3 5]) = 2*n123 - 1;
    eDofs([2 4 6]) = 2*n123;

    ue = U(eDofs);                
    strain = elemB{k} * ue;      
    stress = D * strain;

    sigma_x(k) = stress(1);
    sigma_y(k) = stress(2);
    tau_xy(k)  = stress(3);
    vonMises(k) = sqrt(stress(1)^2 - stress(1)*stress(2) + stress(2)^2 + 3*stress(3)^2);
end

centroids = zeros(numElements,2);
for k = 1:numElements
    centroids(k,:) = mean(nodes(elements(k,:),:), 1);
end

[maxVM, idxMax] = max(vonMises);
fprintf('Global max von Mises stress: %.3f MPa, element %d, at (x=%.3f, y=%.3f)\n', ...
        maxVM, idxMax, centroids(idxMax,1), centroids(idxMax,2));
fprintf(['  NOTE: the global max often sits right at the fixed-end corner. That is a\n' ...
         '  known FEM modeling artifact (idealizing the support as a rigid line of fixed\n' ...
         '  nodes creates a stress singularity there), not a physical result — it would\n' ...
         '  keep growing as the mesh is refined near that corner.\n']);


distFromHole = hypot(centroids(:,1)-hole_cx, centroids(:,2)-hole_cy);
nearHole = distFromHole < 2*hole_r;
[maxVM_hole, idxHole] = max(vonMises(nearHole));
holeElems = find(nearHole);
idxHoleGlobal = holeElems(idxHole);
fprintf('Max von Mises stress NEAR THE HOLE (within 2r): %.3f MPa, at (x=%.3f, y=%.3f)\n', ...
        maxVM_hole, centroids(idxHoleGlobal,1), centroids(idxHoleGlobal,2));

Kt_net   = maxVM_hole / sigma_net;
Kt_gross = maxVM_hole / sigma_gross;
fprintf('Kt (net-section, standard handbook convention):   %.3f   <- compare this to published charts\n', Kt_net);
fprintf('Kt (gross-section, full un-notched width):        %.3f\n', Kt_gross);
fprintf(['(Net-section Kt should be <= 3, approaching 3 only as d/W -> 0, and decreasing\n' ...
         ' toward ~2 as d/W -> 1 — matching standard Peterson/Roark/Shigley charts.)\n\n']);

% Nodal-averaged von Mises
nodeVM = zeros(numNodes,1);
nodeCount = zeros(numNodes,1);
for k = 1:numElements
    nodeVM(elements(k,:))    = nodeVM(elements(k,:)) + vonMises(k);
    nodeCount(elements(k,:)) = nodeCount(elements(k,:)) + 1;
end
nodeVM = nodeVM ./ max(nodeCount,1);


% VISUALIZATION: VON MISES STRESS CONTOUR
figure('Name','Von Mises Stress','Color','w');
patch('Faces', elements, 'Vertices', nodes, ...
      'FaceVertexCData', nodeVM, 'FaceColor', 'interp', 'EdgeColor', 'none');
axis equal tight; colorbar; colormap(jet);
xlabel('X'); ylabel('Y');
title('Von Mises Stress (MPa) — nodal-averaged contour');
hold on;
plot(nodes(fixedNodes,1), nodes(fixedNodes,2), 'k>', 'MarkerFaceColor','k','MarkerSize',3);
plot(nodes(loadNodes,1),  nodes(loadNodes,2),  'kv', 'MarkerFaceColor','y','MarkerSize',3);

%  VISUALIZATION: DISPLACEMENT MAGNITUDE
figure('Name','Displacement Magnitude','Color','w');
patch('Faces', elements, 'Vertices', nodes, ...
      'FaceVertexCData', Umag, 'FaceColor', 'interp', 'EdgeColor', 'none');
axis equal tight; colorbar; colormap(parula);
xlabel('X'); ylabel('Y');
title('Displacement Magnitude (mm)');