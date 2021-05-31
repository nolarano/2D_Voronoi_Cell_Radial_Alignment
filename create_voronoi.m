% Create voronoi data 
 function [Cdat,omega,c2v,vorder,LUT,adj,c2c,pair_index,pair,omega_eqa]= create_voronoi(x,y)
global K A0 width height outside_cells inside_cells perfcirc_innotbound perfcirc_boundary_list rho_g g
%re-create voronoi tessellation 
%Input: (x,y) positions of voronoi centers

% x and y are column vectors
N = length(x);%number of cells
cData = 6;

%replicate the origin points 8 times to surround the original domain;
Nx = x;
Ny = y+height;
Ex = x+width;
Ey = y;
Sx = x;
Sy = y-height;
Wx = x-width;
Wy = y;
NEx = x+width;
NEy = y+height;
SEx = x+width;
SEy = y-height;
SWx = x-width;
SWy = y-height;
NWx = x-width;
NWy = y+height;


%all points
Tx = [x;Nx;NEx;Ex;SEx;Sx;SWx;Wx;NWx];
Ty = [y;Ny;NEy;Ey;SEy;Sy;SWy;Wy;NWy];
TPos = [Tx,Ty];
%voronoi(Tx,Ty)
%axis([0 900 0 900])
[Tv,Tc] = voronoin(TPos);

%Define omega
omega = zeros(2,size(Tv,1));
count_o = 1;
count_c = 1;

%Cut out the original domain
for i = 1:N
    Vcc = Tc{i}';
    %rv_check = (Tv(Vcc,:) <= 900) & (Tv(Vcc,:) >= 0);
    %if sum(rv_check(:,1) & rv_check(:,2)) > 0
    New_Tc{count_c} = count_o:1:count_o+length(Vcc)-1;
    for j = 1:length(Vcc)
        omega(1,count_o) = Tv(Vcc(j),1);
        omega(2,count_o) = Tv(Vcc(j),2);
        
        count_o = count_o+1;
    end
    count_c = count_c+1;
    %end    
end
omega = omega(:,1:count_o-1);
[~,ia,ic] = unique(omega','rows');
omega = omega(:,ia);


%update New_Tc
for i = 1:length(New_Tc)
    New_Tc{i} = ic(New_Tc{i});
end

%Define adj
adj = zeros(size(New_Tc,1),size(New_Tc,1));
for i = 1:length(New_Tc)
    Vcc = New_Tc{i};
    for j = 1:length(Vcc)-1
        adj(Vcc(j),Vcc(j+1)) = 1;
        adj(Vcc(j+1),Vcc(j)) = 1;
    end
    %connect the start and the end point
    adj(Vcc(1),Vcc(end)) = 1;
    adj(Vcc(end),Vcc(1)) = 1;
end
adj = sparse(adj);

%create c2v
c2v = zeros(length(New_Tc),size(omega,2));
for i = 1:length(New_Tc)
    Vcc = New_Tc{i};
    c2v(i,Vcc) = 1;

end

%Create c2c
c2c = zeros(N,N);
for i = 1:N
    verts = New_Tc{i}';
    for j = 1:length(verts)
        c2c(i,find(c2v(:,verts(j)))) = 1;
    end
    c2c(i,i) = 0; %exclude itself
end

%Create vorder
for i = 1:N
    cverts = find(c2v(i,:));
    Z = length(cverts);
    oCverts = zeros(1,Z);
    oCverts(1) = cverts(1);
    for j = 2:Z
        conVerts = find(adj(:,oCverts(j-1)));
        conVerts = conVerts(ismember(conVerts,cverts));
        conVerts = conVerts(~ismember(conVerts,oCverts));
        oCverts(j) = conVerts(1);
    end
    oCverts = [oCverts,oCverts(1)];
    vorder(i).order = oCverts;
    
    G = 0;
    for j = 1:Z
    cG = cross([omega(:,vorder(i).order(j));0],[omega(:,vorder(i).order(j+1));0]);
    G = G+cG(3);
    end
    if G < 0
    vorder(i).order = flip(vorder(i).order);
    end
end



% LUT connection matrix for junctions to vertices: No need for now 
nbonds = sum(adj(:))/2;

LUT = zeros(nbonds,size(adj,1));
ii = 1;
for n = 1:size(adj,1)
    nv = find(adj(n,:));
    nv = nv(nv>n);
    for nn = nv
        LUT(ii,n) = 1; LUT(ii,nn) = -1;
        ii = ii + 1;
    end
end

% Boundary cells
inside = [];
for i = 1:length(inside_cells)
    A = find(c2c(inside_cells(i),:) == 1);
    inside = union(inside,A);
end

outside = [];
for i = 1:length(outside_cells)
    A = find(c2c(outside_cells(i),:) == 1);
    outside = union(outside,A);
end
boundary = intersect(inside, outside);
boundcells = intersect(boundary, inside_cells);

%Build cell data structure.

Cdat = zeros(size(c2v,1),cData);
for c=1:size(Cdat,1)
    ind = find(c2v(c,:));
    Cdat(c,1) = x(c);
    Cdat(c,2) = y(c);
    Cdat(c,3) = polyarea(omega(1,vorder(c).order),omega(2,vorder(c).order));
    Cdat(c,4) = 1*g*A0(c);%0.5*A0(c); %change contractility all cells
    Cdat(c,6) = K*(Cdat(c,3)/Cdat(c,4)-1);
end
    Cdat(boundcells, 4) = rho_g*g*A0(boundcells);%1.2*0.5*A0(c); %change boundary cell contractility
    Cdat(outside_cells,4) = Cdat(outside_cells,3);

%boundary cells and vertices
vert_ext = find(sum(c2v,1)<=2);
b_cells = [];
for ii = 1:length(vert_ext)
b_cells = [b_cells,find(c2v(:,vert_ext(ii)))'];
end
b_cells = unique(b_cells);
Cdat(b_cells',5) = 1;
[pair_index,pair] = boundary_pairing(c2v,Cdat,vorder,vert_ext,omega);


%update c2c
for ii = 1:length(vert_ext)
cells = find(c2v(:,vert_ext(ii)));
for j = 1:length(cells)
    vcon = pair(find(pair_index(:) == vert_ext(ii))).con;
    for i =1:length(vcon)
    ccon = find(c2v(:,vcon(i)));
    c2c(cells(j),ccon) = 1;
    c2c(ccon,cells(j)) = 1;
    end
end
end


m = size(c2v,2); 
omega_eqa = zeros(2,m);
% Find the vertices' positions from the actual equation
% for ii = 1:m    
%         if ismember(ii,pair_index)
%             %t = cputime;
%             pv = find(pair_index(:) == ii);
%             pcells = pair(pv).cells_p;
%             pi = pair(pv).cells_i;
%             for jj = 1:3
%             ri = [pcells(1,jj);pcells(2,jj)];
%             two_cells = setdiff([1,2,3],jj);
%             rj = [pcells(1,two_cells(1));pcells(2,two_cells(1))];
%             rk = [pcells(1,two_cells(2));pcells(2,two_cells(2))];
%             
% 
%             
%             D = 2*((ri(1) - rj(1))*(rj(2) - rk(2)) - (ri(2) - rj(2))*(rj(1) - rk(1)))^2;
%             a = ((rj(1)-rk(1))^2+(rj(2)-rk(2))^2)*[ri(1)-rj(1),ri(2)-rj(2)]*[ri(1)-rk(1);ri(2)-rk(2)]/D;
%             b = ((ri(1)-rk(1))^2+(ri(2)-rk(2))^2)*[rj(1)-ri(1),rj(2)-ri(2)]*[rj(1)-rk(1);rj(2)-rk(2)]/D;
%             c = ((ri(1)-rj(1))^2+(ri(2)-rj(2))^2)*[rk(1)-ri(1),rk(2)-ri(2)]*[rk(1)-rj(1);rk(2)-rj(2)]/D;
%             omega_eqa(1:2,ii) = a*ri+b*rj+c*rk;
% 
%             end
%             
%         else
%         three_cells = find(c2v(:,ii));
%              if length(three_cells) ~= 3
%              error('this vertex connects with more than 3 cells');
%              end
%              
%             %t = cputime; 
%             for jj = 1:3
%                 
%                 ri = [Cdat(three_cells(jj),1);Cdat(three_cells(jj),2)];   
%                 two_cells = setdiff(three_cells,three_cells(jj));
%                 rj = [Cdat(two_cells(1),1);Cdat(two_cells(1),2)];
%                 rk = [Cdat(two_cells(2),1);Cdat(two_cells(2),2)];
%         
%             D = 2*((ri(1) - rj(1))*(rj(2) - rk(2)) - (ri(2) - rj(2))*(rj(1) - rk(1)))^2;
%             a = ((rj(1)-rk(1))^2+(rj(2)-rk(2))^2)*[ri(1)-rj(1),ri(2)-rj(2)]*[ri(1)-rk(1);ri(2)-rk(2)]/D;
%             b = ((ri(1)-rk(1))^2+(ri(2)-rk(2))^2)*[rj(1)-ri(1),rj(2)-ri(2)]*[rj(1)-rk(1);rj(2)-rk(2)]/D;
%             c = ((ri(1)-rj(1))^2+(ri(2)-rj(2))^2)*[rk(1)-ri(1),rk(2)-ri(2)]*[rk(1)-rj(1);rk(2)-rj(2)]/D;
%             omega_eqa(1:2,ii) = a*ri+b*rj+c*rk;
%                 
%             end
%             
%         end        
% end
