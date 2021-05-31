
function [X_dot] = solver(t,X) 

global R_outer R_inner inside_cells outside_cells vertex_list vertex_outside vertex_inside K_b perfcirc_boundary_list perfcirc_innotbound ...
    ring2_boundin ring2_boundout

% Equations to be solved by ODEs functions
% Input: positions of voronoi centers

N = length(X)/2; % number of cells
cri = reshape(X',N,2)'; %cells' centroid positions in the usual format 2 by N
x = cri(1,:)'; %x coordinate of centers
y = cri(2,:)'; %y coordinate of centers

[Cdat,omega,c2v,vorder,~,~,c2c,pair_index,pair,~] = create_voronoi(x,y); %create a new voronoi

[J] = J_function(Cdat,c2v,pair_index,pair);

inside_vertex = [];
for i = 1:length(inside_cells)
    A = find(c2v(inside_cells(i),:) == 1);
    inside_vertex = union(inside_vertex,A);
end
outside_vertex = [];
for i = 1:length(outside_cells)
    A = find(c2v(outside_cells(i),:) == 1);
    outside_vertex = union(outside_vertex,A);
end
vertex_list = intersect(inside_vertex, outside_vertex);
vertex_outside = outside_vertex(~ismember(outside_vertex, vertex_list));
vertex_inside = inside_vertex(~ismember(inside_vertex, vertex_list));

b_cells = [];
for i = 1:length(c2v(:,1))
    for j = 1:length(pair_index)
        if ismember(pair_index(j),find(c2v(i,:) ==1))
            b_cells = union(b_cells, i);
        else
            continue
        end
    end
end
%Find pattern boundary cells
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
boundary_list = intersect(boundary, inside_cells);

% Ring
% boundout = [];
% boundin = [];
% for i = 1:length(boundary_list)
%     if sqrt(((Cdat(boundary_list(i),1)-450)^2 + (Cdat(boundary_list(i),2)-450)^2)) <= ((R_outer-R_inner)/2 + R_inner)
%         boundin = [boundin boundary_list(i)];
%     else
%         boundout = [boundout boundary_list(i)];
%     end
% end

% boundary_list = boundout;  % just for ring IC case rename, comment back for k/A0 simulations


% outin = []; %inner outside cells
% for i = 1:length(outside_cells)
%     if sqrt(((Cdat(outside_cells(i),1)-450)^2 + (Cdat(outside_cells(i),2)-450)^2)) <= ((R_outer-R_inner)/2 + R_inner)
%         outin = [outin outside_cells(i)];
%     end
% end
% 
% outsidein = [];
% for i = 1:length(outside_cells)
%     if sqrt(((Cdat(outside_cells(i),1)-450)^2 + (Cdat(outside_cells(i),2)-450)^2)) <= ((R_outer-R_inner)/2 + R_inner)
%         outsidein = [outsidein outside_cells(i)];
%     end
% end
% Find inside not boundary cells and update lists
% innotbound = inside_cells;
% for i = 1:length(boundary_list)
%     a = find(innotbound == boundary_list(i));
%     innotbound(a) = [];
% end
% perfcirc_boundary_list = boundary_list;
% perfcirc_innotbound = innotbound;

%% Order outside boundary cells
coord = [Cdat(boundary_list,1) Cdat(boundary_list,2)];
A = find(coord(:,2) > 450); 
B = find((coord(A,2)) == min(coord(A,2)));
start = coord(A(B),:);
angle_list = acos(((coord(:,1)-450).*(start(1)-450) + (coord(:,2)-450).*(start(2)-450))./(sqrt((coord(:,1)-450).^2 + (coord(:,2)-450).^2)*sqrt((start(1)-450).^2 + (start(2)-450).^2)));
cell_b = zeros(length(angle_list),1);
for i = 2:length(A)
    cell_b(1) = boundary_list(A(B));
    C = sort(angle_list(A));
    cell_b(i) = boundary_list(find(angle_list == C(i)));
end
A2 = find(coord(:,2) < 450);
for i = 1:length(A2)
    C = sort(angle_list(A2),'descend');
    cell_b(i+length(A)) = boundary_list(find(angle_list == C(i)));
end
if Cdat(cell_b(1),1)<450
    cell_b = cell_b';
    cell_b = fliplr(cell_b);
    cell_b = cell_b';
end
% Check cell order
checkorder = zeros(length(cell_b),1);
for i = 1:length(cell_b)
    if i == length(cell_b)
        cell_b(i + 1) = cell_b(1);
    end
    for j = 1:length(find(c2c(cell_b(i),:) ==1))
        cells = find(c2c(cell_b(i),:) == 1);
        if ismember(cells(j),cell_b(i+1))
            checkorder(i) = 1; % all ones indicates correct cell order
        end
    end
end
cell_b(end) = []; % don't want repeated cell used to check order
% Order boundary vertices and delete repetitions
vert_b = [];
for i = 1:length(cell_b)
    vert_all = vorder(cell_b(i)).order;
    vert_all(end) = []; % delete repeated vertex in vorder
    A = ismember(vert_all, vertex_list);
    vert_b = [vert_b vert_all(A)];
end
for i = 1:length(vert_b)
    if i > length(vert_b)
        break
    end
    A = find(vert_b == vert_b(i));
    if length(A) > 1
        for j = 2
            vert_b(A(j)) = [];
            B = find(vert_b == vert_b(i));
            if length(B) > 1
                for k = 2
                    vert_b(B(k)) = [];
                    C = find(vert_b == vert_b(i));
                    if length(C) > 1
                        for l = 2
                            vert_b(C(l)) = [];
                        end
                    end
                end
            end
        end
    end
end
vert_b = [vert_b vert_b(1)];

%% Order inside boundary cells
% coord = [Cdat(boundin,1) Cdat(boundin,2)];
% A = find(coord(:,2) > 450); 
% B = find((coord(A,2)) == min(coord(A,2)));
% start = coord(A(B),:);
% angle_list = acos(((coord(:,1)-450).*(start(1)-450) + (coord(:,2)-450).*(start(2)-450))./(sqrt((coord(:,1)-450).^2 + (coord(:,2)-450).^2)*sqrt((start(1)-450).^2 + (start(2)-450).^2)));
% cell_b = zeros(length(angle_list),1);
% for i = 2:length(A)
%     cell_b(1) = boundin(A(B));
%     C = sort(angle_list(A));
%     cell_b(i) = boundin(find(angle_list == C(i)));
% end
% A2 = find(coord(:,2) < 450);
% for i = 1:length(A2)
%     C = sort(angle_list(A2),'descend');
%     cell_b(i+length(A)) = boundin(find(angle_list == C(i)));
% end
% if Cdat(cell_b(1),1)<450
%     cell_b = cell_b';
%     cell_b = fliplr(cell_b);
%     cell_b = cell_b';
% end
% % Check cell order
% checkorder = zeros(length(cell_b),1);
% for i = 1:length(cell_b)
%     if i == length(cell_b)
%         cell_b(i + 1) = cell_b(1);
%     end
%     for j = 1:length(find(c2c(cell_b(i),:) ==1))
%         cells = find(c2c(cell_b(i),:) == 1);
%         if ismember(cells(j),cell_b(i+1))
%             checkorder(i) = 1; % all ones indicates correct cell order
%         end
%     end
% end
% cell_b(end) = []; % don't want repeated cell used to check order
% % Order boundary vertices and delete repetitions
% vert_bin = [];
% for i = 1:length(cell_b)
%     vert_all = vorder(cell_b(i)).order;
%     vert_all(end) = []; % delete repeated vertex in vorder
%     A = ismember(vert_all, vertex_list);
%     vert_bin = [vert_bin vert_all(A)];
% end
% for i = 1:length(vert_bin)
%     if i > length(vert_bin)
%         break
%     end
%     A = find(vert_bin == vert_bin(i));
%     if length(A) > 1
%         for j = 2
%             vert_bin(A(j)) = [];
%             B = find(vert_bin == vert_bin(i));
%             if length(B) > 1
%                 for k = 2
%                     vert_bin(B(k)) = [];
%                     C = find(vert_bin == vert_bin(i));
%                     if length(C) > 1
%                         for l = 2
%                             vert_bin(C(l)) = [];
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% vert_bin = [vert_bin vert_bin(1)];

[E] = diff_E_omega(Cdat,c2v,vorder,omega,pair_index,pair, vertex_list, vertex_outside, vertex_inside, vert_b, boundary_list);%boundary_list); %, vert_bin, outsidein, boundout);

% R = sqrt((Cdat(:,1)-450).^2 + (Cdat(:,2)-450).^2);
% 
% for i = 1:length(R)
%     if ~ismember(i, boundary_list)
%         R(i) = R_outer;
%     end
% end
% Rx = zeros(400,1);
% Ry = zeros(400,1);
% for i = 1:length(R)
%     Rx(i) = (Cdat(i,1)-450)/sqrt((Cdat(i,1)-450).^2 + (Cdat(i,2)-450).^2);
%     Ry(i) = (Cdat(i,2)-450)/sqrt((Cdat(i,1)-450).^2 + (Cdat(i,2)-450).^2);
% end

X_dot = -J*E; %+ K_b*[(R_outer-R).*Rx; (R_outer-R).*Ry];
% X_dot(outin) = 0;
% X_dot(outin + N) = 0;
% X_dot(boundary_list) = 0;
% X_dot(boundary_list + N) = 0;
% X_dot(outside_cells) = 0;
% X_dot(outside_cells + N) = 0;


