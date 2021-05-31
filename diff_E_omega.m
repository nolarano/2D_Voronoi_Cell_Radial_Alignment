function [E_omega_dot] = diff_E_omega(Cdat,c2v,vorder,omega_v,pair_index,pair, vertex_list, vertex_outside, vertex_inside, vert_b, boundary_list) %,vert_bin, outsidein, boundout) 

global K Lambda1 Lambda2 Lambda3 Lambda4 B inside_cells K_b counter A0 perfcirc_boundary_list rho

n = size(c2v,1);%number of cells
m = size(c2v,2);%number of vertices

counter = zeros(size(c2v,2),1);

E_omega_dot = zeros(2*m,1);

Lambda_1 = Lambda1;
K_1 = K;

for ii = 1:n % go through all cells
     if ismember(ii, boundary_list) %boundout to change outer boundary cells x*K; inside_cells for ring IC = 4*k
         K_1 = rho*K;%1*K; %change boundary cell stiffness
     else
         K_1 = 1*K;
%     elseif ismember(ii, outsidein)
%         K_1 = 4*K;
%     else
%         K_1 = 4*K;
     end
%      if ismember(ii, boundout)
%          Lambda_1 = 4*Lambda1; % ring case
%      else
%          Lambda_1 = Lambda1;
%      end
    
    verts = vorder(ii).order;
    
    for jj = 1:length(verts)-1 %go through all vertices of cell ii
        
        nu = verts(jj);
        nu_neg = verts(jj+1);
        if jj == 1
            nu_pos = verts(end-1);
        else
            nu_pos = verts(jj-1);
        end
       
        counter(nu) = counter(nu) +1;
        
        A = Cdat(ii,3);
        A_0 = Cdat(ii,4);
        
        if ismember(nu, vertex_list)
            if ismember(nu_pos, vertex_list) && ismember(nu_neg, vertex_list)
                Ex = K_1*(A/A_0-1)*(omega_v(2,nu_neg)-omega_v(2,nu_pos)) + ...
                + Lambda_1/2*((omega_v(1,nu)-omega_v(1,nu_pos))/norm(omega_v(:,nu)-omega_v(:,nu_pos)))-Lambda_1/2*((omega_v(1,nu_neg)-omega_v(1,nu))/norm(omega_v(:,nu_neg)-omega_v(:,nu)));
        
                Ey = K_1*(A/A_0-1)*(omega_v(1,nu_pos)-omega_v(1,nu_neg)) + ...
                + Lambda_1/2*((omega_v(2,nu)-omega_v(2,nu_pos))/norm(omega_v(:,nu)-omega_v(:,nu_pos)))-Lambda_1/2*((omega_v(2,nu_neg)-omega_v(2,nu))/norm(omega_v(:,nu_neg)-omega_v(:,nu)));
            elseif ismember(nu_pos, vertex_list) && ismember(nu_neg, vertex_inside)
                Ex = K_1*(A/A_0-1)*(omega_v(2,nu_neg)-omega_v(2,nu_pos)) + ...
                + Lambda_1/2*((omega_v(1,nu)-omega_v(1,nu_pos))/norm(omega_v(:,nu)-omega_v(:,nu_pos)))-(Lambda2/2+B*norm(omega_v(:,nu_neg)-omega_v(:,nu)))*((omega_v(1,nu_neg)-omega_v(1,nu))/norm(omega_v(:,nu_neg)-omega_v(:,nu)));
        
                Ey = K_1*(A/A_0-1)*(omega_v(1,nu_pos)-omega_v(1,nu_neg))+ ...
                + Lambda_1/2*((omega_v(2,nu)-omega_v(2,nu_pos))/norm(omega_v(:,nu)-omega_v(:,nu_pos)))-(Lambda2/2+B*norm(omega_v(:,nu_neg)-omega_v(:,nu)))*((omega_v(2,nu_neg)-omega_v(2,nu))/norm(omega_v(:,nu_neg)-omega_v(:,nu)));
            elseif ismember(nu_pos, vertex_list) && ismember(nu_neg, vertex_outside)
                Ex = K_1*(A/A_0-1)*(omega_v(2,nu_neg)-omega_v(2,nu_pos)) + ...
                + Lambda_1/2*((omega_v(1,nu)-omega_v(1,nu_pos))/norm(omega_v(:,nu)-omega_v(:,nu_pos)))-Lambda3/2*((omega_v(1,nu_neg)-omega_v(1,nu))/norm(omega_v(:,nu_neg)-omega_v(:,nu)));
        
                Ey = K_1*(A/A_0-1)*(omega_v(1,nu_pos)-omega_v(1,nu_neg)) + ...
                + Lambda_1/2*((omega_v(2,nu)-omega_v(2,nu_pos))/norm(omega_v(:,nu)-omega_v(:,nu_pos)))-Lambda3/2*((omega_v(2,nu_neg)-omega_v(2,nu))/norm(omega_v(:,nu_neg)-omega_v(:,nu)));
            elseif ismember(nu_pos, vertex_outside) && ismember(nu_neg, vertex_list)
                Ex = K_1*(A/A_0-1)*(omega_v(2,nu_neg)-omega_v(2,nu_pos)) + ...
                + Lambda3/2*((omega_v(1,nu)-omega_v(1,nu_pos))/norm(omega_v(:,nu)-omega_v(:,nu_pos)))-Lambda_1/2*((omega_v(1,nu_neg)-omega_v(1,nu))/norm(omega_v(:,nu_neg)-omega_v(:,nu)));
        
                Ey = K_1*(A/A_0-1)*(omega_v(1,nu_pos)-omega_v(1,nu_neg)) + ...
                + Lambda3/2*((omega_v(2,nu)-omega_v(2,nu_pos))/norm(omega_v(:,nu)-omega_v(:,nu_pos)))-Lambda_1/2*((omega_v(2,nu_neg)-omega_v(2,nu))/norm(omega_v(:,nu_neg)-omega_v(:,nu)));
            elseif ismember(nu_pos, vertex_inside) && ismember(nu_neg, vertex_list)
                Ex = K_1*(A/A_0-1)*(omega_v(2,nu_neg)-omega_v(2,nu_pos))+ ...
                + Lambda2/2*((omega_v(1,nu)-omega_v(1,nu_pos))/norm(omega_v(:,nu)-omega_v(:,nu_pos)))-Lambda_1/2*((omega_v(1,nu_neg)-omega_v(1,nu))/norm(omega_v(:,nu_neg)-omega_v(:,nu)));
        
                Ey = K_1*(A/A_0-1)*(omega_v(1,nu_pos)-omega_v(1,nu_neg))+ ...
                + Lambda2/2*((omega_v(2,nu)-omega_v(2,nu_pos))/norm(omega_v(:,nu)-omega_v(:,nu_pos)))-Lambda_1/2*((omega_v(2,nu_neg)-omega_v(2,nu))/norm(omega_v(:,nu_neg)-omega_v(:,nu)));
            end
        elseif ismember(nu, vertex_outside) %outside of pattern
            Ex = K_1*(A/A_0-1)*(omega_v(2,nu_neg)-omega_v(2,nu_pos))...
                + Lambda3/2*((omega_v(1,nu)-omega_v(1,nu_pos))/norm(omega_v(:,nu)-omega_v(:,nu_pos)))-Lambda3/2*((omega_v(1,nu_neg)-omega_v(1,nu))/norm(omega_v(:,nu_neg)-omega_v(:,nu)));
        
            Ey = K_1*(A/A_0-1)*(omega_v(1,nu_pos)-omega_v(1,nu_neg))...
                + Lambda3/2*((omega_v(2,nu)-omega_v(2,nu_pos))/norm(omega_v(:,nu)-omega_v(:,nu_pos)))-Lambda3/2*((omega_v(2,nu_neg)-omega_v(2,nu))/norm(omega_v(:,nu_neg)-omega_v(:,nu)));
        else %nu inside pattern
%             if ismember(nu_neg, vertex_list) % touches boundary
%                 Ex = 0.2*K*(A/A_0-1)*(omega_v(2,nu_neg)-omega_v(2,nu_pos))...
%                 + (Lambda4/2+B*norm(omega_v(:,nu)-omega_v(:,nu_pos)))*((omega_v(1,nu)-omega_v(1,nu_pos))/norm(omega_v(:,nu)-omega_v(:,nu_pos)))-Lambda2/2*((omega_v(1,nu_neg)-omega_v(1,nu))/norm(omega_v(:,nu_neg)-omega_v(:,nu)));
%         
%                 Ey = 0.2*K*(A/A_0-1)*(omega_v(1,nu_pos)-omega_v(1,nu_neg))...
%                 + (Lambda4/2+B*norm(omega_v(:,nu)-omega_v(:,nu_pos)))*((omega_v(2,nu)-omega_v(2,nu_pos))/norm(omega_v(:,nu)-omega_v(:,nu_pos)))-Lambda2/2*((omega_v(2,nu_neg)-omega_v(2,nu))/norm(omega_v(:,nu_neg)-omega_v(:,nu)));
%             elseif ismember(nu_pos, vertex_list) % touches boundary
%                 Ex = 0.2*K*(A/A_0-1)*(omega_v(2,nu_neg)-omega_v(2,nu_pos))...
%                 +(Lambda2/2+B*norm(omega_v(:,nu)-omega_v(:,nu_pos)))*((omega_v(1,nu)-omega_v(1,nu_pos))/norm(omega_v(:,nu)-omega_v(:,nu_pos)))-Lambda4/2*((omega_v(1,nu_neg)-omega_v(1,nu))/norm(omega_v(:,nu_neg)-omega_v(:,nu)));
%         
%                 Ey = 0.2*K*(A/A_0-1)*(omega_v(1,nu_pos)-omega_v(1,nu_neg))...
%                 +(Lambda2/2+B*norm(omega_v(:,nu)-omega_v(:,nu_pos)))*((omega_v(2,nu)-omega_v(2,nu_pos))/norm(omega_v(:,nu)-omega_v(:,nu_pos)))-Lambda4/2*((omega_v(2,nu_neg)-omega_v(2,nu))/norm(omega_v(:,nu_neg)-omega_v(:,nu)));
            % all inside
            Ex = K_1*(A/A_0-1)*(omega_v(2,nu_neg)-omega_v(2,nu_pos))...
            + (Lambda4/2+B*norm(omega_v(:,nu)-omega_v(:,nu_pos)))*((omega_v(1,nu)-omega_v(1,nu_pos))/norm(omega_v(:,nu)-omega_v(:,nu_pos)))-(Lambda4/2+B*norm(omega_v(:,nu_neg)-omega_v(:,nu)))*((omega_v(1,nu_neg)-omega_v(1,nu))/norm(omega_v(:,nu_neg)-omega_v(:,nu)));
        
            Ey = K_1*(A/A_0-1)*(omega_v(1,nu_pos)-omega_v(1,nu_neg))...
            + (Lambda4/2+B*norm(omega_v(:,nu)-omega_v(:,nu_pos)))*((omega_v(2,nu)-omega_v(2,nu_pos))/norm(omega_v(:,nu)-omega_v(:,nu_pos)))-(Lambda4/2+B*norm(omega_v(:,nu_neg)-omega_v(:,nu)))*((omega_v(2,nu_neg)-omega_v(2,nu))/norm(omega_v(:,nu_neg)-omega_v(:,nu)));
            
        end
        
        E_omega_dot(nu) = E_omega_dot(nu) + Ex;
        E_omega_dot(m+nu) = E_omega_dot(m+nu) + Ey;
        
        %boundary vertices
        if ismember(nu,pair_index)
            conV = pair(find(pair_index(:) == nu)).con;
            
            for jjj = 1:length(conV)
                
                E_omega_dot(conV(jjj)) = E_omega_dot(conV(jjj)) + Ex;
                E_omega_dot(m+conV(jjj)) = E_omega_dot(m+conV(jjj)) + Ey;
                
            end
        end
    end
    for i = 1:(length(vert_b)-1) %outer boundary cells
        nu_b = vert_b(i);
        nu_neg_b = vert_b(i+1);
        if i == 1
            nu_pos_b = vert_b(end-1);
        else
            nu_pos_b = vert_b(i-1);
        end 
        A_pattern = sum(Cdat(inside_cells,3));
        A0_pattern = length(inside_cells)*A0(1); %sum(Cdat(inside_cells,4));
        E_omega_dot(nu_b) = E_omega_dot(nu_b) + K_b*(A_pattern/A0_pattern-1)*(omega_v(2,nu_neg_b)-omega_v(2,nu_pos_b));
        E_omega_dot(m+nu_b) = E_omega_dot(m+nu_b) + K_b*(A_pattern/A0_pattern-1)*(omega_v(1,nu_pos_b)-omega_v(1,nu_neg_b));
    end
%     for i = 1:(length(vert_bin)-1) %inner boundary cells
%         nu_b = vert_bin(i);
%         nu_neg_b = vert_bin(i+1);
%         if i == 1
%             nu_pos_b = vert_bin(end-1);
%         else
%             nu_pos_b = vert_bin(i-1);
%         end 
%         E_omega_dot(nu_b) = E_omega_dot(nu_b) + K_b*(A_pattern/A0_pattern-1)*(omega_v(2,nu_neg_b)-omega_v(2,nu_pos_b));
%         E_omega_dot(m+nu_b) = E_omega_dot(m+nu_b) + K_b*(A_pattern/A0_pattern-1)*(omega_v(1,nu_pos_b)-omega_v(1,nu_neg_b));
%     end
end

% % % test projection of Kb energy
% xhalf = E_omega_dot(1:m);
% yhalf = E_omega_dot((m+1):2*m);
% xvert = find(xhalf ~= 0);
% yvert = find(yhalf ~= 0);
% cell = [];
% proj = zeros(1,length(xvert));
% 
% for i = 1:length(xvert)
%     cel = find(c2v(:,xvert(i))==1);
%     cin = ismember(cel,boundary_list);
%     cell = cel(cin);
%     proj(i) = (Cdat(cell(1),1)-450)*E_omega_dot(xvert(i)) + (Cdat(cell(1),2)-450)*E_omega_dot(xvert(i)+m);
end

%end
