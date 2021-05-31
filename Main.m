%Main program: simulating cells motion on an Epithelial tissue
%Author: Nonthakorn Olaranont and Sarah St. Pierre (Last update May 31 2020)
%
%Instruction on how to run simulation can be found in :
%https://github.com/nolarano/2D_Voronoi_Cell_Radial_Alignment.git
%
clear all;
clc;
tic
%% parameters
global K K_b Lambda1 Lambda2 Lambda3 Lambda4 A0 width height distance gamma R_outer R_inner inside_cells outside_cells vertex_list vertex_outside B criteria counter ...
    perfcirc_innotbound perfcirc_boundary_list ring2_boundout ring2_boundin lambda_g rho g rho_g
K = 1;
K_b = 0.1; %0.1
Lambda1 = 30; % interfacial tension 20 circle and ring 
Lambda2 = 0; % 15
Lambda3 = 0; % outside
Lambda4 = 0; % innermost vertices: adhesion 15
R_outer = 223; %400 
R_inner = 200; 
width = 900;
height = 900;
distance = 1; 
gamma = 2;
g = 1;
rho = 1;
rho_g = 1;
B = 0; % contractility
criteria = 1e-16; 
N_input = 1; %50
resume = load('IC_625_3');
%resume = load('Circ2.mat');
filename = 'Circ2_3';
Cdat = resume.Cdat_t(end).dat;
A0 = Cdat(:,4);


% boundout = load('ring2layer_boundout');
% ring2_boundout = boundout.ring2layer_boundout;
% boundin  = load('ring2layer_boundin');
% ring2_boundin = boundin.ring2layer_boundin;
% inside_cells = resume.inside_cells;
% outside_cells = resume.outside_cells;
%inside_cells = load('bigcirc_inside_cells.mat');
%inside_cells = inside_cells.bigcirc_inside_cells;
%outside_cells = load('bigcirc_outside_cells.mat');
%outside_cells = outside_cells.bigcirc_outside_cells;
% pc_innotbound = load('perfcirc_innotbound.mat');
% perfcirc_innotbound = pc_innotbound.perfcirc_innotbound;
% pc_boundary_list = load('perfcirc_boundary_list.mat');
% perfcirc_boundary_list = pc_boundary_list.perfcirc_boundary_list;

%% initial configuration

% N = 400; %number of cells (square number only) %100
% A0 = ones(N,1)*(width*height/N);
% 
% boxh = linspace(0,height,sqrt(N)+1);
% boxw = linspace(0,width,sqrt(N)+1);
% box_diffh = boxh(2)-boxh(1);
% box_diffw = boxh(2)-boxh(1);
% 
% x = zeros(N,1);
% y = zeros(N,1);
% for i=1:sqrt(N)
%     for j=1:sqrt(N)
%         x((i-1)*sqrt(N)+j) = rand()*box_diffw+boxw(j);
%         y((i-1)*sqrt(N)+j) = rand()*box_diffh+boxh(i);
%         % A0((i-1)*sqrt(N)+j) = 9900-400*(j-1); % left to right gradient
%         % A0((i-1)*sqrt(N)+j) = 9900-400*(i-1); % bottom to top gradient
%     end
% %     for j = 6:10
% %         A0((i-1)*sqrt(N)+j) = 10900 - 400*(j-1);
% %     end
% %     for j = 1:5
% %         A0((i-1)*sqrt(N)+j) = 7300 + 400*(j-1);
% %     end
% end
% [Cdat,~,c2v,~,~,~,~,~,~,~] = create_voronoi(x,y);
%% Circle
 inside_cells = find(((Cdat(:,1)-450).^2 + (Cdat(:,2)-450).^2) <= R_outer^2);
 outside_cells = find(((Cdat(:,1)-450).^2 + (Cdat(:,2)-450).^2) > R_outer^2);

%% Ring
% inside_cells = find(((Cdat(:,1)-450).^2 + (Cdat(:,2)-450).^2) <= R_outer^2 & ((Cdat(:,1)-450).^2 + (Cdat(:,2)-450).^2) >= R_inner^2);
% outside_cells = find(((Cdat(:,1)-450).^2 + (Cdat(:,2)-450).^2) > R_outer^2 | ((Cdat(:,1)-450).^2 + (Cdat(:,2)-450).^2) < R_inner^2);
%% Proliferation
% [Cdat] = cellDivision(Cdat);
%%
Cdat_t(1).dat = Cdat;
Err = 1;
i = 2;

Er = 1;
%% iteration

for input = 1:N_input
    
    if Er < criteria %Er is error 1
        break;
    end
    
    i
    t = cputime;
    [x,y] = flow(Cdat);
    [Cdat,omega,~,vorder,~,~,~,~,~,~] = create_voronoi(x,y);
    
    %Energy and Position Difference
    % add in v list
    PE = sum(K/2*Cdat(:,4).*(Cdat(:,4)./Cdat(:,3)-1).^2);
    LE= 0;
    
    for ii = 1:length(vorder)
        for jj = 1:length(vorder(ii).order)-1
            if ismember(vorder(ii).order(jj), vertex_list) && ismember(vorder(ii).order(jj+1), vertex_list)
                LE = LE+Lambda1/2*norm(omega(:,vorder(ii).order(jj)) - omega(:,vorder(ii).order(jj+1)));
            elseif ismember(vorder(ii).order(jj), vertex_outside) || ismember(vorder(ii).order(jj+1), vertex_outside)
                LE = LE+Lambda3/2*norm(omega(:,vorder(ii).order(jj)) - omega(:,vorder(ii).order(jj+1)));
            else
                LE = LE+Lambda2/2*norm(omega(:,vorder(ii).order(jj)) - omega(:,vorder(ii).order(jj+1))); 
            end
        end
    end
    
    TE = PE+LE;   
    Err = sum((Cdat_t(i-1).dat(:,1)-Cdat(:,1)).^2+(Cdat_t(i-1).dat(:,2)-Cdat(:,2)).^2)/size(Cdat,1);
    % Ac = width*height/400; % change to N of cells
    %for k = 2:length(Cdat_t)
    % Er = norm([Cdat_t(i-1).dat(:,1);Cdat_t(i-1).dat(:,2)]-[Cdat(:,1);Cdat(:,2)],inf);
        % Er2(k) = norm([Cdat_t(k).dat(:,1);Cdat_t(k).dat(:,2)]-[Cdat_t(k-1).dat(:,1);Cdat_t(k-1).dat(:,2)],2)/sqrt(Ac)/(2*400);
    Er = max(max(abs(Cdat_t(i-1).dat(:,1)-Cdat(:,1))), max(abs(Cdat_t(i-1).dat(:,2)-Cdat(:,2))));

    
    %Store in Data Structures
    PE_data(i-1) = PE;
    LE_data(i-1) = LE;
    TE_data(i-1) = TE;
    Err_data(i-1) = Err;
    Er1_data(i-1) = Er;
    % Er2_data(i-1) = Er2(i);
    
    
    Cdat_t(i).dat = Cdat;
    
    fprintf('Time per iteration %f\nError = %f\nError 1 = %f\n',cputime-t,Err, Er)
    
    save([filename,'.mat'],'Cdat_t','Err_data','Lambda1','Lambda2','K','PE_data','LE_data','TE_data','A0','width','height','inside_cells','outside_cells')
    
    %if ismember(i,[4,6,8,10])
%         criteria = criteria/10;
    %end
        
    i = i+1;
    
    
end

%% Cell Area
% Cdat_inside = Cdat_t(end).dat(inside_cells,:);
% area_inside = Cdat_inside(:,3);
% distance_inside = sqrt((Cdat_inside(:,1)-450).^2 + (Cdat_inside(:,2)-450).^2);
% norm_distance = distance_inside/max(distance_inside);
% matrix = [area_inside norm_distance];
% %plot(norm_distance,area_inside,'.')
% csvwrite('area_distance.txt',matrix)

%% Error Plot
figure
% Er1_data(1) = [];
plot(Er1_data)
saveas(gcf,['Error1',filename,'.fig'])
close all

figure
plot(linspace(1,i-2, i-2), TE_data, 'r')
hold on
plot(linspace(1,i-2, i-2), PE_data, 'b')
plot(linspace(1,i-2, i-2), LE_data, 'g')
ylabel('Energy')
xlabel('iteration #')
legend('TE', 'PE', 'LE')
hold off
saveas(gcf, ['Energy',filename,'.fig'])
close all


%% Video
%Pmax = 2195; %2120 (ring)
%Pmin = 1006; %928.9 (ring); 1006 (circ)
 Pmax = -Inf;
 Pmin = Inf;
for i = 1:length(Cdat_t)

    Cdat = Cdat_t(i).dat;
    x = Cdat(:,1);
    y = Cdat(:,2);
    [~,~,~,vorder,~,~,~,~,~,~] = create_voronoi(x,y);
    mmax = max(Cdat(:,3));
    mmin = min(Cdat(:,3));
    Pmax = max(Pmax,mmax);
    Pmin = min(Pmin,mmin);
%     for ii = 1:length(vorder)
%         measure = Cdat(ii,3);
%         Pmax = max(Pmax,measure);
%         Pmin = min(Pmin,measure);
%     end
    
end

v = VideoWriter(filename,'MPEG-4');
%v = VideoWriter(['CircIC5_1K']);
open(v)
for i = 1:length(Cdat_t)
    Cdat = Cdat_t(i).dat;
    x = Cdat(:,1);
    y = Cdat(:,2);
    [~,omega,~,vorder,~,~,~,~,~,~] = create_voronoi(x,y);
    
    f = figure('Position',[20 20 900 900],'visible','off');
    for ii = 1:length(vorder)
            measure = Cdat(ii,3);
            red_value = (measure-Pmin)/(Pmax-Pmin);
            if ismember(ii,inside_cells)
                fill(omega(1,vorder(ii).order),omega(2,vorder(ii).order),[0 red_value 0])
                plot(omega(1,vorder(ii).order),omega(2,vorder(ii).order),'w')
            else
                fill(omega(1,vorder(ii).order),omega(2,vorder(ii).order),[1 1 1]) %[red_value 0 1-red_value]
                plot(omega(1,vorder(ii).order),omega(2,vorder(ii).order),'w')
            end
            hold on
            
            hold on

    end
    plot([0,width,width,0,0],[0,0,height,height,0],'k');
    hold on
    %map = [zeros(256,1), linspace(0,1,256)',zeros(256,1)];
    %colormap(map)
    %colorbar('Ticks',[])
    axis([-100 width+100 -100 height+100])
    frame = getframe(gcf);
    hold off
    close all
    writeVideo(v,frame);
end
close(v)


%% Print counter
% counter

