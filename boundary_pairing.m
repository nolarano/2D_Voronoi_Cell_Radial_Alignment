function [pair_index,pair] =boundary_pairing(c2v,Cdat,vorder,b,rv)
%classify and pair boundary vetices
global width height delete_cells
brv = b;
et = 10^-6;
for i = 1:length(b)
    pair_index(i) = brv(i);%brv(i,1);
    tp1 = find(abs(rv(1,:) - (rv(1,brv(i))-width)) < et  & abs(rv(2,:) - rv(2,brv(i))) < et );
    tp2 = find(abs(rv(1,:) - (rv(1,brv(i))+width)) < et  & abs(rv(2,:) - rv(2,brv(i))) < et );
    tp3 = find(abs(rv(1,:) - (rv(1,brv(i)))) < et  & abs(rv(2,:) - (rv(2,brv(i))-height)) < et );
    tp4 = find(abs(rv(1,:) - (rv(1,brv(i)))) < et  & abs(rv(2,:) - (rv(2,brv(i))+height)) < et );
    tp5 = find(abs(rv(1,:) - (rv(1,brv(i))-width)) < et  & abs(rv(2,:) - (rv(2,brv(i))-height)) < et );
    tp6 = find(abs(rv(1,:) - (rv(1,brv(i))+width)) < et  & abs(rv(2,:) - (rv(2,brv(i))-height)) < et );
    tp7 = find(abs(rv(1,:) - (rv(1,brv(i))-width)) < et  & abs(rv(2,:) - (rv(2,brv(i))+height)) < et );
    tp8 = find(abs(rv(1,:) - (rv(1,brv(i))+width)) < et  & abs(rv(2,:) - (rv(2,brv(i))+height)) < et );
    
    %% adaptive comparision
    if length(tp1)>1
        [min_d,min_i] = min((rv(1,tp1)+width-rv(1,brv(i))).^2+(rv(2,tp1)-rv(2,brv(i))).^2);
        tp1 = tp1(min_i);
    end
    
    if length(tp2)>1
        [min_d,min_i] = min((rv(1,tp2)-width-rv(1,brv(i))).^2+(rv(2,tp2)-rv(2,brv(i))).^2);
        tp2 = tp2(min_i);
    end
    
    if length(tp3)>1
        [min_d,min_i] = min((rv(1,tp3)-rv(1,brv(i))).^2+(rv(2,tp3)+height-rv(2,brv(i))).^2);
        tp3 = tp3(min_i);
    end
    
    if length(tp4)>1
        [min_d,min_i] = min((rv(1,tp4)-rv(1,brv(i))).^2+(rv(2,tp4)-height-rv(2,brv(i))).^2);
        tp4 = tp4(min_i);
    end
    if length(tp5)>1
        [min_d,min_i] = min((rv(1,tp5)+width-rv(1,brv(i))).^2+(rv(2,tp5)+height-rv(2,brv(i))).^2);
        tp5 = tp5(min_i);
    end
    if length(tp6)>1
        [min_d,min_i] = min((rv(1,tp6)-width-rv(1,brv(i))).^2+(rv(2,tp6)+height-rv(2,brv(i))).^2);
        tp6 = tp6(min_i);
    end
    if length(tp7)>1
        [min_d,min_i] = min((rv(1,tp7)+width-rv(1,brv(i))).^2+(rv(2,tp7)-height-rv(2,brv(i))).^2);
        tp7 = tp7(min_i);
    end
    if length(tp8)>1
        [min_d,min_i] = min((rv(1,tp8)-width-rv(1,brv(i))).^2+(rv(2,tp8)-height-rv(2,brv(i))).^2);
        tp8 = tp8(min_i);
    end
    
    %%
    tp = unique([tp1,tp2,tp3,tp4,tp5,tp6,tp7,tp8]);
    
    %Catch a bug
    if isempty(tp)
        for iiii = 1:size(c2v,1)
            plot(rv(1,vorder(iiii).order),rv(2,vorder(iiii).order),'b')
            hold on
            if ismember(iiii,delete_cells)
                fill(rv(1,vorder(iiii).order),rv(2,vorder(iiii).order),'r')
                hold on
            end
        end
        scatter(rv(1,brv(i)),rv(2,brv(i)))
        hold on
        tp1
        tp2
        tp3
        tp4
        tp5
        tp6
        tp7
        tp8
        tp
        scatter(rv(1,tp),rv(2,tp),'y','filled')
        
        error('matching incorrect')
    end
    
    %Create pairing data
    pair(i).con = tp;
    bi_cells = find(c2v(:,b(i)))';
    cell_pair = Cdat(bi_cells,1:2)';
    cell_index = bi_cells;
    for ii = 1:length(tp)
        cells = find(c2v(:,tp(ii)))';
        xytemp = bsxfun(@minus,Cdat(cells,1:2)',(rv(:,tp(ii))-rv(:,b(i))));
        cell_pair = [cell_pair,xytemp];
        cell_index = [cell_index,cells];
    end
    
    %Catch a bug
    if length(cell_index) > 3
        tri = delaunay(cell_pair(1,:),cell_pair(2,:));
        
        for kk = 1:length(tri)
         tri2 = mod(cell_index(tri(kk,:)),size(Cdat,1));
         cell_index2(kk,:)= cell_index2(tri2);
         cell_pair2(2*(kk-1)+1:2*(kk-1)+2,:) = cell_pair(:,inc);
        end
        %inc = randperm(length(cell_index),3);
        %cell_index= cell_index(inc);
        %cell_pair = cell_pair(:,inc);
            
    pair(i).cells_p = cell_pair2;
    pair(i).cells_i = cell_index2;
    end
    
    if size(cell_pair,2) < 3
        tp1
        tp2
        tp3
        tp4
        tp5
        tp6
        tp7
        tp8
        tp
        bi_cells = find(c2v(:,b(i)))'
        cell_pair = Cdat(bi_cells,1:2)'
        cell_index = bi_cells
        figure
        for ii = 1:length(tp)
            scatter(rv(1,tp(ii)),rv(2,tp(ii)),'r')
            hold on
            cells = find(c2v(:,tp(ii)))'
            xytemp = bsxfun(@minus,Cdat(cells,1:2)',(rv(:,tp(ii))-rv(:,b(i))));
            cell_pair = [cell_pair,xytemp];
            cell_index = [cell_index,cells]
        end
        
        for ii = 1:length(vorder)
            plot(rv(1,vorder(ii).order),rv(2,vorder(ii).order),'k')
            hold on
            if ismember(ii,delete_cells)
                fill(rv(1,vorder(ii).order),rv(2,vorder(ii).order),'r')
                hold on
            end
        end
        scatter(rv(1,b(i)),rv(2,b(i)),'filled')
        hold on
        scatter(Cdat(cell_index,1),Cdat(cell_index,2))
        hold off
        
        error('cell_pair')
    end
    
    %Catch a bug
    if length(cell_index) < 3
        error('cell_index')
    end
    
    pair(i).cells_p = cell_pair;
    pair(i).cells_i = cell_index;
end
end

