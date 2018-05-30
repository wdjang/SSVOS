%% k-Ring Graph Construction
function edge_mat = k_ring_graph_construction(sp_img, param_k)

    sp_num = max(sp_img(:));
    [h_size, w_size] = size(sp_img);    
    edge_mat = zeros(sp_num,sp_num);
    
    for y_id = 1:h_size-1
        for x_id = 1:w_size-1
            if sp_img(y_id,x_id) ~= sp_img(y_id,x_id+1)
                edge_mat(sp_img(y_id,x_id),sp_img(y_id,x_id+1)) = 1;
                edge_mat(sp_img(y_id,x_id+1),sp_img(y_id,x_id)) = 1;
            end
            if sp_img(y_id,x_id) ~= sp_img(y_id+1,x_id)
                edge_mat(sp_img(y_id,x_id),sp_img(y_id+1,x_id)) = 1;
                edge_mat(sp_img(y_id+1,x_id),sp_img(y_id,x_id)) = 1;
            end
            if sp_img(y_id,x_id) ~= sp_img(y_id+1,x_id+1)
                edge_mat(sp_img(y_id,x_id),sp_img(y_id+1,x_id+1)) = 1;
                edge_mat(sp_img(y_id+1,x_id+1),sp_img(y_id,x_id)) = 1;
            end
            if sp_img(y_id+1,x_id) ~= sp_img(y_id,x_id+1)
                edge_mat(sp_img(y_id+1,x_id),sp_img(y_id,x_id+1)) = 1;
                edge_mat(sp_img(y_id,x_id+1),sp_img(y_id+1,x_id)) = 1;
            end
        end
    end
    
    % Connect neighbor's neighbors
    iter_mat = edge_mat;
    for iter_id = 1:param_k-1
        iter_mat = iter_mat + iter_mat*edge_mat;
    end
    edge_mat = iter_mat > 0;
%     edge_mat = edge_mat.*(1-eye(size(edge_mat)));
    
    % Connect boundary superpixels
%     spbdry_all = [sp_img(end,:)'; sp_img(1,:)'; sp_img(:,1); sp_img(:,end)];
%     spbdry_list = unique(spbdry_all);
%     for spi_id = 1:length(spbdry_list)
%         for spj_id = 1:length(spbdry_list)
%             if spi_id == spj_id
%                 continue;
%             end
%             conn_edge(spbdry_list(spi_id),spbdry_list(spj_id)) = 1;
%             conn_edge(spbdry_list(spj_id),spbdry_list(spi_id)) = 1;
%         end
%     end
    
end