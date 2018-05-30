function sp_map = superpixel_generation(in_img, frame_id, result_path, sp_param)

SLIC_path = fullfile(result_path, 'data', 'superpixels/');
if ~exist(SLIC_path,'dir')
    mkdir(SLIC_path);
end

bmp_name = fullfile(SLIC_path,sprintf('%04d.bmp',frame_id));

if exist([bmp_name(1:end-4), '.dat'])
    sp_map = ReadDAT([size(in_img,1),size(in_img,2)],...
        [bmp_name(1:end-4), '.dat']);  % Read per-pixel superpixel index
else
    imwrite(in_img,bmp_name);
    slic_com = ['SLICSuperpixelSegmentation', ' ', bmp_name, ' ', int2str(sp_param.edge), ...
        ' ', int2str(sp_param.area), ' ', SLIC_path];  
    system(slic_com);   fprintf('\n');
    slic_img = imread([bmp_name(1:end-4), '_SLIC.bmp']);
    slic_img = flip(slic_img,2);
    imwrite(slic_img, [bmp_name(1:end-4), '_SLIC.png']);
    delete([bmp_name(1:end-4), '_SLIC.bmp']);
    delete(bmp_name);
    sp_map = ReadDAT([size(in_img,1),size(in_img,2)],...
        [bmp_name(1:end-4), '.dat']);  % Read per-pixel superpixel index
end

end


