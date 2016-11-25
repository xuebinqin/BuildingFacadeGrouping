   function points = doh_Blob(img,num_blobs)
    img = double(img);
       
    if nargin==1
        nb_blobs = 120;
    else
        nb_blobs = num_blobs;
    end
   
    sigma_begin = 2;
    sigma_end  = 15;
    sigma_step = 1;
    sigma_array =sigma_begin:sigma_step:sigma_end;
    sigma_nb   = numel(sigma_array);
       
    [img_height img_width]  = size(img);
       
    snlo =zeros(img_height,img_width,sigma_nb);
    for i=1:sigma_nb
        sigma = sigma_array(i);
        w = 3*sigma;
        [x y]=meshgrid(-w:w,-w:w);
        sigma2 = sigma^2;
        temp =exp(-(x.^2+y.^2)/(2*sigma2))/(2*pi*sigma2);
        mxx = temp .* (x.^2/sigma2-1);
        myy = temp .* (y.^2/sigma2-1);
        mxy = temp .* (x.*y/sigma2);
 
        dxx=imfilter(img,mxx,'replicate');
        dyy=imfilter(img,myy,'replicate');
        dxy=imfilter(img,mxy,'replicate');
        snlo(:,:,i) = dxx.*dyy-dxy.^2;
    end
       
    snlo_dil = imdilate(snlo,ones(3,3,3));
    blob_candidate_index =find(snlo==snlo_dil);
    blob_candidate_value =snlo(blob_candidate_index);
    [~,index] = sort(blob_candidate_value,'descend');
    blob_index = blob_candidate_index(index(1:min(nb_blobs,numel(index))));
    [lig,col,sca] =ind2sub([img_height,img_width,sigma_nb],blob_index);
    points = [lig,col,3*sigma_array(sca)'];
   
end