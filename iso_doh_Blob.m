   function response = iso_doh_Blob(img,regions)
   %input:regions =[x,y,a,b.theta];
  
    img = double(img);
       
    [img_height img_width]  = size(img);

    snlo = zeros(length(regions),1);
    steps1 = length(regions);
    hwait1 = waitbar(0,'Waiting>>>>>>>>');
    for i=1:length(regions)
        
        px=floor(regions(i,1));
        py=floor(regions(i,2));
        sigmax=floor(regions(i,3)/2);
        sigmay=floor(regions(i,4)/2);
%         theta=regions(i,5);
        
% % %         sigma = sigma_array(i);
        sigma = 0;
        if sigmax>sigmay
            sigma = sigmay;
        else
            sigma = sigmax;
        end
        w = 3*sigma;
        
        [x y]=meshgrid(-w:w,-w:w);
        sigma2 = sigma^2;
        
%         u = x.*cos(theta) + y.*sin(theta);
%         v = y.*cos(theta) - x.*sin(theta);
        
        temp =exp(-(x.^2+y.^2)/(2*sigma2))/(2*pi*sigma2);
        mxx = temp .* (x.^2/sigma2-1);
        myy = temp .* (y.^2/sigma2-1);
        mxy = temp .* (x.*y/sigma2);
        
        %pad image
        
        Img = img;
        if px-w<1
            Img = padarray(Img,[0, w-px+1], 'replicate', 'pre');
            px = px+(w-px+1);
        end
        if py-w<1
            Img = padarray(Img,[w-py+1 0], 'replicate', 'pre');
            py = py+(w-py+1);
        end
        if px+w>img_width
            Img = padarray(Img,[0 px+w-img_width], 'replicate', 'post');
        end
        if py+w>img_height
            Img = padarray(Img,[py+w-img_height 0], 'replicate', 'post');
        end
        
        localImg = Img(uint32(py-w):uint32(py+w), uint32(px-w):uint32(px+w));
        
        dxx = localImg.*mxx;
        dyy = localImg.*myy;
        dxy = localImg.*mxy;
        snlo(i) = sum(sum(dxx.*dyy-dxy.^2));

        str=['Compute iso DoH ',num2str(i/steps1*100),'%'];
        waitbar(i/steps1,hwait1,str);
    end
    close(hwait1); %
    response = snlo;   
   
end