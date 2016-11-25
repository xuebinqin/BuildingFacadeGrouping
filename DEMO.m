close all
clear all
clc
tic

thre_qm = 5;%set the quadralateral threshold
% =====load image======
imgname = 'image1';
imgnameHZ=[imgname,'.jpg'];


I=imread(imgnameHZ);
imsize=[size(I,1) size(I,2)];
Igray=rgb2gray(I);
% Igray = I;
%Igray=uint8(filter2(fspecial('average',5),Igray));%/255;%%mean filtering
imshow(Igray);
hold on 
fprintf('detecting MSER regions...\n');


%=================================A extract MSERS features===============================
regions=detectMSERFeatures(Igray,'ThresholdDelta',3,'RegionAreaRange',[4, 10000],'MaxAreaVariation',0.2);
plot(regions);
title('Results of MSERs');

%remove the regions which have extreme b/a
ab_idx = [];
for i = 1:length(regions)
    ab = regions(i).Axes;
    if ab(1)/ab(2)>5 || ab(1)*ab(2)>10000 || ab(1)/ab(2)<0.2
        ab_idx =[ab_idx, i];
    end
end
if ~isempty(ab_idx)
    regions(ab_idx) = [];
end

clear I
%==================================B MSERs Filtering by Steered-DoH===========================

%--------------Compute Steered DoH respnse-----------
fprintf('computing DoH response...\n');
points = zeros(length(regions),5);
for i = 1:length(regions)
    x = regions(i).Location;
    ax = regions(i).Axes;
    theta = regions(i).Orientation;
    points(i,:) = [x, ax, theta];
end
% clear regions

response = Steered_doh_Blob(Igray,points);
response_iso = iso_doh_Blob(Igray,points);

regions = regions(find(response_iso<0));
response = response(find(response_iso<0));

clear Igray
clear points

%-----------------Filtering by Steered-DoH Response----------------
fprintf('filter MSER by Steered DoH response...\n');
% sort descend regions as area
pxlNum=zeros(length(regions),1);
for k=1:length(regions)
    ab = regions(k).Axes;
    pxlNum(k,1)=ab(1);%*ab(2);
end
% regions sorted by area descend
[spxlNum,index] = sort(pxlNum, 'descend');
clear pxlNum
clear spxlNum
% store PixelList of each regions
PL=regions.PixelList;
save('regions.mat','regions');
clear regions

% initialized overlap image for count overlapped pixels
Iovlp=zeros(imsize);
pidx=[];
% figure(6);
for i=1:length(index)
    x=PL{index(i)};

    % get majority pixel and its number of temp
    %[id num_p] = mode(Iovlp(x));
    uval = unique(Iovlp(min(x(:,2)):max(x(:,2)),min(x(:,1)):max(x(:,1))));
    [m,n] = histc(Iovlp(min(x(:,2)):max(x(:,2)),min(x(:,1)):max(x(:,1))),uval);
    num_vec = sum(m,2);
    [~,pos] = max(num_vec);
    id = uval(pos);
    
    if id == 0% no overlap
          Iovlp(min(x(:,2)):max(x(:,2)),min(x(:,1)):max(x(:,1)))= i;
    else % overlap over half of current regions area
        % compare DoH respose of overlapped two regions
        if response(index(i)) < response(index(id)) % current region has large response
            % remove former regions
            y = PL{index(id)};

            Iovlp(min(y(:,2)):max(y(:,2)),min(y(:,1)):max(y(:,1)))= 0;
            index(id) = 0;
            Iovlp(min(x(:,2)):max(x(:,2)),min(x(:,1)):max(x(:,1)))= i;
        else
            index(i) = 0;
        end
    end
end

clear PL
clear Iovlp

load 'regions'
I=imread(imgnameHZ);

figure(2)
imshow(I);
hold on
index(find(index==0))=[];
regions_fil = regions(index);
response_fil = response(index);
plot(regions_fil);
title('Results of Filtering by Steered-DoH');
clear regions


x = regions_fil.Location;
hold on

clear response
clear response_iso
clear response_fil
clear regions
% clear regions_fil
clear index
clear y

%=================C Delaunay Triangulation induced Feature Perceptural Grouping=================
fprintf('construct TIN...\n');
X = double(x(:,1));
Y = double(x(:,2));
% plot(X,Y,'g+');
% Z = zeros(size(x,1));
% figure(3)
tin=delaunay(X,Y);
% % % % % % % % % % % % % % % % triplot(tin,X,Y);
% trisurf(tin,X,Y,Z)
hold on

%-----------------Search for Adjacent Triangulation-----------------
% find tri related to each piont
pt2tin = cell(length(regions_fil),1);
for t = 1:length(tin)
    if length(find(pt2tin{tin(t,1)}==t)) == 0
        pt2tin{tin(t,1)}(length(pt2tin{tin(t,1)})+1)=t;
    end
    if length(find(pt2tin{tin(t,2)}==t)) == 0
        pt2tin{tin(t,2)}(length(pt2tin{tin(t,2)})+1)=t;
    end
    if length(find(pt2tin{tin(t,3)}==t)) == 0
        pt2tin{tin(t,3)}(length(pt2tin{tin(t,3)})+1)=t;
    end
end
% find tri related to each other
fprintf('Find Parallelogram...\n');
quadra_par = zeros(length(tin),3);%store Parallelogram parameters
quadra_idx = zeros(length(tin),3);%store adjacent tri
tin_Facade = zeros(length(tin),1); % denote whether the triangulation is facade
stepsp = length(tin);
hwaitp = waitbar(0,'Waiting>>>>>>>>');
for k = 1:length(tin)
    temp_tri = [pt2tin{tin(k,1)}(:);pt2tin{tin(k,2)}(:);pt2tin{tin(k,3)}(:)];
    tmp_bin = unique(temp_tri);
    [tri_cnt,tri_id] = hist(temp_tri,tmp_bin);
    tmp_ad = tri_id(find(tri_cnt==2));
    
    for p = 1:length(tmp_ad)
        if length(find(tin(tmp_ad(p),:) == tin(k,3))) == 0
            p4 = setdiff(tin(tmp_ad(p),:),intersect(tin(tmp_ad(p),:),tin(k,:)));
            quadra_idx(k,1) = tmp_ad(p);
            quadra_par(k,1) = quadralateral_measure(x(tin(k,3),:),x(tin(k,1),:),x(tin(k,2),:),x(p4,:));
            if quadra_par(k,1) < thre_qm
                tin_Facade(k) = 1;
                tin_Facade(tmp_ad(p)) = 1;
            end
        end
    end
    
    for p = 1:length(tmp_ad)
        if length(find(tin(tmp_ad(p),:) == tin(k,1))) == 0
            p4 = setdiff(tin(tmp_ad(p),:),intersect(tin(tmp_ad(p),:),tin(k,:)));
            quadra_idx(k,2) = tmp_ad(p);
            quadra_par(k,2) = quadralateral_measure(x(tin(k,1),:),x(tin(k,2),:),x(tin(k,3),:),x(p4,:));
            if quadra_par(k,2) < thre_qm
                tin_Facade(k) = 1;
                tin_Facade(tmp_ad(p)) = 1;
            end
        end
    end
    
    for p = 1:length(tmp_ad)
        if length(find(tin(tmp_ad(p),:) == tin(k,2))) == 0
            p4 = setdiff(tin(tmp_ad(p),:),intersect(tin(tmp_ad(p),:),tin(k,:)));
            quadra_idx(k,3) = tmp_ad(p);
            quadra_par(k,3) = quadralateral_measure(x(tin(k,2),:),x(tin(k,3),:),x(tin(k,1),:),x(p4,:));
            if quadra_par(k,3) < thre_qm
                tin_Facade(k) = 1;
                tin_Facade(tmp_ad(p)) = 1;
            end
        end
    end
    str=['Compute Parallelogram ',num2str(k/stepsp*100),'%'];
    waitbar(k/stepsp,hwaitp,str);
end
close(hwaitp); %



%============================================D Noise Removal===========================================
fprintf('Post-processing...\n');

Facades = {}; %define a cell for storing the facade tin
state_tin = zeros(length(tin),1);
Facad_num = 0;
for i = 1:length(tin)
    if state_tin(i) == 0 && tin_Facade(i) == 1
        temp_facade(1) = i;
        state_tin(i) = 1;
        j=1;
        
        while j <= length(temp_facade) || j == 1
            if quadra_idx(temp_facade(j),1) ~= 0 && tin_Facade(quadra_idx(temp_facade(j),1)) ==1 && state_tin(quadra_idx(temp_facade(j),1)) == 0
                temp_facade(length(temp_facade)+1) = quadra_idx(temp_facade(j),1);
                state_tin(quadra_idx(temp_facade(j),1)) = 1;
            end
            if quadra_idx(temp_facade(j),2) ~= 0 && tin_Facade(quadra_idx(temp_facade(j),2)) ==1 && state_tin(quadra_idx(temp_facade(j),2)) == 0
                temp_facade(length(temp_facade)+1) = quadra_idx(temp_facade(j),2);
                state_tin(quadra_idx(temp_facade(j),2)) = 1;
            end
            if quadra_idx(temp_facade(j),3) ~= 0 && tin_Facade(quadra_idx(temp_facade(j),3)) ==1 && state_tin(quadra_idx(temp_facade(j),3)) == 0
                temp_facade(length(temp_facade)+1) = quadra_idx(temp_facade(j),3);
                state_tin(quadra_idx(temp_facade(j),3)) = 1;
            end
            j = j+1;
        end
%         if length(temp_facade) > 2
            Facades{i} = temp_facade;
%         end
        temp_facade = [];
    end
end

for i = 1:length(Facades)
    if length(Facades{i}) == 2
        for j = 1:length(Facades{i})
            tin_Facade(Facades{i}(j))=0;
        end
    end
end

% dispaly results
figure(3)
imshow(I);
hold on
plot(regions_fil);
hold on
% % % % % % triplot(tin,X,Y);

hold on
triplot(tin(find(tin_Facade==1),:),X,Y,'r','LineWidth',1.5);
title('Results of DT induced grouping');

fprintf('process over...\n');
toc