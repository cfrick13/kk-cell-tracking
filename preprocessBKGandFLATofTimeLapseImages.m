
%% fullbkg
function preprocessBKGandFLATofTimeLapseImages(A,B)
% A = 'D:\Frick\';
% B = '2015_12_15 smad3g smFISH';
C = B;
E={B}; %beginning of file name for saving
% channelstoinput = {'_mKate','_EGFP','_CFP','DIC'};
channelstoinput = {'_mKate','_EGFP','_CFP','DIC','_Hoechst'};
% channelstoinput = {'mKate','_EGFP','_CFP','_DIC'};
channelinputs = '(';
for i=1:length(channelstoinput) % creates a string of from '(c1|c2|c3|c4)' for regexp functions
    if i ==1
    channelinputs = strcat(channelinputs,channelstoinput{i});
    elseif i < length(channelstoinput)
        channelinputs = strcat(channelinputs,'|',channelstoinput{i});
    else
        channelinputs = strcat(channelinputs,'|',channelstoinput{i},')');
    end
end

BACKGROUND = [56:60];

DIRone = strcat(A,B,'\',C);
cd(DIRone);
dirlist = dir('*.tif');
[~,~,~,ScenesListed] = regexp([dirlist.name],'s[0-9]+');
sceneList = unique(ScenesListed);
[~,~,~,timeFrameList] = regexp([dirlist.name],'t[0-9]+');
timeList = unique(timeFrameList);
[~,~,~,channelsListed] = regexp([dirlist.name],channelinputs);
channelList = unique(channelsListed);

cd .. 
mkdir(strcat(A,B,'\flatfield_corrected'));
cd (strcat(A,B,'\flatfield_corrected'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%make folders for each scene
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for scenedir = sceneList
    scene = scenedir{1};
%     sceneN = str2num(scene(2:end));
%     scenenum = scenestring(sceneN); 
scenefile{1} = strcat(E{1},'_scene_',scene);
mkdir(scenefile{1});
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for DIRflats = {DIRone} %choose directory (this is a remant of splitz)
    DIRflat = char(DIRflats);
  
    cd (DIRflat)

for channel = channelList %cycle through one channel at a time
filelist = dir(strcat('*',char(channel),'*.tif'));
cfile = {filelist.name};

% for j = 1:length(cfile) %load image and determine scene and timepoint and channel

parfor j = 1:length(cfile) %load image and determine scene and timepoint and channel
    filename = char(cfile(j));
    img = imread(filename);
    
    [a,b] = regexp(filename,'s[0-9]+'); %determine scene
    scene = filename(a:b);
    
    [c,d] = regexp(filename,'t[0-9]+'); %determine timepoint
    if ~isempty(c)
        tpoint = filename(c:d);
        tpoints = tpoint;
    else
        tpoint='t00';
        tpoints = '*';
    end
    
    [e,f] = regexp(filename,channelinputs); %determine channel
    chan = filename(e:f);
    
    bkimg = medianBKG(BACKGROUND,tpoints,chan); %load background images and compile into one median image
    bkShape = reshape(bkimg,[1 size(bkimg,2).*size(bkimg,1)]);
    bkSort = sort(bkShape(~isnan(bkShape)));
    normalizedBkgimg = (double(bkimg)./mean([bkSort(round(length(bkSort).*0.9999)) bkSort(round(length(bkSort).*0.99999))]));
    
    flat = uint16(double(img)./normalizedBkgimg); %flatten the experimental image with the background image
    savename = strcat(E,'_flat_',scene,'_',tpoint,'_',chan,'.tif');
    SAVdir = strcat(A,B,'\flatfield_corrected\',E{1},'_scene_',scene);
    savethatimage(savename,SAVdir,flat,DIRflat,j);
end
end
end


% lineageyo(A,B);
% cd('D:\Users\zeiss\Documents\MATLAB')
% % bleachcorrection(A,B)
% bleachcorrectionNew(A,B)
% cd('D:\Users\zeiss\Documents\MATLAB')
% toughsegmentationforstacks(A,B)
end


function bkimg = medianBKG(BACKGROUND,tpoints,chan)
bkimgs = zeros(512,512,length(BACKGROUND));
i = 1;
for scenenumber = BACKGROUND
scene = scenestring(scenenumber);
% files = strcat('*',scene,'*',tpoints,'*',chan,'*');
files = strcat('*',chan,'*',scene,'*',tpoints,'*');
filelist = dir(files); cfile = {filelist.name}; filename = char(cfile);
bkk = imread(filename);
 
% %  if ~strcmp(chan,'c3')
% if ~strcmp(chan,'DIC')
%      bkklog = segmentationBKGsecond(bkk,[],scenenumber,[],[]);
%      bkkg = regionfill(bkk,bkklog);
%  else
%      bkkg=bkk;
% end
%  
bkkg=bkk;
 
 bkimgs(:,:,i) = bkkg;
i=i+1;
end
bkimg = median(bkimgs,3);
end


function If = segmentationBKGsecond(FinalImage,subdirname,scenename,filename,channel)
fig=1;
% mkdir(strcat('c8_flat'));
% mkdir(strcat('c9_flat'));
% parameters
dimdiff = 2048./size(FinalImage(:,:,1),1);

zerostrel = 2;
% firststrel = round(30./dimdiff);
firststrel = round(20./dimdiff);
sigmafirst = firststrel.*3;
kernelgsizefirst = firststrel.*6;
fracsmoothing = 0.5.*dimdiff;

%set the time


% start
img = FinalImage;
imgorig = img;

img = wiener2(img,[5 5]);
se =strel('disk',zerostrel);
Ie = imerode(imgorig,se);
Iobr = imreconstruct(Ie,img);
Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
IobrcbrF = imcomplement(Iobrcbr);
gaus = double(IobrcbrF);

se =strel('disk',firststrel);
Ie = imerode(gaus,se);
Iobr = imreconstruct(Ie,gaus);
Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
gaus = double(Iobrcbr);


sigma = sigmafirst;
kernelgsize = kernelgsizefirst;
gaustwo = gaussianBlurz(double(gaus),sigma,kernelgsize);

sub = double(gaus) -double(gaustwo);%%%%%%% key step!
b = find(sub == min(min(sub)));
rattio = gaustwo(b)./gaus(b);
gaustwocorr = gaustwo./rattio;
sub_scale_corr = double(gaus) - double(gaustwocorr);




subtractionref = sub_scale_corr;
vec = reshape(subtractionref,size(subtractionref,1)^2,1);
[numbers,bincenters] = hist(double(vec),0:fracsmoothing/10:10000);
numbers = medfilt1(numbers, 10); %smooths curve
fraction = numbers./sum(numbers);

mf = max(fraction);
%%%%%%%%%%%%%%%%%%%%
left=0.5*mf;
slopedown=0.4*mf;
%%%%%%%%%%%%%%%%%%%%%

leftedge = find(fraction > left,1,'first');
insideslopedown = find(fraction(leftedge:end) < slopedown,1,'first');
threshlocation = bincenters(leftedge+insideslopedown);

% figure(22)
%     bar(bincenters,fraction./mf);hold on
%     xlim([-500 1000])
%     ylim([0 0.8])
% stem ([threshlocation threshlocation],[0 10],'g');hold off
% drawnow

subtractionthreshold = threshlocation;

if size(subtractionthreshold,1)==size(subtractionthreshold,2)
else
     subtractionthreshold = mean(threshlocation);
end

% subtractionthreshold = graythresh(subtractionref);

subtracted = sub_scale_corr-subtractionthreshold;
subzero = (subtracted<0);
subtractedzero = subtracted.*(~subzero);


Ie = subtractedzero;
a = find(Ie>0);
submax = zeros(size(Ie));
Ie(a)=50;





Ih = Ie>0;
Igclose = imclose(Ih,strel('disk',round(30./dimdiff)));
Igclosemax = imclose(Ih,strel('disk',round(80./dimdiff)));
Igcopenmax = imopen(Igclosemax,strel('disk',round(10./dimdiff)));
Igcopen = imopen(Igclose,strel('disk',2));
Igcofill = imfill(Igcopen,'holes');
Igcfopen = bwareaopen(Igcofill,round(5000./(dimdiff.^2)));
Igcfopendil = imerode(Igcfopen,strel('disk',round(5)));

% finalsigma=20;
% finalkernelgsize=40;
% gaus = gaussianBlurz(IobrcbrF,round(sigmafirst./dimdiff),(kernelgsizefirst./dimdiff));
% 
% % sigma=40;
% % kernelgsize=80;
% % gaus = gaussianBlurz(gaus,sigmafirst.*2,kernelgsizefirst.*2);
% 
% imgt = -double(gaus);
% % % imgt(~(Igcfopen>0)) = -Inf;
% % imgt(~(Igcopenmax>0)) = -Inf;
% imgt(~(Igcopenmax>0)) = 0;
% 
% % L=watershed(imgt);
% L=imgt;
L = Igcopenmax;

% L(Igcfopendil<1) = 0;
% imagesc(L)
% colormap parula
If = L>0;
% se = strel('disk',16);
se = strel('disk',12);
ifd = imdilate(If,se);
ifdc = imclose(ifd,se);
    ifdc(1:512,1)=255;
    ifdc(1:512,512)=255;
    ifdc(1,1:512)=255;
    ifdc(512,1:512)=255;
    

If = ifdc;
prim = bwperim(If);
se = strel('disk',2);
primd = imdilate(prim,se);
imgorig(primd) = max(max(imgorig));
% figure(99),imagesc(imgorig)
% drawnow
stophere=1;



end


function scene = scenestring(scenenumber)
    scene = char('s00');
        reg = num2str(scenenumber);
        c = size(reg,2);
        for b = 1:c
            scene(3+1-b) = reg(c+1-b);
        end
end

function tpoint = timestring(time)
    tpoint = char('t00');
        reg = num2str(time);
        c = size(reg,2);
        for b = 1:c
            tpoint(3+1-b) = reg(c+1-b);
        end
end
 
function savethatimage(savename,SAVdir,flat,DIRone,j)
disp(strcat(savename,'...',num2str(j)));
cd (SAVdir);
imwrite(flat,char(savename),'tiff');
cd (DIRone);
end



function bw = gaussianBlurz(im,sigma,kernelgsize,varargin)

filtersize = [kernelgsize kernelgsize];
kernelg = fspecial('gaussian',filtersize,sigma);

%% image filtering
gFrame = imfilter(im,kernelg,'repl');

if ~isempty(varargin)
    bw=gFrame.*uint16(varargin{1}>0);
else
    bw=gFrame;
end
end


function bw = logMasked(im,ksize,varargin)
% Discrete Laplacian
kernel = chooseKernel(ksize);
%% image filtering
lapFrame = imfilter(im,kernel,'repl');
if ~isempty(varargin)
    bw=lapFrame.*uint16(varargin{1}>0);
else
    bw=lapFrame;
end
end



function kernel = chooseKernel(ksize)
if ksize ==5
kernel = [-4 -1  0 -1 -4;...
     -1  2  3  2 -1;...
      0  3  4  3  0;...
     -1  2  3  2 -1;...
     -4 -1  0 -1 -4];


% % % -4 -1  0 -1 -4
% % % -1  2  3  2 -1
% % % 0  3  4  3  0
% % % -1  2  3  2 -1
% % % -4 -1  0 -1 -4

elseif ksize == 7
kernel =[-10 -5 -2 -1 -2 -5 -10;... 
    -5  0  3  4  3  0  -5;... 
    -2  3  6  7  6  3  -2;... 
    -1  4  7  8  7  4  -1;... 
    -2  3  6  7  6  3  -2;... 
    -5  0  3  4  3  0  -5;... 
    -10 -5 -2 -1 -2 -5 -10];... 
    
% % % -10 -5 -2 -1 -2 -5 -10 
% % % -5  0  3  4  3  0  -5 
% % % -2  3  6  7  6  3  -2 
% % % -1  4  7  8  7  4  -1 
% % % -2  3  6  7  6  3  -2 
% % % -5  0  3  4  3  0  -5
% % % -10 -5 -2 -1 -2 -5 -10
end
end



function  LoGstack = LaplacianOfGaussianStack(imgstack,dims,ksize)
        LoGstack = zeros(dims(1),dims(2),dims(3));
        for i = 1:size(imgstack,3)
        LoGstack(:,:,i) = logMasked(imgstack(:,:,i),ksize);
        end
        
end


