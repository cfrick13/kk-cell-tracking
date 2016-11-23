function  SegmentationOfTimeLapseImages(A,B)

choosefiles(A,B)
end

function choosefiles(A,B)
global nucleus_seg 
cd (strcat(A,B,'\flatfield_corrected'));

nucleus_seg = '_Hoechst_flat';
% nucleus_seg = 'CFP';



primarylist = dir('*_s*');
% primarylist = dir('*_s02');




subd = {primarylist.name};
% parfor i=1:length(subd)
for i=1:length(subd)
    subdir = subd{i};
    subdirname = char(subdir);
        sceneinfo = regexp(subdirname,'s[0-9]+');
        scenename = subdirname(sceneinfo:sceneinfo+2);
        cd(subdirname)
    finaldir = dir('*tiffs*');
        finaldirname = char({finaldir.name});
        cd(finaldirname)
    % file = dir('*EGFP_flat.tif*');

cd .. 
dirlist = dir(nucleus_seg);
if isempty(dirlist)
    dirlist = dir('mKate_flat');
    foldername = 'mKate_flat';
else
    foldername = nucleus_seg;
end
cd(finaldirname)
    
    file = dir(strcat('*',foldername,'*tif'));
    filename = char({file.name});
    FinalImage = loadStack(filename);
    cd ..
    cd (foldername)
    filelist = dir('*.tif');
    fname = filelist(1).name;
    cd ..
    disp(scenename)
    segmentationNucleus(FinalImage,subdirname,scenename,fname,'NucleusBinary_flat');

dirlist = dir('mKate_flat');
if isempty(dirlist)
    dirlist = dir(nucleus_seg);
    foldername = nucleus_seg;
else
    foldername = 'mKate_flat';
end   
cd(finaldirname)
    file = dir(strcat('*',foldername,'*tif'));
    filename = char({file.name});
    FinalImage = loadStack(filename);
        cd ..
        cd (foldername)
        filelist = dir('*.tif');
        fname = filelist(1).name;
        cd ..
%     segmentationRFP(FinalImage,subdirname,scenename,fname,'mKatebinary_flat');
    

dirlist = dir('_EGFP_flat');
if isempty(dirlist)
    dirlist = dir('mKate_flat');
    foldername = 'mKate_flat';
else
    foldername = '_EGFP_flat';
end
cd(finaldirname)
    
    file = dir(strcat('*',foldername,'*tif'));
    filename = char({file.name});
    FinalImage = loadStack(filename);
    cd ..
    cd (foldername)
    filelist = dir('*.tif');
    fname = filelist(1).name;
    cd ..
    segmentationMNG(FinalImage,subdirname,scenename,fname,'EGFPbinary_flat');
    stophere=1;
    segmentationBKGsecond(FinalImage,subdirname,scenename,fname,'BKGbinary_flat');
    cd ..
end

end

function FinalImage=loadStack(FileTif)
% [a,b] = uigetfile;
% FileTif = a;
% cd (b)
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
FinalImage=zeros(nImage,mImage,NumberImages,'uint16');
 
TifLink = Tiff(FileTif, 'r');
for i=1:NumberImages
   TifLink.setDirectory(i);
   FinalImage(:,:,i)=TifLink.read();
end
TifLink.close();
end







function If = segmentationNucleus(FinalImage,subdirname,scenename,filename,channel)
global nucleus_seg
fig=1;
mkdir(strcat('NucleusBinary_flat'));
mkdir(strcat('c5_flat'));

% parameters
dimdiff = 2048./size(FinalImage(:,:,1),1);
zerostrel = 5;
firststrel = round(50./(dimdiff.^2));
sigmafirst = firststrel.*5;
kernelgsizefirst = firststrel.*10;
% fracsmoothing = 0.5.*dimdiff;
fracsmoothing = 0.5;
weiner2p = 20;

dirlist = dir(nucleus_seg);
if isempty(dirlist)
    dirlist = dir('mKate_flat');
    foldername = 'mKate_flat';
else
    foldername = nucleus_seg;
end

tsn = determineTimeFrame(foldername);

% start
for frames = 1:size(FinalImage,3)
img = FinalImage(:,:,frames); 
imgorig = img;


img = wiener2(img,[weiner2p weiner2p]);
se =strel('disk',zerostrel);
Ie = imerode(imgorig,se);
Iobr = imreconstruct(Ie,img);
Iobrone = Iobr;
Iobrd = imdilate(Iobr,se);
Iobrdone = Iobrd;
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
IobrcbrF = imcomplement(Iobrcbr);
gaus = double(IobrcbrF);

% figure(99)
% subplot(3,3,1);
% imagesc(Iobr);
% subplot(3,3,2);
% imagesc(Iobrd);
% subplot(3,3,3);
% imagesc(Iobrcbr);
% subplot(3,3,4);
% imagesc(Iobrcbr);
% subplot(3,3,5);
% imagesc(IobrcbrF);
% subplot(3,3,6);
% imagesc(gaus);

se =strel('disk',firststrel);
Ie = imerode(gaus,se);
Iobr = imreconstruct(Ie,gaus);
Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
gaus = double(Iobrcbr);
% 
% figure(999)
% subplot(3,3,1);
% imagesc(Iobr);
% subplot(3,3,2);
% imagesc(Iobrd);
% subplot(3,3,3);
% imagesc(Iobrcbr);
% subplot(3,3,4);
% imagesc(Iobrcbr);
% subplot(3,3,5);
% imagesc(IobrcbrF);
% subplot(3,3,6);
% imagesc(gaus);
% subplot(3,3,7)
% imagesc(imgorig)


sigma = sigmafirst;
kernelgsize = kernelgsizefirst;
gaustwo = gaussianBlurz(double(gaus),sigma,kernelgsize);

sub = double(gaus) -double(gaustwo);%%%%%%% key step!
b = find(sub == min(min(sub)),1,'first');
rattio = gaustwo(b)./gaus(b);
gaustwocorr = gaustwo./rattio;
sub_scale_corr = double(gaus) - double(gaustwocorr);




subtractionref = sub_scale_corr;
vec = reshape(subtractionref,size(subtractionref,1)^2,1);
[numbers,bincenters] = hist(double(vec),0:fracsmoothing:10000);



numbers = medfilt1(numbers, 10); %smooths curve
fraction = numbers./sum(numbers);

mf = max(fraction);

%%%%%%%%%%%%%%%%%%%%
left=0.5*mf;
slopedown=0.4*mf;
%%%%%%%%%%%%%%%%%%%%%

leftedge = find(fraction > left,1,'first');
insideslopedown = find(fraction(leftedge:end) < slopedown,1,'first');
insideslopeup = find(fraction(leftedge+insideslopedown:end) >0.0012,1,'first');
trough = min(fraction(leftedge+insideslopedown:leftedge+insideslopedown+insideslopeup));
troughindex = find(fraction(leftedge+insideslopedown:leftedge+insideslopedown+insideslopeup) == trough);
troughindexrounded = round(median(troughindex));



threshlocation = bincenters(leftedge+insideslopedown);
% 
% figure(22)
%     bar(bincenters,fraction);hold on
%     xlim([-500 1000])
%     ylim([0 0.1])
% stem ([threshlocation threshlocation],[0 1]);hold off
% drawnow

subtractionthreshold = threshlocation;

if size(subtractionthreshold,1)==size(subtractionthreshold,2)
else
     subtractionthreshold = mean(threshlocation);
end

% subtractionthreshold = graythresh(subtractionref);

subtracted = sub_scale_corr-subtractionthreshold.*1.1;
% subtracted = sub_scale_corr-subtractionthreshold.*1.1;
subzero = (subtracted<0);
subtractedzero = subtracted.*(~subzero);


Ie = subtractedzero;
a = find(Ie>0);
submax = zeros(size(Ie));
Ie(a)=50;





Ih = Ie>0;
Ihd = imdilate(Ih,strel('disk',1));
Ihdc = imclose(Ihd,strel('disk',4));
Ihdcf = imfill(Ihdc,'holes');
Im = Ihdcf;



%%%%% this is the ultimate addition for watershed segmentation!!!
see = strel('disk',1);
seo = strel('disk',8);
Ier = Im;
Isum = Ier;
for i=1:30
    Ier = imerode(Ier,see);
    Iero = imopen(Ier,seo);
%     Isum = Isum+(Iero.*i);
    Isum = Isum+(Iero);
    Ier=Iero;
end
Isum(Isum>130) = 130;
figure(66)
imagesc(Isum);
gausshed = gaussianBlurz(Isum,round(sigmafirst./dimdiff),round(kernelgsizefirst./dimdiff));
imgt = -double(gausshed);
waterBoundary = imerode(Im,strel('disk',1));
imgt(~(waterBoundary>0)) = -Inf;
L=watershed(imgt);

L(waterBoundary<1) = 0;
If = L>1;
If = imerode(If,strel('disk',2));




time = tsn{frames};
tim = time(2:end);
if frames==10
figure(1)
imagesc(If)
stophere=1;
end
savethatimage(scenename,time,If.*255,frames,filename,channel)
end


end

function If = segmentationRFP(FinalImage,subdirname,scenename,filename,channel)
fig=1;
mkdir(strcat(channel));

% parameters
left = 0.004;
slopedown = 0.003;


dimdiff = 2048./size(FinalImage(:,:,1),1);

zerostrel = 2;
firststrel = round(30./(dimdiff.^2));
sigmafirst = firststrel.*3;
kernelgsizefirst = firststrel.*6;
% fracsmoothing = 0.5.*dimdiff;
fracsmoothing = 0.5;

foldername = 'mKate_flat';
tsn = determineTimeFrame(foldername);



% start
for frames = 1:size(FinalImage,3)
img = FinalImage(:,:,frames); 
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
b = find(sub == min(min(sub)),1,'first');
rattio = gaustwo(b)./gaus(b);
gaustwocorr = gaustwo./rattio;
sub_scale_corr = double(gaus) - double(gaustwocorr);




subtractionref = sub_scale_corr;
vec = reshape(subtractionref,size(subtractionref,1)^2,1);
[numbers,bincenters] = hist(double(vec),0:fracsmoothing:10000);



numbers = medfilt1(numbers, 10); %smooths curve
fraction = numbers./sum(numbers);

mf = max(fraction);

%%%%%%%%%%%%%%%%%%%%
left=0.5*mf;
slopedown=0.4*mf;
%%%%%%%%%%%%%%%%%%%%%

leftedge = find(fraction > left,1,'first');
insideslopedown = find(fraction(leftedge:end) < slopedown,1,'first');
insideslopeup = find(fraction(leftedge+insideslopedown:end) >0.0012,1,'first');
trough = min(fraction(leftedge+insideslopedown:leftedge+insideslopedown+insideslopeup));
troughindex = find(fraction(leftedge+insideslopedown:leftedge+insideslopedown+insideslopeup) == trough);
troughindexrounded = round(median(troughindex));



threshlocation = bincenters(leftedge+insideslopedown);
% 
% figure(22)
%     bar(bincenters,fraction);hold on
%     xlim([-500 1000])
%     ylim([0 0.1])
% stem ([threshlocation threshlocation],[0 1]);hold off
% drawnow

subtractionthreshold = threshlocation;

if size(subtractionthreshold,1)==size(subtractionthreshold,2)
else
     subtractionthreshold = mean(threshlocation);
end

% subtractionthreshold = graythresh(subtractionref);

subtracted = sub_scale_corr-subtractionthreshold.*1.04;
% subtracted = sub_scale_corr-subtractionthreshold.*1.1;
subzero = (subtracted<0);
subtractedzero = subtracted.*(~subzero);


Ie = subtractedzero;
a = find(Ie>0);
submax = zeros(size(Ie));
Ie(a)=50;


if frames==10
figure(2)
imagesc(Ie)
stophere=1;
end

%%%%%%%%%%%%%%%
Ih = Ie>0;
se = strel('disk',5);
Ihe = imerode(Ih,se);
Ihed = imdilate(Ihe,se);
If = Ihed;



stophere=1;
% time = settimecharacter(frames);
time = tsn{frames};
tim = time(2:end);
savethatimage(scenename,time,If.*255,frames,filename,channel)
end


end

function If = segmentationMNG(FinalImage,subdirname,scenename,filename,channel)
global nucleus_seg
fig=1;
mkdir(strcat(channel));

% parameters
left = 0.004;
slopedown = 0.003;

dimdiff = 2048./size(FinalImage(:,:,1),1);

zerostrel = 2;
firststrel = round(30./dimdiff);
sigmafirst = firststrel.*3;
kernelgsizefirst = firststrel.*6;
fracsmoothing = 0.5.*dimdiff;

dirlist = dir(nucleus_seg);
if isempty(dirlist)
    dirlist = dir('mKate_flat');
    foldername = 'mKate_flat';
else
    foldername = nucleus_seg;
end

tsn = determineTimeFrame(foldername);


% start
for frames = 1:size(FinalImage,3)
img = FinalImage(:,:,frames); 
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
b = find(sub == min(min(sub)),1,'first');
rattio = gaustwo(b)./gaus(b);
gaustwocorr = gaustwo./rattio;
sub_scale_corr = double(gaus) - double(gaustwocorr);




subtractionref = sub_scale_corr;
vec = reshape(subtractionref,size(subtractionref,1)^2,1);
[numbers,bincenters] = hist(double(vec),0:fracsmoothing:10000);



numbers = medfilt1(numbers, 10); %smooths curve
fraction = numbers./sum(numbers);

mf = max(fraction);

%%%%%%%%%%%%%%%%%%%%
left=0.3*mf;
slopedown=0.2*mf;
%%%%%%%%%%%%%%%%%%%%%

leftedge = find(fraction > left,1,'first');
insideslopedown = find(fraction(leftedge:end) < slopedown,1,'first');
insideslopeup = find(fraction(leftedge+insideslopedown:end) >0.0012,1,'first');
trough = min(fraction(leftedge+insideslopedown:leftedge+insideslopedown+insideslopeup));
troughindex = find(fraction(leftedge+insideslopedown:leftedge+insideslopedown+insideslopeup) == trough);
troughindexrounded = round(median(troughindex));



threshlocation = bincenters(leftedge+insideslopedown);
% 
% figure(22)
%     bar(bincenters,fraction);hold on
%     xlim([-500 1000])
%     ylim([0 0.1])
% stem ([threshlocation threshlocation],[0 1]);hold off
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
gaus = gaussianBlurz(IobrcbrF,round(sigmafirst./dimdiff),(kernelgsizefirst./dimdiff));

% sigma=40;
% kernelgsize=80;
% gaus = gaussianBlurz(gaus,sigmafirst.*2,kernelgsizefirst.*2);

imgt = -double(gaus);
% imgt(~(Igcfopen>0)) = -Inf;
imgt(~(Igcopenmax>0)) = -Inf;

L=watershed(imgt);

L(Igcfopendil<1) = 0;
% imagesc(L)
% colormap parula
If = L>1;

stophere=1;
% time = settimecharacter(frames);
time = tsn{frames};
tim = time(2:end);
savethatimage(scenename,time,If.*255,frames,filename,channel)
end


end

function If = segmentationBKGsecond(FinalImage,subdirname,scenename,filename,channel)
global nucleus_seg
fig=1;
mkdir(strcat(channel));
% parameters
dimdiff = 2048./size(FinalImage(:,:,1),1);

zerostrel = 2;
% firststrel = round(30./dimdiff);
firststrel = round(20./dimdiff);
sigmafirst = firststrel.*3;
kernelgsizefirst = firststrel.*6;
fracsmoothing = 0.5.*dimdiff;

%set the time
dirlist = dir(nucleus_seg);
if isempty(dirlist)
    dirlist = dir('mKate_flat');
    foldername = 'mKate_flat';
else
    foldername = nucleus_seg;
end

tsn = determineTimeFrame(foldername);

% start
for frames = 1:size(FinalImage,3)
img = FinalImage(:,:,frames); 
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
b = find(sub == min(min(sub)),1,'first');
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
cd('NucleusBinary_flat')
time = tsn{frames};
tim = time(2:end);


% time = settimecharacter(frames);
% tim = time(1:end);
filelist = dir(strcat('*t',tim,'_*'));
% if isempty(filelist)
%     stophere=1;
%     filelist = dir(strcat('*t',tim(2:end),'_*'));
% end
clog = imread(char(filelist.name));

% se = strel('disk',4);
se = strel('disk',8);
clogd = imdilate(clog,se);
If(logical(clogd))=1;
cd .. 

prim = bwperim(If);
se = strel('disk',2);
primd = imdilate(prim,se);
imgorig(primd) = max(max(imgorig));
% figure(99),imagesc(imgorig)
% drawnow
stophere=1;
savethatimage(scenename,time,If.*255,frames,filename,channel)
end


end


function savethatimage(scenename,time,Ie,frames,filename,channel)
cd(channel)

[a,b] = regexp(filename,'(_mKate|CFP|EGFP|DIC)');
fname = strcat(filename(1:a-1),channel,filename(b+1:end));
fname = filename;

[a,b] = regexp(fname,'_t[0-9]+');
fname(a:b) = strcat('_',time);
% % % if isempty(a)
% % %     [a,b] = regexp(fname,'t[0-9][0-9]');
% % %     tnum = str2double(fname(a+1:b))-1;
% % %     tm = round(tnum+str2double(time));
% % %     time = settimecharacter(tm);
% % %     if (length(time)>length(fname(a+1:b)))
% % %     time(1) = 't';
% % %     else
% % %     time = horzcat('t',time);
% % % %     disp(strcat('time',time));
% % %     end
% % %     fname(a:b) = time;
% % % else
% % %     time = horzcat('t',time);
% % % %     disp(strcat('2time',time));
% % %     fname(a:b) = time;
% % % end
% if length(time)==length(fname(a+1:b))
% time = horzcat('t',time);
% fname(a:b) = time;
% elseif length(fname(a+1:b))<length(time)
% time = horzcat('t',time(2:end));
% fname(a:b) = time;
% end

imwrite(uint8(Ie),fname,'tiff','WriteMode','overwrite');
cd .. 

% cd('tiffs')
% % imwrite(uint8(Ie),strcat(scenename,'_','t',time,'_NucleusBinary_flat.tif'));
% if ~isempty(dir(strcat(channel(1:2),'*'))) && frames == 1
% imwrite(uint8(Ie),strcat(channel,'.tif'),'tiff','WriteMode','overwrite');
% else
% imwrite(uint8(Ie),strcat(channel,'.tif'),'tiff','WriteMode','append');
% end
% cd ..

end

function displaynormalizationfactor(b,c,frames,scenename)
disp(strcat(num2str(b),'_',num2str(c),'_',num2str(frames),'_',scenename))
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
function segmentation(FinalImage,subdirname,scenename,filename,channel)
fig=1;
mkdir(strcat('NucleusBinary_flat'));
mkdir(strcat('c5_flat'));
% mkdir(strcat('NucleusBinary_flat'));
for frames = 1:size(FinalImage,3)
% for frames = 1:size(FinalImage,3)
img = FinalImage(:,:,frames); 

%% parameters
strelsize           =    3;      %3
peakthresh          =    200;    %set low...will be corrected by the segthresh.
sigma               =    20;          %80
kernelgsize         =    15;   %15 for 512x512
meanregion          =    2;
segthresh           =    2200;
strelsizegaus       =    3;  %10
strelsizesub        =    2;
subtractionthreshold =   1400;
Upz = 0.07;
Downz = 0.06;

mask_em         = zeros(size(img));

%% start

%smooth initial image
se = strel('disk',strelsize);
Ie = imerode(img,se);
Iobr = imreconstruct(Ie,img);
Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
img = Iobrcbr;


%     figure(fig)
%     fig = fig+1;
%     imagesc(img);
%     title('smoothed input image');


            %ONLY NECESSARY IF TRYING TO FIND THE PEAKS IN THE SMOOTHED GAUS 
            
            gaus = gaussianblur(img);
           
    % figure(fig)
    % fig = fig+1;
    % hold off
    % bar(bincenters,fraction)
    % ylim([0 0.01])
    % ylim([0 0.003])
    % xlim([0 6000])
    % hold on
    % stem(segthresh,1,'g');
    % hold off


% findpeaksgaus(gaus) % homemade function to find peaks after applying gaussian blur 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub = double(img) -double(gaus);%%%%%%% key step!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
div = double(gaus)./double(img);
maxlocale = max(max(div(10:500,10:500)));
d = find(div == maxlocale);
div(d);
positiond = d(1);
[e,f] = ind2sub(size(img),positiond);


 minlocale = min(min(sub(10:500,10:500)));
a = find(sub == minlocale);
sub(a);

position = a(1);
[b,c] = ind2sub(size(img),position);

%     figure(fig)
%     fig = fig+1;
%     imagesc(sub);
%     title('subtracted');
%     hold on
%     plot(c,b,'y+');


displaynormalizationfactor(e,f,frames,scenename)
%normalizationfactors = double(gaus(b,c))./double(img(b,c));
normalizationfactors = double(gaus(e,f))./double(img(e,f));
normalizationfactor = mean(mean(normalizationfactors));

% figure(fig)
% fig = fig+1;
% imagesc(subgauss);hold on
% plot(c,b,'y+')

scaledgaus = double(gaus)./normalizationfactor;
sub_scale_corr = double(img)-double(scaledgaus);

 
%     figure(fig)
%     fig = fig+1;
%     imagesc(sub_scale_corr);
%     title('subtracted_scaled');


subtractionref = sub_scale_corr;
vec = reshape(subtractionref,size(subtractionref,1)^2,1);
[numbers,bincenters] = hist(double(vec),0:2:10000);
numbers = medfilt1(numbers, 15); %smooths curve
fraction = numbers./sum(numbers);
% rightedge = find(fraction > 0.003,1,'last');
% slopedown = find(fraction(rightedge:end) <0.002,1,'first');
mf = max(fraction);
Upz = mf.*0.8;
Downz = mf.*0.7;
rightedge = find(fraction > Upz,1,'last');
slopedown = find(fraction(rightedge:end) <Downz,1,'first');
threshlocation=[];
if ~isempty(slopedown)
thresh = bincenters(rightedge+slopedown);
threshlocation = thresh - 0.05*thresh;
else
    if ~isempty(threshlocation)
    threshlocation = threshlocation;
    else
        threshlocation = 300;
    end
end

%  figure(fig)
% fig = fig+1;
% bar(bincenters,fraction);hold on
% ylim([0 0.01])
% xlim([0 2000]);hold on
% stem(threshlocation,1,'r');
%     hold off
%     
    subtractionthreshold = threshlocation;


subtracted = sub_scale_corr-ones(size(sub_scale_corr)).*subtractionthreshold;
subzero = (subtracted<0);
subtractedzero = subtracted.*(~subzero);

%     figure(fig)
%     fig = fig+1;
%     imagesc(subtractedzero,[0 20]);
%     title('subtracted_scaled_zeroed');


se = strel('disk',strelsizesub);
Ie = imerode(subtractedzero,se);
a = find(Ie>1);
submax = zeros(size(Ie));
Ie(a)=50;
% 
%     figure(fig)
%     fig=fig+1;
%     imagesc(Ie);
%     title('subtractedzeromax');
%     fig =1;


time = settimecharacter(frames);
        
stophere=1;
savethatimage(scenename,time,Ie,frames,filename,channel)
% imwrite(uint8(Ie),strcat(scenename,'_','t',time,'_NucleusBinary_flat.tif'));
end
stophere=1;
end


function tsn = determineTimeFrame(foldername)
cd(foldername)
fflist = dir(strcat('*.tif'));
ffnames = {fflist.name};
[a,b,c,d] = regexp(ffnames,'_t[0-9]+');
[a,b,c,d] = regexp(ffnames,'_t[0-9]++');
tnames = cellfun(@(x) x{1},d,'UniformOutput',0);
[a,b,c,d] = regexp(tnames,'t[0-9]++');%added extra step because sometimes the time frame parsin was mistaken
tnames = cellfun(@(x) x{1},d,'UniformOutput',0);
tsn = sort(tnames);
cd ..
end

