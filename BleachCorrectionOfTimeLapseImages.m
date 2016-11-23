function BleachCorrectionOfTimeLapseImages(A,B,channelstoinput)
primarydir = strcat(A,B,'\','flatfield_corrected');
cd(primarydir)
primarylist = dir('*_s*');
primarynames = {primarylist.name};
primaryone = char(primarynames{1});
    cd(strcat(A,B,'\flatfield_corrected\',primaryone))
    finaldir = dir('*flat\');
    finaldirnames = {finaldir.name};
    
    channelinputs =channelregexpmaker(channelstoinput);
   

%     p = regexp(finaldirnames,'c[0-3]');
%     p = regexp(finaldirnames,'(mKate|EGFP|CFP|DIC)');
    p = regexp(finaldirnames,channelinputs);
    px = cellfun(@isempty,p,'UniformOutput',1);
    finaldirnames = finaldirnames(~px);

    
for i = 1:length(finaldirnames);    
    cycle=1;
    for subdir=primarynames
    subdirname = char(subdir);
    cd(strcat(A,B,'\flatfield_corrected\',subdirname))
    SAVdir = strcat(A,B,'\flatfield_corrected\',subdirname,'\tiffs\');
    finaldirname = char(finaldirnames{i});
    cd (strcat(A,B,'\flatfield_corrected\',subdirname,'\',finaldirname));
    filelist = dir('*.tif*');
        if ~isempty(filelist)
            filenames = {filelist.name};
            medz = loadintostackanddeterminemedian(filenames,SAVdir,subdirname,channelinputs);
            eachmedian(cycle,:) = medz;
            cycle=cycle+1;
        end
    cd ..
    end
    
    medmedian = median(eachmedian,1);
    normmedmedian = medmedian./median(medmedian);
    
    
    for subdir=primarynames
    subdirname = char(subdir);
    cd(strcat(A,B,'\flatfield_corrected\',subdirname))
    SAVdir = strcat(A,B,'\flatfield_corrected\',subdirname,'\tiffs\');
    finaldirname = char(finaldirnames{i});
    cd (strcat(A,B,'\flatfield_corrected\',subdirname,'\',finaldirname));
    filelist = dir('*.tif*');
        if ~isempty(filelist)
            filenames = {filelist.name};
            loadintostackandbleachcorrect(filenames,SAVdir,normmedmedian,channelinputs);
        end
    cd ..
    end
    
    
end 
end



function med = loadintostackanddeterminemedian(filenames,SAVdir,subdirname,channelinputs)
med = zeros(1,length(filenames));
FinalImage=zeros(512,512,length(filenames),'double');
for j = 1:length(filenames)
    filename = char(filenames{j});
    img = imread(filename);
    FinalImage(:,:,j) = double(img);
    med(j) = median(median(double(img)));
end

cd (SAVdir)
[a,b] = regexp(filename,channelinputs);
channel = filename(a:b);
savename = strcat(channel,'_flat_bleach_corr.tif');
disp(subdirname)
%     for j = 1:length(filenames)
%     medcorr = med(j)./medmed;
%     img_bleach_corr = FinalImage(:,:,j)./medcorr;
% %     imwrite(uint16(img_bleach_corr),char(savename),'WriteMode','append');
%     end
end


function loadintostackandbleachcorrect(filenames,SAVdir,normmedian,channelinputs)

FinalImage=zeros(512,512,length(filenames),'double');
for j = 1:length(filenames)
    filename = char(filenames{j});
    img = imread(filename);
    FinalImage(:,:,j) = double(img);
end

cd (SAVdir)
[a,b] = regexp(filename,channelinputs);
channel = filename(a:b);
savename = strcat(channel,'_flat_bleach_corr.tif');
disp(savename)

filepresent = dir(savename);
if ~isempty(filepresent)
    delete(savename)
end

    for j = 1:length(filenames)
    img_bleach_corr = FinalImage(:,:,j)./normmedian(j);
    imwrite(uint16(img_bleach_corr),char(savename),'WriteMode','append');
    end
end


function channelinputs =channelregexpmaker(channelstoinput)
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
end
