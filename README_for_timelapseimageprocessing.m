
%*files must be exported as TIFFs from Zen using channel names
%*the export folder must have the same name as the parent folder containing the zen experiment files
% such as parent folder = 'D:/kkim/2016_10_25 bcat cells' then export folder should be 'D:/kkim/2016_10_25 bcat cells/2016_10_25 bcat cells'

% parentFolderName = 'FrickPaperData';
% parentFolderName = 'kkim';
% mfile = mfilename('fullpath');
% [a,b] = regexp(mfile,parentFolderName);
% mfiledir = mfile(1:b+1);
% parentdir = mfiledir;
% A = strcat(parentdir);

%channels to be processed and analyzed
channelstoinput = {'_mKate','_EGFP','_CFP','DIC','_Hoechst'};

%set parent directory
A = 'C:\Users\zeiss\Pictures\KKIM\Bcat-Citrine Homozygote Experiments\';

%for loop iterating through subdirectories
for BB = {'2016-11-18 bcat homozygotes'};
B = BB{1};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% renamemCherrytoMkate(A,B)
% renamemWRONGtoRIGHT(A,B)
% cd('D:\Users\zeiss\Documents\MATLAB')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BackgroundAndFlatfieldCorrectionOfTimeLapseImages(A,B,channelstoinput);
% cd('D:\Users\zeiss\Documents\MATLAB')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SortBKGandFLATcorrectedImagesIntoFolders(A,B,channelstoinput);
% cd('D:\Users\zeiss\Documents\MATLAB')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BleachCorrectionOfTimeLapseImages(A,B,channelstoinput);
% cd('D:\Users\zeiss\Documents\MATLAB')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SegmentationOfTimeLapseImages(A,B);
% renameDirectoryItems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end



% A = 'D:\Frick\';
% B = '2016_01_25 smad3g smFISH';
% % 
% % preprocessBKGandFLATnew(A,B)
% % % lineageyo(A,B);
% % lineageyoz(A,B);
% % cd('D:\Users\zeiss\Documents\MATLAB')
% % % bleachcorrection(A,B)
% % bleachcorrectionNew(A,B)
% cd('C:\Users\zeiss\Documents\MATLAB')
% toughsegmentationforstacks(A,B)

