%% STEP 1: (If need be) Uncompress Data.tar.gz if need be
[file,path]= uigetfile({'*.tar.gz'},'Select compressed Data file');%,'');
          if isnumeric(file) || isnumeric(path)
              errordlg('No valid file selected','ERROR!!!')
              return 
          end
          try  
            DataFile=[path file]; cd(path)
            h=waitbar(0.2,'Uncompressing .gz','Name',['Uncompressing '  file]);
            gunzip(DataFile)
            waitbar(0.5,h,'Uncompressing .tar')
            untar(DataFile(1:end-3))
            waitbar(1,h,['Done Uncompressing ' file]), pause(1), close(h)
          catch
            errordlg('File selected invalid','ERROR!!!')
          end

%% STEP 2: Run G59 analysis
DataDir = uigetdir(pwd,'Select uncompressed Data directory');

%select dataset to analyze from drop-down options
%DON'T close this uifigure
fig = uifigure('Name','Select dataset'); fig.Position=[fig.Position(1),fig.Position(2),200 100]; 
dd = uidropdown(fig,'Position',[5,5,195 22],'Items',{'GSE45581';'GSE111477';'GSE5847';'TCGA'}); 

%% STEP 3: Run this AFTER selecting dataset to analyze above, run the G59 analysis code
Dataset=dd.Value;disp([Dataset 'analysis...'])%selected dataset
IBCsignature(DataDir,Dataset)% IBC G59 analysis
