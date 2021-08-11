function T=IBCsignature(DataDir,Dataset)
h=waitbar(0,'progress','Name','IBC signature analysis progress');
waitbar(0.01,h,'Loading G59 genes')
%import signature genes and arrange them into a vector of strings.
%previous gene symbols are included for backward compatibility.
waitbar(.01,h,'Reading G59 genes file') 
G59 = readtable([DataDir filesep 'IBC59signatureList.txt'],'ReadVariableNames',false); G59=table2cell(G59);
G59=G59(:); G59=unique(G59(~cellfun(@isempty,G59)));

%% Determine Dataset loaded
%%
if strcmp(Dataset,'GSE45581')

%load GSE45581 data and its genes annotation file
waitbar(.02,h,['Reading ' Dataset ' dataset file']) 
GSE45581 = table2cell(readtable([DataDir filesep 'GSE45581-GPL6480_series_matrix.txt']));%Data
waitbar(.05,h,['Reading ' Dataset ' annotation file']) 
GPL6480 = table2cell(readtable([DataDir filesep 'GPL6480-9577.txt']));%annotation file

Probes=GSE45581(37:end-1,1); %get probes
Data=cellfun(@str2double,GSE45581(37:end-1,2:end)); %get expression
Samples=GSE45581([9:11,36],2:end); %get samples names and include ER and HER2 status
Samples=strrep(Samples,' (FISH+)',''); Samples=strrep(Samples,' (FISH-)','');Samples=strrep(Samples,'sample type: ','');
Samples=Samples'; for i=1:size(Samples,1), Samples(i,1)={strjoin(Samples(i,:))}; end; Samples(:,2:end)=[];Samples=Samples';

%align probes with gene symbols
[Lia,Locb]=ismember(Probes,GPL6480(:,1)); Locb(Locb==0)=[];
Data(~Lia,:)=[]; Genes=GPL6480(Locb,7);
idx=cellfun(@isempty,Genes); %probes without gene annotation
Data(idx,:)=[]; Genes(idx)=[];

%for gene with multiple probes, use one with highest variance
G=unique(Genes);G(1)=[];
data=nan(numel(G),size(Data,2));% new data matrix
for i=1:size(G,1)
waitbar(((i/size(G,1))*0.7)+.05,h,[Dataset ': Processing gene symbols: ' num2str(i) '/' num2str(size(G,1))])    
Lia=ismember(Genes,G{i});

if sum(Lia)>1
idx=find(Lia);        
Vary=nanstd(Data(Lia,:),0,2);
idx=idx(Vary==max(Vary));
    data(i,:)=Data(idx,:);
else
    data(i,:)=Data(Lia,:);
end
end
waitbar(.75,h,['Done annotating ' Dataset ' data']) 
Data=data; Genes=G; %final expression variables
clearvars -except Data Genes Samples G59 h DataDir Dataset

%%%PLOT HEATMAP
waitbar(.8,h,['Plotting ' Dataset ' heatmap']) 
Lia=ismember(Genes,G59);
cmap=redbluecmap(11);
cgo = clustergram(log2(Data(Lia,:))','Standardize',1,...
    'RowLabels',Samples','Colormap',cmap);
addTitle(cgo, [Dataset ' G59 unsupervised hierarchical clustering heatmap'])

%%%PLOT PCA scatter
waitbar(.8,h,['Plotting ' Dataset ' heatmap']) 
Lia=ismember(Genes,G59);
[coeff, zscores, pcvars] = pca(log2(Data(Lia,:))');
pcvars=pcvars./sum(pcvars) * 100;
idx=nan(45,1); idx(21:40)=1;idx(1:20)=2; idx(41:45)=3; 
figure('Renderer','Painters')
plot(zscores(idx==1,1),zscores(idx==1,2),'om'), hold on
plot(zscores(idx==2,1),zscores(idx==2,2),'og'), hold on
plot(zscores(idx==3,1),zscores(idx==3,2),'ob'), box off
legend('IBC samples','non-IBC samples','Normal samples','Location','northwest')
xlabel(['1st PC (' num2str(round(pcvars(1),1)) '%)']);
ylabel(['2nd PC (' num2str(round(pcvars(2),1)) '%)']);
title([Dataset ' PCA plot']);

%%%% Random forest training and IBC scoring
Lia=ismember(Genes,G59);
Res=[21,22,23,25,27,28,29,31,32,35,38,26];% 12 IBC samples for training (Fig 1a Zare et al)
nonRes=[9,3,7,12,2,13,1,6,5,4,15,17];% 12 nonIBC samples for training (Fig 1a Zare et al)

% SET features and train
waitbar(.9,h,{'Selecting samples and trainning a ';['G59 model in ' Dataset ' dataset']}) 
rng(1); %for reproducibility
nTrees=5000; %number of random forest trees
trainData=[[Data(Lia,Res)',ones(length(Res),1)];[Data(Lia,nonRes)',zeros(length(nonRes),1)]];
features = trainData(:,(1:end-1));
classLabels = trainData(:,end);
B = TreeBagger(nTrees,features,classLabels,'OOBPrediction','On','Method','classification',...
    'OOBPredictorImportance','on');% trained model

waitbar(.95,h,'Predicting IBC score for all samples') 
[~,scores] = predict(B,Data(Lia,:)');%predict all samples

Summary=[Samples', num2cell(scores), num2cell(1:numel(Samples))'];
Summary=sortrows(Summary,3,'descend');%sort based on IBC probability score column
figure('Renderer','Painters');
bar(1:size(Summary,1),[Summary{:,3}]), ylabel 'IBC probability score';
box off, H=gca; H.XTick=1:length(Summary(:,1)); H.XTickLabel=Summary(:,1);H.XTickLabelRotation=90; ylim([0,1])
X=xlim; line([min(X),max(X)],[0.5,0.5],'Color','magenta','LineStyle','--'),H=gca; H.YTick=0:0.5:1;
title([Dataset ' IBC probability scores']);

T = table(Summary(:,1),Summary(:,3));
figure('Name','IBC probability scores')
uitable('Data',T{:,:},'ColumnName',{'Samples','IBCscores'},'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
waitbar(1,h,['Done with ' Dataset ' analysis']), pause(2), close (h) 

%%
elseif strcmp(Dataset,'GSE111477')

%load GSE111477 data and its genes annotation file
waitbar(.02,h,['Reading ' Dataset ' dataset file']) 
GSE111477 = table2cell(readtable([DataDir filesep 'GSE111477-GPL5175_series_matrix.txt']));%Data
waitbar(.05,h,['Reading ' Dataset ' annotation file']) 
GPL5175 = table2cell(readtable([DataDir filesep 'GPL5175-3188.txt']));%annotation file

Probes=GSE111477(35:end-1,1); %get probes
Data=cellfun(@str2double,GSE111477(35:end-1,2:end)); %get expressions
Samples=GSE111477([9,34],2:end); %get samples names
Samples=strrep(Samples,'condition: Non-Inflammatory','non-IBC');
Samples=strrep(Samples,'condition: Inflammatory','IBC');
Samples=Samples'; for i=1:size(Samples,1), Samples(i,1)={strjoin(Samples(i,:))}; end; Samples(:,2:end)=[];Samples=Samples';

%align probes with gene symbols
[Lia,Locb]=ismember(cellfun(@num2str,Probes,'UniformOutput',0),cellfun(@num2str,GPL5175(:,1),'UniformOutput',0)); 
Locb(Locb==0)=[];
Data(~Lia,:)=[]; 

% extract gene symbols
Genes=GPL5175(Locb,10);
idx=strfind(Genes,'//'); idx(cellfun(@isempty,idx))={[1,1]};
idx=cellfun(@(x) x(1:2),idx,'UniformOutput',0);
for i=1:numel(idx); Genes(i)={Genes{i}(idx{i,1}(1)+2:idx{i,1}(2)-1)};end
Genes=strrep(Genes,' ','');
idx=cellfun(@isempty,Genes); %probes without gene annotation
Data(idx,:)=[]; Genes(idx)=[];

%for gene with multiple probes, use one with highest variance
G=unique(Genes);
data=nan(numel(G),size(Data,2));% new data matrix
for i=1:size(G,1)
waitbar(((i/size(G,1))*0.7)+.05,h,[Dataset ': Processing gene symbols: ' num2str(i) '/' num2str(size(G,1))])    
Lia=ismember(Genes,G{i});

if sum(Lia)>1
idx=find(Lia);        
Vary=nanstd(Data(Lia,:),0,2);
idx=idx(Vary==max(Vary));
    data(i,:)=Data(idx,:);
else
    data(i,:)=Data(Lia,:);
end
end
waitbar(.75,h,['Done annotating ' Dataset ' data']) 
Data=data; Genes=G; %final expression variables
clearvars -except Data Genes Samples G59 h DataDir Dataset

%%%PLOT HEATMAP
waitbar(.8,h,['Plotting ' Dataset ' heatmap']) 
Lia=ismember(Genes,G59);
cmap=redbluecmap(11);
cgo = clustergram(Data(Lia,:)','Standardize',1,...
    'RowLabels',Samples','Colormap',cmap,'RowPDist','spearman');
addTitle(cgo, [Dataset ' G59 unsupervised hierarchical clustering heatmap'])

%%%PLOT PCA scatter
waitbar(.8,h,['Plotting ' Dataset ' heatmap']) 
Lia=ismember(Genes,G59);
[coeff, zscores, pcvars] = pca(Data(Lia,:)');
pcvars=pcvars./sum(pcvars) * 100;%pcvars(1)
idx=cellfun(@isempty,strfind(Samples,'non-IBC'));%IBC
figure('Renderer','Painters')
plot(zscores(idx,1),zscores(idx,2),'or'), hold on
plot(zscores(~idx,1),zscores(~idx,2),'ob'), box off
legend('IBC samples','non-IBC samples','Location','northeast')
xlabel(['1st PC (' num2str(round(pcvars(1),1)) '%)']);
ylabel(['2nd PC (' num2str(round(pcvars(2),1)) '%)']);
title([Dataset ' PCA plot']);Hfig=gcf; Hfig.Renderer='Painters';


%%%% Random forest training and IBC scoring
Lia=ismember(Genes,G59);
Res=[14:37,53:61];Res=Res(round(1:numel(Res)/16:numel(Res)));% 16 IBC samples for training 
nonRes=[1:13,38:52];nonRes=nonRes(round(1:numel(nonRes)/14:numel(nonRes)));% 14 nonIBC samples for training 

% SET features and train
waitbar(.9,h,{'Selecting samples and trainning a ';['G59 model in ' Dataset ' dataset']}) 
rng(1); %for reproducibility
nTrees=5000; %number of random forest trees
trainData=[[Data(Lia,Res)',ones(length(Res),1)];[Data(Lia,nonRes)',zeros(length(nonRes),1)]];
features = trainData(:,(1:end-1));
classLabels = trainData(:,end);
B = TreeBagger(nTrees,features,classLabels,'OOBPrediction','On','Method','classification',...
    'OOBPredictorImportance','on');% trained model

waitbar(.95,h,'Predicting IBC score for all samples') 
[~,scores] = predict(B,Data(Lia,:)');%predict all samples

Summary=[Samples', num2cell(scores), num2cell(1:numel(Samples))'];
Summary=sortrows(Summary,3,'descend');%sort based on IBC probability score column
figure('Renderer','Painters');
bar(1:size(Summary,1),[Summary{:,3}]), ylabel 'IBC probability score';
box off, H=gca; H.XTick=1:length(Summary(:,1)); H.XTickLabel=Summary(:,1);H.XTickLabelRotation=90; ylim([0,1])
X=xlim; line([min(X),max(X)],[0.5,0.5],'Color','magenta','LineStyle','--'),H=gca; H.YTick=0:0.5:1;
title([Dataset ' IBC probability scores']);

T = table(Summary(:,1),Summary(:,3));
figure('Name','IBC probability scores')
uitable('Data',T{:,:},'ColumnName',{'Samples','IBCscores'},'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
waitbar(1,h,['Done with ' Dataset ' analysis']), pause(2), close (h) 
        
%%
elseif strcmp(Dataset,'GSE5847')
%load GSE5847 data and its genes annotation file
waitbar(.02,h,['Reading ' Dataset ' dataset file']) 
GSE5847 = table2cell(readtable([DataDir filesep 'GSE5847-GPL96_series_matrix.txt']));%Data
waitbar(.05,h,['Reading ' Dataset ' annotation file']) 
GPL96 = table2cell(readtable([DataDir filesep 'GPL96-57554.txt']));%annotation file

Probes=GSE5847(39:end-1,1); %get probes
Data=cellfun(@str2double,GSE5847(39:end-1,2:end)); %get expressions
Samples=GSE5847([7,9,38],2:end); %get samples names
Samples=strrep(Samples,'human breast cancer ',''); Samples=strrep(Samples,'diagnosis: ','');
Samples=Samples'; for i=1:size(Samples,1), Samples(i,1)={strjoin(Samples(i,:))}; end; Samples(:,2:end)=[];Samples=Samples';
% Get tumor epithelium
Lia=cellfun(@isempty,strfind(Samples,'tumor'));
Data(:,Lia)=[]; Samples(Lia)=[];
Samples=strrep(Samples,'tumor epithelium ',''); 

%align probes with gene symbols
[Lia,Locb]=ismember(Probes,GPL96(:,1)); Locb(Locb==0)=[];
Data(~Lia,:)=[]; Genes=GPL96(Locb,11);

idx=strfind(Genes,'///'); 
for i=1:numel(idx); try Genes(i)={Genes{i}(1:idx{i,1}(1)-1)}; catch, end;end
Genes=strrep(Genes,' ','');
idx=cellfun(@isempty,Genes); %probes without gene annotation
Data(idx,:)=[]; Genes(idx)=[];

%for gene with multiple probes, use one with highest variance
G=unique(Genes);
data=nan(numel(G),size(Data,2));% new data matrix
for i=1:size(G,1)
waitbar(((i/size(G,1))*0.7)+.05,h,[Dataset ': Processing gene symbols: ' num2str(i) '/' num2str(size(G,1))])    
Lia=ismember(Genes,G{i});

if sum(Lia)>1
idx=find(Lia);        
Vary=nanstd(Data(Lia,:),0,2);
idx=idx(Vary==max(Vary));
    data(i,:)=Data(idx,:);
else
    data(i,:)=Data(Lia,:);
end
end
waitbar(.75,h,['Done annotating ' Dataset ' data']) 
Data=data; Genes=G; %final expression variables

clearvars -except Data Genes Samples G59 h DataDir Dataset

 %%%PLOT HEATMAP
waitbar(.8,h,['Plotting ' Dataset ' heatmap']) 
Lia=ismember(Genes,G59);
cmap=redbluecmap(11);
cgo = clustergram(Data(Lia,:)','Standardize',1,...
    'RowLabels',Samples','Colormap',cmap,...
'RowPDist','spearman');
addTitle(cgo, [Dataset ' G59 unsupervised hierarchical clustering heatmap'])

%%%PLOT PCA scatter
waitbar(.8,h,['Plotting ' Dataset ' heatmap']) 
[coeff, zscores, pcvars] = pca(Data(Lia,:)');
pcvars=pcvars./sum(pcvars) * 100;%pcvars(1)
idx=cellfun(@isempty,strfind(Samples,'non-IBC'));%IBC
figure('Renderer','Painters')
plot(zscores(idx,1),zscores(idx,2),'or'), hold on
plot(zscores(~idx,1),zscores(~idx,2),'ob'), box off

legend('IBC samples','non-IBC samples','Location','southeast')
xlabel(['1st PC (' num2str(round(pcvars(1),1)) '%)']);
ylabel(['2nd PC (' num2str(round(pcvars(2),1)) '%)']);
title([Dataset ' PCA plot'])


%%%% Random forest training and IBC scoring
Lia=ismember(Genes,G59);
Res=1:2:13;% 7 IBC samples for training
nonRes=14:2:48;% 18 nonIBC samples for training

% SET features and train
waitbar(.9,h,{'Selecting samples and trainning a ';['G59 model in ' Dataset ' dataset']}) 
rng(1); %for reproducibility
nTrees=5000; %number of random forest trees
trainData=[[Data(Lia,Res)',ones(length(Res),1)];[Data(Lia,nonRes)',zeros(length(nonRes),1)]];
features = trainData(:,(1:end-1));
classLabels = trainData(:,end);
B = TreeBagger(nTrees,features,classLabels,'OOBPrediction','On','Method','classification',...
    'OOBPredictorImportance','on');% trained model

waitbar(.95,h,'Predicting IBC score for all samples') 
[~,scores] = predict(B,Data(Lia,:)');%predict all samples

Summary=[Samples', num2cell(scores), num2cell(1:numel(Samples))'];
Summary=sortrows(Summary,3,'descend');%sort based on IBC probability score column
figure('Renderer','Painters');
bar(1:size(Summary,1),[Summary{:,3}]), ylabel 'IBC probability score';
box off, H=gca; H.XTick=1:length(Summary(:,1)); H.XTickLabel=Summary(:,1);H.XTickLabelRotation=90; ylim([0,1])
X=xlim; line([min(X),max(X)],[0.5,0.5],'Color','magenta','LineStyle','--'),H=gca; H.YTick=0:0.5:1;
title([Dataset ' IBC probability scores']);

T = table(Summary(:,1),Summary(:,3));
figure('Name','IBC probability scores')
uitable('Data',T{:,:},'ColumnName',{'Samples','IBCscores'},'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
waitbar(1,h,['Done with ' Dataset ' analysis']), pause(2), close (h) 
    
%%
elseif strcmp(Dataset,'TCGA')
%First, load GSE45581 data and its genes annotation file
waitbar(.02,h,['Reading GSE45581 dataset before ' Dataset]) 
GSE45581 = table2cell(readtable([DataDir filesep 'GSE45581-GPL6480_series_matrix.txt']));%Data
waitbar(.05,h,['Reading GSE45581 annotation before ' Dataset]) 
GPL6480 = table2cell(readtable([DataDir filesep 'GPL6480-9577.txt']));%annotation file

Probes=GSE45581(37:end-1,1); %get probes
Data=cellfun(@str2double,GSE45581(37:end-1,2:end)); %get expression
Samples=GSE45581([9:11,36],2:end); %get samples names and include ER and HER2 status
Samples=strrep(Samples,' (FISH+)',''); Samples=strrep(Samples,' (FISH-)','');Samples=strrep(Samples,'sample type: ','');
Samples=Samples'; for i=1:size(Samples,1), Samples(i,1)={strjoin(Samples(i,:))}; end; Samples(:,2:end)=[];Samples=Samples';

%align probes with gene symbols
[Lia,Locb]=ismember(Probes,GPL6480(:,1)); Locb(Locb==0)=[];
Data(~Lia,:)=[]; Genes=GPL6480(Locb,7);
idx=cellfun(@isempty,Genes); %probes without gene annotation
Data(idx,:)=[]; Genes(idx)=[];

%for gene with multiple probes, use one with highest variance
G=unique(Genes);G(1)=[];
data=nan(numel(G),size(Data,2));% new data matrix
for i=1:size(G,1)
waitbar(((i/size(G,1))*0.5)+.05,h,['GSE45581: Processing gene symbols: ' num2str(i) '/' num2str(size(G,1))])    
Lia=ismember(Genes,G{i});

if sum(Lia)>1
idx=find(Lia);        
Vary=nanstd(Data(Lia,:),0,2);
idx=idx(Vary==max(Vary));
    data(i,:)=Data(idx,:);
else
    data(i,:)=Data(Lia,:);
end
end
waitbar(.55,h,'Done annotating GSE45581 data') 
Data=data; Genes=G; %final expression variables
clearvars -except Data Genes Samples G59 h DataDir Dataset

%%% LOAD %% Test TCGA FireHouse Legacy

waitbar(.6,h,{['Reading and extracting ' Dataset ' dataset file....'];...
    'might take a while...'}) 
Raw=table2cell(readtable([DataDir filesep 'TCGA_FirehoseLegacy_clinical_RNA_Seq.txt'],'ReadVariableNames',false));

Dtcga=Raw(2:end,3:end);Dtcga=cellfun(@str2double,Dtcga,'UniformOutput',0);Dtcga=cell2mat(Dtcga);%TCGA RNAseq
Stcga=Raw(1,3:end);%TCGA samples
Gtcga=Raw(2:end,1);%TCGA Genes

waitbar(.65,h,['Reading ' Dataset ' clinical file']) 
ClinicalRaw=table2cell(readtable([DataDir filesep 'TCGA_FirehoseLegacy_clinical.txt'],'ReadVariableNames',false));

[Soverlap,Locb]=ismember(strrep(Stcga,'-01',''),ClinicalRaw(:,2)); Locb(Locb==0)=[];%find TCGA samples overlap with clinical samples

ClinicalDF=ClinicalRaw(Locb,[2,108:109]); % test=[Stcga(Soverlap)',ClinicalDF];%get overall survival
ClinicalDF=strrep(ClinicalDF,':LIVING',''); ClinicalDF=strrep(ClinicalDF,':DECEASED','');
ClinicalDF(cellfun(@isempty,ClinicalDF))={nan};
ClinicalDF(:,2:end)=cellfun(@str2double,ClinicalDF(:,2:end),'UniformOutput',0);
%Combine & Qnormalize
waitbar(.7,h,['Quantile normalizing ' Dataset ' with GSE45581']) 
[Lia,Locb]=ismember(Genes,Gtcga); Locb(Locb==0)=[];
DataNorm=quantilenorm([Data(Lia,:),Dtcga(Locb,:)],'median', 1);
% 
Data=DataNorm(:,1:45);
Dtcga=DataNorm(:,46:end);
Genes=Genes(Lia);%common in all

%%%% Random forest training and IBC scoring
Lia=ismember(Genes,G59);
Res=21:40;% Use all GSE45581 IBC samples for a more robust model
nonRes=1:20;% Use all GSE45581 nonIBC samples for a more robust model

% SET features and train
waitbar(.8,h,{'Selecting samples and trainning a G59 model in';...
    ['GSE45581 samples to predict IBC-like in ' Dataset ' dataset']}) 
rng(1); %for reproducibility
nTrees=5000; %number of random forest trees
trainData=[[Data(Lia,Res)',ones(length(Res),1)];[Data(Lia,nonRes)',zeros(length(nonRes),1)]];
features = trainData(:,(1:end-1));
classLabels = trainData(:,end);
B = TreeBagger(nTrees,features,classLabels,'OOBPrediction','On','Method','classification',...
    'OOBPredictorImportance','on');% trained model

waitbar(.9,h,{'Using GSE45581 model to predict';['IBC score for ' Dataset ' samples']}) 
% score TCGA
[~,scores] = predict(B,Dtcga(Lia,:)');
waitbar(.95,h,[Dataset ' KM plot for IBC-like samples']) 
scores=scores(Soverlap,:);%sample scores with clinical data
idx=scores(:,2)>=0.5; %index IBC-like samples

OSyrs=cell2mat(ClinicalDF(:,3))./12; %Overall survival in years; ClinicalDF(idx,:)
Event=cell2mat(ClinicalDF(:,2));

%KM plot
KM=[[OSyrs(idx),Event(idx)];...
    [OSyrs(~idx),Event(~idx)]];
KMgroups=[repmat({'IBC-like'},sum(idx),1);...
    repmat({'nonIBC-like'},sum(~idx),1)];

 [p, fh, stats] = MatSurv(KM(:,1),KM(:,2),KMgroups,'Xstep',2,'LineColor',[1,0,0;0,0,1],...'DrawMSL',1,...
     'Xlabel',{'Time (years)'},'Ylabel',{'Proportion of overall survival'},'BaseFontSize',10, 'Title', Dataset);

T = table(ClinicalDF,num2cell(scores(:,2)));
figure('Name',['IBC probability scores in ' Dataset])
uitable('Data',T{:,:},'ColumnName',{'Patient','OSstatus','OSmonths','IBCscores'},'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
waitbar(1,h,['Done with ' Dataset ' analysis']), pause(2), close (h)     
end
