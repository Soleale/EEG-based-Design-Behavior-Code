%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Shumin Data Processing Pipeline for all datasets%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Baseline Event: Eyes Open + Fix Cross %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% To be done

% Move the following raw data into a new folder
% E:\OneDrive - Politecnico di Milano\201920\PhD\Shumin_Research\Protocol_shumin\Data\RawData
%
eeglab
close all
clear
clc


% Data Processing
% Step 1 - Info Prepare
ParaEM.fs=128;ParaEM.chanN=14;ParaEM.ChanLoc='emotiv14loc.ced';
TaskMarkerInfo=load("MarkerCSVRow.txt");
taskN=1:13; baselineevent='eyesopen';
reliablerange=16; % if the median is 16 larger/smaller than ref, than it should be removed
% Step 2 - Get Threshold for data cleaning and POW of Baseline
% Input: 
% TaskMarkerInfo: loaded at step 1
% ParaEM:         parameter relates to the headset, loaded at step 1
% baselineevet:   set in step 1 -- 'eyesopen'
% Output:
% normalEEGthreshold.eyesopen{subjecti,1}.FreqPOW.bandpow(:,bandi)-for step3
% GroupEEGthreshold.eyesopen(bandi).band(subjecti,:) - for step4 
[normalEEGthreshold.eyesopen,GroupEEGthreshold.eyesopen]=getthreshold(TaskMarkerInfo,ParaEM,'eyesopen',reliablerange);
% [normalEEGthreshold.eyesclosed,GroupEEGthreshold.eyesclosed]=getthreshold(TaskMarkerInfo,ParaEM,'eyesclosed',reliablerange);
% Step 3 - Get Processed POW of Games and Tasks/Stages
% Input: 
% TaskMarkerInfo: loaded at step 1
% TaskN:          1:13 - N. tasks to be processed, loaded at step 1
% Reliablerange:  % if the median is this range times larger/smaller than ref, than it should be removed, loaded at step 1
% ParaEM:         parameter relates to the headset, loaded at step 1
% baselineevet:   set in step 1 -- 'eyesopen'
% normalEEGthreshold.eyesopen{subjecti,1}.FreqPOW.bandpow(:,bandi) - get
% from step 2
% Output:
% ProcessedData(subjecti).Task{taski,1}.FreqPOW.bandpow(:,bandi)
% ProcessedTask4Stage{subjecti,1}(stagei).FreqPOW.bandpow(:,bandi)
% GroupEEGpow(taski).bandN(bandi).pow(subjecti,:)
% GroupEEGpow(stagei+13).bandN(bandi).powstage(subjecti,:)
[ProcessedData,GroupEEGpow]=...
ProcessWholeGroup4POW(TaskMarkerInfo,taskN,reliablerange,ParaEM,...
normalEEGthreshold,baselineevent);
% Step 4, get TRP and Remove Outlier
GroupTRP=GetTRPfromPOW(GroupEEGpow,GroupEEGthreshold,1:16);
%
% RemainData
PercentageRemain=getremainpercentage(ProcessedData);
save(['ProcessedGroup2DwMTsketch120s_reliablerange' num2str(reliablerange) 'subbandsab.mat'],'ProcessedData','normalEEGthreshold','reliablerange','PercentageRemain','GroupEEGpow','GroupTRP','GroupEEGthreshold','-v7.3')
%%
% save(['ProcessedDataGroup2_reliablerange' num2str(reliablerange) '.mat'],'GroupEEGthreshold','-append')
%% cluster
Group=ReadGroups();

%% Statistical Tests
% load('ProcessedDataGroup2_reliablerange16.mat', 'GroupTRP')%, 'GroupEEGthreshold','ParaEM')
bandi=5;%1-all, 2-theta, 3-alpha, 4-beta, 5-gamma
logscale=0;
% Between Groups
[H_BG,P_BG,TableinExcel_BG]=RunStatitis4TasksBetweenGroups(GroupTRP,Group,bandi,logscale);
% Between Tasks
[H_BT,P_BT,TableinExcel_BT]=RunStatitis4TasksBetweenTasks(GroupTRP,Group,bandi,logscale);

%% radarplotTRPMean(GroupCluster,GroupTRP,bandN,13,[0 4])
close all
load('ProcessedDataGroup2_reliablerange16.mat', 'GroupTRP')%, 'GroupEEGthreshold','ParaEM')
bandi=[3 4];%1-all, 2-theta, 3-alpha, 4-beta, 5-gamma
MakeRadarPlots(Group,GroupTRP,bandi,[0 3])
% MakeRadarPlots(Group,GroupTRP,bandN,RLim)
% MakeRadarPlotsBetweenTasks(Group,GroupTRP,bandi,[0 3])
% MakeRadarPlotsBetweenTasks(Group,GroupTRP,bandN,RLim)
%% boxplotTRPDistribution(GroupCluster,GroupTRP,bandN,13,[0 4])
close all
% load('ProcessedDataGroup2_reliablerange16.mat', 'GroupTRP')%, 'GroupEEGthreshold','ParaEM')
bandN=4;%1-all, 2-theta, 3-alpha, 4-beta, 5-gamma
Boxlimit=[0 4];
taskN=1:13;
% boxplotTRPDistribution(Group.All,GroupTRP,bandN,taskN,Boxlimit,' Whole Group')
% %% Profile
% boxplotTRPDistribution(Group.Cluster_Age,GroupTRP,bandN,taskN,Boxlimit,' Clustered by Age')
% boxplotTRPDistribution(Group.Cluster_Gender,GroupTRP,bandN,taskN,Boxlimit,' Clustered by Gender')
% boxplotTRPDistribution(Group.Cluster_Language,GroupTRP,bandN,taskN,Boxlimit,' Clustered by Language')
% % boxplotTRPDistribution(Group.Cluster_Hand,GroupTRP,bandN,taskN,RLim,' Clustered by Hand')
% 
% %% Task 3 - TTCT
% boxplotTRPDistribution(Group.Cluster_3V,GroupTRP,bandN,11,Boxlimit,' Clustered by Variety')
% boxplotTRPDistribution(Group.Cluster_3F,GroupTRP,bandN,11,Boxlimit,' Clustered by Fluency')
% boxplotTRPDistribution(Group.Cluster_3N,GroupTRP,bandN,11,Boxlimit,' Clustered by Novelty')
% %% Task 4 - Design with morphological table
% txt=input("plot the task in stages? Y/N [Y]","s");%taskN=14:16;
% if isempty(txt)
%     txt = 'Y';
% end
% switch txt
%     case "Y"
%         taskN=14:16;
%     case "N"
%         taskN=12;
% end
boxplotTRPDistribution(Group.Cluster_totalDwMT,GroupTRP,bandN,taskN,Boxlimit,' Clustered by Total Score')
boxplotTRPDistribution(Group.Cluster_LoD,GroupTRP,bandN,taskN,Boxlimit,' Clustered by Level of Detail')
boxplotTRPDistribution(Group.Cluster_4V,GroupTRP,bandN,taskN,Boxlimit,' Clustered by Viability')
% %% Task 5 - Empty the Glass
% boxplotTRPDistribution(Group.Cluster_5V,GroupTRP,bandN,13,Boxlimit,' Clustered by Variety')
% boxplotTRPDistribution(Group.Cluster_5F,GroupTRP,bandN,13,Boxlimit,' Clustered by Fluency')
% boxplotTRPDistribution(Group.Cluster_5N,GroupTRP,bandN,13,Boxlimit,' Clustered by Novelty')
%% Functions
radarplotTRPMeanbetweenGroups(GroupCluster,GroupTRP,bandN,taskN,RLim,AddTitle)
radarplotTRPMeanbetweenTasks(GroupCluster,GroupTRP,bandN,RLim,AddTitle)
function [H,P,TableinExcel] = RunStatitis4TasksBetweenTasks(GroupTRP,Group,bandN,logscale)
clusternames=fieldnames(Group);
TableinExcel.P.Tasks=zeros(1,14);
for clusteri=1:length(clusternames)
    for subclusteri=1:length(Group.(clusternames{clusteri,1}))
        Cluster=Group.(clusternames{clusteri,1}){1,subclusteri};
        [H.task2new_5,P.task2new_5]=statisticaltestsbetweentasks(GroupTRP,Cluster,bandN,[10 13],logscale);
        [H.task3_4s3,P.task3_4s3]=statisticaltestsbetweentasks(GroupTRP,Cluster,bandN,[11 16],logscale);
        [H.task2know_2new,P.task2know_2new]=statisticaltestsbetweentasks(GroupTRP,Cluster,bandN,[9 10],logscale);
        [H.task4s1_4s2,P.task4s1_4s2]=statisticaltestsbetweentasks(GroupTRP,Cluster,bandN,[14 15],logscale);
        [H.task4s1_4s3,P.task4s1_4s3]=statisticaltestsbetweentasks(GroupTRP,Cluster,bandN,[14 16],logscale);
        [H.task4s2_4s3,P.task4s2_4s3]=statisticaltestsbetweentasks(GroupTRP,Cluster,bandN,[15 16],logscale);
        TableinExcel.P.Tasks(end+1,:)=P.task2new_5.p_ks';
        TableinExcel.P.Tasks(end+1,:)=P.task3_4s3.p_ks';
        TableinExcel.P.Tasks(end+1,:)=P.task2know_2new.p_ks';
        TableinExcel.P.Tasks(end+1,:)=P.task4s1_4s2.p_ks';
        TableinExcel.P.Tasks(end+1,:)=P.task4s1_4s3.p_ks';
        TableinExcel.P.Tasks(end+1,:)=P.task4s2_4s3.p_ks';
    end
    TableinExcel.P.Tasks(end+1,:)=zeros(1,14);
end
end
function [H,P,TableinExcel]=RunStatitis4TasksBetweenGroups(GroupTRP,Group,bandi,logscale)
GroupN=1:2;
% %Task 3 - TTCT
[H.task3_V,P.task3_V]=statisticaltestsbetweengroups(GroupTRP,Group.Cluster_3V,bandi,11,GroupN,logscale);
[H.task3_N,P.task3_N]=statisticaltestsbetweengroups(GroupTRP,Group.Cluster_3N,bandi,11,GroupN,logscale);
[H.task3_F,P.task3_F]=statisticaltestsbetweengroups(GroupTRP,Group.Cluster_3F,bandi,11,GroupN,logscale);
%Task 4 - Design with morphological table
[H.task4_LoD,P.task4_LoD]=statisticaltestsbetweengroups(GroupTRP,Group.Cluster_LoD,bandi,[12 14:16],GroupN,logscale);
[H.task4_V,P.task4_V]=statisticaltestsbetweengroups(GroupTRP,Group.Cluster_4V,bandi,[12 14:16],GroupN,logscale);
[H.task4_totalDwMT1,P.task4_totalDwMT1]=statisticaltestsbetweengroups(GroupTRP,Group.Cluster_totalDwMT,bandi,[12 14:16],GroupN,logscale);%group A vs group B
% [H.task4_totalDwMT2,P.task4_totalDwMT2]=statisticaltestsbetweengroups(GroupTRP,Group.Cluster_totalDwMT,bandi,[12 14:16],[1 3],logscale);%group A vs group C
% [H.task4_totalDwMT3,P.task4_totalDwMT3]=statisticaltestsbetweengroups(GroupTRP,Group.Cluster_totalDwMT,bandi,[12 14:16],[2 3],logscale);%group B vs group C
%Task 5 - Empty the Glass
[H.task5_F,P.task5_F]=statisticaltestsbetweengroups(GroupTRP,Group.Cluster_5F,bandi,13,GroupN,logscale);
[H.task5_N,P.task5_N]=statisticaltestsbetweengroups(GroupTRP,Group.Cluster_5N,bandi,13,GroupN,logscale);
[H.task5_V1,P.task5_V1]=statisticaltestsbetweengroups(GroupTRP,Group.Cluster_5V,bandi,13,GroupN,logscale);
[H.task5_V2,P.task5_V2]=statisticaltestsbetweengroups(GroupTRP,Group.Cluster_5V,bandi,13,[1 3],logscale);
[H.task5_V3,P.task5_V3]=statisticaltestsbetweengroups(GroupTRP,Group.Cluster_5V,bandi,13,[2 3],logscale);
% %Profile
[H.tasks_Age,P.tasks_Age]=statisticaltestsbetweengroups(GroupTRP,Group.Cluster_Age,bandi,1:16,GroupN,logscale);
[H.tasks_Gender,P.tasks_Gender]=statisticaltestsbetweengroups(GroupTRP,Group.Cluster_Gender,bandi,1:16,GroupN,logscale);
[H.tasks_Language,P.tasks_Language]=statisticaltestsbetweengroups(GroupTRP,Group.Cluster_Language,bandi,1:16,GroupN,logscale);
% [H.tasks_Hand,P.tasks_Hand]=statisticaltestsbetweengroups(GroupTRP,Group.Cluster_Hand,bandi,1:16,GroupN,logscale);
%Arrange to excel
TableinExcel.P.Task3(1,:)=P.task3_V.p_ks';
TableinExcel.P.Task3(2,:)=P.task3_N.p_ks';
TableinExcel.P.Task3(3,:)=P.task3_F.p_ks';

TableinExcel.P.Task4(1:4,:)=P.task4_V.p_ks';
TableinExcel.P.Task4(5:8,:)=P.task4_LoD.p_ks';
TableinExcel.P.Task4(9:12,:)=P.task4_totalDwMT1.p_ks';
% TableinExcel.P.Task4(13:16,:)=P.task4_totalDwMT2.p_ks';
% TableinExcel.P.Task4(17:20,:)=P.task4_totalDwMT3.p_ks';

TableinExcel.P.Task5(1,:)=P.task5_F.p_ks';
TableinExcel.P.Task5(2,:)=P.task5_N.p_ks';
TableinExcel.P.Task5(3,:)=P.task5_V1.p_ks';
TableinExcel.P.Task5(4,:)=P.task5_V2.p_ks';
TableinExcel.P.Task5(5,:)=P.task5_V3.p_ks';

TableinExcel.P.Profile(1:14,:)=P.tasks_Age.p_ks;
TableinExcel.P.Profile((1:14)+15,:)=P.tasks_Gender.p_ks;
TableinExcel.P.Profile((1:14)+15*2,:)=P.tasks_Language.p_ks;
% TableinExcel.P.Profile((1:14)+15*3,:)=P.tasks_Hand.p_ks;

end
function Group=ReadGroups()
T = readtable('GroupAssignment _All.xlsx',"Range",'A1:Q34');
Group.All={1:33};
Group.Cluster_Age = {T.Age>=24,T.Age<24};
Group.Cluster_Gender = {strcmp(T.Gender,'M'),strcmp(T.Gender,'F')};
% Group.Cluster_Hand = {strcmp(T.Hand,'R'),strcmp(T.Hand,'E')};
Group.Cluster_Language = {strcmp(T.Language,'ITA'),strcmp(T.Language,'CHN')};
% GroupCluster_Hand = {strcmp(T.Hand,'R'),strcmp(T.Hand,'L'),strcmp(T.Hand,'E')};
%Task 3 - TTCT
Group.Cluster_3V = {strcmp(T.Task3Sketch_Variety,'A'),strcmp(T.Task3Sketch_Variety,'B')};
Group.Cluster_3F = {strcmp(T.Task3Sketch_Fluency,'A'),strcmp(T.Task3Sketch_Fluency,'B')};
Group.Cluster_3N = {strcmp(T.Task3Sketch_Novelty,'A'),strcmp(T.Task3Sketch_Novelty,'B')};
%Task 4 - Design with morphological table
Group.Cluster_LoD = {strcmp(T.Task4DwMT_LevelofDetail,'A'),strcmp(T.Task4DwMT_LevelofDetail,'B')};
Group.Cluster_4V = {strcmp(T.Task4DwMT_Viability,'A'),strcmp(T.Task4DwMT_Viability,'B')};
Group.Cluster_totalDwMT =  {strcmp(T.Task4DwMT_Total,'A'),strcmp(T.Task4DwMT_Total,'B')};%,strcmp(T.Task4DwMT_Total,'C')};
%Task 5 - Empty the Glass
Group.Cluster_5V = {strcmp(T.Task5PS_Variety,'A'),strcmp(T.Task5PS_Variety,'B'),strcmp(T.Task5PS_Variety,'C')};
Group.Cluster_5F = {strcmp(T.Task5PS_Fluency,'A'),strcmp(T.Task5PS_Fluency,'B')};
Group.Cluster_5N = {strcmp(T.Task5PS_Novelty,'A'),strcmp(T.Task5PS_Novelty,'B')};
 %Task 2 - AU
 Group.Cluster_2F = {strcmp(T.Task2AU_Fluency,'A'),strcmp(T.Task2AU_Fluency,'B')};
end
function PercentageRemain=getremainpercentage(ProcessedData)
for subjecti=1:33
    for taski=1:16        
            PercentageRemain.TotalRawDataPoints(subjecti,taski)=size(ProcessedData(subjecti).Task{taski,1}.rawEEG,2);
        for bandi=1:5
            PercentageRemain.Band(bandi).BandLeftDataPoints(subjecti,taski)=size(ProcessedData(subjecti).Task{taski,1}.SubbandSignal(bandi).s,2);            
            PercentageRemain.Band(bandi).PercentageLeft(subjecti,taski)=PercentageRemain.Band(bandi).BandLeftDataPoints(subjecti,taski)/PercentageRemain.TotalRawDataPoints(subjecti,taski);            
        end
    end
end
end
function boxplotTRPDistribution(GroupCluster,GroupTRP,bandN,taskN,Boxlimit,AddTitle)
taskname={'eyesopen','eyesclosed','game1','game2','game3','game4','game5',...
    'warm-up task','AU-Known','AU-New','Complete Figure','DwMT','Empty the Glass',...
    'Reading Stage','Selecting Stage','Sketching Stage'};
BandTag={'All','Theta','Alpha','Beta','Gamma'};
GroupTag={'Group A','Group B','Group C'};
GroupN = size(GroupCluster,2);
ColorTag={[0.8500, 0.3250, 0.0980],[0, 0.4470, 0.7410],[0.4660, 0.6740, 0.1880]};%["r","b","g"];


for bandi=1:length(bandN)     
    for taski=1:length(taskN)
        groupMember=GroupCluster{1,1};
        OutlierRemoved2plot=nan(min(sum(groupMember),length(groupMember)),14*GroupN);
        OutlierRemoved2plot(:,1:GroupN:GroupN*14) = GroupTRP(taskN(taski)).bandN(bandN(bandi)).trp_OutlierRemoved(groupMember,:);
        figure
        boxplot(OutlierRemoved2plot,'PlotStyle','compact','whisker', inf,"Colors",ColorTag{1,1});
        ax=gca;
        hold on
        plot(1:GroupN:GroupN*14,nanmedian(OutlierRemoved2plot(:,1:GroupN:GroupN*14)),'Color',ColorTag{1,1})
        if GroupN>1
            hold on
            for groupi=2:GroupN
                groupMember=GroupCluster{1,groupi};
                OutlierRemoved2plot=nan(min(sum(groupMember),length(groupMember)),14*GroupN);
                OutlierRemoved2plot(:,groupi:GroupN:GroupN*14) = GroupTRP(taskN(taski)).bandN(bandN(bandi)).trp_OutlierRemoved(groupMember,:);
                boxplot(ax,OutlierRemoved2plot,'PlotStyle','compact','whisker', inf,"Colors",ColorTag{1,groupi});
                hold on
                plot(groupi:GroupN:GroupN*14,nanmedian(OutlierRemoved2plot(:,groupi:GroupN:GroupN*14)),'Color',ColorTag{1,groupi})
                hold on
            end            
        end
        grid on        
        ax.YLabel.String="TRP";
        ax.YLim=Boxlimit;
        ax.XTick=(GroupN/2:GroupN:GroupN*14)-0.5;
        ax.XTickLabelRotation=0;
        ax.XTickLabel=(pad({' AF3',' F7',' F3',' FC5',' T7',' P7',' O1',' O2',' P8',' T8',' FC6',' F4',' F8',' AF4'},12,'left'));        
%         set(ax,'XTickLabel',{});
        ax.Title.String=([taskname{1,taskN(taski)} ' TRP - ' BandTag{1,bandN(bandi)} ' Band' AddTitle]);
        legend(GroupTag{1,1:GroupN});
    end
end
end
function GroupTRP=GetTRPfromPOW(GroupEEGpow,GroupEEGthreshold,tasks)
for taski=tasks
    for bandi=1:5
        % GroupEEGpow(taski).bandN(bandi).pow  33*14
        GroupTRP(taski).bandN(bandi).trp=GroupEEGpow(taski).bandN(bandi).pow./GroupEEGthreshold.eyesopen(bandi).band;%#ok
        GroupTRP(taski).bandN(bandi).trp_OutlierRemoved = removeOutlierbytrim(GroupTRP(taski).bandN(bandi).trp);%#ok
    end
end
end
function datatable_OutlierRemoved = removeOutlierbytrim(datatable)
[subjectN,channelN]=size(datatable);
N_subject_2getThreshold=ceil(0.1*subjectN/2);
for chani=1:channelN
    sortedtrp=sort(datatable(:,chani));
    sortedtrp=sortedtrp(N_subject_2getThreshold+1:end+1-N_subject_2getThreshold);
    thresholdh=mean(sortedtrp)+std(sortedtrp)*1.96;
    thresholdl=mean(sortedtrp)-std(sortedtrp)*1.96;
    datatable((datatable(:,chani)>thresholdh),chani)=nan;
    datatable((datatable(:,chani)<thresholdl),chani)=nan;
end
datatable_OutlierRemoved=datatable;
end
function [normalEEGthreshold,GroupEEGthreshold]=getthreshold(TaskMarkerInfo,ParaEM,baselineevent,reliablerange)
subjectN=TaskMarkerInfo(:,1);%:128;
TaskOrder=TaskMarkerInfo(:,2);
TaskMarkers=TaskMarkerInfo(:,3:end);

for subjecti=1:length(subjectN)
    ParaSub(subjecti).subjectistr=num2str(subjectN(subjecti));%#ok
    ParaSub(subjecti).taskorder=TaskOrder(subjecti);%#ok
    ParaSub(subjecti).taskmarkers=TaskMarkers(subjecti,:);%#ok
    ParaSub(subjecti).MarkerEventsRow=FromTaskOrder2MarkerRow(ParaSub(subjecti).taskorder);%#ok

    EEGdata = csvread(['D:\2022 Shumin_Thesis Prep\RawData\EmotivCSV\' ParaSub(subjecti).subjectistr '.csv'],2,0);% This excludes the first row of csv file(since they are info, not eeg data)
    EEGdata(:,[1:3 18:end])=[];

    normalEEGthreshold{subjecti,1} = outlierthreshold(EEGdata,ParaSub(subjecti),ParaEM,baselineevent,reliablerange);%#ok
    for bandi=1:5
        GroupEEGthreshold(bandi).band(subjecti,:)=normalEEGthreshold{subjecti,1}.FreqPOW.bandpow(:,bandi)';%#ok
    end
end
end
function [ProcessedData,GroupEEGpow]=ProcessWholeGroup4POW(TaskMarkerInfo,taskN,reliablerange,ParaEM,normalEEGthreshold,baselineevent)
subjectN=TaskMarkerInfo(:,1);%:128;
TaskOrder=TaskMarkerInfo(:,2);
TaskMarkers=TaskMarkerInfo(:,3:end);
for subjecti=1:length(subjectN)
    ParaSub(subjecti).subjectistr=num2str(subjectN(subjecti));%#ok
    ParaSub(subjecti).taskorder=TaskOrder(subjecti);%#ok
    ParaSub(subjecti).taskmarkers=TaskMarkers(subjecti,:);%#ok
    ParaSub(subjecti).MarkerEventsRow=FromTaskOrder2MarkerRow(ParaSub(subjecti).taskorder);%#ok
    EEGdata = csvread(['D:\2022 Shumin_Thesis Prep\RawData\EmotivCSV\' ParaSub(subjecti).subjectistr '.csv'],2,0);% This excludes the first row of csv file(since they are info, not eeg data)
    EEGdata(:,[1:3 18:end])=[];

    for taski=taskN
        % 1 eyes open, 2 eyes closed
        % 3-7 Game in time order(First Game-Last Game), 8-13 Task in task number order(Task1-Task5)
        % Game are not aligned with the task, Last Game might/not be the one before
        % Task 5
        [ProcessedData(subjecti).Task{taski,1},TimeStartEnd(subjecti,taski*2-1),TimeStartEnd(subjecti,taski*2)]=DataProcessingPipeline(EEGdata,ParaSub(subjecti),ParaEM,normalEEGthreshold.(baselineevent){subjecti,1}.FreqPOW.bandpow,reliablerange,taski);%#ok
        for bandi=1:5
            GroupEEGpow(taski).bandN(bandi).pow(subjecti,:)=ProcessedData(subjecti).Task{taski, 1}.FreqPOW.bandpow(:,bandi)';%#ok
        end

    end
    if find(taskN==12)
        ProcessedTask4Stage{subjecti,1}=DataProcessingPipelineTask4(ProcessedData(subjecti).Task{12,1}.rawEEG,ParaSub(subjecti),ParaEM,normalEEGthreshold.(baselineevent){subjecti,1}.FreqPOW.bandpow,reliablerange);%#ok
        for stagei=1:3 %%modified, now 3 stages in DwMT
            ProcessedData(subjecti).Task{stagei+13,1}=ProcessedTask4Stage{subjecti,1}(stagei);
            for bandi=1:5
                GroupEEGpow(stagei+13).bandN(bandi).pow(subjecti,:)=ProcessedTask4Stage{subjecti,1}(stagei).FreqPOW.bandpow(:,bandi)';
            end
        end

    end
end
end
function [ProcessedData,TimeStart,TimeEnd]=DataProcessingPipeline(eegdataraw,ParaSub,ParaDevice,normalEEGthreshold,reliablerange,taski)
TimeStart=ParaSub.taskmarkers(ParaSub.MarkerEventsRow(taski));
TimeEnd=ParaSub.taskmarkers(ParaSub.MarkerEventsRow(taski)+1);

if size(eegdataraw,1)~=14
    eegdataraw=eegdataraw';
end
eegdata=eegdataraw(:,TimeStart:TimeEnd);

ProcessedData.rawEEG=eegdata;

EEGraw = pop_importdata('dataformat','array','nbchan',ParaDevice.chanN,'data',eegdata,'setname','eegdata','srate',ParaDevice.fs,'pnts',0,'xmin',0,'chanlocs',ParaDevice.ChanLoc);
EEGdcrm = dcremove(ParaDevice,EEGraw);
ProcessedData.dcrmEEG = EEGdcrm.data;

ThinkTS=LoadThinkingInterval(ParaSub,ParaDevice,taski);
EEGdcrm=removeSignalByAudio(ThinkTS,EEGdcrm);

[ProcessedData.FreqPOW,ProcessedData.SubbandSignal]=getPowerFreq(EEGdcrm,ParaDevice,normalEEGthreshold,reliablerange);


end
function ProcessedTask4Stage=DataProcessingPipelineTask4(Task4RawData,ParaSub,ParaDevice,normalEEGthreshold,reliablerange)
marker = readtable('Task4Segments_Group2.txt');
subjectid=find(marker.SubjectID==str2double(ParaSub.subjectistr));

TimeStart=floor([marker.S1_S(subjectid) marker.S2_S(subjectid) marker.S3_S(subjectid)]*ParaDevice.fs)+1;
TimeEnd=floor([marker.S1_E(subjectid) marker.S2_E(subjectid) min(marker.S3_E(subjectid),marker.S3_S(subjectid)+120)]*ParaDevice.fs)+1;%MODIFIED shumin - 2min of skeching
TimeEnd(end)=min(TimeEnd(end),length(Task4RawData));

if size(Task4RawData,1)~=14
    Task4RawData=Task4RawData';
end
%% Join 2STAGES 
% data.stage1=Task4RawData(:,[TimeStart(1):TimeEnd(1) TimeStart(2):TimeEnd(2)]);
% data.stage2 =Task4RawData(:,TimeStart(3):TimeEnd(3));

data.stage1 =Task4RawData(:,TimeStart(1):TimeEnd(1));
data.stage2 =Task4RawData(:,TimeStart(2):TimeEnd(2));
data.stage3 =Task4RawData(:,TimeStart(3):TimeEnd(3));

for stagei=1:3%MODIFIED,NOW THREE STAGES IN DwMT

    eegdatastagei= data.(['stage' num2str(stagei)]);

    ProcessedTask4Stage(stagei).rawEEG=eegdatastagei;%#ok

    EEGraw = pop_importdata('dataformat','array','nbchan',ParaDevice.chanN,'data',eegdatastagei,'setname','eegdata','srate',ParaDevice.fs,'pnts',0,'xmin',0,'chanlocs',ParaDevice.ChanLoc);
    EEGdcrm = dcremove(ParaDevice,EEGraw);
    ProcessedTask4Stage(stagei).dcrmEEG = EEGdcrm.data;%#ok

    [ProcessedTask4Stage(stagei).FreqPOW,ProcessedTask4Stage(stagei).SubbandSignal]=getPowerFreq(EEGdcrm,ParaDevice,normalEEGthreshold,reliablerange);%#ok
end
end
function EEGdcrm=dcremove(Para,EEG) %eegdata size (14*length)
EEGdcrm=EEG;
eeg=EEG.data;
IIR_TC = Para.fs*2; %2 second time constant- adjust as required
[rows, cols]=size(eeg);
if cols == Para.chanN
    eeg=eeg';
    [rows, cols]=size(eeg);
end
EEGdcrm.data = zeros(rows,cols);
back= eeg(:,1);
for ii = 2 : cols
    back = (back*(IIR_TC-1) + eeg(:,ii)) / IIR_TC;
    EEGdcrm.data(:,ii)=eeg(:,ii)- back;
end
end
function [FreqPOW,SubbandSignal]=getPowerFreq(EEGdcrm,ParaDevice,normalEEGthreshold,reliablerange)
% BandN=[4 4 7 13 30; 45 7 13 30 45];

BandN=[7 10 12 15 20; 10 12 15 20 30];


if ~isempty(EEGdcrm.data)
    for bandi=1:length(BandN)
        EEGfilter = pop_eegfiltnew(EEGdcrm, 'locutoff',BandN(1,bandi),'hicutoff',BandN(2,bandi),'plotfreqz',0);
        SubbandSignal(bandi).s_br=EEGfilter.data;%#ok
        if ~isnan(normalEEGthreshold)
            EEGfilter=removeoutlier(EEGfilter,ParaDevice,normalEEGthreshold(:,bandi)*reliablerange);
            FreqPOW.median(bandi).win = EEGfilter.medianineachwin;
        end
        %             FreqPOW.bandpow(:,bandi) = mean(abs(EEGfilter.data).^2,2);
        %         else
        FreqPOW.bandpow(:,bandi) = median(EEGfilter.data.^2,2);       
        %         end
        SubbandSignal(bandi).s = EEGfilter.data;%#ok
    end
else
    FreqPOW.bandpow(:,:)=nan;
    SubbandSignal(1).s = nan;
end
end
function MarkerEventsRow=FromTaskOrder2MarkerRow(TaskOrder)
TaskOrderstr=num2str(TaskOrder);
Task2Position=find(TaskOrderstr=='2');
switch Task2Position
    case 2
        MarkerEventsRow=[1 3 5  9 15 19 23 7 11 13 17 21 25];
    case 3
        MarkerEventsRow=[1 3 5  9 13 19 23 7 15 17 11 21 25];
    case 4
        MarkerEventsRow=[1 3 5  9 13 17 23 7 19 21 11 15 25];
    case 5
        MarkerEventsRow=[1 3 5  9 13 17 21 7 23 25 11 15 19];
end
TaskOrderstr_n2=TaskOrderstr(2:end);
TaskOrderstr_n2(TaskOrderstr_n2=='2')=[];

[~,I]=sort([str2double(TaskOrderstr_n2(1)) str2double(TaskOrderstr_n2(2)) str2double(TaskOrderstr_n2(3))]);
MarkerEventsRow(end-2:end)=MarkerEventsRow(I+10);

end
function ThinkTS=LoadThinkingInterval(ParaSub,ParaDevice,taski)
TaskOrderstr=num2str(ParaSub.taskorder);
TAGi=taski;
if taski>=8
    if taski==10
        TAGi=find(TaskOrderstr=='2')*2+3;
    else
        TAGi=find(TaskOrderstr==num2str(max(taski-7-floor(taski/11),1)))*2+2;
        if TAGi>find(TaskOrderstr=='2')*2+2
            TAGi=TAGi+1;
        end
    end
elseif taski>=3
    TAGi=2*taski-3;
    if taski>find(TaskOrderstr=='2')+2
        TAGi=TAGi+1;
    end
end
ThinkTS=floor((load(['D:\2022 Shumin_Thesis Prep\ProcessedData\Video\' ParaSub.subjectistr '_ThinkTS.mat']).ThinkTS.(['TaG' num2str(TAGi)])(1:end-1)-48000)/48000*ParaDevice.fs)+1;
end
function EEGdcrm=removeSignalByAudio(ThinkTStask,EEGdcrm)
removeSpeakdata=[];
ThinkTStask(end)=min(ThinkTStask(end),length(EEGdcrm.data));
for intervali=1:length(ThinkTStask)/2
    if ThinkTStask(intervali*2)-ThinkTStask(intervali*2-1)<EEGdcrm.srate*0.5
        ThinkTStask(intervali*2)=nan;
        ThinkTStask(intervali*2-1)=nan;
    end
end
ThinkTStask(isnan(ThinkTStask))=[];

for intervali=1:length(ThinkTStask)/2
    removeSpeakdata=[removeSpeakdata EEGdcrm.data(:,ThinkTStask(intervali*2-1):ThinkTStask(intervali*2))];%#ok
end

EEGdcrm.data=removeSpeakdata;
EEGdcrm.pnts=size(EEGdcrm.data,2);
EEGdcrm.times(size(EEGdcrm.data,2)+1:end)=[];
end
function EEGfilter=removeoutlier(EEGfilter,ParaDevice,normalEEGthreshold)
filtereddata=EEGfilter.data;
filtereddataupdate=filtereddata;
winLSpercentage=0.25;
%filtereddata to find interval
fs=ParaDevice.fs;
winlength=floor(winLSpercentage*fs);
winshift=floor(winLSpercentage*winlength);
winstart=1:winshift:(length(filtereddata)-winlength+winshift-1);
winend=winstart+winlength;
winend(end)=min([winend(end) length(filtereddata)]);
cutinterval=zeros(size(winstart));
for chani=1:14
    for wini=1:length(winstart)
        EEGfilter.medianineachwin(chani,wini)=median(filtereddata(chani,winstart(wini):winend(wini)).^2);
        if median(filtereddata(chani,winstart(wini):winend(wini)).^2)>normalEEGthreshold(chani)
            filtereddataupdate(:,winstart(wini):winend(wini))=nan;
            cutinterval(wini)=1;
        end
    end
    cutintervali=[find(cutinterval==1) length(winstart)+1];
    if cutintervali(1)<=1/winLSpercentage+1
        filtereddataupdate(:,winstart(1):winend(cutintervali(1)))=nan;
    end
    for ii=1:length(cutintervali)-1
        if cutintervali(ii+1)-cutintervali(ii)-1<1/winLSpercentage+1
            filtereddataupdate(:,winstart(cutintervali(ii)):winend(min(cutintervali(ii+1),length(winstart))))=nan;
        end
    end
end

filtereddataupdate(:,isnan(filtereddataupdate(1,:)))=[];
EEGfilter.data=filtereddataupdate;
EEGfilter.pnts=size(filtereddataupdate,2);
EEGfilter.times(size(filtereddataupdate,2)+1:end)=[];
end
function normalEEGthreshold2 = outlierthreshold(eegdataraw,ParaSub,ParaDevice,baselineevent,reliablerange)
switch baselineevent
    case 'eyesopen'
        taski=1;threshold=nan;
    case 'eyesclosed'
        taski=2;threshold=nan;
end
normalEEGthreshold1=DataProcessingPipeline(eegdataraw,ParaSub,ParaDevice,threshold,nan,taski);
normalEEGthreshold2=DataProcessingPipeline(eegdataraw,ParaSub,ParaDevice,normalEEGthreshold1.FreqPOW.bandpow,reliablerange,taski);
end
function radarplot(data,RLim)
theta = deg2rad(107.5:25:432.5);
figure
polarplot(theta,data,'-*','LineWidth',2.5);%app.Pax,
ax=gca;
ax.RLim=RLim;
ax.LineWidth = 0.5;
ax.ThetaTick =([22.5 47.5 72.5 107.5:25:360 ]);
ax.ThetaTickLabel=({'F4','F8','AF4','AF3','F7','F3','FC5','T7','P7','O1','O2','P8','T8','FC6'});
end
function radarplotTRPMeanbetweenGroups(GroupCluster,GroupTRP,bandN,taskN,RLim,AddTitle)
taskname={'eyesopen','eyesclosed','game1','game2','game3','game4','game5',...
    'warm-up task','AU-Known','AU-New','Complete Figure','DwMT','Empty the Glass',...
    'Reading Stage','Selecting Stage','Sketching Stage'};
BandTag={'All','Theta','Alpha','Beta','Gamma'};
GroupTag={'A','B','C'};
GroupN = length(GroupCluster);

if GroupN>1
    for bandi=1:length(bandN)
        for taski=1:length(taskN)
            for groupi=1:GroupN
                groupMember=GroupCluster{1,groupi};
                Group(bandi).OutlierRemoved2plot(groupi,:)= nanmean(GroupTRP(taskN(taski)).bandN(bandN(bandi)).trp_OutlierRemoved(groupMember,:));%#ok
            end
            radarplot(Group(bandi).OutlierRemoved2plot,RLim)
            title([taskname{1,taskN(taski)} ' TRP - ' BandTag{1,bandN(bandi)} ' Band' AddTitle]);
            legend(GroupTag{1,1:GroupN});
        end
    end
else
    groupMember=GroupCluster{1,1};
    for bandi=1:length(bandN)
        for taski=1:length(taskN)
            Task(bandi).OutlierRemoved2plot(taski,:)= nanmean(GroupTRP(taskN(taski)).bandN(bandN(bandi)).trp_OutlierRemoved(groupMember,:));%#ok
        end
        radarplot(Task(bandi).OutlierRemoved2plot,RLim)
        title(['TRP - ' BandTag{1,bandN(bandi)} ' Band' AddTitle]);
        legend(taskname{1,taskN});
    end
end
end
function radarplotTRPMeanbetweenTasks(GroupCluster,GroupTRP,bandN,RLim,AddTitle)
taskname={'eyesopen','eyesclosed','game1','game2','game3','game4','game5',...
    'warm-up task','AU-Known','AU-New','Complete Figure','DwMT','Empty the Glass',...
    'Reading Stage','Selecting Stage','Sketching Stage'};
BandTag={'All','Theta','Alpha','Beta','Gamma'};
GroupTag={'A','B','C'};
GroupN = length(GroupCluster);
taskparis=[10 16 10 14 14 15;
           13 11 9 15 16 16];
for groupi=1:GroupN
    if GroupN>1
        addTitle=[addTitle ' ' GroupTag{1,groupi}];%#ok
    end
    groupMember=GroupCluster{1,groupi};
    for taskpairi=1:length(taskparis)
        for bandi=1:length(bandN)
            for taski=1:2
                Tasks(bandi).OutlierRemoved2plot(taski,:)= nanmean(GroupTRP(taskparis(taski,taskpairi)).bandN(bandN(bandi)).trp_OutlierRemoved(groupMember,:));%#ok
            end
            radarplot(Tasks(bandi).OutlierRemoved2plot,RLim)
            title([AddTitle 'TRP - ' BandTag{1,bandN(bandi)} ' Band']);
            legend(taskname{1,taskparis(:,taskpairi)});
        end
    end
end
end
% function MakeRadarPlots(Group,GroupTRP,bandN,RLim)
%% Profile
% taskN=8:16;
% radarplotTRPMeanbetweenGroups(Group.Cluster_Age,GroupTRP,bandN,taskN,RLim,' Clustered by Age')
% radarplotTRPMeanbetweenGroups(Group.Cluster_Gender,GroupTRP,bandN,taskN,RLim,' Clustered by Gender')
% % radarplotTRPMean(Group.Cluster_Hand,GroupTRP,bandN,taskN,RLim,' Clustered by Hand')
% radarplotTRPMeanbetweenGroups(Group.Cluster_Language,GroupTRP,bandN,taskN,RLim,' Clustered by Language')
%% Task 3 - TTCT
% radarplotTRPMeanbetweenGroups(Group.Cluster_3V,GroupTRP,bandN,11,RLim,' Clustered by Variety')
% radarplotTRPMeanbetweenGroups(Group.Cluster_3F,GroupTRP,bandN,11,RLim,' Clustered by Fluency')
% radarplotTRPMeanbetweenGroups(Group.Cluster_3N,GroupTRP,bandN,11,RLim,' Clustered by Novelty')
% %% Task 4 - Design with morphological table
% radarplotTRPMeanbetweenGroups(Group.Cluster_totalDwMT,GroupTRP,bandN,[12 14:16],RLim,' Clustered by Total Score') 
% radarplotTRPMeanbetweenGroups(Group.Cluster_LoD,GroupTRP,bandN,[12 14:16],RLim,' Clustered by Level of Detail')
% radarplotTRPMeanbetweenGroups(Group.Cluster_4V,GroupTRP,bandN,[12 14:16],RLim,' Clustered by Viability')
% %% Task 5 - Empty the Glass
% radarplotTRPMeanbetweenGroups(Group.Cluster_5V,GroupTRP,bandN,13,RLim,' Clustered by Variety')
% radarplotTRPMeanbetweenGroups(Group.Cluster_5F,GroupTRP,bandN,13,RLim,' Clustered by Fluency')
% radarplotTRPMeanbetweenGroups(Group.Cluster_5N,GroupTRP,bandN,13,RLim,' Clustered by Novelty')
% end
% function MakeRadarPlotsBetweenTasks(Group,GroupTRP,bandN,RLim)
% 
% radarplotTRPMeanbetweenTasks(Group.All,GroupTRP,bandN,RLim, 'Whole Group');
% 
% end
function [H,P]=statisticaltestsbetweengroups(GroupTRP,GroupCluster,bandN,taskN,GroupN,logscale)
% GroupTRP(14).bandN(3).trp_OutlierRemoved  
% taskname={'eyesopen','eyesclosed','game1','game2','game3','game4','game5',...
%     'warm-up task','AU-Known','AU-New','Complete Figure','DwMT','Empty the Glass',...
%     'Reading Stage','Selecting Stage','Sketching Stage'};
% BandTag={'All','Theta','Alpha','Beta','Gamma'};
% GroupTag={'A','B','C'};

% STEP 1 - ks test -> Normal Distribution?
% Kolmogorov-Smirnov test%0-yes normally distributed, 1-not normally distributed->Non parametric tests
% determine whether to check varience
for taski=1:length(taskN)
    thistask=taskN(taski);
    for chani=1:14
        for groupi=1:2
            thisgroup=GroupN(groupi);
            groupMember=GroupCluster{1,thisgroup};
            if logscale==1
                Group(groupi).Task(taski).Chan(chani).OutlierRemoved2plot= log(GroupTRP(thistask).bandN(bandN).trp_OutlierRemoved(groupMember,chani));%#ok
            else
                Group(groupi).Task(taski).Chan(chani).OutlierRemoved2plot= GroupTRP(thistask).bandN(bandN).trp_OutlierRemoved(groupMember,chani);%#ok
            end
            x=Group(groupi).Task(taski).Chan(chani).OutlierRemoved2plot;            
            x(isnan(x))=[];
            if ~isempty(x) && length(x)>1
                if std(x)~=0
                    x=(x-mean(x))/std(x);
                end
                [H.h(taski).chan(chani,groupi),P.p(taski).chan(chani,groupi)]=kstest(x);
            else
                 H.h(taski).chan(chani,groupi)=1;
                 P.p(taski).chan(chani,groupi)=nan;
            end
            % h = kstest(x) returns a test decision for the null hypothesis that the data in vector x comes from a standard normal distribution,
            % against the alternative that it does not come from such a distribution, using the one-sample Kolmogorov-Smirnov test.
            % The result h is 1 if the test rejects the null hypothesis at the 5% significance level, or 0 otherwis
        end
    end
end
%% STEP 2 (IF STEP 1 GETS 0 FOR BOTH GROUPS | normally distributed) - two-sample F-test
% to check if two groups come from normal distributions with same variances
% 0-same vaience 1-diff. vaience
% HOWEVER, WE GOT ALL 1 FROM STEP 1 -> Should use Heteroskedastic t-test
% NO NEED TO CHECK VARIENCE
%% STEP 3 (IF STEP 2 GETS 0 - Two-sample t-test, IF STEP 2 GETS 1)
%If Normally & same varience -> two-sample t-test
%0-same mean, no sta.sig.diff. between groups, 1-diff.mean, sta.sig.diff. between groups

%If not (Normally && same varience)
%Method 1 -> two-sample Kolmogorov-Smirnov test
%0-from same distribution, no sta.sig.diff. between groups, 1-diff. distribution, sta.sig.diff. between groups
%Method 2 -> two-sided Wilcoxon rank sum test
%0-from same distribution, no sta.sig.diff. between groups, 1-diff. distribution, sta.sig.diff. between groups
%We Hope to get 1 !!!
for taski= 1:length(taskN)    
    for chani=1:14
        x=Group(1).Task(taski).Chan(chani).OutlierRemoved2plot;
        y=Group(2).Task(taski).Chan(chani).OutlierRemoved2plot;
        x(isnan(x))=[];y(isnan(y))=[];
        if isempty(x)||isempty(y)
            H.h_var(chani,taski)=0; P.p_var(chani,taski)=nan;
            H.h_tt2(chani,taski)=0; P.p_tt2(chani,taski)=nan; 
            H.h_ks(chani,taski)=0; P.p_ks(chani,taski)=nan; 
            H.h_rk(chani,taski)=0; P.p_rk(chani,taski)=nan; 
        else
            [H.h_var(chani,taski),P.p_var(chani,taski),H.h_tt2(chani,taski),P.p_tt2(chani,taski),H.h_ks(chani,taski),P.p_ks(chani,taski),H.h_rk(chani,taski),P.p_rk(chani,taski)]=MeanEqual(x,y,H.h(taski).chan(chani,:));
        end
    end
end

end
function [H,P]=statisticaltestsbetweentasks(GroupTRP,Cluster,bandN,taskN,logscale)
% GroupTRP(14).bandN(3).trp_OutlierRemoved  
% taskname={'eyesopen','eyesclosed','game1','game2','game3','game4','game5',...
%     'warm-up task','AU-Known','AU-New','Complete Figure','DwMT','Empty the Glass',...
%     'Reading Stage','Selecting Stage','Sketching Stage'};
% BandTag={'All','Theta','Alpha','Beta','Gamma'};
% GroupTag={'A','B','C'};

% STEP 1 - ks test -> Normal Distribution?
% Kolmogorov-Smirnov test%0-yes normally distributed, 1-not normally distributed->Non parametric tests
% determine whether to check varience

for chani=1:14
    for taski=1:2       
        if logscale==1
            Task(taski).Chan(chani).OutlierRemoved2plot= log(GroupTRP(taskN(taski)).bandN(bandN).trp_OutlierRemoved(Cluster,chani));%#ok
        else
            Task(taski).Chan(chani).OutlierRemoved2plot= GroupTRP(taskN(taski)).bandN(bandN).trp_OutlierRemoved(Cluster,chani);%#ok
        end
        x=Task(taski).Chan(chani).OutlierRemoved2plot;
        x(isnan(x))=[];
        if ~isempty(x) && length(x)>1
            if mean(x)~=1 || std(x)~=0
                x=(x-mean(x))/std(x);
            end
            [H.h.chan(chani,taski),P.p.chan(chani,taski)]=kstest(x);
        else
            H.h.chan(chani,taski)=1;
            P.p.chan(chani,taski)=nan;
        end
        % h = kstest(x) returns a test decision for the null hypothesis that the data in vector x comes from a standard normal distribution,
        % against the alternative that it does not come from such a distribution, using the one-sample Kolmogorov-Smirnov test.
        % The result h is 1 if the test rejects the null hypothesis at the 5% significance level, or 0 otherwis
    end
end

%% STEP 2 (IF STEP 1 GETS 0 FOR BOTH GROUPS | normally distributed) - two-sample F-test
% to check if two groups come from normal distributions with same variances
% 0-same vaience 1-diff. vaience 
%% STEP 3 (IF STEP 2 GETS 0 - Two-sample t-test, IF STEP 2 GETS 1)
%If Normally & same varience -> two-sample t-test
%0-same mean, no sta.sig.diff. between groups, 1-diff.mean, sta.sig.diff. between groups

%If not (Normally && same varience)
%Method 1 -> two-sample Kolmogorov-Smirnov test
%0-from same distribution, no sta.sig.diff. between groups, 1-diff. distribution, sta.sig.diff. between groups
%Method 2 -> two-sided Wilcoxon rank sum test
%0-from same distribution, no sta.sig.diff. between groups, 1-diff. distribution, sta.sig.diff. between groups
%We Hope to get 1 !!!

for chani=1:14
    x=Task(1).Chan(chani).OutlierRemoved2plot;
    y=Task(2).Chan(chani).OutlierRemoved2plot;
    x(isnan(x))=[];y(isnan(y))=[];
    [H.h_var(chani),P.p_var(chani),H.h_tt2(chani),P.p_tt2(chani),H.h_ks(chani),P.p_ks(chani),H.h_rk(chani),P.p_rk(chani)]=MeanEqual(x,y,H.h.chan(chani,:));
end

end
function [h_var,p_var,h_tt2,p_tt2,h_ks,p_ks,h_rk,p_rk]=MeanEqual(x,y,Hnorm)% distinguish if the two groups are sta.sig.diff.

if sum(Hnorm)==0 && min(length(x),length(y))>1 % if ks test for both groups are 0 (both dist. normaly)
    [h_var,p_var]= vartest2(x,y);
    % returns a test decision for the null hypothesis that the data in vectors x and y comes from normal distributions with the same variance,
    % using the two-sample F-test. The alternative hypothesis is that they come from normal distributions with different variances.
    % result h is 1 if the test rejects the null hypothesis at the 5% significance level, and 0 otherwise.
    if h_var==0 % equal variance homoskedastic test
        [h_tt2,p_tt2] = ttest2(x,y,'Vartype','equal');
    else
        [h_tt2,p_tt2] = ttest2(x,y,'Vartype','unequal');
    end
    h_ks=h_tt2;p_ks=p_tt2;p_rk=p_ks;h_rk=h_ks;
    % h = ttest2(x,y) returns a test decision for the null hypothesis that the data in vectors x and y comes from independent random samples
    % from normal distributions with equal means and equal but unknown variances, using the two-sample t-test.
    % The alternative hypothesis is that the data in x and y comes from populations with unequal means.
    % The result h is 1 if the test rejects the null hypothesis at the 5% significance level, and 0 otherwise.
elseif sum(Hnorm)>0 && min(length(x),length(y))>1
    
    h_var=nan;p_var=nan;h_tt2=nan;p_tt2=nan;
    [h_ks,p_ks] = kstest2(x,y);
    % h = kstest2(x1,x2) returns a test decision for the null hypothesis that the data in vectors x1 and x2 are from the same continuous distribution,
    % using the two-sample Kolmogorov-Smirnov test. The alternative hypothesis is that x1 and x2 are from different continuous distributions.
    % The result h is 1 if the test rejects the null hypothesis at the 5% significance level, and 0 otherwise.
    [p_rk,h_rk] = ranksum(x,y);
    % p = ranksum(x,y) returns the p-value of a two-sided Wilcoxon rank sum test. ranksum tests the null hypothesis that data in x and y
    % are samples from continuous distributions with equal medians, against the alternative that they are not.
    % The test assumes that the two samples are independent. x and y can have different lengths. %This test is equivalent to a Mann-Whitney U-test.
else
    h_var=nan;p_var=nan;h_tt2=nan;p_tt2=nan;h_ks=nan;p_ks=nan;p_rk=nan;h_rk=nan;
end

end
% function radarplotEyesPOWMean(GroupEEGthreshold,bandN,RLim)
% BandTag={'All','Theta','Alpha','Beta','Gamma'};
% 
% for bandi=1:length(bandN)
%     EyeMove(bandi).OutlierRemoved2plot(1,:)= nanmean(removeOutlierbytrim(GroupEEGthreshold.eyesopen(bandN(bandi)).band));
%     EyeMove(bandi).OutlierRemoved2plot(2,:)=  nanmean(removeOutlierbytrim(GroupEEGthreshold.eyesclosed(bandN(bandi)).band));    
%     radarplot(EyeMove(bandi).OutlierRemoved2plot,RLim)
%     title(['POW - ' BandTag{1,bandN(bandi)} ' Band']);
%     legend({'Eyes Open','Eyes Closed'})
% end
% end
% function topoplotTask4StageMean(GroupAssignment,GroupTRP,bandN,colorbarlimit,Para)
% Stagetag={'Read','Select','Sketch'}; BandTag={'All','Theta','Alpha','Beta','Gamma'};
% for groupi=1:length(GroupAssignment)
%     plottitle=char('A'+groupi-1);
%     groupmember=GroupAssignment{1,groupi};
%     for bandi=bandN
%         figure
%         for stagei=1:3
%             subplot(1,3,stagei)
%             topoplot(nanmean(GroupTRP(stagei+3).bandN(bandi).trp_OutlierRemoved(groupmember,:)),Para.ChanLoc,'maplimits',colorbarlimit,'electrodes','labels')
%             title(['Stage ' Stagetag{1,stagei}]);
%         end        
%        annotation('textbox', [0.7, 0.2, 0.1, 0.1], 'String', ['Group ' plottitle '  ' BandTag{1,bandi} ' Band'])
%     end    
% end
% end
% function topoplotTaskMean(GroupAssignment,GroupTRP,taskN,bandN,colorbarlimit,Para)
% taskname={'eyesopen','eyesclosed','game1','game2','game3','game4','game5','task1','task2a','task2b','task3','task4','task5'};
% GroupN = length(GroupAssignment);
% for bandi=bandN
%     figure
%     for groupi=1:GroupN
%         plottitle=char('A'+groupi-1);
%         groupmember=GroupAssignment{1,groupi};
%         for taski=taskN
%             subplot(GroupN,length(taskN),taski+(GroupN-1)*length(taskN))
%             topoplot(nanmean(GroupTRP(taski).bandN(bandi).trp_OutlierRemoved(groupmember,:)),Para.ChanLoc,'maplimits',colorbarlimit,'electrodes','labels')
%             title(['Task ' taskname{1,taski} ' Group ' plottitle]);
%         end
%     end
% end
% end