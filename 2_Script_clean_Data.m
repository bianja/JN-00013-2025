% Script that loads recording data and temperature data and gives
% temperature informaation

clearvars
close all
cd \\132.187.28.171\home\rest\data\Analysis\matlab\scripts\temperature
% load table and sort animal numbers
load('Temperature_Data.mat')
load('PhysitempData.mat')
[~,index] = sortrows([AllAni.Animal].'); AllAni = AllAni(index); clear index

% set stimuli values
stim1.xvelo.w2 = [5.625 11.25 56.25 112.5 225 450 675 900 1120 1350];
stim1.xfreq.w2 = [.5 1 5 10 20 40 60 80 100 120];
stim2.xvelo.w2 = [5.625 11.25 56.25 112.5 225 450 675 900 1125 1350];
stim2.xfreq.w2 = [.5 1 5 10 20 40 60 80 100 120];

stim1.xvelo.w4 = [8.4375 16.875 84.375 168.75 337.5 675 1012.5 1350 1687.5 2025];
stim1.xfreq.w4 = [.38 .75 3.75 7.5 15 30 45 60 75 90];
stim2.xvelo.w4 = [11.25 22.5 112.5 225 450 900 1350 1800 2250 2700];
stim2.xfreq.w4 = [.5 1 5 10 20 40 60 80 100 120];

stim1.xvelo.w8 = [14.0625 28.125 140.625 281.25 562.5 1125 1687.5 2250 2812.5];
stim1.xfreq.w8 = [0.31 0.63 3.1 6.2 12.5 25 37.5 50 62.5];
stim2.xvelo.w8 = [22.5 45 225 450 900 1800 2700];
stim2.xfreq.w8 = [.5 1 5 10 20 40 60];

%% update stimuli information and remove not presented stimuli (120 Hz, 8x8)
for i = 1 : size(AllAni,2)
    if AllAni(i).Animal <= 136 % old stimuli
        AllAni(i).xfreq = stim1.xfreq;
        AllAni(i).xvelo = stim1.xvelo;
        AllAni(i).R01.yfreq.mean.translation.bw.w8(end) = [];
    else % new stimuli
        AllAni(i).xfreq = stim2.xfreq;
        AllAni(i).xvelo = stim2.xvelo;
    end
end

%% add temperature information
structNbr = [131 133];
Ani = 0;
for i = 55 : 59%62 %62 %size(AllAni,2)
    if AllAni(i).Animal == Ani % Ani for current unit has not changed to previous one
        if ChE == 1
            for k = 1 : num % assign saved temperature information
%                 clear(char(AllAni(i).(['R0',num2str(k)]).Temp))
                AllAni(i).(['R0',num2str(k)]).Temp = AllAni(i-1).(['R0',num2str(k)]).Temp;
                AllAni(i).StimTimes.timesStart = timesStart';
                AllAni(i).StimTimes.timesEnd = timesEnd';
                AllAni(i).StimTimes.Type = stimFilesOrder'; % 1 = RF, 2 = OF, 3 = CS
            end
        else
            structNbr(end) = i;
        end
    else % Ani has changed; load new animal data
        clear filesMat files OF stim1 stim2 TData temp Tmax Tmin Tmed folders stimFilesOrder timesStart timesEnd
        ChE = 1; % set to 0 if Channel with temperature information is unavailable to reduce computation time
        if AllAni(i).Animal <= 118
            cd \\132.187.28.171\home\data118\
        else
            cd \\132.187.28.171\home\rest\data\Ephys\
        end
        folders = dir; % list all folders and change path to folder of current animal
        for k = 1 : size(folders,1)
            if contains(folders(k).name,['A',num2str(AllAni(i).Animal)])
                cd(folders(k).name)
                cd(AllAni(i).File(strfind(AllAni(i).File,'Rec'):strfind(AllAni(i).File,'Rec')+11))
            end
        end
        files = dir; % list all files
        for k = 1 : 9
            if ~isempty(AllAni(i).(['R0',num2str(k)]))
                num = k; % tells about number of recordings of the current unit
            end
        end
        filesMat = dir('*mat*');
        for k = 1 : size(filesMat,1) % correct german month labeling to english
            if ~isempty(strfind(filesMat(k).date,'Mrz'))
                filesMat(k).date(strfind(filesMat(k).date,'Mrz'):strfind(filesMat(k).date,'Mrz')+2) = 'Mar';
            elseif ~isempty(strfind(filesMat(k).date,'Mai'))
                filesMat(k).date(strfind(filesMat(k).date,'Mai'):strfind(filesMat(k).date,'Mai')+2) = 'May';
            elseif ~isempty(strfind(filesMat(k).date,'Okt'))
                filesMat(k).date(strfind(filesMat(k).date,'Okt'):strfind(filesMat(k).date,'Okt')+2) = 'Oct';
            elseif ~isempty(strfind(filesMat(k).date,'Dez'))
                filesMat(k).date(strfind(filesMat(k).date,'Dez'):strfind(filesMat(k).date,'Dez')+2) = 'Dec';
            end
        end
        c = 1;
        for k = 1 : size(filesMat,1)
            if contains(filesMat(k).name,'test_Arena')
                timesEnd(c) = datetime(filesMat(k).date);
                timesStart(c) = timesEnd(c) - minutes(4);
                stimFilesOrder(c,1) = 1;
                c = c + 1; ho = k;
            elseif contains(filesMat(k).name,'Arena') && ~contains(filesMat(k).name,'test')
                timesEnd(c) = datetime(filesMat(k).date);
                timesStart(c) = timesEnd(c) - minutes(28);
                stimFilesOrder(c,1) = 2;
                c = c + 1; ho = k;
            elseif contains(filesMat(k).name,'Continuous')
                if contains(filesMat(ho+2).name,'test_Arena')
                    timesEnd(c) = datetime(filesMat(ho+2).date) - minutes(4);
                elseif contains(filesMat(ho+2).name,'Arena') && ~contains(filesMat(ho+2).name,'test')
                    timesEnd(c) = datetime(filesMat(ho+2).date) - minutes(28);
                end
                timesStart(c) = datetime(filesMat(ho).date);
                stimFilesOrder(c,1) = 3;
                c = c + 1; ho = k;
            end
        end
        if ~isempty(dir('*xls*')) % temperature tracked by multimeter, excel file exits
            % load excel file
            temp = dir('*xls*');
            TData = readtable(temp.name);
            % save stim start and end times and order of stimulus type in struct
            AllAni(i).StimTimes.timesStart = timesStart';
            AllAni(i).StimTimes.timesEnd = timesEnd';
            AllAni(i).StimTimes.Type = stimFilesOrder'; % 1 = RF, 2 = OF, 3 = CS
            for k = 1 : num
                OF = find(stimFilesOrder == 2);
                index = [0 0];
                [~,index(1)] = ismember(timesStart(OF(k)),TData.Date_Time);
                [~,index(2)] = ismember(timesEnd(OF(k)),TData.Date_Time);
                if index(1) == 0
                    in = 1;
                    while index(1) == 0
                        [~,index(1)] = ismember(timesStart(OF(k))+seconds(in),TData.Date_Time);
                        in = in + 1;
                    end
                end
                if index(2) == 0
                    in = 1;
                    while index(2) == 0
                        [~,index(2)] = ismember(timesEnd(OF(k))+seconds(1),TData.Date_Time);
                        if Ani == 147 % Batterie leer, Temperature zum Ende nicht mehr aufgenommen
                            index(2) = length(TData.Date_Time);
                        end
                        in = in + 1;
                    end
                end
                Tmed(k) = median(str2num(cell2mat(TData.Value(index(1):index(2)))))+1.0222;
                Tmin(k) = min(str2num(cell2mat(TData.Value(index(1):index(2)))))+1.0222;
                Tmax(k) = max(str2num(cell2mat(TData.Value(index(1):index(2)))))+1.0222;
                clear(char(AllAni(i).(['R0',num2str(k)]).Temp))
                AllAni(i).(['R0',num2str(k)]).Temp = [Tmin(k) Tmed(k) Tmax(k)];
            end
        elseif ~isnan(AllAni(i).R01.Temp) % temperature tracked by Physitemp
            clearvars -except AllAni ChE structNbr i filesMat num stimFilesOrder timesStart timesEnd fitT fitmV
            % temperature calibration
            % load files and extract min, max and median temperature data
            % separate stim and spike files
            stimnum = 1; spikenum = 1;
            for k = 1 : size(filesMat,1)
                if contains(filesMat(k).name,'Arena') || contains(filesMat(k).name,'Continuous') % stim files
                    stimC (stimnum) = k;
                    stimnum = stimnum + 1;
                else % spike files
                    spikeC(spikenum) = k;
                    spikenum = spikenum + 1;
                end
            end
            % find start and end times of stimulus sets
            for k = 1 : length(spikeC)
                if contains(filesMat(spikeC(k)).name,'R0') && ChE == 1
                    load(filesMat(spikeC(k)).name);
                    if exist('Ch17','var')
                        Tmed(k) = round(fitT(find(round(fitmV) == round(median(Ch17.values)*10000),1,'first')),1);
                        Tmin(k) = round(fitT(find(round(fitmV) == round(min(Ch17.values)*10000),1,'first')),1);
                        Tmax(k) = round(fitT(find(round(fitmV) == round(max(Ch17.values)*10000),1,'first')),1);
                        AllAni(i).(['R0',num2str(k)]).Temp = [Tmin(k) Tmed(k) Tmax(k)];
                    else
                        ChE = 0;
                        structNbr(end) = i; % contains line number of units without temperature information
                    end
                end
            end
            AllAni(i).StimTimes.timesStart = timesStart';
            AllAni(i).StimTimes.timesEnd = timesEnd';
            AllAni(i).StimTimes.Type = stimFilesOrder'; % 1 = RF, 2 = OF, 3 = CS
        else % no temperature was tracked
            structNbr(end) = i; % contains line number of units without temperature information
        end
        Ani = AllAni(i).Animal;
    end
end
if exist('structNbr') % display all units that needs to get checked on
    disp(['Check on temperature data of units: ',num2str(structNbr)])
end

%% calc mean in case temperatures were recorded with repetitions
% clearvars
% cd \\132.187.28.171\home\rest\data\Analysis\matlab\scripts\temperature
% load('Temperature_Data_Temperature_Zwischenspeicher_backup7.mat')
for i = 5 : 7 %1 : size(AllAni,2)
    if ~isnan(AllAni(i).R01.Temp)
        clear tempC
        add = 1; c = 1;
        for k = 1 : 9
            if ~isempty(AllAni(i).(['R0',num2str(k)]))
                posOF = k;
            end
        end
%         posOF = find(AllAni(i).StimTimes.Type == 2);
        for k = 1 : posOF-1
            if abs(AllAni(i).(['R0',num2str(k)]).Temp(2) - AllAni(i).(['R0',num2str(k+1)]).Temp(2)) < 3
                tempC.(['var',num2str(add)])(c) = k;
                c = c + 1;
            else
                tempC.(['var',num2str(add)])(c) = k;
                add = add + 1;
                c = 1;
            end
%             if posOF(k)+1 == posOF(k+1)
%                 tempC.(['var',num2str(add)])(c) = k;
%                 c = c + 1;
%             elseif k == 1
%                 tempC.(['var',num2str(add)])(c) = k;
%                 add = add + 1;
%                 c = 1;
%             elseif posOF(k)-1 == posOF(k-1)
%                 tempC.(['var',num2str(add)])(c) = k;
%                 add = add + 1;
%                 c = 1;
%             else
%                 tempC.(['var',num2str(add)])(c) = k;
%                 add = add + 1;
%                 c = 1;
%             end
        end
        tempC.(['var',num2str(add)])(c) = k+1;
        tempAni = AllAni(i);
        for k = 1 : posOF
            tempAni.(['R0',num2str(k)]) = [];
        end
        clear meanR
        for k = 1 : size(struct2table(tempC),2)
            pos = tempC.(['var',num2str(k)]);
            if length(pos) > 1
                for j = 1 : length(pos)
                    meanR(j) = AllAni(i).(['R0',num2str(pos(j))]);
                end
                tempAni.(['R0',num2str(k)]) = mean_temperature_repetitions(meanR);
                tempAni.(['R0',num2str(k)]).rep = j;
                tempAni.rep = 0;
            else
                tempAni.(['R0',num2str(k)]) = AllAni(i).(['R0',num2str(k)]);
            end
        end
        AllAni(i) = tempAni;
    end
end
%%
for i = 1 : size(AllAni,2)
    r1(i) = AllAni(i).R01.Temp(2);
    r2(i) = AllAni(i).R02.Temp(2);
    if ~isempty(AllAni(i).R03)
        r3(i) = AllAni(i).R03.Temp(2);
    else
        r3(i) = NaN;
    end
    if ~isempty(AllAni(i).R04)
        r4(i) = AllAni(i).R04.Temp(2);
    else
        r4(i) = NaN;
    end
    if ~isempty(AllAni(i).R05)
        r5(i) = AllAni(i).R05.Temp(2);
    else
        r5(i) = NaN;
    end
    if ~isempty(AllAni(i).R06)
        r6(i) = AllAni(i).R06.Temp(2);
    else
        r6(i) = NaN;
    end
end
table = [r1' r2' r3' r4' r5' r6'];

%% delete data points that were recorded with wrong stim and add a column with visual resolution values
for i = 1 : size(AllAni,2) % start value = 46
    if AllAni(i).Animal < 122
        AllAni(i).VR.fw = 2.8125;
        AllAni(i).VR.bw = 2.8125;
    else
        for j = 1 : 9
            if ~isempty(AllAni(i).(['R0',num2str(j)]))
                AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.w2 = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
            end
        end
        AllAni(i).VR.fw = 2.8125*2;
        AllAni(i).VR.bw = 2.8125;
    end
end

for i = 1 : size(AllAni,2)
    for j = 1 : 9
        if ~isempty(AllAni(i).(['R0',num2str(j)]))
            if AllAni(i).Animal < 136
                if AllAni(i).Animal == 80
                    AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.w8 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.w8(1:8) NaN];
                    AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.w8 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.w8(1:8) NaN];
                    AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.cw.w8 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.cw.w8(1:8) NaN];
                    AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.ccw.w8 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.ccw.w8(1:8) NaN];
                    AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.w4 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.w4(1:9)];
                    AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.w4 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.w4(1:9)];
                    AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.cw.w4 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.cw.w4(1:9)];
                    AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.ccw.w4 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.ccw.w4(1:9)];
                    AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.w2 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.w2(1:9)];
                    AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.w2 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.w2(1:9)];
                    AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.cw.w2 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.cw.w2(1:9)];
                    AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.ccw.w2 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.ccw.w2(1:9)];
        display('test')
                else
                    AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.w8(10) = NaN;
                    AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.w8(10) = NaN;
                    AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.cw.w8(10) = NaN;
                    AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.ccw.w8(10) = NaN;
                end
            else
                AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.w8(8:10) = NaN;
                AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.w8(8:10) = NaN;
                AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.cw.w8(8:10) = NaN;
                AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.ccw.w8(8:10) = NaN;
            end
        end
    end
end

%% correct x vector (add NaN)
for i = 1 : size(AllAni,2)
    if AllAni(i).Animal < 137
        AllAni(i).xfreq.w8 = [AllAni(i).xfreq.w8 NaN];
        AllAni(i).xvelo.w8 = [AllAni(i).xvelo.w8 NaN];
    else
        AllAni(i).xfreq.w8 = [AllAni(i).xfreq.w8 NaN NaN NaN];
        AllAni(i).velo.w8 = [AllAni(i).xvelo.w8 NaN NaN NaN];
    end
end

%% correct A80
for i = 4 : 7
    for j = 1 : 6
        AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.w2 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.w2];
        AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.w2 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.w2];
        AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.cw.w2 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.cw.w2];
        AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.ccw.w2 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.ccw.w2];
        
        AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.w4 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.w4];
        AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.w4 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.w4];
        AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.cw.w4 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.cw.w4];
        AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.ccw.w4 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.ccw.w4];
        
        AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.w8 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.w8];
        AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.w8 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.w8];
        AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.cw.w8 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.cw.w8];
        AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.ccw.w8 = [NaN AllAni(i).(['R0',num2str(j)]).yfreq.mean.rotation.ccw.w8];
    end
end

%% fit curve on data points and save new x and y data
i = 129;
close all
figure
x = [.5 1 5 10 20 40 60 80 100 120];
for i = 1 : 101
    y = Data(i).R02.yfreq.mean.translation.bw.w8;
    f = createFit(x,y); % Interpolant - Shape preservin (PCHIP)
    X = linspace(min(x), max(x),1000); % get x values of fit
    Y = f(X); % get y values of fit
    subplot(8,14,i)
%     plot(f,x,y)
    plot(X,Y,'-r')
    hold on
    plot(x,y,'ro')
% plot(x,y,'--k')
    set(gca,'xscale','log','xticklabels',{})
    set(gcf,'position',[210 300 370 330])
end

%%
cd \\132.187.28.171\home\rest\data\Analysis\matlab\scripts\temperature
