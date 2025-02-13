%% Plot data
% 
% line 85: change value for type for changing the response type
% 2 = bidir shift, 15 = bidir no shift, 8 = unidir shift btf, 9 = unidir shift ftb
%

cd \\132.187.28.171\home\rest\data\Analysis\matlab
clearvars
close all

load('Temperature_Data_correct_backup4.mat')

% to plot bidir shift units that have a non-detectable peak for at least 
% one condition: 1. set type = 2, 2. run the following two lines of code 
% (typeC can be assigned to the values 7, 8, 9):
% x = [AllAni.type];
% AllAni = AllAni(find(x == typeC));

Data = AllAni;
DataBU = Data;

%% differentiate types
% peak is shifted for both directions, peak is shifted for only on
% direction, peak is not shifted for both direcctions
load('XvalRange.mat')
Xval = X;
Data = DataBU;
for k = 1 : 2
%     clear TCmin TCmax 
    if k == 1
        dir = 'ftb';
        dir2 = 'fw';
    elseif k == 2
        dir = 'btf';
        dir2 = 'bw';
    end
    for  i = 1 : size(Data,2)
        Temp(1:6,1:2) = NaN;
        for j = Data(i).first : Data(i).last
            Temp(j,1) = Data(i).(['R0',num2str(j)]).Temp(2);
            Temp(j,2) = Data(i).peak.(char(dir)).loc(j);
            Temp(j,3) = Data(i).peak.(char(dir)).pk(j);
            Temp(j,4) = Data(i).(['R0',num2str(j)]).background.sum(1);
        end
        mit = min(Temp(:,1));
        mat = max(Temp(:,1));
        Tval(i,1) = mit(1);
        Tval(i,2) = mat(1);
        smit = Temp(find(Temp(:,1) == mit(1)),2);
        smat = Temp(find(Temp(:,1) == mat(1)),2);
        Spikeval.(char(dir))(i,1) = smit(1);
        Spikeval.(char(dir))(i,2) = smat(1);
        smit = Temp(find(Temp(:,1) == mit(1)),3);
        smat = Temp(find(Temp(:,1) == mat(1)),3);
        Spikestreng.(char(dir))(i,1) = smit(1);
        Spikestreng.(char(dir))(i,2) = smat(1);
        smit = Temp(find(Temp(:,1) == mit(1)),3);
        smat = Temp(find(Temp(:,1) == mat(1)),3);
        Baseline(i,1) = smit(1);
        Baseline(i,2) = smat(1);
        Tmi = find(Temp(:,1) == mit(1));
        Tma = find(Temp(:,1) == mat(1));
        ymin(i,:) = Data(i).(['R0',num2str(Tmi(1))]).yfreq.mean.translation.(char(dir2)).w8;
        ymax(i,:) = Data(i).(['R0',num2str(Tma(1))]).yfreq.mean.translation.(char(dir2)).w8;
        x = Data(i).xvelo.w8;
        x = x((~isnan(ymax(i,:))));
        f = createFit(x,ymin(i,(~isnan(ymin(i,:))))); % Interpolant - shape preserving (PCHIP)
        X = linspace(22.5, 2700, 1000); % get x values
        Ymin = f(X); % get y values
        f = createFit(x,ymax(i,(~isnan(ymax(i,:))))); % Interpolant - shape preserving (PCHIP)
        X = linspace(22.5, 2700, 1000); % get x values
        Ymax= f(X); % get y values
        maxval = max([Ymin; Ymax]);
        TCmin.(char(dir))(i,:) = Ymin/maxval(1);
        TCmax.(char(dir))(i,:) = Ymax/maxval(1);
    end
    for i = 1 : size(Spikeval.(char(dir)),1)
        if Spikeval.(char(dir))(i,1) < Spikeval.(char(dir))(i,2) % für warme T ist das Sensitivitätsmaximum höher
            sensRes(i,k) = 1;
        elseif Spikeval.(char(dir))(i,1) == Spikeval.(char(dir))(i,2) % Sensitivitätsmaximum ändert sich nicht
            sensRes(i,k) = 6+k;%2
        elseif Spikeval.(char(dir))(i,1) > Spikeval.(char(dir))(i,2) % table1(i,1) > table1(i,2) % für warme T ist das Sensitivit#tsmaximum niedriger
            sensRes(i,k) = 4;%0;
        else
            sensRes(i,k) = 1; % these are larger for high T as well, but dont have peak
        end
    end
end
sensRes = sum(sensRes,2);

% find units
close all
clear Q10
% 2 = bidir shift, 15 = bidir no shift, 8 = unidir shift btf, 9 = unidir shift ftb
type = 2;
if type == 2
    savename = char([cd,'\fig_manuscriptOF\ExcBi_BiDirShift']);
elseif type == 15
    savename = char([cd,'\fig_manuscriptOF\ExcBi_BiDirNoShift']);
elseif type == 8
    savename = char([cd,'\fig_manuscriptOF\ExcBi_UniDirShiftBTF']);
elseif type == 9
    savename = char([cd,'\fig_manuscriptOF\ExcBi_UniDirShiftFTB']);
elseif type == 0
    savename = char([cd,'\fig_manuscriptOF\ExcBi_BiDirShift_High']);
else
    error('wrong type')
end
if type == 0
    plotpos = 1:size(Data,2);
else
    plotpos = find(sensRes == type);
end
Data = DataBU(plotpos);
%% plot
for k = 1 : 2
    if k == 1
        dir = 'ftb';
        dir2 = 'fw';
        if type == 0 && typeC == 7
            bpt = Spikeval.(char(dir))(plotpos,:);
            savebpt = bpt/0.0222;
            save highPeakftb_sensMax savebpt
        elseif type == 0 && typeC == 9
            bpt = Spikeval.(char(dir))(plotpos,:);
        else
            bpt = Spikeval.(char(dir))(plotpos,:);
            bpt(find(isnan(bpt(:,2)) == 1),1) = NaN;
            bpt(find(isnan(bpt(:,1)) == 1),2) = NaN;
            savebpt = bpt/0.0222;
            save lowPeakftb_sensMax savebpt
        end
    elseif k == 2
        dir = 'btf';
        dir2 = 'bw';
        if type == 0 && typeC == 8
            bpt = Spikeval.(char(dir))(plotpos,:);
            savebpt = bpt/0.0222;
            save highPeakbtf_sensMax savebpt
        elseif type == 0 && typeC == 9
            bpt = Spikeval.(char(dir))(plotpos,:);
        else
            bpt = Spikeval.(char(dir))(plotpos,:);
            bpt(find(isnan(bpt(:,2)) == 1),1) = NaN;
            bpt(find(isnan(bpt(:,1)) == 1),2) = NaN;
            savebpt = bpt/0.0222;
            save lowPeakbtf_sensMax savebpt
        end
    end
    % bpt = Spikeval.(char(dir))(plotpos,:);
    % bpt(find(isnan(bpt(:,2)) == 1),1) = NaN;
    % bpt(find(isnan(bpt(:,1)) == 1),2) = NaN;
    Boxplot_B(bpt/0.0222,2,15,[18/255 136/255 215/255; 170/255 32/255 55/255],{[num2str(mean(Ttemp(plotpos,1))),'+-',num2str(std(Ttemp(plotpos,1)))],[num2str(mean(Ttemp(plotpos,2))),'+-',num2str(std(Ttemp(plotpos,2)))]},[1 2])
    writematrix(bpt(:,1)/0.0222,['SensMax_',char(dir),'_L.txt']) % for stats with R
    writematrix(bpt(:,2)/0.0222,['SensMax_',char(dir),'_H.txt']) % for stats with R
    xlabel('temperature (°C)')
    title(dir)
    ylabel('sensitivity maximum (Hz)')
    ylim([0 3000])
    hold on
%     plot([0 7],[60 60],'--','Color',[.7 .7 .7 .5],'LineWidth',2)
%     plot([0 7],[90 90],'--','Color',[.7 .7 .7 .5],'LineWidth',1.5)
%     plot([0 7],[120 120],'--','Color',[.7 .7 .7 .5],'LineWidth',1)
    for j = 1 : 2
        x = ones(size(bpt(:,j)/0.0222)).*(1+(rand(size(bpt(:,j)/0.0222))*2-1)/10)+j-1;
        scatter(x,bpt(:,j)/0.0222,10,'filled','MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',1)
    end
    set(gcf,'Position',[820 450 320 330])
    
    % statistics
    [pSM{k},hSM{k},statsSM{k}] = signrank(bpt(:,1)/0.0222,bpt(:,2)/0.0222);
    text(1,2750,['W-test:',num2str(pSM{k})])
    print([char(savename),'_bp_sensMax_',dir],'-dsvg','-r300','-painters')
    savefig([char(savename),'_bp_sensMax_',dir,'.fig'])

    % Line plot
    figure
    hold on
    for i = 1 : length(plotpos)
        plot(1:2,[Spikeval.(char(dir))(plotpos(i),1)/0.0222,Spikeval.(char(dir))(plotpos(i),2)/0.0222],'-','LineWidth',1.5,'Color',[.7 .7 .7])
    end
    tempplot = Spikeval.(char(dir))(plotpos,1:2);
    tempplot([find(isnan(Spikeval.(char(dir))(plotpos,1))); find(isnan(Spikeval.(char(dir))(plotpos,2)))],1:2) = NaN;
    plot(1:2,[nanmedian(tempplot(:,1)/0.0222),nanmedian(tempplot(:,2)/0.0222)],'-','LineWidth',2.5,'Color',[.3 .3 .3])
    title(dir)
    xlim([0.75 2.25])
    set(gca,'XTick',[1 2],'XTickLabels',{'23','32'})
    xlabel('temperature (°C)')
    ylabel('preferred temporal frequency (Hz)')
    box on
    axis square
    set(gcf,'position',[830 50 370 295])
    ylim([0 3000])
    print([char(savename),'_lp_sensMax_',dir],'-dsvg','-r300','-painters')
    savefig([char(savename),'_lp_sensMax_',dir,'.fig'])


    % RespStreng
    if k == 1
        dir = 'ftb';
        dir2 = 'fw';
        if type == 0 && typeC == 7
            bpt = Spikestreng.(char(dir))(plotpos,:);
            savebpt = bpt;
            save highPeakftb_respStreng savebpt
        elseif type == 0 && typeC == 9
            bpt = Spikestreng.(char(dir))(plotpos,:);
        else
            bpt = Spikestreng.(char(dir))(plotpos,:);
            bpt(find(isnan(bpt(:,2)) == 1),1) = NaN;
            bpt(find(isnan(bpt(:,1)) == 1),2) = NaN;
            savebpt = bpt;
            save lowPeakftb_respStreng savebpt
        end
    elseif k == 2
        dir = 'btf';
        dir2 = 'bw';
        if type == 0 && typeC == 8
            bpt = Spikestreng.(char(dir))(plotpos,:);
            savebpt = bpt;
            save highPeakbtf_respStreng savebpt
        elseif type == 0 && typeC == 9
            bpt = Spikestreng.(char(dir))(plotpos,:);
        else
            bpt = Spikestreng.(char(dir))(plotpos,:);
            bpt(find(isnan(bpt(:,2)) == 1),1) = NaN;
            bpt(find(isnan(bpt(:,1)) == 1),2) = NaN;
            savebpt = bpt;
            save lowPeakbtf_respStreng savebpt
        end
    end
    % bpt = Spikestreng.(char(dir))(plotpos,:);
    % bpt(find(isnan(bpt(:,2)) == 1),1) = NaN;
    % bpt(find(isnan(bpt(:,1)) == 1),2) = NaN;
    Boxplot_B(bpt,2,15,[18/255 136/255 215/255; 170/255 32/255 55/255],{[num2str(mean(Ttemp(plotpos,1))),'+-',num2str(std(Ttemp(plotpos,1)))],[num2str(mean(Ttemp(plotpos,2))),'+-',num2str(std(Ttemp(plotpos,2)))]},[1 2])
    writematrix(bpt(:,1),['RespAmpl_',char(dir),'_L.txt']) % for stats with R
    writematrix(bpt(:,2),['RespAmpl_',char(dir),'_H.txt']) % for stats with R
    xlabel('temperature (°C)')
    title(dir)
    ylabel('response strength (s-1)')
    ylim([0 500])
    hold on
    for j = 1 : 2
        x = ones(size(bpt(:,j))).*(1+(rand(size(bpt(:,j)))*2-1)/10)+j-1;
        scatter(x,bpt(:,j),10,'filled','MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',1)
    end
    set(gcf,'Position',[820 450 320 330])
    
    % statistics
    [pRS{k},hRS{k},statsRS{k}] = signrank(bpt(:,1),bpt(:,2));
    text(1,480,['W-test:',num2str(pRS{k})])
    print([char(savename),'_bp_respStren_',dir],'-dsvg','-r300','-painters')
    savefig([char(savename),'_bp_respStren_',dir,'.fig'])

    % Line plot
    figure
    hold on
    for i = 1 : length(plotpos)
        plot(1:2,[Spikestreng.(char(dir))(plotpos(i),1),Spikestreng.(char(dir))(plotpos(i),2)],'-','LineWidth',1.5,'Color',[.7 .7 .7])
    end
    tempplot = Spikestreng.(char(dir))(plotpos,1:2);
    tempplot([find(isnan(Spikestreng.(char(dir))(plotpos,1))); find(isnan(Spikestreng.(char(dir))(plotpos,2)))],1:2) = NaN;
    plot(1:2,[nanmedian(tempplot(:,1)),nanmedian(tempplot(:,2))],'-','LineWidth',2.5,'Color',[.3 .3 .3])
    title(dir)
    xlim([0.75 2.25])
    set(gca,'XTick',[1 2],'XTickLabels',{'23','32'})
    xlabel('temperature (°C)')
    ylabel('response strength (s-1)')
    box on
    axis square
    set(gcf,'position',[830 50 370 295])
    ylim([0 500])
    print([char(savename),'_lp_respStren_',dir],'-dsvg','-r300','-painters')
    savefig([char(savename),'_lp_respStren_',dir,'.fig'])

    % Line plot
    figure
    hold on
    for i = 1 : length(plotpos)
        plot(1:2,[Baseline(plotpos(i),1),Baseline(plotpos(i),2)],'-','LineWidth',1.5,'Color',[.7 .7 .7])
    end
    tempplot = Baseline(plotpos,1:2);
    tempplot([find(isnan(Baseline(plotpos,1))); find(isnan(Baseline(plotpos,2)))],1:2) = NaN;
    plot(1:2,[nanmedian(tempplot(:,1)),nanmedian(tempplot(:,2))],'-','LineWidth',2.5,'Color',[.3 .3 .3])
    xlim([0.75 2.25])
    set(gca,'XTick',[1 2],'XTickLabels',{'23','32'})
    xlabel('temperature (°C)')
    ylabel('baseline (s-1)')
    box on
    axis square
    set(gcf,'position',[830 50 370 295])
    ylim([0 500])
    print([char(savename),'_lp_baseline_',dir],'-dsvg','-r300','-painters')
    savefig([char(savename),'_lp_baseline_',dir,'.fig'])

    
    % Q10 val
    for i = 1 : length(plotpos)
        Q10.SM(i,k) = ((Spikeval.(char(dir))(plotpos(i),2))/(Spikeval.(char(dir))(plotpos(i),1)))^(10/(Ttemp(plotpos(i),2)-Ttemp(plotpos(i),1)));
        Q10.RS(i,k) = (Spikestreng.(char(dir))(plotpos(i),2)/Spikestreng.(char(dir))(plotpos(i),1))^(10/(Ttemp(plotpos(i),2)-Ttemp(plotpos(i),1)));
        Q10.B(i,k) = (Baseline(plotpos(i),2)/Baseline(plotpos(i),1))^(10/(Ttemp(plotpos(i),2)-Ttemp(plotpos(i),1)));
    end

    % plot hot and cold
    figure
    hold on
    if k == 1
        dir = 'ftb';
        dir2 = 'fw';
        if type == 0 && typeC == 7
            bpt = Spikestreng.(char(dir))(plotpos,:);
        elseif type == 0 && typeC == 9
            bpt = Spikestreng.(char(dir))(plotpos,:);
        else
            bpt = Spikestreng.(char(dir))(plotpos,:);
            bpt(find(isnan(bpt(:,2)) == 1),1) = NaN;
            bpt(find(isnan(bpt(:,1)) == 1),2) = NaN;
            TCmin.(char(dir))(plotpos(find(isnan(bpt(:,1)) == 1)),:) = NaN;
            TCmax.(char(dir))(plotpos(find(isnan(bpt(:,1)) == 1)),:) = NaN;
        end
    elseif k == 2
        dir = 'btf';
        dir2 = 'bw';
        if type == 0 && typeC == 8
            bpt = Spikestreng.(char(dir))(plotpos,:);
        elseif type == 0 && typeC == 9
            bpt = Spikestreng.(char(dir))(plotpos,:);
        else
            bpt = Spikestreng.(char(dir))(plotpos,:);
            bpt(find(isnan(bpt(:,2)) == 1),1) = NaN;
            bpt(find(isnan(bpt(:,1)) == 1),2) = NaN;
            TCmin.(char(dir))(plotpos(find(isnan(bpt(:,1)) == 1)),:) = NaN;
            TCmax.(char(dir))(plotpos(find(isnan(bpt(:,1)) == 1)),:) = NaN;
        end
    end
    % bpt = Spikestreng.(char(dir))(plotpos,:);
    % bpt(find(isnan(bpt(:,2)) == 1),1) = NaN;
    % bpt(find(isnan(bpt(:,1)) == 1),2) = NaN;
    % TCmin.(char(dir))(plotpos(find(isnan(bpt(:,1)) == 1)),:) = NaN;
    % TCmax.(char(dir))(plotpos(find(isnan(bpt(:,1)) == 1)),:) = NaN;
    plot_distribution_prctile(X,TCmin.(char(dir))(plotpos,:),'Color',[18/255 136/255 215/255],'LineWidth',1.5,'Alpha',0.15,'Prctile',[50])
    plot_distribution_prctile(X,TCmax.(char(dir))(plotpos,:),'Color',[170/255 32/255 55/255],'LineWidth',1.5,'Alpha',0.15,'Prctile',[50])
    ylabel('spike rate (s^-^1)')
    xlabel('velocity(°/s)')
    title(dir)
    axis square
    set(gcf,'Position',[820 400 320 330])
    set(gca,'xscale','log','XTick',[10,100,1000],'XTickLabels',{'10','100','1000'})
    ylim([0 1])
%     xlim([10^1.1 10^3.6])
    xlim([X(1) X(end)])
    box on
    % some colors (b,b,r) [103/255 152/255 192/255]; [37/255 123/255 195/255]; [163/255 0 2/255]
    % print([pwd '/figsave/',num2str(respType),'_meanTC_',char(dir),'_w',num2str(wid)],'-dsvg','-r300','-painters')   
    print([char(savename),'_TC_Med25P_',dir],'-dsvg','-r300','-painters')
    savefig([char(savename),'_TC_Med25P_',dir,'.fig'])

    
    figure
    hold on
    % plot([10 10],[0 1],'--','Color',[.7 .7 .7 .5],'LineWidth',1.5)
    % plot([20 20],[0 1],'--','Color',[.7 .7 .7 .5],'LineWidth',1.5)
    % plot([40 40],[0 1],'--','Color',[.7 .7 .7 .5],'LineWidth',1.5)
    for i = 1 : length(plotpos)
        plot(X,TCmin.(char(dir))(plotpos(i),:),'Color',[18/255 136/255 215/255 .2],'LineWidth',1)
        plot(X,TCmax.(char(dir))(plotpos(i),:),'Color',[170/255 32/255 55/255 .2],'LineWidth',1)
    end
    plot(X,nanmedian(TCmin.(char(dir))(plotpos,:)),'Color',[18/255 136/255 215/255],'LineWidth',2)
    plot(X,nanmedian(TCmax.(char(dir))(plotpos,:)),'Color',[170/255 32/255 55/255],'LineWidth',2)
    ylabel('spike rate (s^-^1)')
    xlabel('velocity (°/s)')
    title(dir)
    axis square
    set(gcf,'Position',[820 400 320 330])
    set(gca,'xscale','log','XTick',[10,100,1000],'XTickLabels',{'10','100','1000'})
    ylim([0 1])
%     xlim([10^1.1 10^3.6])
    xlim([X(1) X(end)])
    box on
    print([char(savename),'_TC_MedIndi_',dir],'-dsvg','-r300','-painters')
    savefig([char(savename),'_TC_MedIndi_',dir,'.fig'])
end

% Q10 bp
Boxplot_B(Q10.RS,2,15,[150/255 150/255 150/255; 150/255 150/255 150/255],{'ftb','btf'},[1 2])
writematrix(Q10.RS(:,1),'Q10_RespAmpl_ftb.txt') % for stats with R
writematrix(Q10.RS(:,2),'Q10_RespAmpl_btf.txt') % for stats with R
xlabel('movement direction')
ylabel('Q10 response strength')
ylim([0 8])
hold on
for j = 1 : 2
    x = ones(size(Q10.RS(:,j))).*(1+(rand(size(Q10.RS(:,j)))*2-1)/10)+j-1;
    scatter(x,Q10.RS(:,j),10,'filled','MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',1)
end
set(gcf,'Position',[820 450 320 330])
text(1,16,'ttest(1)')
[~,a] = ttest(Q10.RS(:,1),ones(length(Q10.RS(:,1)),1));
text(.5,14,num2str(a));
[~,a] = ttest(Q10.RS(:,2),ones(length(Q10.RS(:,2)),1));
text(1.5,14,num2str(a));
print([char(savename),'_bp_Q10_respStren_'],'-dsvg','-r300','-painters')
savefig([char(savename),'_bp_Q10_respStren_.fig'])

Boxplot_B(Q10.SM,2,15,[150/255 150/255 150/255; 150/255 150/255 150/255],{'ftb','btf'},[1 2])
writematrix(Q10.SM(:,1),'Q10_SensMax_ftb.txt') % for stats with R
writematrix(Q10.SM(:,2),'Q10_SensMax_btf.txt') % for stats with R
xlabel('movement direction')
ylabel('Q10 sensitivity maximum')
ylim([0 5])
hold on
for j = 1 : 2
    x = ones(size(Q10.SM(:,j))).*(1+(rand(size(Q10.SM(:,j)))*2-1)/10)+j-1;
    scatter(x,Q10.SM(:,j),10,'filled','MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',1)
end
set(gcf,'Position',[820 450 320 330])
text(1,5,'ttest(1)')
[~,a] = ttest(Q10.SM(:,1),ones(length(Q10.SM(:,1)),1));
text(.5,3,num2str(a));
[~,a] = ttest(Q10.SM(:,2),ones(length(Q10.SM(:,2)),1));
text(1.5,3,num2str(a));
print([char(savename),'_bp_Q10_sensMax_'],'-dsvg','-r300','-painters')
savefig([char(savename),'_bp_Q10_sensMax_.fig'])

Boxplot_B(Q10.B,2,15,[150/255 150/255 150/255; 150/255 150/255 150/255],{'ftb','btf'},[1 2])
% writematrix(Q10.B(:,1),'Q10_SensMax_B.txt') % for stats with R
xlabel('movement direction')
ylabel('Q10 baseline')
ylim([0 5])
hold on
for j = 1 : 2
    x = ones(size(Q10.B(:,j))).*(1+(rand(size(Q10.B(:,j)))*2-1)/10)+j-1;
    scatter(x,Q10.B(:,j),10,'filled','MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',1)
end
set(gcf,'Position',[820 450 320 330])
print([char(savename),'_bp_Q10_baseline_'],'-dsvg','-r300','-painters')
savefig([char(savename),'_bp_Q10_baseline_.fig'])

%% corr - peak - T
% Data = DataBU(plotpos);
for k = 1 : 2
    figure(50+k)
    axis square
    hold on
    if k == 1
        dir = 'ftb';
        dir2 = 'fw';
    else 
        dir = 'btf';
        dir2 = 'bw';
    end
    for  i = 1 : size(Data,2)
        for j = Data(i).first : Data(i).last
            if ~isnan(Data(i).peak.(char(dir)).pk(j))
                plot(Data(i).(['R0',num2str(j)]).Temp(2),Data(i).peak.(char(dir)).loc(j),'.','Color',[0.3 0.3 0.3])
            end
        end
    end
    box on
    axis square
    xlim([20 36])
    ylim([0 80])
    xlabel('Temperature (°C)')
    ylabel('Preferred temporal frequency (Hz)')
    set(gcf,'position',[700 500 400 350])
    title(dir)
    set(gcf,'Position',[420 250 320 330])
end
figure(51)
print([char(savename),'_scatter_sensMax_ftb'],'-dsvg','-r300','-painters')
savefig([char(savename),'_scatter_sensMax_ftb.fig'])
figure(52)
print([char(savename),'_scatter_sensMax_btf'],'-dsvg','-r300','-painters')
savefig([char(savename),'_scatter_sensMax_btf.fig'])

%% baseline
if type == 2
    fv = 12;
    a = [1 1 1 1 1 1 1 1 1 1 1 1];
    b = [2 5 5 6 3 3 2 2 6 6 2 2];
    b = [5 4 6 3 2 2 3 4 5 7 3 2];
    b = [3 2 2 6 3 2 2 3 5 4 3 2];
    b = b-1;
elseif type == 15
    fv = 2;
    a = [1 1];
    b = [1 2];
else
    if type == 8 % btf dir (8)
        fv = 3;
        a = [2 1];
        b = [2 3];
    elseif type == 9 % ftb dir (9)
        fv = 2;
        a = [1 1 3 3 1];
        b = [2 2 4 5 6];
    end
end
close(figure(80))
close(figure(85))
% close(figure(84))

for k = fv : size(Data,2)
    tempT = []; tempS = []; tempRSbtf = []; tempRSbtfT = []; tempRSftb = []; tempRSftbT = [];
    for i = a(k-fv+1) : b(k-fv+1) %Data(k).first : Data(k).last-1
        tempT = [tempT Data(k).(char(['C0',num2str(i)])).bg.T];
        tempS = [tempS Data(k).(char(['C0',num2str(i)])).bg.spikes];
        tempRSbtf = [tempRSbtf Data(k).(char(['C0',num2str(i)])).spikes.btf];
        tempRSbtfT = [tempRSbtfT Data(k).(char(['C0',num2str(i)])).T.btf];
        tempRSftb = [tempRSftb Data(k).(char(['C0',num2str(i)])).spikes.ftb];
        tempRSftbT = [tempRSftbT Data(k).(char(['C0',num2str(i)])).T.ftb];
    end

    table2(k,1) = mean(tempS(find(tempT == min(tempT)))); % min bg
    table2(k,2) = mean(tempS(find(tempT == max(tempT)))); % max bg

    figure(80)
    subplot(3,4,k-fv+1)
    hold on
    axis square
    box on
    plot(tempT,tempS,'.k')
    [P,S] = polyfit(tempT,tempS,1);
    rsq(k) = S.rsquared;
    % [y_fit,delta] = polyval(P,tempT,S);
    % figure
    % plot(tempT,tempS/max(tempS),'bo')
    % hold on
    % plot(tempT,y_fit,'r-')
    % plot(tempT,y_fit+2*delta,'m--',tempT,y_fit-2*delta,'m--')
    % title('Linear Fit of Data with 95% Prediction Interval')
    % legend('Data','Linear Fit','95% Prediction Interval')
    yfit = P(1)*tempT+P(2);
    slope(k) = P(1);
    ycut(k) = P(2);
    plot(tempT,yfit,'r-.','LineWidth',1,'Color','k')
    xlim([22 34])

    % Steigungsgerade für ftb und btf stimuli
    [Pftb,Sftb] = polyfit(tempRSftbT,tempRSftb,1);
    B26ftb(k) = Pftb(1)*26+Pftb(2);
    B28ftb(k) = Pftb(1)*28+Pftb(2);
    B30ftb(k) = Pftb(1)*30+Pftb(2);
    B32ftb(k) = Pftb(1)*32+Pftb(2);
    B34ftb(k) = Pftb(1)*34+Pftb(2);
    [Pbtf,Sbtf] = polyfit(tempRSbtfT,tempRSbtf,1);
    B26btf(k) = Pbtf(1)*26+Pbtf(2);
    B28btf(k) = Pbtf(1)*28+Pbtf(2);
    B30btf(k) = Pbtf(1)*30+Pbtf(2);
    B32btf(k) = Pbtf(1)*32+Pbtf(2);
    B34btf(k) = Pbtf(1)*34+Pbtf(2);
    slopeftb(k) = Pftb(1);
    slopebtf(k) = Pbtf(1);
    rsqftb(k) = Sftb.rsquared;
    rsqbtf(k) = Sbtf.rsquared;
    
    % plot die Steigungsgeraden für ftb und btf
    yfitftb = Pftb(1)*tempRSftbT+Pftb(2);
    yfitbtf = Pbtf(1)*tempRSbtfT+Pbtf(2);
    plot(tempRSftbT,tempRSftb,'.c')
    plot(tempRSftbT,yfitftb,'r-.','LineWidth',1,'Color','c')
    plot(tempRSbtfT,tempRSbtf,'.m')
    plot(tempRSbtfT,yfitbtf,'r-.','LineWidth',1,'Color','m')
    if k == 23
        legend('baseline','','ftb','','btf')
    end

    figure(85)
    hold on
    axis square
    box on
    % plot(tempT,tempS,'.')%,'Color',[.7 .7 .7 .5])
    plot(tempT,yfit,'-','LineWidth',1.5,'Color',[.3 .3 .3])
    tempc = corrcoef(tempT,tempS);
    corres(k) = tempc(1,2);
    corrp(k) = tempc(1,1);
    B26(k) = P(1)*26+P(2);
    B28(k) = P(1)*28+P(2);
    B30(k) = P(1)*30+P(2);
    B32(k) = P(1)*32+P(2);
    B34(k) = P(1)*34+P(2);


    % Q10.B(k) = (B34/B26)^(10/(34-26));
    % TBaselineftb(k,1) = tempRSftbT(find(tempRSftbT == min(tempRSftbT),1));
    % TBaselineftb(k,2) = tempRSftbT(find(tempRSftbT == max(tempRSftbT),1));
    % SBaselineftb(k,1) = tempRSftb(find(tempRSftbT == min(tempRSftbT),1));
    % SBaselineftb(k,2) = tempRSftb(find(tempRSftbT == max(tempRSftbT),1));
    ylim([0 300])
    xlim([22 34])
end
% plot(pX,pY,'.','Color',[.75 .45 .65 .5]) % die Daten sind komisch und
% falsch in die Datentabelle geschrieben, ich weiß nicht wieso und weil mir
% auch irgendwie niemand helfen will, werde ich die Daten wohl rausnehmen
% müssen. Ich liebs hier.
% plot(pX,pF,'-','LineWidth',1.5,'Color',[.75 .45 .65])

% repetitive
% close(figure(81))
for k = fv : size(Data,2)-1
    tempTf = []; tempSf = []; tempTb = []; tempSb = [];
    for i = a(k-fv+1) : b(k-fv+1) % Data(k).first : Data(k).last-1
        tempTf = [tempTf Data(k).(char(['C0',num2str(Data(i).first)])).T.ftb];
        tempSf = [tempSf Data(k).(char(['C0',num2str(Data(i).first)])).spikes.ftb];
        tempTb = [tempTb Data(k).(char(['C0',num2str(Data(i).first)])).T.btf];
        tempSb = [tempSb Data(k).(char(['C0',num2str(Data(i).first)])).spikes.btf];
    end

    % figure(81)
    % subplot(5,5,k)
    % hold on
    % axis square
    % box on
    % plot(tempTf,tempSf/max(tempSf),'.m')
    % plot(tempTb,tempSb/max(tempSb),'.b')
    P = polyfit(tempTf,tempSf,1);
    yfitf = P(1)*tempTf+P(2);
    % plot(tempTf,yfitf,'r-.','LineWidth',1)
    P = polyfit(tempTb,tempSb,1);
    yfitb = P(1)*tempTb+P(2);
    % plot(tempTb,yfitb,'r-.','LineWidth',1)

    table1(k,1) = mean(tempSf(find(tempTf == min(tempTf)))); % minftb
    table1(k,2) = mean(tempSf(find(tempTf == max(tempTf)))); % maxftb
    table1(k,3) = mean(tempSb(find(tempTb == min(tempTb)))); % minbtf
    table1(k,4) = mean(tempSb(find(tempTb == max(tempTb)))); % maxbtf

    % figure(84)
    % hold on
    % axis square
    % box on
    % plot(tempTf,tempSf/max(tempSf),'.','Color',[.3 .3 .7 .3])
    % plot(tempTb,tempSb/max(tempSb),'.','Color',[.3 .7 .7 .3])
    % plot(tempTf,yfitf,'-','LineWidth',1.5,'Color',[.3 .3 .7])
    % plot(tempTb,yfitb,'-','LineWidth',1.5,'Color',[.3 .7 .7])
end
table1(find(table1 == 0)) = NaN;
% Boxplot_B(table1,4,15,[.75 .75 .75; .75 .75 .75; .75 .75 .75; .75 .75 .75],{'min f','max f','min b','max b'},[1,2,3,4])
% for j = 1 : 4
%     x = ones(size(table1(:,j),1),1).*(1+(rand(size(table1(:,j),1),1)*2-1)/10)+j-1;
%     scatter(x,table1(:,j),10,'filled','MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',1)
% end
% title('repetitive')
% set(gcf,'Position',[720 250 320 330])
figure(80)
set(gcf,'Position',[10 70 900 900])
print([char(savename),'_scatter_slope_baseline_individual'],'-dsvg','-r300','-painters')
savefig([char(savename),'_scatter_slope_baseline_individual.fig'])
% figure(81)
% set(gcf,'Position',[420 250 320 330])
figure(85)
set(gcf,'Position',[1320 250 320 330])
print([char(savename),'_scatter_slope_baseline'],'-dsvg','-r300','-painters')
savefig([char(savename),'_scatter_slope_baseline.fig'])
% figure(84)
% set(gcf,'Position',[1620 250 320 330])

% table2(find(table2 == 0)) = NaN;
% Boxplot_B(table2,2,15,[.75 .75 .75; .75 .75 .75],{'min','max'},[1,2])
% for j = 1 : 2
%     x = ones(size(table2(:,j),1),1).*(1+(rand(size(table2(:,j),1),1)*2-1)/10)+j-1;
%     scatter(x,table2(:,j),10,'filled','MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',1)
% end
% title('baseline')
% set(gcf,'Position',[1020 250 320 330])
% signrank(table2(:,1),table2(:,2))
% signrank(table1(:,1),table1(:,2))
% signrank(table1(:,3),table1(:,4))

% für Steigung m ist 1 der höchstmögliche Wert (unter der Annahme, dass wir
% von 24°C auf 34°C eine maximale Steigung von 0 auf 1 haben). In dem Fall
% ist 0.1 der höchstmögliche Wert, um das zu vereinfachen rechne ich alle
% Werte mal 10, sodass der mögliche Wertebereich von 0 bis 1 ist
slope(find(slope == 0)) = NaN;
slopeftb(find(slopeftb == 0)) = NaN;
slopebtf(find(slopebtf == 0)) = NaN;
corres(find(corres == 0)) = NaN;
rsq(find(rsq == 0)) = NaN;
rsqftb(find(rsqftb == 0)) = NaN;
rsqbtf(find(rsqbtf == 0)) = NaN;
% slope = slope*10;
Boxplot_B([slope' slopeftb' slopebtf'],3,15,[150/255 150/255 150/255; 150/255 150/255 150/255; 150/255 150/255 150/255],{'m','corr','R^2'},[1 2 3])
x = ones(size(slope)).*(1+(rand(size(slope))*2-1)/25);
scatter(x,slope,10,'filled','MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',1)
x = ones(size(slopeftb)).*(1+(rand(size(slopeftb))*2-1)/25)+1;
scatter(x,slopeftb,10,'filled','MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',1)
x = ones(size(slopebtf)).*(1+(rand(size(slopebtf))*2-1)/25)+2;
scatter(x,slopebtf,10,'filled','MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',1)
xlim([.5 3.5])
ylim([0 20])
ylabel('slope')
set(gca,'XTick',[1 2 3],'XTickLabel',{'baseline','ftb','btf'})
set(gcf,'Position',[1020 250 320 330])
print([char(savename),'_bp_slope_baseline'],'-dsvg','-r300','-painters')
savefig([char(savename),'_bp_slope_baseline.fig'])

Boxplot_B([rsq' rsqftb' rsqbtf'],3,15,[150/255 150/255 150/255; 150/255 150/255 150/255; 150/255 150/255 150/255],{'m','corr','R^2'},[1 2 3])
x = ones(size(rsq)).*(1+(rand(size(rsq))*2-1)/25)+0;
scatter(x,rsq,10,'filled','MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',1)
x = ones(size(rsqftb)).*(1+(rand(size(rsqftb))*2-1)/25)+1;
scatter(x,rsqftb,10,'filled','MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',1)
x = ones(size(rsqbtf)).*(1+(rand(size(rsqbtf))*2-1)/25)+2;
scatter(x,rsqbtf,10,'filled','MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',1)
xlim([.5 3.5])
ylim([0 1])
ylabel('R^2')
set(gca,'XTick',[1 2 3],'XTickLabel',{'baseline','ftb','btf'})
set(gcf,'Position',[1020 250 320 330])
print([char(savename),'_bp_rsq_baseline'],'-dsvg','-r300','-painters')
savefig([char(savename),'_bp_rsq_baseline.fig'])

% LP für slope und r^2
figure
hold on
plot([1 2 3],[rsq; rsqftb; rsqbtf],'LineWidth',1.2,'Color',[.7 .7 .7])
plot(1:3,nanmedian([rsq' rsqftb' rsqbtf']),'LineWidth',1.7,'Color',[.3 .3 .3])
set(gcf,'Position',[820 450 320 330])
ylabel('R^2')
set(gca,'XTick',[1 2 3],'XTickLabel',{'baseline','ftb','btf'})
xlim([.5 3.5])
ylim([0 1])
box on
print([char(savename),'_lp_rsq_baseline'],'-dsvg','-r300','-painters')
savefig([char(savename),'_lp_rsq_baseline.fig'])
figure
hold on
plot([1 2 3],[slope; slopeftb; slopebtf],'LineWidth',1.2,'Color',[.7 .7 .7])
plot(1:3,nanmedian([slope' slopeftb' slopebtf']),'LineWidth',1.7,'Color',[.3 .3 .3])
set(gcf,'Position',[820 450 320 330])
ylabel('m')
set(gca,'XTick',[1 2 3],'XTickLabel',{'baseline','ftb','btf'})
xlim([.5 3.5])
ylim([0 20])
box on
print([char(savename),'_lp_slope_baseline'],'-dsvg','-r300','-painters')
savefig([char(savename),'_lp_slope_baseline.fig'])


%% Q10 für baseline
% TBaseline und SBaseline
% Q10 val

if type == 2 
    writematrix(Q10.B(12:end)',['Q10_B.txt']) % for stats with R
    writematrix(Q10.Bftb(12:end)',['Q10_Bftb.txt']) % for stats with R
    writematrix(Q10.Bbtf(12:end)',['Q10_Bbtf.txt']) % for stats with R
    Q10.B = (B34./B26).^(10/(34-26));
    Q10.Bftb = (B34ftb./B26ftb).^(10/(34-26));
    Q10.Bbtf = (B34btf./B26btf).^(10/(34-26));
elseif type == 15
    % Q10.B = ((B34+1)./(B26+1)).^(10/(34-26));
    % Q10.Bftb = ((B34ftb+1)./(B26ftb+1)).^(10/(34-26));
    % Q10.Bbtf = ((B34btf+1)./(B26btf+1)).^(10/(34-26));
    Q10.B = ((B32)./(B28)).^(10/(32-28));
    Q10.Bftb = ((B32ftb)./(B28ftb)).^(10/(32-28));
    Q10.Bbtf = ((B32btf)./(B28btf)).^(10/(32-28));
    Q10.B(1) = NaN;
    Q10.Bftb(1) = NaN;
    Q10.Bbtf(1) = NaN;    
else
    Q10.B = (B34./B26).^(10/(34-26));
    Q10.Bftb = (B34ftb./B26ftb).^(10/(34-26));
    Q10.Bbtf = (B34btf./B26btf).^(10/(34-26));
end

% SBaseline(find(SBaseline(:,1) == 0),:) = NaN;
% TBaseline(find(SBaseline(:,1) == 0),:) = NaN;
% SBaseline(11,:) = NaN;
% TBaseline(11,:) = NaN;
% for i = 1 : length(TBaseline)
%     Q10.B(i) = (SBaseline(i,2)/SBaseline(i,1))^(10/(TBaseline(i,2)-TBaseline(i,1)));
% end
% Q10.SM(i,k) = (a/b)^(10/(c-d));

Boxplot_B([Q10.B' Q10.Bftb' Q10.Bbtf'],3,15,[150/255 150/255 150/255; 150/255 150/255 150/255; 150/255 150/255 150/255],{'baseline','ftb','btf'},[1 2 3])
ylabel('Q10 repetitive stim')
ylim([0 10])
xlim([.5 3.5])
hold on
x = ones(size(Q10.B(1,:))).*(1+(rand(size(Q10.B(1,:)))*2-1)/10)+1-1;
scatter(x,Q10.B(1,:),10,'filled','MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',1)
x = ones(size(Q10.Bftb(1,:))).*(1+(rand(size(Q10.Bftb(1,:)))*2-1)/10)+2-1;
scatter(x,Q10.Bftb(1,:),10,'filled','MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',1)
x = ones(size(Q10.Bbtf(1,:))).*(1+(rand(size(Q10.Bbtf(1,:)))*2-1)/10)+3-1;
scatter(x,Q10.Bbtf(1,:),10,'filled','MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',1)
set(gcf,'Position',[820 450 320 330])
text(1,3.5,'ttest(1)')
[~,ab] = ttest(Q10.B',ones(length(Q10.B),1));
[~,aftb] = ttest(Q10.Bftb',ones(length(Q10.Bftb),1));
[~,abtf] = ttest(Q10.Bbtf',ones(length(Q10.Bbtf),1));
text(.5,3,[num2str(ab),', ',num2str(aftb),', ',num2str(abtf)]);
print([char(savename),'_bp_Q10_repetitiveStim_'],'-dsvg','-r300','-painters')
savefig([char(savename),'_bp_Q10_repetitiveStim_.fig'])

figure
hold on
plot(1:3,[Q10.B;Q10.Bftb; Q10.Bbtf],'LineWidth',1.2,'Color',[.7 .7 .7])
plot(1:3,nanmedian([Q10.B' Q10.Bftb' Q10.Bbtf']),'LineWidth',1.7,'Color',[.3 .3 .3])
set(gcf,'Position',[820 450 320 330])
ylabel('Q10 repetitive stim')
set(gca,'XTick',[1 2 3],'XTickLabel',{'baseline','ftb','btf'})
xlim([.5 3.5])
ylim([0 10])
box on
ft = friedman([Q10.B(12:end)' Q10.Bftb(12:end)' Q10.Bbtf(12:end)']);
close
abftb = signrank(Q10.B(12:end), Q10.Bftb(12:end));
abbtf = signrank(Q10.B(12:end), Q10.Bbtf(12:end));
aftbbtf = signrank(Q10.Bftb(12:end), Q10.Bbtf(12:end));
text(.5,4.5,['friedman(0.01 0.003 0.0003): ',num2str(ft)]);
text(.5,3.5,'wtest(a-ftb;a-btf;ftb-btf');
text(.5,3,[num2str(abftb),', ',num2str(abbtf),', ',num2str(aftbbtf)]);
print([char(savename),'_lp_Q10_repetitiveStim_'],'-dsvg','-r300','-painters')
savefig([char(savename),'_lp_Q10_repetitiveStim_.fig'])

%% hist sens max difference of T conditions
hdataftb = abs(Spikeval.ftb(plotpos,1)/0.0222-Spikeval.ftb(plotpos,2)/0.0222);
hdataftb = abs(Spikeval.ftb(plotpos,1)/0.0222).*Q10.SM(:,1);
hdatabtf = abs(Spikeval.btf(plotpos,1)/0.0222-Spikeval.btf(plotpos,2)/0.0222);
hdatabtf = abs(Spikeval.btf(plotpos,1)/0.0222).*Q10.SM(:,2);
hdataftb(find(isnan(hdataftb))) = 3000; 
hdatabtf(find(isnan(hdatabtf))) = 3000; 
figure; hold on;
histogram(hdataftb,'BinEdges',0:200:3000)
histogram(hdatabtf,'BinEdges',0:200:3000)
xlim([0 3500])
ylim([0 12])
set(gcf,'Position',[850 250 500 330])
box on
legend({'ftb','btf'},'Location','north')
xlabel('sensitivity maximum (°/s)')
ylabel('count')
print([char(savename),'_hist_sensMax'],'-dsvg','-r300','-painters')
savefig([char(savename),'_hist_sensMax.fig'])

figure
hold on
axis square
plot(abs(Ttemp(plotpos,1)-Ttemp(plotpos,2)),hdataftb,'.b')
plot(abs(Ttemp(plotpos,1)-Ttemp(plotpos,2)),hdatabtf,'.r')
set(gcf,'Position',[1020 250 320 330])
box on
xlim([0 12])
ylim([0 2250])
xlabel('diff(temperature (°C))')
ylabel('diff(sensitivity maximum (°/s))')
legend({'ftb','btf'},'Location','north')
print([char(savename),'_corr_sensMax_T_diff'],'-dsvg','-r300','-painters')
savefig([char(savename),'_corr_sensMax_T_diff.fig'])

%% hist sens max at T = 33 °C
hdataftb = (abs(Spikeval.ftb(plotpos,1)/0.0222).*Q10.SM(:,1)).*((Ttemp(plotpos,2)-Ttemp(plotpos,1))/10);
hdatabtf = (abs(Spikeval.btf(plotpos,1)/0.0222).*Q10.SM(:,2)).*((Ttemp(plotpos,2)-Ttemp(plotpos,1))/10);

Boxplot_B([hdataftb hdatabtf],2,15,[150/255 150/255 150/255; 150/255 150/255 150/255],{'ftb','btf'},[1 2])
ylabel('expected sensitivity maximum at T = 33 °C (° s^-^1)')
ylim([0 2500])
xlim([.5 2.5])
hold on
x = ones(size(hdataftb)).*(1+(rand(size(hdataftb))*2-1)/10)+1-1;
scatter(x,hdataftb,10,'filled','MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',1)
x = ones(size(hdataftb)).*(1+(rand(size(hdataftb))*2-1)/10)+2-1;
scatter(x,hdataftb,10,'filled','MarkerFaceColor',[.3 .3 .3],'MarkerFaceAlpha',1)
set(gcf,'Position',[820 450 320 330])
exppval = signrank(hdataftb, hdatabtf);
text(1,100,['wtest, p-val: ',num2str(exppval)]);
print([char(savename),'_bp_expected_SensMax_'],'-dsvg','-r300','-painters')
savefig([char(savename),'_bp_expected_SensMax_.fig'])

hdataftb(find(isnan(hdataftb))) = 3000; 
hdatabtf(find(isnan(hdatabtf))) = 3000; 
figure; hold on;
histogram(hdataftb,'BinEdges',0:200:3000)
histogram(hdatabtf,'BinEdges',0:200:3000)
xlim([0 3500])
ylim([0 12])
set(gcf,'Position',[850 250 500 330])
box on
legend({'ftb','btf'},'Location','north')
xlabel('sensitivity maximum (°/s)')
ylabel('count')
print([char(savename),'_hist_sensMax_T33'],'-dsvg','-r300','-painters')
savefig([char(savename),'_hist_sensMax_T33.fig'])

figure
hold on
axis square
plot(abs(Ttemp(plotpos,1)-Ttemp(plotpos,2)),hdataftb,'.b')
plot(abs(Ttemp(plotpos,1)-Ttemp(plotpos,2)),hdatabtf,'.r')
set(gcf,'Position',[1020 250 320 330])
box on
xlim([0 12])
ylim([0 2250])
xlabel('diff(temperature (°C))')
ylabel('diff(sensitivity maximum (°/s))')
legend({'ftb','btf'},'Location','north')
print([char(savename),'_corr_sensMax_T33'],'-dsvg','-r300','-painters')
savefig([char(savename),'_corr_sensMax_T33.fig'])
abftb = signrank(Q10.B(12:end), Q10.Bftb(12:end));
abbtf = signrank(Q10.B(12:end), Q10.Bbtf(12:end));
aftbbtf = signrank(Q10.Bftb(12:end), Q10.Bbtf(12:end));
text(.5,4.5,['friedman(0.01 0.003 0.0003): ',num2str(ft)]);
text(.5,3.5,'wtest(a-ftb;a-btf;ftb-btf');
text(.5,3,[num2str(abftb),', ',num2str(abbtf),', ',num2str(aftbbtf)]);

%% lineplot to compare sens max of ftb and btf
figure; subplot(1,2,1); hold on; axis square
temp = [Spikeval.ftb(plotpos,1) Spikeval.btf(plotpos,1)];
for i = 1 : length(plotpos)
    plot(temp(i,:),'b-','LineWidth',1.2)
end
xlim([.5 2.5])
ylim([0 70])
set(gca,'XTick',[1,2],'XTickLabels',{'ftb','btf'})
box on

signrank(temp(:,1),temp(:,2))


temp = [Spikeval.ftb(plotpos,2) Spikeval.btf(plotpos,2)];
temp(find(isnan(temp) == 1)) = 65;
subplot(1,2,2); hold on; axis square
for i = 1 : length(plotpos)
    plot(temp(i,:),'r-','LineWidth',1.2)
end
xlim([.5 2.5])
ylim([0 70])
set(gca,'XTick',[1,2],'XTickLabels',{'ftb','btf'})
box on
set(gcf,'Position',[680 250 700 330])

signrank(temp(:,1),temp(:,2))


%% plot individual units
close(figure(80))
for i = 1:23
    test = Data(i);
    T1 = test.(char(['R0',num2str(test.first)])).yfreq.mean.translation.fw.w8(~isnan(test.R01.yfreq.mean.translation.bw.w8));
    T2 = test.(char(['R0',num2str(test.first+1)])).yfreq.mean.translation.fw.w8(~isnan(test.R02.yfreq.mean.translation.bw.w8));
    % maxval = max([T1a,T2a]);
    % T1 = T1a/maxval;
    % T2 = T2a/maxval;
    xvar = test.xvelo.w8(~isnan(test.xvelo.w8));
    xvar = xvar(1:length(T1));

    %% plot unit

    figure(80)
    subplot(6,4,i)
    axis square
    hold on
    plot(xvar,T1,'b-')
    plot(xvar,T2,'m-')
    set(gcf,'position',[785 200 315 315])
    box on
    text(2250,0.9,[num2str(test.R01.Temp(1)),' °C'],'Color','b')
    text(2250,0.8,[num2str(test.R02.Temp(1)),' °C'],'Color','m')
end

%% close all
close all
    
figure
hold on
for i = 1 : 23 
    plot([Ttemp(plotpos(i),1) Ttemp(plotpos(i),2)],[Spikeval.ftb(plotpos(i),1) Spikeval.ftb(plotpos(i),2)],'b-')
    plot([Ttemp(plotpos(i),1) Ttemp(plotpos(i),2)],[Spikeval.btf(plotpos(i),1) Spikeval.btf(plotpos(i),2)],'r-')
end
% plot(Ttemp(plotpos,1),Spikeval.ftb(plotpos,1),'ob')
% plot(Ttemp(plotpos,2),Spikeval.ftb(plotpos,2),'ob')
% plot(Ttemp(plotpos,1),Spikeval.btf(plotpos,1),'or')
% plot(Ttemp(plotpos,2),Spikeval.btf(plotpos,2),'or')
xlim([20 40])
ylim([0 60])
set(gcf,'position',[785 200 315 315])
box on
title('ftb b - btf - r')
xlabel('Temperature')
ylabel('sensitivity maximum (° s^-^1)')

m = (Spikeval.ftb(plotpos,1).*Q10.SM(:,1)-Spikeval.ftb(plotpos,1))/10;
b = Spikeval.ftb(plotpos,1)-(m.*Ttemp(plotpos,1));
y = m.*x+b;
x = 20:40;
figure
hold on
for i = 1 :23
    % m(i) = (Spikeval.ftb(plotpos(i),1)*Q10.SM(i,1)-Spikeval.ftb(plotpos(i),1))/10;
    % b(i) = Spikeval.ftb(plotpos(i),1)-(m(i)*Ttemp(plotpos(i),1));
    % y(i,:) = m(i)*x+b(i);
    plot(x,y(i,:),'-b')
end

m = (Spikeval.btf(plotpos,1).*Q10.SM(:,2)-Spikeval.btf(plotpos,1))/10;
b = Spikeval.btf(plotpos,1)-(m.*Ttemp(plotpos,2));
y = m.*x+b;
x = 20:40;
for i = 1 : 23
    % m(i) = (Spikeval.ftb(plotpos(i),1)*Q10.SM(i,1)-Spikeval.ftb(plotpos(i),1))/10;
    % b(i) = Spikeval.ftb(plotpos(i),1)-(m(i)*Ttemp(plotpos(i),1));
    % y(i,:) = m(i)*x+b(i);
    plot(x,y(i,:),'-r')
end
xlim([15 45])
ylim([-50 120])
set(gcf,'position',[785 200 315 315])
box on
title('ftb b - btf - r')
xlabel('Temperature')
ylabel('expected sensitivity maximum (° s^-^1)')