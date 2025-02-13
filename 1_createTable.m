% to add individual units to large data-table

%%
% k = 1; % row in table
clear AllAni
% load AllResults_OF_T.mat
% save AllResults_OF_T_backup.mat AllAni
load Results_all_manuskript
save Results_all_manuskript_backup 'AllAni'

%% find correct line in table to save data 
Rec = 'R01';
unit = 2;
AniNum = str2double(file.name(strfind(file.name,'Ephys')+14:strfind(file.name,'Rec')-2)); %str2double(file.name(65:66)); %str2double(file.name(42:44));
switch Rec
    case 'R01'
        k = size(AllAni,2)+1;
        % save basic information
        AllAni(k).File = file.name;
        AllAni(k).Animal = AniNum;
        AllAni(k).UnitNbr = unit;
        AllAni(k).rep = rep;
        AllAni(k).xfreq.w2 = x_t_p2_plot;
        AllAni(k).xfreq.w4 = x_t_p4_plot; %/2;
        AllAni(k).xfreq.w8 = x_t_p8_plot; %/4;
        
        AllAni(k).xvelo.w2 = x2;
        AllAni(k).xvelo.w4 = x4;
        AllAni(k).xvelo.w8 = x8;
    otherwise
        AniF = find([AllAni.Animal] == AniNum);
        UnitF = find([AllAni.UnitNbr] == unit);
        k = AniF(find(AniF == UnitF));
end

%% save data to table
AllAni(k).(Rec).yfreq.mean.translation.bw.w2 = y_t_p2_plot.(['unit0',num2str(unit)]);
AllAni(k).(Rec).yfreq.mean.translation.bw.w4 = y_t_p4_plot.(['unit0',num2str(unit)]);
AllAni(k).(Rec).yfreq.mean.translation.bw.w8 = y_t_p8_plot.(['unit0',num2str(unit)]);
AllAni(k).(Rec).yfreq.mean.translation.fw.w2 = y_t_n2_plot.(['unit0',num2str(unit)]);
AllAni(k).(Rec).yfreq.mean.translation.fw.w4 = y_t_n4_plot.(['unit0',num2str(unit)]);
AllAni(k).(Rec).yfreq.mean.translation.fw.w8 = y_t_n8_plot.(['unit0',num2str(unit)]);
AllAni(k).(Rec).yfreq.mean.rotation.cw.w2 = y_r_p2_plot.(['unit0',num2str(unit)]);
AllAni(k).(Rec).yfreq.mean.rotation.cw.w4 = y_r_p4_plot.(['unit0',num2str(unit)]);
AllAni(k).(Rec).yfreq.mean.rotation.cw.w8 = y_r_p8_plot.(['unit0',num2str(unit)]);
AllAni(k).(Rec).yfreq.mean.rotation.ccw.w2 = y_r_n2_plot.(['unit0',num2str(unit)]);
AllAni(k).(Rec).yfreq.mean.rotation.ccw.w4 = y_r_n4_plot.(['unit0',num2str(unit)]);
AllAni(k).(Rec).yfreq.mean.rotation.ccw.w8 = y_r_n8_plot.(['unit0',num2str(unit)]);

AllAni(k).(Rec).yfreq.sd.translation.bw.w2 = y_t_p2_sd.(['unit0',num2str(unit)]);
AllAni(k).(Rec).yfreq.sd.translation.bw.w4 = y_t_p4_sd.(['unit0',num2str(unit)]);
AllAni(k).(Rec).yfreq.sd.translation.bw.w8 = y_t_p8_sd.(['unit0',num2str(unit)]);
AllAni(k).(Rec).yfreq.sd.translation.fw.w2 = y_t_n2_sd.(['unit0',num2str(unit)]);
AllAni(k).(Rec).yfreq.sd.translation.fw.w4 = y_t_n4_sd.(['unit0',num2str(unit)]);
AllAni(k).(Rec).yfreq.sd.translation.fw.w8 = y_t_n8_sd.(['unit0',num2str(unit)]);
AllAni(k).(Rec).yfreq.sd.rotation.cw.w2 = y_r_p2_sd.(['unit0',num2str(unit)]);
AllAni(k).(Rec).yfreq.sd.rotation.cw.w4 = y_r_p4_sd.(['unit0',num2str(unit)]);
AllAni(k).(Rec).yfreq.sd.rotation.cw.w8 = y_r_p8_sd.(['unit0',num2str(unit)]);
AllAni(k).(Rec).yfreq.sd.rotation.ccw.w2 = y_r_n2_sd.(['unit0',num2str(unit)]);
AllAni(k).(Rec).yfreq.sd.rotation.ccw.w4 = y_r_n4_sd.(['unit0',num2str(unit)]);
AllAni(k).(Rec).yfreq.sd.rotation.ccw.w8 = y_r_n8_sd.(['unit0',num2str(unit)]);

AllAni(k).(Rec).times.translation.bw.w2 = stimTimes.(['unit0',num2str(unit)]).translation.bw.w2;
AllAni(k).(Rec).times.translation.bw.w4 = stimTimes.(['unit0',num2str(unit)]).translation.bw.w4;
AllAni(k).(Rec).times.translation.bw.w8 = stimTimes.(['unit0',num2str(unit)]).translation.bw.w8;
AllAni(k).(Rec).times.translation.fw.w2 = stimTimes.(['unit0',num2str(unit)]).translation.fw.w2;
AllAni(k).(Rec).times.translation.fw.w4 = stimTimes.(['unit0',num2str(unit)]).translation.fw.w4;
AllAni(k).(Rec).times.translation.fw.w8 = stimTimes.(['unit0',num2str(unit)]).translation.fw.w8;
AllAni(k).(Rec).times.rotation.cw.w2 = stimTimes.(['unit0',num2str(unit)]).rotation.cw.w2;
AllAni(k).(Rec).times.rotation.cw.w4 = stimTimes.(['unit0',num2str(unit)]).rotation.cw.w4;
AllAni(k).(Rec).times.rotation.cw.w8 = stimTimes.(['unit0',num2str(unit)]).rotation.cw.w8;
AllAni(k).(Rec).times.rotation.ccw.w2 = stimTimes.(['unit0',num2str(unit)]).rotation.ccw.w2;
AllAni(k).(Rec).times.rotation.ccw.w4 = stimTimes.(['unit0',num2str(unit)]).rotation.ccw.w4;
AllAni(k).(Rec).times.rotation.ccw.w8 = stimTimes.(['unit0',num2str(unit)]).rotation.ccw.w8;

AllAni(k).(Rec).delay.translation.bw.w2 = delpos.y_t_p2_plot.(['unit0',num2str(unit)]);
AllAni(k).(Rec).delay.translation.bw.w4 = delpos.y_t_p4_plot.(['unit0',num2str(unit)]);
AllAni(k).(Rec).delay.translation.bw.w8 = delpos.y_t_p8_plot.(['unit0',num2str(unit)]);
AllAni(k).(Rec).delay.translation.fw.w2 = delpos.y_t_n2_plot.(['unit0',num2str(unit)]);
AllAni(k).(Rec).delay.translation.fw.w4 = delpos.y_t_n4_plot.(['unit0',num2str(unit)]);
AllAni(k).(Rec).delay.translation.fw.w8 = delpos.y_t_n8_plot.(['unit0',num2str(unit)]);
AllAni(k).(Rec).delay.rotation.cw.w2 = delpos.y_r_p2_plot.(['unit0',num2str(unit)]);
AllAni(k).(Rec).delay.rotation.cw.w4 = delpos.y_r_p4_plot.(['unit0',num2str(unit)]);
AllAni(k).(Rec).delay.rotation.cw.w8 = delpos.y_r_p8_plot.(['unit0',num2str(unit)]);
AllAni(k).(Rec).delay.rotation.ccw.w2 = delpos.y_r_n2_plot.(['unit0',num2str(unit)]);
AllAni(k).(Rec).delay.rotation.ccw.w4 = delpos.y_r_n4_plot.(['unit0',num2str(unit)]);
AllAni(k).(Rec).delay.rotation.ccw.w8 = delpos.y_r_n8_plot.(['unit0',num2str(unit)]);

AllAni(k).(Rec).background.sum = background.sum.(['unit0',num2str(unit)]);
AllAni(k).(Rec).background.rawmean = background.rawmean.(['unit0',num2str(unit)]);
AllAni(k).(Rec).background.rawsd = background.rawsd.(['unit0',num2str(unit)]);
AllAni(k).(Rec).peak.translation.bw = half.tp(:,unit);
AllAni(k).(Rec).peak.translation.fw = half.tn(:,unit);
AllAni(k).(Rec).peak.rotation.cw = half.rp(:,unit);
AllAni(k).(Rec).peak.rotation.ccw = half.rn(:,unit);
AllAni(k).(Rec).width.translation.bw = width.tp(:,unit);
AllAni(k).(Rec).width.translation.fw = width.tn(:,unit);
AllAni(k).(Rec).width.rotation.cw = width.rp(:,unit);
AllAni(k).(Rec).width.rotation.ccw = width.rn(:,unit);

AllAni(k).(Rec).RF = RF.(['unit0',num2str(unit)]);

if exist('Ch190') && floor(sum(Ch19.title == 'Temp')/length(Ch19.title))
    AllAni(k).(Rec).Temp.mean = mean(Ch19.values(find(Ch19.times > stim_times_start(1) & Ch19.times < stim_times_end(end))));
    AllAni(k).(Rec).Temp.sd = std(Ch19.values(find(Ch19.times > stim_times_start(1) & Ch19.times < stim_times_end(end))));
else
    AllAni(k).(Rec).Temp.mean = NaN;
    AllAni(k).(Rec).Temp.sd = NaN;
end

% save AllResultsRF.mat AllAni
% save 'AllResults_OF_T.mat' 'AllAni'
save Results_all_manuskript 'AllAni'



