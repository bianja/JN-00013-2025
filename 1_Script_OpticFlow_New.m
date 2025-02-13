% analyse Arena data
% script to analyse optic flow responses for a quick overview. Load files 
% exported from spike2 before running. Gives velocity and temporal 
% frequency plot for translation and rotation for three spatial 
% frequencies.
%
% 1. load mat file containing spike information (exported from Spike2)
% 2. load mat file containing stimulus information
% 3. run script 
%
% Version 03-03-2023
%% set filename to save figures
savename = 'Fig1_Example_TC_T32_2';

%% load data
% close all
clc
clearvars -except Ch1 Ch2 Ch8 Ch9 Ch10 Ch11 Ch12 Ch13 Ch14 Ch15 Ch16 Ch17 Ch18 Ch19 Ch20 Ch21 Ch22 Ch23 Ch24 Ch25 Ch26 Ch27 Ch28 Ch29 Ch30 Ch31 Ch32 Ch33 Ch34 Ch35 file rec savename
vars = inputdlg({'Number of units','Channel number of first unit','Number of stimuli repetitions'},...
              'Customer', [1 30; 1 30; 1 30]); 
unitnbr = str2double(vars{1}); % number of channels containing spike information (unit)
firstunit = str2double(vars{2})-1; % channel containing spike information of the first unit
rep = str2double(vars{3}); % number of repetitions of each stimuli
len = 108; % length of one repetition (results of number of temporal and spatial frequencies)
num = 10; % number of temporal frequencies used (for w8 individual - less number of stimuli)

rec_old = rec;
[stim_times_start, stim_times_end, rec, delete] = findTimesCorrected_new(Ch1.times, Ch1.values, -0.02, rec);

%% calculate spike count (per s), delay, and spike timings
for k = 1 : unitnbr
    count.(['unit0',num2str(k)]) = zeros(1,len*rep);
    f = 1;
        for i = 1 : len*rep
            temp = eval(['Ch',num2str(firstunit+k)]);
            count.(['unit0',num2str(k)])(i) = (sum((temp.times > stim_times_start(f)...
                & temp.times < stim_times_end(f)) == 1))/(stim_times_end(f) - stim_times_start(f));
            temp2 = find((temp.times > stim_times_start(f) & temp.times < stim_times_end(f)));
            temp3 = find((temp.times > stim_times_start(f)-5 & temp.times < stim_times_end(f)+5));
            if isempty(temp2) == 1
                temp2 = NaN;            
                delay.(['unit0',num2str(k)])(i) = NaN;
            else
                delay.(['unit0',num2str(k)])(i) = (temp.times(temp2(1)) - stim_times_start(f));
            end
            if isempty(temp3) == 1
                temp3 = NaN;            
                timings.(['unit0',num2str(k)]){i} = NaN;
            else
                timings.(['unit0',num2str(k)]){i} = (temp.times(temp3) - stim_times_start(f));
            end
            f = f + 1;
        end
end

if stim_times_end(1)-stim_times_start(1) < 6
else
    error('sth. wrong with start end times')
end

%% find positions and extract counts
pos.tp2 = find([rec(1:rep*len).vel] > 0 & [rec(1:rep*len).width] == 2 & contains({rec(1:rep*len).type},'translational'));
pos.tn2 = find([rec(1:rep*len).vel] < 0 & [rec(1:rep*len).width] == 2 & contains({rec(1:rep*len).type},'translational'));
pos.rp2 = find([rec(1:rep*len).vel] > 0 & [rec(1:rep*len).width] == 2 & contains({rec(1:rep*len).type},'rotational'));
pos.rn2 = find([rec(1:rep*len).vel] < 0 & [rec(1:rep*len).width] == 2 & contains({rec(1:rep*len).type},'rotational'));
pos.tp4 = find([rec(1:rep*len).vel] > 0 & [rec(1:rep*len).width] == 4 & contains({rec(1:rep*len).type},'translational'));
pos.tn4 = find([rec(1:rep*len).vel] < 0 & [rec(1:rep*len).width] == 4 & contains({rec(1:rep*len).type},'translational'));
pos.rp4 = find([rec(1:rep*len).vel] > 0 & [rec(1:rep*len).width] == 4 & contains({rec(1:rep*len).type},'rotational'));
pos.rn4 = find([rec(1:rep*len).vel] < 0 & [rec(1:rep*len).width] == 4 & contains({rec(1:rep*len).type},'rotational'));
pos.tp8 = find([rec(1:rep*len).vel] > 0 & [rec(1:rep*len).width] == 8 & contains({rec(1:rep*len).type},'translational'));
pos.tn8 = find([rec(1:rep*len).vel] < 0 & [rec(1:rep*len).width] == 8 & contains({rec(1:rep*len).type},'translational'));
pos.rp8 = find([rec(1:rep*len).vel] > 0 & [rec(1:rep*len).width] == 8 & contains({rec(1:rep*len).type},'rotational'));
pos.rn8 = find([rec(1:rep*len).vel] < 0 & [rec(1:rep*len).width] == 8 & contains({rec(1:rep*len).type},'rotational'));

% calculate background activity 
for i = 1 : unitnbr
    [background.sum.(['unit0',num2str(i)]), background.raw.(['unit0',num2str(i)]), spikes.(['unit0',num2str(i)]), exclude.(['unit0',num2str(i)])] = testData(stim_times_start, stim_times_end, eval(['Ch',num2str(i+firstunit)]), rep,len);
    [~, posexclude.tp2.(['unit0',num2str(i)])] = intersect(pos.tp2,exclude.(['unit0',num2str(i)]));  
    [~, posexclude.tp4.(['unit0',num2str(i)])] = intersect(pos.tp4,exclude.(['unit0',num2str(i)]));
    [~, posexclude.tp8.(['unit0',num2str(i)])] = intersect(pos.tp8,exclude.(['unit0',num2str(i)]));
end

for k = 1 : unitnbr % this loop is for the units
    tp2.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tp2);
    delpos.tp2.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.tp2);
    ttpb2.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tp2);
    
    tn2.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tn2);
    delpos.tn2.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.tn2);
    ttnb2.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tn2);
    
    rp2.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tp2);
    delpos.rp2.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.rp2);
    trpb2.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tp2);
    
    rn2.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tn2);
    delpos.rn2.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.rn2);
    trnb2.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tn2);
    
    tp4.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tp4);
    delpos.tp4.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.tp4);
    ttpb4.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tp4);
    
    tn4.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tn4);
    delpos.tn4.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.tn4);
    ttnb4.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tn4);
    
    rp4.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tp4);
    delpos.rp4.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.rp4);
    trpb4.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tp4);
    
    rn4.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tn4);
    delpos.rn4.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.rn4);
    trnb4.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tn4);
    
    tp8.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tp8);
    delpos.tp8.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.tp8);
    ttpb8.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tp8);
    
    tn8.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tn8);
    delpos.tn8.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.tn8);
    ttnb8.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tn8);
    
    rp8.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tp8);
    delpos.rp8.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.rp8);
    trpb8.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tp8);
    
    rn8.(['unit0',num2str(k)]) = count.(['unit0',num2str(k)])(pos.tn8);
    delpos.rn8.(['unit0',num2str(k)]) = delay.(['unit0',num2str(k)])(pos.rn8);
    trnb8.(['unit0',num2str(k)]) = background.raw.(['unit0',num2str(k)])(pos.tn8);
    
    j = 1;
    for i = 1 : rep*10
        tpST2.(['unit0',num2str(k)]){j} = timings.(['unit0',num2str(k)]){pos.tp2(j)};
        tnST2.(['unit0',num2str(k)]){j} = timings.(['unit0',num2str(k)]){pos.tn2(j)};
        rpST2.(['unit0',num2str(k)]){j} = timings.(['unit0',num2str(k)]){pos.tp2(j)};
        rnST2.(['unit0',num2str(k)]){j} = timings.(['unit0',num2str(k)]){pos.tn2(j)};
        tpST4.(['unit0',num2str(k)]){j} = timings.(['unit0',num2str(k)]){pos.tp4(j)};
        tnST4.(['unit0',num2str(k)]){j} = timings.(['unit0',num2str(k)]){pos.tn4(j)};
        rpST4.(['unit0',num2str(k)]){j} = timings.(['unit0',num2str(k)]){pos.tp4(j)};
        rnST4.(['unit0',num2str(k)]){j} = timings.(['unit0',num2str(k)]){pos.tn4(j)};
        if j < 8
            tpST8.(['unit0',num2str(k)]){j} = timings.(['unit0',num2str(k)]){pos.tp8(j)};
            tnST8.(['unit0',num2str(k)]){j} = timings.(['unit0',num2str(k)]){pos.tn8(j)};
            rpST8.(['unit0',num2str(k)]){j} = timings.(['unit0',num2str(k)]){pos.tp8(j)};
            rnST8.(['unit0',num2str(k)]){j} = timings.(['unit0',num2str(k)]){pos.tn8(j)};
        end
        j = j + 1;
    end
end

xtp2 = [rec(pos.tp2).vel];
xtn2 = [rec(pos.tn2).vel]*(-1);
xrp2 = [rec(pos.rp2).vel];
xrn2 = [rec(pos.rn2).vel]*(-1);
xtp4 = [rec(pos.tp4).vel];
xtn4 = [rec(pos.tn4).vel]*(-1);
xrp4 = [rec(pos.rp4).vel];
xrn4 = [rec(pos.rn4).vel]*(-1);
xtp8 = [rec(pos.tp8).vel];
xtn8 = [rec(pos.tn8).vel]*(-1);
xrp8 = [rec(pos.rp8).vel];
xrn8 = [rec(pos.rn8).vel]*(-1);
    
%% sort stimuli and plot
% sort stimuli
for k = 1 : unitnbr
    [x_t_p2,temp] = sort(xtp2);
    y_t_p2.(['unit0',num2str(k)]) = tp2.(['unit0',num2str(k)])(temp);
    y_t_p2_delay.(['unit0',num2str(k)]) = delpos.tp2.(['unit0',num2str(k)])(temp);
    tpb2.(['unit0',num2str(k)]) = ttpb2.(['unit0',num2str(k)])(temp);
    times.(['unit0',num2str(k)]).translation.bw.w2 = tpST2.(['unit0',num2str(k)])(temp);
    
    [x_t_n2,temp] = sort(xtn2);
    y_t_n2.(['unit0',num2str(k)]) = tn2.(['unit0',num2str(k)])(temp);
    y_t_n2_delay.(['unit0',num2str(k)]) = delpos.tn2.(['unit0',num2str(k)])(temp);
    tnb2.(['unit0',num2str(k)]) = ttnb2.(['unit0',num2str(k)])(temp);
    times.(['unit0',num2str(k)]).translation.fw.w2 = tnST2.(['unit0',num2str(k)])(temp);
    
    [x_r_p2,temp] = sort(xrp2);
    y_r_p2.(['unit0',num2str(k)]) = rp2.(['unit0',num2str(k)])(temp);
    y_r_p2_delay.(['unit0',num2str(k)]) = delpos.rp2.(['unit0',num2str(k)])(temp);
    rpb2.(['unit0',num2str(k)]) = trpb2.(['unit0',num2str(k)])(temp);
    times.(['unit0',num2str(k)]).rotation.cw.w2 = rpST2.(['unit0',num2str(k)])(temp);
    
    [x_r_n2,temp] = sort(xrn2);
    y_r_n2.(['unit0',num2str(k)]) = rn2.(['unit0',num2str(k)])(temp);
    y_r_n2_delay.(['unit0',num2str(k)]) = delpos.rn2.(['unit0',num2str(k)])(temp);
    rnb2.(['unit0',num2str(k)]) = trnb2.(['unit0',num2str(k)])(temp);
    times.(['unit0',num2str(k)]).rotation.ccw.w2 = rnST2.(['unit0',num2str(k)])(temp);
    
    [x_t_p4,temp] = sort(xtp4);
    y_t_p4.(['unit0',num2str(k)]) = tp4.(['unit0',num2str(k)])(temp);
    y_t_p4_delay.(['unit0',num2str(k)]) = delpos.tp4.(['unit0',num2str(k)])(temp);
    tpb4.(['unit0',num2str(k)]) = ttpb4.(['unit0',num2str(k)])(temp);
    times.(['unit0',num2str(k)]).translation.bw.w4 = tpST4.(['unit0',num2str(k)])(temp);
    
    [x_t_n4,temp] = sort(xtn4);
    y_t_n4.(['unit0',num2str(k)]) = tn4.(['unit0',num2str(k)])(temp);
    y_t_n4_delay.(['unit0',num2str(k)]) = delpos.tn4.(['unit0',num2str(k)])(temp);
    tnb4.(['unit0',num2str(k)]) = ttnb4.(['unit0',num2str(k)])(temp);
    times.(['unit0',num2str(k)]).translation.fw.w4 = tnST4.(['unit0',num2str(k)])(temp);
    
    [x_r_p4,temp] = sort(xrp4);
    y_r_p4.(['unit0',num2str(k)]) = rp4.(['unit0',num2str(k)])(temp);
    y_r_p4_delay.(['unit0',num2str(k)]) = delpos.rp4.(['unit0',num2str(k)])(temp);
    rpb4.(['unit0',num2str(k)]) = trpb4.(['unit0',num2str(k)])(temp);
    times.(['unit0',num2str(k)]).rotation.cw.w4 = rpST4.(['unit0',num2str(k)])(temp);
    
    [x_r_n4,temp] = sort(xrn4);
    y_r_n4.(['unit0',num2str(k)]) = rn4.(['unit0',num2str(k)])(temp);
    y_r_n4_delay.(['unit0',num2str(k)]) = delpos.rn4.(['unit0',num2str(k)])(temp);
    rnb4.(['unit0',num2str(k)]) = trnb4.(['unit0',num2str(k)])(temp);
    times.(['unit0',num2str(k)]).rotation.ccw.w4 = rnST4.(['unit0',num2str(k)])(temp);

    [x_t_p8,temp] = sort(xtp8);
    y_t_p8.(['unit0',num2str(k)]) = tp8.(['unit0',num2str(k)])(temp);
    y_t_p8_delay.(['unit0',num2str(k)]) = delpos.tp8.(['unit0',num2str(k)])(temp);
    tpb8.(['unit0',num2str(k)]) = ttpb8.(['unit0',num2str(k)])(temp);
    times.(['unit0',num2str(k)]).translation.bw.w8 = tpST8.(['unit0',num2str(k)])(temp);
    
    [x_t_n8,temp] = sort(xtn8);
    y_t_n8.(['unit0',num2str(k)]) = tn8.(['unit0',num2str(k)])(temp);
    y_t_n8_delay.(['unit0',num2str(k)]) = delpos.tn8.(['unit0',num2str(k)])(temp);
    tnb8.(['unit0',num2str(k)]) = ttnb8.(['unit0',num2str(k)])(temp);
    times.(['unit0',num2str(k)]).translation.fw.w8 = tnST8.(['unit0',num2str(k)])(temp);
    
    [x_r_p8,temp] = sort(xrp8);
    y_r_p8.(['unit0',num2str(k)]) = rp8.(['unit0',num2str(k)])(temp);
    y_r_p8_delay.(['unit0',num2str(k)]) = delpos.rp8.(['unit0',num2str(k)])(temp);
    rpb8.(['unit0',num2str(k)]) = trpb8.(['unit0',num2str(k)])(temp);
    times.(['unit0',num2str(k)]).rotation.cw.w8 = rpST8.(['unit0',num2str(k)])(temp);
    
    [x_r_n8,temp] = sort(xrn8);
    y_r_n8.(['unit0',num2str(k)]) = rn8.(['unit0',num2str(k)])(temp);
    y_r_n8_delay.(['unit0',num2str(k)]) = delpos.rn8.(['unit0',num2str(k)])(temp);
    rnb8.(['unit0',num2str(k)]) = trnb8.(['unit0',num2str(k)])(temp);
    times.(['unit0',num2str(k)]).rotation.ccw.w8 = rnST8.(['unit0',num2str(k)])(temp);
end

for k = 1 : unitnbr
    c = 1;
    for i = 1 : 10
        for j = 1 : rep
            stimTimes.(['unit0',num2str(k)]).translation.bw.w2(j,i) = times.(['unit0',num2str(k)]).translation.bw.w2(c);
            stimTimes.(['unit0',num2str(k)]).translation.fw.w2(j,i) = times.(['unit0',num2str(k)]).translation.fw.w2(c);
            stimTimes.(['unit0',num2str(k)]).rotation.cw.w2(j,i) = times.(['unit0',num2str(k)]).rotation.cw.w2(c);
            stimTimes.(['unit0',num2str(k)]).rotation.ccw.w2(j,i) = times.(['unit0',num2str(k)]).rotation.ccw.w2(c);
            stimTimes.(['unit0',num2str(k)]).translation.bw.w4(j,i) = times.(['unit0',num2str(k)]).translation.bw.w4(c);
            stimTimes.(['unit0',num2str(k)]).translation.fw.w4(j,i) = times.(['unit0',num2str(k)]).translation.fw.w4(c);
            stimTimes.(['unit0',num2str(k)]).rotation.cw.w4(j,i) = times.(['unit0',num2str(k)]).rotation.cw.w4(c);
            stimTimes.(['unit0',num2str(k)]).rotation.ccw.w4(j,i) = times.(['unit0',num2str(k)]).rotation.ccw.w4(c);
            if i < 8
                stimTimes.(['unit0',num2str(k)]).translation.bw.w8(j,i) = times.(['unit0',num2str(k)]).translation.bw.w8(c);
                stimTimes.(['unit0',num2str(k)]).translation.fw.w8(j,i) = times.(['unit0',num2str(k)]).translation.fw.w8(c);
                stimTimes.(['unit0',num2str(k)]).rotation.cw.w8(j,i) = times.(['unit0',num2str(k)]).rotation.cw.w8(c);
                stimTimes.(['unit0',num2str(k)]).rotation.ccw.w8(j,i) = times.(['unit0',num2str(k)]).rotation.ccw.w8(c);
            end
            c = c + 1;
        end
    end
end

x_t_p2_plot = zeros(length(1 : rep : num*rep)); x_t_n2_plot = zeros(length(1 : rep : num*rep)); x_r_p2_plot = zeros(length(1 : rep : num*rep)); x_r_n2_plot = zeros(length(1 : rep : num*rep));
x_t_p4_plot = zeros(length(1 : rep : num*rep)); x_t_n4_plot = zeros(length(1 : rep : num*rep)); x_r_p4_plot = zeros(length(1 : rep : num*rep)); x_r_n4_plot = zeros(length(1 : rep : num*rep));
x_t_p8_plot = zeros(length(1 : rep : num*rep)); x_t_n8_plot = zeros(length(1 : rep : num*rep)); x_r_p8_plot = zeros(length(1 : rep : num*rep)); x_r_n8_plot = zeros(length(1 : rep : num*rep));
j = 1;
for i = 1 : rep : num*rep % size(rec,2)/(3*rep)
    x_t_p2_plot(j) = x_t_p2(i);
    x_t_n2_plot(j) = x_t_n2(i);
    x_r_p2_plot(j) = x_r_p2(i);
    x_r_n2_plot(j) = x_r_n2(i);  
    x_t_p4_plot(j) = x_t_p4(i);
    x_t_n4_plot(j) = x_t_n4(i);
    x_r_p4_plot(j) = x_r_p4(i);
    x_r_n4_plot(j) = x_r_n4(i);
    if j < 8
        x_t_p8_plot(j) = x_t_p8(i);
        x_t_n8_plot(j) = x_t_n8(i);
        x_r_p8_plot(j) = x_r_p8(i);
        x_r_n8_plot(j) = x_r_n8(i);
    end
   
    for k = 1 : unitnbr
        for c = 1 : rep
            raw.y_t_p2.(['unit0',num2str(k)])(c,j) = y_t_p2.(['unit0',num2str(k)])(i+c-1);
            raw.y_t_n2.(['unit0',num2str(k)])(c,j) = y_t_n2.(['unit0',num2str(k)])(i+c-1);
            raw.y_r_p2.(['unit0',num2str(k)])(c,j) = y_r_p2.(['unit0',num2str(k)])(i+c-1);
            raw.y_r_n2.(['unit0',num2str(k)])(c,j) = y_r_n2.(['unit0',num2str(k)])(i+c-1);
            raw.y_t_p4.(['unit0',num2str(k)])(c,j) = y_t_p4.(['unit0',num2str(k)])(i+c-1);
            raw.y_t_n4.(['unit0',num2str(k)])(c,j) = y_t_n4.(['unit0',num2str(k)])(i+c-1);
            raw.y_r_p4.(['unit0',num2str(k)])(c,j) = y_r_p4.(['unit0',num2str(k)])(i+c-1);
            raw.y_r_n4.(['unit0',num2str(k)])(c,j) = y_r_n4.(['unit0',num2str(k)])(i+c-1);
            if j < 8
                raw.y_t_p8.(['unit0',num2str(k)])(c,j) = y_t_p8.(['unit0',num2str(k)])(i+c-1);
                raw.y_t_n8.(['unit0',num2str(k)])(c,j) = y_t_n8.(['unit0',num2str(k)])(i+c-1);
                raw.y_r_p8.(['unit0',num2str(k)])(c,j) = y_r_p8.(['unit0',num2str(k)])(i+c-1);
                raw.y_r_n8.(['unit0',num2str(k)])(c,j) = y_r_n8.(['unit0',num2str(k)])(i+c-1);
            end
        end
        y_t_p2_plot.(['unit0',num2str(k)])(j) = nanmedian(y_t_p2.(['unit0',num2str(k)])(i:i+rep-1)); %#ok<*NANMEDIAN>
        y_t_n2_plot.(['unit0',num2str(k)])(j) = nanmedian(y_t_n2.(['unit0',num2str(k)])(i:i+rep-1));
        y_r_p2_plot.(['unit0',num2str(k)])(j) = nanmedian(y_r_p2.(['unit0',num2str(k)])(i:i+rep-1));
        y_r_n2_plot.(['unit0',num2str(k)])(j) = nanmedian(y_r_n2.(['unit0',num2str(k)])(i:i+rep-1));
        y_t_p4_plot.(['unit0',num2str(k)])(j) = nanmedian(y_t_p4.(['unit0',num2str(k)])(i:i+rep-1)); %#ok<*NANMEDIAN>
        y_t_n4_plot.(['unit0',num2str(k)])(j) = nanmedian(y_t_n4.(['unit0',num2str(k)])(i:i+rep-1));
        y_r_p4_plot.(['unit0',num2str(k)])(j) = nanmedian(y_r_p4.(['unit0',num2str(k)])(i:i+rep-1));
        y_r_n4_plot.(['unit0',num2str(k)])(j) = nanmedian(y_r_n4.(['unit0',num2str(k)])(i:i+rep-1));        
        
        y_t_p2_sd.(['unit0',num2str(k)])(j) = nanstd(y_t_p2.(['unit0',num2str(k)])(i:i+rep-1)); %#ok<*NANSTD>
        y_t_n2_sd.(['unit0',num2str(k)])(j) = nanstd(y_t_n2.(['unit0',num2str(k)])(i:i+rep-1));
        y_r_p2_sd.(['unit0',num2str(k)])(j) = nanstd(y_r_p2.(['unit0',num2str(k)])(i:i+rep-1));
        y_r_n2_sd.(['unit0',num2str(k)])(j) = nanstd(y_r_n2.(['unit0',num2str(k)])(i:i+rep-1));
        y_t_p4_sd.(['unit0',num2str(k)])(j) = nanstd(y_t_p4.(['unit0',num2str(k)])(i:i+rep-1)); %#ok<*NANSTD>
        y_t_n4_sd.(['unit0',num2str(k)])(j) = nanstd(y_t_n4.(['unit0',num2str(k)])(i:i+rep-1));
        y_r_p4_sd.(['unit0',num2str(k)])(j) = nanstd(y_r_p4.(['unit0',num2str(k)])(i:i+rep-1));
        y_r_n4_sd.(['unit0',num2str(k)])(j) = nanstd(y_r_n4.(['unit0',num2str(k)])(i:i+rep-1));
        
        delpos.y_t_p2_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.tp2.(['unit0',num2str(k)])(i:i+rep-1));
        delpos.y_t_n2_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.tn2.(['unit0',num2str(k)])(i:i+rep-1));
        delpos.y_r_p2_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.rp2.(['unit0',num2str(k)])(i:i+rep-1));
        delpos.y_r_n2_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.rn2.(['unit0',num2str(k)])(i:i+rep-1));      
        delpos.y_t_p4_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.tp4.(['unit0',num2str(k)])(i:i+rep-1));
        delpos.y_t_n4_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.tn4.(['unit0',num2str(k)])(i:i+rep-1));
        delpos.y_r_p4_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.rp4.(['unit0',num2str(k)])(i:i+rep-1));
        delpos.y_r_n4_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.rn4.(['unit0',num2str(k)])(i:i+rep-1));           
        
        background.rawmean.(['unit0',num2str(k)]).translation.bw.w2(j) = nanmean(tpb2.(['unit0',num2str(k)])(i:i+rep-1));
        background.rawmean.(['unit0',num2str(k)]).translation.fw.w2(j) = nanmean(tnb2.(['unit0',num2str(k)])(i:i+rep-1));
        background.rawmean.(['unit0',num2str(k)]).rotation.cw.w2(j) = nanmean(rpb2.(['unit0',num2str(k)])(i:i+rep-1));
        background.rawmean.(['unit0',num2str(k)]).rotation.ccw.w2(j) = nanmean(rnb2.(['unit0',num2str(k)])(i:i+rep-1));
        background.rawmean.(['unit0',num2str(k)]).translation.bw.w4(j) = nanmean(tpb4.(['unit0',num2str(k)])(i:i+rep-1));
        background.rawmean.(['unit0',num2str(k)]).translation.fw.w4(j) = nanmean(tnb4.(['unit0',num2str(k)])(i:i+rep-1));
        background.rawmean.(['unit0',num2str(k)]).rotation.cw.w4(j) = nanmean(rpb4.(['unit0',num2str(k)])(i:i+rep-1));
        background.rawmean.(['unit0',num2str(k)]).rotation.ccw.w4(j) = nanmean(rnb4.(['unit0',num2str(k)])(i:i+rep-1));
        
        background.rawsd.(['unit0',num2str(k)]).translation.bw.w2(j) = nanstd(tpb2.(['unit0',num2str(k)])(i:i+rep-1));
        background.rawsd.(['unit0',num2str(k)]).translation.fw.w2(j) = nanstd(tnb2.(['unit0',num2str(k)])(i:i+rep-1));
        background.rawsd.(['unit0',num2str(k)]).rotation.cw.w2(j) = nanstd(rpb2.(['unit0',num2str(k)])(i:i+rep-1));
        background.rawsd.(['unit0',num2str(k)]).rotation.ccw.w2(j) = nanstd(rnb2.(['unit0',num2str(k)])(i:i+rep-1));
        background.rawsd.(['unit0',num2str(k)]).translation.bw.w4(j) = nanstd(tpb4.(['unit0',num2str(k)])(i:i+rep-1));
        background.rawsd.(['unit0',num2str(k)]).translation.fw.w4(j) = nanstd(tnb4.(['unit0',num2str(k)])(i:i+rep-1));
        background.rawsd.(['unit0',num2str(k)]).rotation.cw.w4(j) = nanstd(rpb4.(['unit0',num2str(k)])(i:i+rep-1));
        background.rawsd.(['unit0',num2str(k)]).rotation.ccw.w4(j) = nanstd(rnb4.(['unit0',num2str(k)])(i:i+rep-1));
        
        if j < 8
            y_t_p8_plot.(['unit0',num2str(k)])(j) = nanmedian(y_t_p8.(['unit0',num2str(k)])(i:i+rep-1)); %#ok<*NANMEDIAN>
            y_t_n8_plot.(['unit0',num2str(k)])(j) = nanmedian(y_t_n8.(['unit0',num2str(k)])(i:i+rep-1));
            y_r_p8_plot.(['unit0',num2str(k)])(j) = nanmedian(y_r_p8.(['unit0',num2str(k)])(i:i+rep-1));
            y_r_n8_plot.(['unit0',num2str(k)])(j) = nanmedian(y_r_n8.(['unit0',num2str(k)])(i:i+rep-1));
            y_t_p8_sd.(['unit0',num2str(k)])(j) = nanstd(y_t_p8.(['unit0',num2str(k)])(i:i+rep-1)); %#ok<*NANSTD>
            y_t_n8_sd.(['unit0',num2str(k)])(j) = nanstd(y_t_n8.(['unit0',num2str(k)])(i:i+rep-1));
            y_r_p8_sd.(['unit0',num2str(k)])(j) = nanstd(y_r_p8.(['unit0',num2str(k)])(i:i+rep-1));
            y_r_n8_sd.(['unit0',num2str(k)])(j) = nanstd(y_r_n8.(['unit0',num2str(k)])(i:i+rep-1));
            delpos.y_t_p8_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.tp8.(['unit0',num2str(k)])(i:i+rep-1));
            delpos.y_t_n8_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.tn8.(['unit0',num2str(k)])(i:i+rep-1));
            delpos.y_r_p8_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.rp8.(['unit0',num2str(k)])(i:i+rep-1));
            delpos.y_r_n8_plot.(['unit0',num2str(k)])(j) = nanmedian(delpos.rn8.(['unit0',num2str(k)])(i:i+rep-1));
            background.rawmean.(['unit0',num2str(k)]).translation.bw.w8(j) = nanmean(tpb8.(['unit0',num2str(k)])(i:i+rep-1));
            background.rawmean.(['unit0',num2str(k)]).translation.fw.w8(j) = nanmean(tnb8.(['unit0',num2str(k)])(i:i+rep-1));
            background.rawmean.(['unit0',num2str(k)]).rotation.cw.w8(j) = nanmean(rpb8.(['unit0',num2str(k)])(i:i+rep-1));
            background.rawmean.(['unit0',num2str(k)]).rotation.ccw.w8(j) = nanmean(rnb8.(['unit0',num2str(k)])(i:i+rep-1));
            background.rawsd.(['unit0',num2str(k)]).translation.bw.w8(j) = nanstd(tpb8.(['unit0',num2str(k)])(i:i+rep-1));
            background.rawsd.(['unit0',num2str(k)]).translation.fw.w8(j) = nanstd(tnb8.(['unit0',num2str(k)])(i:i+rep-1));
            background.rawsd.(['unit0',num2str(k)]).rotation.cw.w8(j) = nanstd(rpb8.(['unit0',num2str(k)])(i:i+rep-1));
            background.rawsd.(['unit0',num2str(k)]).rotation.ccw.w8(j) = nanstd(rnb8.(['unit0',num2str(k)])(i:i+rep-1));
        end
    end
    j = j + 1;
end

%% calculate receptive field data
for i = 1 : unitnbr
    RF.(['unit0',num2str(i)]) = receptiveFieldData(file.name(1:strfind(file.name, 'Recording01')+10),i,firstunit,unitnbr);
end

%% translation
% calculate values for x axis
x = [0.5, 1, 5, 10, 20, 40, 60, 80, 100, 120];
x2 = calcVelocity(x,2);
x4 = calcVelocity(x,4);
x8 = calcVelocity(x(1:7),8);

figure('Name','btf','NumberTitle','off'); hold on % btf
plotFreqTemp(x2, y_t_p2_plot,[1/255 102/255 94/255],'x:',background.sum);
plotFreqTemp(x4, y_t_p4_plot,[1/255 102/255 94/255],'x--',background.sum);
plotFreqTemp(x8, y_t_p8_plot,[1/255 102/255 94/255],'x-',background.sum);
% createLegend({'bw2','bw4','bw8'},[1/255 102/255 94/255; 1/255 102/255 94/255; 1/255 102/255 94/255],{':','--','-'},'best',1.5)
% ylim([0 80])
set(gcf,'position',[100 40 800 150])
axis square
% print([savename,'_btf'],'-depsc','-r300','-tiff','-painters')
% savefig([savename,'_btf.fig'])
    
figure('Name','ftb','NumberTitle','off'); hold on % ftb
plotFreqTemp(x2, y_t_n2_plot,[140/255 81/255 10/255],'x:',background.sum);
plotFreqTemp(x4, y_t_n4_plot,[140/255 81/255 10/255],'x--',background.sum);
plotFreqTemp(x8, y_t_n8_plot,[140/255 81/255 10/255],'x-',background.sum);
% createLegend({'fw2','fw4','fw8'},[140/255 81/255 10/255; 140/255 81/255 10/255; 140/255 81/255 10/255],{':','--','-'},'best',1.5)
% ylim([0 80])
set(gcf,'position',[100 270 800 150])
axis square
% print([savename,'_ftb'],'-depsc','-r300','-tiff','-painters')
% savefig([savename,'_ftb.fig'])
    
%% rotation
figure('Name','cw','NumberTitle','off'); hold on % cw
plotFreqTemp(x2, y_r_p2_plot,[1/255 102/255 94/255],'x:',background.sum);
plotFreqTemp(x4, y_r_p4_plot,[1/255 102/255 94/255],'x--',background.sum);
plotFreqTemp(x8, y_r_p8_plot,[1/255 102/255 94/255],'x-',background.sum);
% createLegend({'cw2','cw4','cw8'},[1/255 102/255 94/255; 1/255 102/255 94/255; 1/255 102/255 94/255],{':','--','-'},'best',1.5)
% ylim([0 450])
set(gcf,'position',[100 520 800 100])
axis square
% print([savename,'_cw'],'-depsc','-r300','-tiff','-painters')
% savefig([savename,'_cw.fig'])

figure('Name','ccw','NumberTitle','off'); hold on % ccw
plotFreqTemp(x2, y_r_n2_plot,[140/255 81/255 10/255],'x:',background.sum);
plotFreqTemp(x4, y_r_n4_plot,[140/255 81/255 10/255],'x--',background.sum);
plotFreqTemp(x8, y_r_n8_plot,[140/255 81/255 10/255],'x-',background.sum);
% createLegend({'ccw2','ccw4','ccw8'},[140/255 81/255 10/255; 140/255 81/255 10/255; 140/255 81/255 10/255],{':','--','-'},'best',1.5)
% ylim([0 450])
set(gcf,'position',[100 650 800 100])
axis square
% print([savename,'_ccw'],'-depsc','-r300','-tiff','-painters')
% savefig([savename,'_ccw.fig'])
