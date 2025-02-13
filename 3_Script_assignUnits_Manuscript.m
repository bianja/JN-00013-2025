cd \\132.187.28.171\home\rest\data\Analysis\matlab\scripts\temperature
clearvars
close all
% load('Temperature_Data_Temperature_Zwischenspeicher_backup3.mat')
load('Temperature_Data_Temperature_Zwischenspeicher_backup4n.mat')
% load('\\132.187.28.171\home\rest\data\Analysis\matlab\scripts\temperature\Data_Files\Temperature_Data_Temperature_Zwischenspeicher_backup8.mat')
width = [8 4 2];
pksDiffLim = .75; % diff between peaks needs to be at least 25 % diff of maximal acitivty
pksValDiffLim = .5; % if not the case: valley between peaks (defined as minimum value) needs to be smaller 50 % of maximal value to the smaller peak
highPeakLim = 0.75;
cur = 1;

%% baseline correction
dir2 = 'bw';
for i = cur : size(AllAni,2)
    for j = 1 : 6
        if ~isempty(AllAni(i).(['R0',num2str(j)]))
            for k = width
                if AllAni(i).Animal == 80
                    AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)])(2:find(isnan(AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)])) == 1,1)-1) = ...
                        AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)])(2:find(isnan(AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)])) == 1,1)-1) ...
                        - AllAni(i).(['R0',num2str(j)]).background.rawmean.translation.(char(dir2)).(['w',num2str(k)])(2:find(isnan(AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)])) == 1,1));
                else
                    if ~isempty(find(isnan(AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)])) == 1))
                        AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)])(1:find(isnan(AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)])) == 1,1)-1) = ...
                            AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)])(1:find(isnan(AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)])) == 1,1)-1) ...
                            - AllAni(i).(['R0',num2str(j)]).background.rawmean.translation.(char(dir2)).(['w',num2str(k)])(1:find(isnan(AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)])) == 1,1)-1);
                    else
                        AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)]) = ...
                            AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)]) ...
                            - AllAni(i).(['R0',num2str(j)]).background.rawmean.translation.(char(dir2)).(['w',num2str(k)]);
                    end
                end
            end
        end
    end
end

dir2 = 'fw';
for i = cur : size(AllAni,2)
    for j = 1 : 6
        if ~isempty(AllAni(i).(['R0',num2str(j)]))
            for k = width
                if AllAni(i).Animal == 80
                    AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)])(2:find(isnan(AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)])) == 1,1)-1) = ...
                        AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)])(2:find(isnan(AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)])) == 1,1)-1) ...
                        - AllAni(i).(['R0',num2str(j)]).background.rawmean.translation.(char(dir2)).(['w',num2str(k)])(2:find(isnan(AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)])) == 1,1));
                else
                    if ~isempty(find(isnan(AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)])) == 1))
                        AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)])(1:find(isnan(AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)])) == 1,1)-1) = ...
                            AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)])(1:find(isnan(AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)])) == 1,1)-1) ...
                            - AllAni(i).(['R0',num2str(j)]).background.rawmean.translation.(char(dir2)).(['w',num2str(k)])(1:find(isnan(AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)])) == 1,1)-1);
                    else
                        AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)]) = ...
                            AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.(char(dir2)).(['w',num2str(k)]) ...
                            - AllAni(i).(['R0',num2str(j)]).background.rawmean.translation.(char(dir2)).(['w',num2str(k)]);
                    end
                end
            end
        end
    end
end

%% assign response types
outfw = zeros(size(AllAni,2),6);
outbw = zeros(size(AllAni,2),6);
for i = cur : size(AllAni,2) % unit
    clear type vallim resplim bgmedian bgsd resposfw resposbw consecfw consecbw
    for j = 1 : 6 % temperature conditions/ repetitions
        if ~isempty(AllAni(i).(['R0',num2str(j)]))
            %             bgmedian = median([AllAni(i).(['R0',num2str(j)]).background.rawmean.translation.bw.w2 AllAni(i).(['R0',num2str(j)]).background.rawmean.translation.bw.w4 ...
            %                 AllAni(i).(['R0',num2str(j)]).background.rawmean.translation.bw.w8 AllAni(i).(['R0',num2str(j)]).background.rawmean.translation.fw.w2 ...
            %                 AllAni(i).(['R0',num2str(j)]).background.rawmean.translation.fw.w4 AllAni(i).(['R0',num2str(j)]).background.rawmean.translation.fw.w8]);
            %             bgsd = std([AllAni(i).(['R0',num2str(j)]).background.rawmean.translation.bw.w2 AllAni(i).(['R0',num2str(j)]).background.rawmean.translation.bw.w4 ...
            %                 AllAni(i).(['R0',num2str(j)]).background.rawmean.translation.bw.w8 AllAni(i).(['R0',num2str(j)]).background.rawmean.translation.fw.w2 ...
            %                 AllAni(i).(['R0',num2str(j)]).background.rawmean.translation.fw.w4 AllAni(i).(['R0',num2str(j)]).background.rawmean.translation.fw.w8]);
            
            % in some cases stimuli were not presented, i.e. there is no
            % response. These cases are included as NaN values for the x
            % vector (stimuli) but as real spike rate values in the y
            % vector (because we don't want to rewrite everything).
            % Therefore, use the x vector to get all stimuli that were
            % presented to exclude the not presented ones afterwards
            vallim = [1 find(isnan(AllAni(i).xfreq.w8),1,'first')-1];
            AniY = [AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.w8(vallim(1):vallim(2)-0) AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.w8(vallim(1):vallim(2)-1)];
            diffBW = AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.w8(vallim(1):vallim(2)-0);
            diffFW = AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.w8(vallim(1):vallim(2)-0);
            BWEXC = find(diff(find([true, diff(find(diffBW > max([diffBW diffFW])*0.1)) ~= 1, true])) > 3);  % if ~isempty -> 4 successive TF are greater median+10%of max [btf ftb]
            BWINH = find(diff(find([true, diff(find(diffBW < min([diffBW diffFW])*0.1)) ~= 1, true])) > 3);
            FWEXC = find(diff(find([true, diff(find(diffFW > max([diffBW diffFW])*0.1)) ~= 1, true])) > 3);
            FWINH = find(diff(find([true, diff(find(diffFW < min([diffBW diffFW])*0.1)) ~= 1, true])) > 3);
            
            
            
            %             BWEXC = find(diff(find([true, diff(find(AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.w8(vallim(1):vallim(2)) > max(AniY)*0.1)) ~= 1, true])) >= 2);  % if ~isempty -> 4 successive TF are greater median+10%of max [btf ftb]
            %             BWINH = find(diff(find([true, diff(find(AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.w8(vallim(1):vallim(2)) < min(AniY)*0.1)) ~= 1, true])) >= 2);
            %             FWEXC = find(diff(find([true, diff(find(AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.w8(vallim(1):vallim(2)) > max(AniY)*0.1)) ~= 1, true])) >= 2);
            %             FWINH = find(diff(find([true, diff(find(AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.w8(vallim(1):vallim(2)) < min(AniY)*0.1)) ~= 1, true])) >= 2);
            if ~isempty(BWEXC) && ~isempty(FWEXC) % exc exc
                AllAni(i).type(j) = 3;
            elseif ~isempty(BWEXC) && ~isempty(FWINH) % exc inh - TN1
                AllAni(i).type(j) = 4;
            elseif ~isempty(FWEXC) && ~isempty(BWINH) % exc inh - TN2
                AllAni(i).type(j) = 5;
            elseif ~isempty(BWINH) && ~isempty(FWINH) % inh inh
                AllAni(i).type(j) = 6;
            elseif ~isempty(BWEXC) % exc uni btf
                AllAni(i).type(j) = 1;
            elseif ~isempty(FWEXC) % exc uni ftb
                AllAni(i).type(j) = 2;
            elseif ~isempty(BWINH) % inh uni
                AllAni(i).type(j) = 6;
            elseif ~isempty(FWINH) % inh uni
                AllAni(i).type(j) = 6;
            end
        end
    end
    
    % if for each repetition the same response type was assigned, response
    % type is true
    if length(find(histc(AllAni(i).type,0:6) > 0)) == 1
        AllAni(i).first = 1;
        AllAni(i).last = length(AllAni(i).type);
        AllAni(i).type = AllAni(i).type(1); % for all repetitions same response type
    end
    
    % is type 3 is true for consecutive recordings cut rec
    if histc(AllAni(i).type,3) > 1 % response type 3 is at least two times detected within the whole vector
        pos = find(diff(find(AllAni(i).type == 3)) == 1); % find the position of the consecutive type 3 repetitions
        if ~isempty(pos) % at least two consecutive type 3 exist
%             if length(pos) == 1 % two consecutive rep type 3
% %                 AllAni(i).first = 1;
% %                 AllAni(i).last = length(AllAni(i).type);
% %                 AllAni(i).type = AllAni(i).type(1); % for all repetitions same response type
%             elseif length(pos) > 1  % more than two consecutive type 3 
                seq = find(AllAni(i).type == 3,1,'first');
                vec = find(AllAni(i).type == 3);
                diffpos = diff(vec);
                seqA = {};
                for j = 1 : length(vec) - 1
                    if diffpos(j) == 1
                        seq = [seq, vec(j+1)];
                    else
                        if length(seq) > 1
                            seqA{end+1} = seq;
                        end
                        seq = vec(j + 1);
                    end
                end
                if length(seq) > 1
                    seqA{end+1} = seq;
                end
                % differentiate the ones that are equal length (take first
                % one) from the ones that have different length (take
                % longest one)
                for j = 1 : size(seqA,2)
                    len(j) = length(seqA{j});
                end
                [~,maxpos] = max(len);
                AllAni(i).first = seqA{maxpos}(1);
                AllAni(i).last = seqA{maxpos}(end);
                AllAni(i).type = 3; % for all repetitions same response type
%             end
        end
    end
    
    % if a repetition at the beginning or at the end of all repetitions was
    % assigned to a different response type than the rest of the
    % repetitions, exclude that repetition and assign response type of the
    % remaining repetitions as true (but ensure that at least two
    % repetitions are available, otherwise exclude unit)
    if length(find(histc(AllAni(i).type,0:3) > 0)) > 1
        if length(find(histc(AllAni(i).type,1:3) > 0)) == 1 % at some point response potentially bad, find cutoff repetition
            if length(find(AllAni(i).type == find(histc(AllAni(i).type,1:3) > 0))) > 1 % only units with at least two repetitions of one response times are included
                if find(AllAni(i).type ~= 0,1,'first') < 2 % case, the beginning of recording was okay
                    AllAni(i).first = find(AllAni(i).type ~= 0,1,'first');
                    AllAni(i).last = find(AllAni(i).type ~= 0,1,'last');
                    AllAni(i).type = AllAni(i).type(1);
                else %if find(AllAni(i).type == 0) == 1 && find(AllAni(i).type == 0) ~= 2 % case, end of recording was okay
                    AllAni(i).first = find(AllAni(i).type ~= 0,1,'first');
                    AllAni(i).last = find(AllAni(i).type ~= 0,1,'last');
                    AllAni(i).type = AllAni(i).type(2);
                end
            end
        end
    end
    
%     % if 3 AND 1 OR 2 -> response type to 
%     if find(histc(AllAni(i).type,3) > 0) % at least one rep is type 3
%         if find(histc(AllAni(i).type,2) > 0) && find(histc(AllAni(i).type,1) > 0) % exclude because no consistent response type
%         elseif find(histc(AllAni(i).type,2) > 0) % set to response type 2
%         elseif find(histc(AllAni(i).type,3) > 0) % set to response type 1
%         else % do nothing
%         end
%     end

    
    % if a unit is assigned to response type 1 OR 2 (unidir), but for some of the
    % repetitions to response type 3, assign 1/2 as response type
    if length(AllAni(i).type) > 1 && length(find(histc(AllAni(i).type,2:3) > 0)) == 2
        if length(find(AllAni(i).type == 2)) > 1
            if ~isempty(find(histc(AllAni(i).type,0:1) > 0))
                if ~isempty(find(diff(find(AllAni(i).type == 2)) == 1))
                    allPos = find(AllAni(i).type == 2);
                    AllAni(i).first = allPos(find(diff(find(AllAni(i).type == 2)) == 1,1,'first'));
                    AllAni(i).last = allPos(find(diff(find(AllAni(i).type == 2)) == 1,1,'last'))+1;
                    AllAni(i).type = 2;
                end
            else
                AllAni(i).first = 1;
                AllAni(i).last = length(AllAni(i).type);
                AllAni(i).type = 2;
            end
        elseif histc(AllAni(i).type,3) > 1
            if ~isempty(find(diff(find(AllAni(i).type == 3)) == 1))
                allPos = find(AllAni(i).type == 3);
                AllAni(i).first = allPos(find(diff(find(AllAni(i).type == 3)) == 1,1,'first'));
                AllAni(i).last = allPos(find(diff(find(AllAni(i).type == 3)) == 1,1,'last'))+1;
                AllAni(i).type = 3;
            end % no else because it can still be response type 2
        end
    end
    if length(AllAni(i).type) > 1 && length(find(histc(AllAni(i).type,[1 3]) > 0)) == 2
        if length(find(AllAni(i).type == 1)) >= 1
            if ~isempty(find(histc(AllAni(i).type,0) > 0)) && ~isempty(find(histc(AllAni(i).type,2) > 0))
                if ~isempty(find(diff(find(AllAni(i).type == 1)) == 1))
                    allPos = find(AllAni(i).type == 1);
                    AllAni(i).first = allPos(find(diff(find(AllAni(i).type == 1)) == 1,1,'first'));
                    AllAni(i).last = allPos(find(diff(find(AllAni(i).type == 1)) == 1,1,'last'))+1;
                    AllAni(i).type = 1;
                else
                    AllAni(i).first = NaN;
                    AllAni(i).last = NaN;
                    AllAni(i).type = 0;
                end
            else
                AllAni(i).first = 1;
                AllAni(i).last = length(AllAni(i).type);
                AllAni(i).type = 1;
            end
        elseif histc(AllAni(i).type,3) > 1
            if ~isempty(find(diff(find(AllAni(i).type == 3)) == 1))
                allPos = find(AllAni(i).type == 3);
                AllAni(i).first = allPos(find(diff(find(AllAni(i).type == 3)) == 1,1,'first'));
                AllAni(i).last = allPos(find(diff(find(AllAni(i).type == 3)) == 1,1,'last'))+1;
                AllAni(i).type = 3;
            else
                AllAni(i).first = NaN;
                AllAni(i).last = NaN;
                AllAni(i).type = 0;
            end
        end
    else
        %         AllAni(i).first = NaN;
        %         AllAni(i).last = NaN;
        %         AllAni(i).type = 0;
    end
    
        
    % exclude all units that have a different response type for every
    % repetition
    if length(AllAni(i).type) > 1
        if histc(AllAni(i).type,1:3) <= 1
            AllAni(i).type = 0;
        end
    end
    
    % if a unit is assigned to response type 1 AND 2, exclude
    if length(AllAni(i).type) > 1 && length(find(histc(AllAni(i).type,1:2) > 0)) == 2
        AllAni(i).type = 0;
    end
    
    % all 'not responsding units' to NaN
    if length(AllAni(i).type) == 1 && AllAni(i).type == 0
        AllAni(i).first = NaN;
        AllAni(i).last = NaN;
    end
    
    for j = 1 : 6 % temperature conditions/ repetitions (individual for each animal)
        if ~isempty(AllAni(i).(['R0',num2str(j)]))
            field = {'peak', 'width'}; % clear because this is old code that potentially is wrong and bad code
            if isfield(AllAni(i).(['R0',num2str(j)]),field)
                AllAni(i).(['R0',num2str(j)]) = rmfield(AllAni(i).(['R0',num2str(j)]),field);
            end
        end
    end
    
    %% update first and last and type based on rec quality
    AllAni(72).type = 7;
    AllAni(35).type = 0;
    AllAni(35).last = NaN;
    AllAni(35).first = NaN;
    AllAni(34).type = 0;
    AllAni(34).last = NaN;
    AllAni(34).first = NaN;
    AllAni(69).type = 7;
    AllAni(71).last = 2;
    AllAni(90).first = 2;
    AllAni(98).first = 4;
    AllAni(113).last = 3;
%     AllAni(11).type = 0;
%     AllAni(11).first = NaN;
%     AllAni(11).last = NaN;
if exist('but')
else
    but = load('\\132.187.28.171\home\rest\data\Analysis\matlab\scripts\temperature\Data_Files\Temperature_Data_correct_backup2.mat');
end
for q = 1 : size(AllAni,2)
    AllAni(q).first = but.AllAni(q).first;
    AllAni(q).last = but.AllAni(q).last;
end

    %% find peaks and fwhm for all stim types using the function findpeaks
    % in diesem und dem darauffolgenden part werden die Kurvenverläufe
    % geprüft und ggf aussortiert. Wenn eine Kurve unschön ist, wird das in
    % die Matrizen outfw oder outbw hinterlegt (Wert == 1).
    for j = AllAni(i).first : AllAni(i).last % temperature conditions/ repetitions (individual for each animal)
        if ~isnan(AllAni(i).first)
            if ~isempty(AllAni(i).(['R0',num2str(j)]))
                width = [8,4,2];
                for k = width % width
                    % find peaks - fw
                    clear x y pks locs w p
                    y = AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.(['w',num2str(k)]);
                    x = AllAni(i).xfreq.(['w',num2str(k)]);
                    if ~isempty(y(~(isnan(y))))
                        [pks,locs,w,p] = findpeaks(y(~isnan(y)),x(~isnan(y)),'WidthReference','halfheight');
                        [pksmin,locsmin] = findpeaks(-y(~isnan(-y)),x(~isnan(-y)));
                    end
                    % if one peak perfectly fine
                    if exist('pks','var')
                        if length(pks) == 1
                            % if last spike rate value higher than peak and
                            % peak smaller than 0.5*max spike rate of this 
                            % unit: define as high peak response 
                            if pks < y(vallim(2)) && pks < max(y)*highPeakLim
                                AllAni(i).peak.ftb.pk(j).(['w',num2str(k)]) = max(y);
                                AllAni(i).peak.ftb.loc(j).(['w',num2str(k)]) = x(find(y == (max(y))));
                                AllAni(i).peak.ftb.wid(j).(['w',num2str(k)]) = NaN;
                                AllAni(i).type = 7;
                            else 
                                AllAni(i).peak.ftb.pk(j).(['w',num2str(k)]) = pks;
                                AllAni(i).peak.ftb.loc(j).(['w',num2str(k)]) = locs;
                                AllAni(i).peak.ftb.wid(j).(['w',num2str(k)]) = w;
                            end
                        elseif length(locs) >= 3 || length(pksmin) >= 3 % if 4 or more peaks || if three or more valleys, but less peaks: exclude
                            if k == 8 % response type depends on width == 8
                                outfw(i,j) = 1;
                            end
                            AllAni(i).peak.ftb.pk(j).(['w',num2str(k)]) = NaN;
                            AllAni(i).peak.ftb.loc(j).(['w',num2str(k)]) = NaN;
                            AllAni(i).peak.ftb.wid(j).(['w',num2str(k)]) = NaN;
                        elseif isempty(pks) % no peak deetcted -> peak likely to be at highest presented temporal frequency
                            if x(find(y == (max(y)))) == max(x) % highest ft has highest spike rate
                                AllAni(i).peak.ftb.pk(j).(['w',num2str(k)]) = max(y);
                                AllAni(i).peak.ftb.loc(j).(['w',num2str(k)]) = x(find(y == (max(y))));
                                AllAni(i).peak.ftb.wid(j).(['w',num2str(k)]) = NaN;
                                AllAni(i).type = 7;
                            end
                        else % if not, check for some criteria
                            % diff between peaks needs to be at least 25 % of
                            % maximal acitivty
                            if(max(pks((pks ~= max(pks)))))/(max(pks)) < pksDiffLim
                                AllAni(i).peak.ftb.pk(j).(['w',num2str(k)]) = max(pks);
                                AllAni(i).peak.ftb.loc(j).(['w',num2str(k)]) = locs(pks==max(pks));
                                AllAni(i).peak.ftb.wid(j).(['w',num2str(k)]) = w(pks==max(pks));
                            else
                                % if not the case: valley between peaks
                                % (defined as minimum value) needs to be
                                % smaller 25 % of maximal value to the smaller
                                % peak
                                if (max(abs(pksmin)))/(min(pks)) > pksValDiffLim
                                    AllAni(i).peak.ftb.pk(j).(['w',num2str(k)]) = max(pks);
                                    AllAni(i).peak.ftb.loc(j).(['w',num2str(k)]) = locs(pks==max(pks));
                                    AllAni(i).peak.ftb.wid(j).(['w',num2str(k)]) = w(pks==max(pks));
                                else % exclude
                                    if k == 8 % response type depends on width == 8
                                        outfw(i,j) = 1;
                                    end
                                    AllAni(i).peak.ftb.pk(j).(['w',num2str(k)]) = NaN;
                                    AllAni(i).peak.ftb.loc(j).(['w',num2str(k)]) = NaN;
                                    AllAni(i).peak.ftb.wid(j).(['w',num2str(k)]) = NaN;
                                end
                            end
                            
                            % if peaks with exact same amplitude and valley
                            % between peaks smaller 25 % of maximum: real peak
                            % is the one that has the smaller amplitude
                            % difference to one neighbour point (nur einmal
                            % der Fall? Lohnt sich der Aufwand?) neeeee
                            
                        end
                    else % in these cases stimulus was wrong (update rate) - no data available
                        AllAni(i).peak.ftb.pk(j).(['w',num2str(k)]) = NaN;
                        AllAni(i).peak.ftb.loc(j).(['w',num2str(k)]) = NaN;
                        AllAni(i).peak.ftb.wid(j).(['w',num2str(k)]) = NaN;
                    end
                end
            end
        end
    end
    
    %% find peaks and fwhm for all stim types using the function findpeaks
    % I am sorry, das wäre bestimmt auch eine gute Funktion gewesen oder
    % weniger code durch einen for-loop. Dieser Block ist exakt der gleich
    % code wie der Block oben drüber, nur anstatt für ftb für btf.
    for j = AllAni(i).first : AllAni(i).last % temperature conditions/ repetitions (individual for each animal)
        if ~isnan(AllAni(i).first)
            if ~isempty(AllAni(i).(['R0',num2str(j)]))
                width = [8,4,2];
                for k = width % width
                    % find peaks - bw
                    clear x y pks locs w p
                    y = AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.(['w',num2str(k)]);
                    x = AllAni(i).xfreq.(['w',num2str(k)]);
                    if ~isempty(y(~(isnan(y))))
                        [pks,locs,w,p] = findpeaks(y(~isnan(y)),x(~isnan(y)),'WidthReference','halfheight');
                        [pksmin,locsmin] = findpeaks(-y(~isnan(-y)),x(~isnan(-y)));
                    end
                    % if one peak perfectly fine
                    if exist('pks','var')
                        if length(pks) == 1
                            % if last spike rate value higher than peak and
                            % peak smaller than 0.5*max spike rate of this 
                            % unit: define as high peak response 
                            if pks < y(vallim(2)) && pks < max(y)*highPeakLim
                                AllAni(i).peak.btf.pk(j).(['w',num2str(k)]) = max(y);
                                AllAni(i).peak.btf.loc(j).(['w',num2str(k)]) = x(find(y == (max(y))));
                                AllAni(i).peak.btf.wid(j).(['w',num2str(k)]) = NaN;
                                if AllAni(i).type == 7
                                    AllAni(i).type = 9;
                                else
                                    AllAni(i).type = 8;
                                end
                            else
                                AllAni(i).peak.btf.pk(j).(['w',num2str(k)]) = pks;
                                AllAni(i).peak.btf.loc(j).(['w',num2str(k)]) = locs;
                                AllAni(i).peak.btf.wid(j).(['w',num2str(k)]) = w;
                            end
                        elseif length(locs) >= 3 || length(pksmin) >= 3 % length(find(locs >= 10)) >= 3 % % if 3 or more peaks || if 3 or more valleys, but less peaks: exclude
                            if k == 8 % response type depends on width == 8
                                outbw(i,j) = 1;
                            end
                            AllAni(i).peak.btf.pk(j).(['w',num2str(k)]) = NaN;
                            AllAni(i).peak.btf.loc(j).(['w',num2str(k)]) = NaN;
                            AllAni(i).peak.btf.wid(j).(['w',num2str(k)]) = NaN;
                        elseif isempty(pks) % no peak deetcted -> peak likely to be at highest presented temporal frequency
                            if x(find(y == (max(y)))) == max(x) % highest ft has highest spike rate
                                AllAni(i).peak.btf.pk(j).(['w',num2str(k)]) = max(y);
                                AllAni(i).peak.btf.loc(j).(['w',num2str(k)]) = x(find(y == (max(y))));
                                AllAni(i).peak.btf.wid(j).(['w',num2str(k)]) = NaN;
                                if AllAni(i).type == 7
                                    AllAni(i).type = 9;
                                else
                                    AllAni(i).type = 8;
                                end
                            end
                        else % if not, check for some criteria
                            % diff between peaks needs to be at least 25 % of
                            % maximal acitivty
                            if(max(pks((pks ~= max(pks)))))/(max(pks)) < pksDiffLim
                                AllAni(i).peak.btf.pk(j).(['w',num2str(k)]) = max(pks);
                                AllAni(i).peak.btf.loc(j).(['w',num2str(k)]) = locs(pks==max(pks));
                                AllAni(i).peak.btf.wid(j).(['w',num2str(k)]) = w(pks==max(pks));
                            else
                                % if not the case: valley between peaks
                                % (defined as minimum value) needs to be
                                % smaller 25 % of maximal value to the smaller
                                % peak
                                if (max(abs(pksmin)))/(min(pks)) > pksValDiffLim
                                    AllAni(i).peak.btf.pk(j).(['w',num2str(k)]) = max(pks);
                                    AllAni(i).peak.btf.loc(j).(['w',num2str(k)]) = locs(pks==max(pks));
                                    AllAni(i).peak.btf.wid(j).(['w',num2str(k)]) = w(pks==max(pks));
                                else % exclude
                                    if k == 8 % response type depends on width == 8
                                        outbw(i,j) = 1;
                                    end
                                    AllAni(i).peak.btf.pk(j).(['w',num2str(k)]) = NaN;
                                    AllAni(i).peak.btf.loc(j).(['w',num2str(k)]) = NaN;
                                    AllAni(i).peak.btf.wid(j).(['w',num2str(k)]) = NaN;
                                end
                            end
                            
                            % if peaks with exact same amplitude and valley
                            % between peaks smaller 25 % of maximum: real peak
                            % is the one that has the smaller amplitude
                            % difference to one neighbour point (nur einmal
                            % der Fall? Lohnt sich der Aufwand?) neeeee
                            
                        end
                    else % in these cases stimulus was wrong (update rate) - no data available
                        AllAni(i).peak.btf.pk(j).(['w',num2str(k)]) = NaN;
                        AllAni(i).peak.btf.loc(j).(['w',num2str(k)]) = NaN;
                        AllAni(i).peak.btf.wid(j).(['w',num2str(k)]) = NaN;
                    end
                end
            end
        end
    end
    
    %%
    % Bisher basiert die Einteilung der units auf der Abweichung zur
    % Hintergrundaktivität. In den beiden Blöcken über diesen hier hab ich
    % für jede Kurve geschaut ob was komisches passiert, also diese
    % mehrfachen peaks oder Kurven, die eigentlich gar keine Kurven sind.
    % In diesem part schaue ich dann jetzt, wo diese 'komischen' Kurven
    % liegen und überarbeite ggf die Einteilung in die Gruppen bzw. ändere
    % die Temperaturkondition (also dass z.B. am Anfang oder Ende Kurven
    % rausgenommen werden, weil das recording nicht stabil war)
    clear len
    if isnan(AllAni(i).first)
        AllAni(i).type = 0;
    end
    if AllAni(i).type ~= 0 % for responding units check curves and correct repetitions, first, and last for each unit
        % define start and end repetition for each recording in outfw
        if AllAni(i).first ~= 1
            outfw(i,1:AllAni(i).first-1) = NaN;
        end
        if AllAni(i).last ~= 6
            outfw(i,AllAni(i).last+1:6) = NaN;
        end
        
        % in case of one or several curves for one unit were
        % bad, cut recording length. Information in
        % outfw; if outfw has value == 1, repetition should get
        % excluded. Either cut recording before or after.
        if find(outfw(i,:) == 1)
            c = find(outfw(i,:) == 1);
            for p =  1 : length(c) %length(find(outfw(i,:) == 1))
                % find gaps between ones and count numbers of
                % repetitions. At least two consecutive
                % repetitions are neccessary to include a unit.
                % In case of several consecutive repetitions
                % divided by a one) choose the one with more
                % repetitions.
                if p == 1
                    len(p) = c(p) - AllAni(i).first;
                    first(p) = AllAni(i).first;
                    last(p) = c(p);
                elseif p == length(c)
                    len(p) = c(p) - 1 - c(p-1);
                    first(p) = c(p-1) + 1;
                    last(p) = c(p);
                    len(p+1) = AllAni(i).last - (c(p) + 1);
                    first(p+1) = c(p) + 1;
                    last(p+1) = AllAni(i).last;
                else
                    len(p) = c(p) - c(p-1);
                    first(p) = c(p-1) + 1;
                    last(p) = c(p);
                end
            end
            lenmax = find(len == max(len)); % position of longest sequence(s)
            if lenmax >= 2 % at least two consecutive repetitions required
                if length(len) > 2 % no consecutive repetitions
                    if length(lenmax) >= 2 % equal long sequences, choose first one
                        AllAni(i).first = first(lenmax(1));
                        AllAni(i).last = last(lenmax(1));
                    else % only one sequence with more than two consecutive repetitions
                        AllAni(i).first = first(lenmax); % start
                        AllAni(i).last = last(lenmax); % end
                    end
                else % exclude unit for current movement direction
                    outfw(i,1:6) = 1;
                end
            else
                outfw(i,1:6) = 1;
            end
        end
    else % NaN
        outfw(i,1:6) = 1;
    end
    
    %% Sorry again. Hier ist das auch nur wieder der exakt gleiche code wie
    % oben, nur für btf anstatt ftb
    clear len
    if AllAni(i).type ~= 0 % for responding units check curves and correct repetitions, first, and last for each unit
        % define start and end repetition for each recording in outfw
        if AllAni(i).first ~= 1
            outbw(i,1:AllAni(i).first-1) = NaN;
        end
        if AllAni(i).last ~= 6
            outbw(i,AllAni(i).last+1:6) = NaN;
        end
        
        % in case of one or several curves for one unit were
        % bad, cut recording length. Information in
        % outbw; if outbw has value == 1, repetition should get
        % excluded. Either cut recording before or after.
        if find(outbw(i,:) == 1)
            c = find(outbw(i,:) == 1);
            for p =  1 : length(c) %length(find(outbw(i,:) == 1))
                % find gaps between ones and count numbers of
                % repetitions. At least two consecutive
                % repetitions are neccessary to include a unit.
                % In case of several consecutive repetitions
                % divided by a one) choose the one with more
                % repetitions.
                if p == 1
                    len(p) = c(p) - AllAni(i).first;
                    first(p) = AllAni(i).first;
                    last(p) = c(p);
                elseif p == length(c)
                    len(p) = c(p) - 1 - c(p-1);
                    first(p) = c(p-1) + 1;
                    last(p) = c(p);
                    len(p+1) = AllAni(i).last - (c(p) + 1);
                    first(p+1) = c(p) + 1;
                    last(p+1) = AllAni(i).last;
                else
                    len(p) = c(p) - c(p-1);
                    first(p) = c(p-1) + 1;
                    last(p) = c(p);
                end
            end
            lenmax = find(len == max(len)); % position of longest sequence(s)
            if lenmax >= 2 % at least two condecutive repetitions required
                if length(len) > 2 % no consecutive repetitions
                    if length(lenmax) >= 2 % equal long sequences, choose first one
                        AllAni(i).first = first(lenmax(1));
                        AllAni(i).last = last(lenmax(1));
                    else % only one sequence with more than two consecutive repetitions
                        AllAni(i).first = first(lenmax); % start
                        AllAni(i).last = last(lenmax); % end
                    end
                else % exclude unit for current movement direction
                    outbw(i,1:6) = 1;
                end
            else
                outbw(i,1:6) = 1;
            end
        end
    else % NaN
        outbw(i,1:6) = 1;
    end
    
    %% if spike rate at smallest temporal frequency higher than spike rate at
    % peak: 2 peaks -> if valley between to peaks greater .5 to peak -> exclude
    field = {'w8'};
    if ~isnan(AllAni(i).first)
        for j = AllAni(i).first : AllAni(i).last
            if isfield(AllAni(i).peak.ftb.pk,field)
                if ~isnan(AllAni(i).peak.ftb.pk(j).w8)
                    if AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.w8(1) > AllAni(i).peak.ftb.pk(j).w8
                        y = AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.fw.w8;
                        x = AllAni(i).xfreq.w8;
                        [pkstemp, locstemp] = findpeaks(-y(~isnan(-y)),x(~isnan(-y)));
                        if (max(abs(pkstemp)))/AllAni(i).peak.ftb.pk(j).w8 > pksValDiffLim % tiefpunkt zwischen peak (erster Messwert) und peak (realpeak) darf nur um 0.5 abweichen vom peak
%                             AllAni(i).peak.ftb.pk(j).(['w',num2str(k)]) = max(pks);
%                             AllAni(i).peak.ftb.loc(j).(['w',num2str(k)]) = locs(pks==max(pks));
%                             AllAni(i).peak.ftb.wid(j).(['w',num2str(k)]) = w(pks==max(pks));
                        else % exclude
                            outfw(i,j) = 1;
                            AllAni(i).peak.ftb.pk(j).(['w',num2str(k)]) = NaN;
                            AllAni(i).peak.ftb.loc(j).(['w',num2str(k)]) = NaN;
                            AllAni(i).peak.ftb.wid(j).(['w',num2str(k)]) = NaN;
                        end
                    end
                end
            end
            
            if isfield(AllAni(i).peak.btf.pk,field)
                if isfield(AllAni(i).peak.btf.pk,field)
                    if ~isnan(AllAni(i).peak.btf.pk(j).w8)
                        if AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.w8(1) > AllAni(i).peak.btf.pk(j).w8
                            y = AllAni(i).(['R0',num2str(j)]).yfreq.mean.translation.bw.w8;
                            x = AllAni(i).xfreq.w8;
                            [pkstemp, locstemp] = findpeaks(-y(~isnan(-y)),x(~isnan(-y)));
                            if (max(abs(pkstemp)))/AllAni(i).peak.btf.pk(j).w8 > pksValDiffLim
%                                 AllAni(i).peak.btf.pk(j).(['w',num2str(k)]) = max(pks);
%                                 AllAni(i).peak.btf.loc(j).(['w',num2str(k)]) = locs(pks==max(pks));
%                                 AllAni(i).peak.btf.wid(j).(['w',num2str(k)]) = w(pks==max(pks));
                            else % exclude
                                outbw(i,j) = 1;
                                AllAni(i).peak.btf.pk(j).(['w',num2str(k)]) = NaN;
                                AllAni(i).peak.btf.loc(j).(['w',num2str(k)]) = NaN;
                                AllAni(i).peak.btf.wid(j).(['w',num2str(k)]) = NaN;
                            end
                        end
                    end
                end
            end
        end
    end
end

%% update all units based on outfw and outbw
for i = cur : size(AllAni,2)
    if nanmean(outfw(i,:)) == 1 && nanmean(outbw(i,:)) == 1 % all curves for both movement directions do not meet the criteria, no response -> response type 0
        AllAni(i).type = 0;
        AllAni(i).first = NaN;
        AllAni(i).last = NaN;
    elseif nanmean(outfw(i,:)) == 1 % only curves for btf movement meet criteria, unidir response to btf -> response type 1
        AllAni(i).type = 1;
    elseif nanmean(outbw(i,:)) == 1 % only curves for ftb movement meet criteria, unidir response to ftb -> response type 2
        AllAni(i).type = 2;
    end
    if ~isnan(AllAni(i).first)
        if AllAni(i).first == AllAni(i).last % only one good curve, update response type to 0
            AllAni(i).type = 0;
            AllAni(i).first = NaN;
            AllAni(i).last = NaN;
        end
    end
end

% case: only one curve meet criteria
for i = cur : size(AllAni,2)
    if length(find(outfw(i,:) == 0) == 1) == 1 && length(find(outbw(i,:) == 0) == 1) == 1 % all curves for both movement directions do not meet the criteria, no response -> response type 0
        AllAni(i).type = 0;
        AllAni(i).first = NaN;
        AllAni(i).last = NaN;
    elseif length(find(outfw(i,:) == 0) == 1) == 1 % only curves for btf movement meet criteria, unidir response to btf -> response type 1
        AllAni(i).type = 1;
    elseif length(find(outbw(i,:) == 0) == 1) == 1 % only curves for ftb movement meet criteria, unidir response to ftb -> response type 2
        AllAni(i).type = 2;
    end
    if ~isnan(AllAni(i).first)
        if AllAni(i).first == AllAni(i).last % only one good curve, update response type to 0
            AllAni(i).type = 0;
            AllAni(i).first = NaN;
            AllAni(i).last = NaN;
        end
    end
end

% set all units that have only one curve to response type 0
for i = cur : size(AllAni,2)
    if ~isnan(AllAni(i).first)
        if AllAni(i).first == AllAni(i).last
            AllAni(i).type = 0;
            AllAni(i).first = NaN;
            AllAni(i).last = NaN;
        end
    else
        AllAni(i).type = 0;
    end
end


% set units with inhibitory responses in any way to NaN (for later analysis
% as potential TN cells)
% vecTN = [13 15 27 28 33 42 53 57 59 62 87 116 118 119 120 122 123 124];
% for i = vecTN
%     AllAni(i).type = 4;
%     AllAni(i).first = 1;
%     AllAni(i).last = 3;
% end

AllAni(108).type = 2;
