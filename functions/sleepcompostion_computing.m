function [Sleep_Composition, wake_binary_vector, sws_binary_vector, REM_binary_vector, NREMinclMA_binary_vector,MA_binary_vector,wake_woMA_binary_vector] = sleepcompostion_computing(FullHypno)

wake_binary_vector = FullHypno==1;
sws_binary_vector = FullHypno==2;
REM_binary_vector = FullHypno==3;

[wake_onset, wake_offset] = binary_to_OnOff(wake_binary_vector);
[sws_onset, sws_offset] = binary_to_OnOff(sws_binary_vector);
[REM_onset, REM_offset] = binary_to_OnOff(REM_binary_vector);

wake_duration = wake_offset-wake_onset;
sws_duration = sws_offset-sws_onset;
REM_duration = REM_offset-REM_onset;

MA_maxdur = 15; % maximum duration of microarrousal
MA_idx = find(wake_duration < MA_maxdur);
MA_onset = wake_onset(MA_idx);
MA_duration = wake_duration(MA_idx);
MA_binary_vector =zeros([1,length(wake_binary_vector)]);
for i=1:length(MA_onset) % making time vector for EEG scoring (frequency = 1Hz)
    t = MA_onset(i)+1;
    d = MA_duration(i)-1;
    MA_binary_vector(t:t+d) = 1;
end

% remove micrarrousal from wake vectors
wake_woMA_onset = wake_onset;
wake_woMA_onset(MA_idx) = [];
wake_woMA_duration = wake_duration;
wake_woMA_duration(MA_idx) = [];
wake_woMA_binary_vector = zeros([1,length(wake_binary_vector)]);
for i=1:length(wake_woMA_onset) % making time vector for EEG scoring (frequency = 1Hz)
    t = wake_woMA_onset(i)+1;
    d = wake_woMA_duration(i)-1;
    wake_woMA_binary_vector(t:t+d) = 1;
end

%State transitions (uncut vectors)
% Creating one vector with different behaviors represented by unique
% numbers (1=wake, 4=sws, 9=REM, 15=MA) at frequency 1Hz
boutscore_vector = zeros([1,length(wake_binary_vector)]);

%Here using the unaligned "uncut" vectors
for i=1:length(wake_woMA_onset)
    t = wake_woMA_onset(i)+1;
    d = wake_woMA_duration(i)-1;
    boutscore_vector(t:t+d) = 1; % wake=1
end

for i=1:length(sws_onset)
    t = sws_onset(i)+1;
    d = sws_duration(i)-1;
    boutscore_vector(t:t+d) = 4; % sws=4
end

if ~isnan(REM_onset)
    for i=1:length(REM_onset)
        t = REM_onset(i)+1;
        d = REM_duration(i)-1;
        boutscore_vector(t:t+d) = 9; %REM=9
    end
end

for i=1:length(MA_onset)
    t = MA_onset(i)+1;
    d = MA_duration(i)-1;
    boutscore_vector(t:t+d) = 15; %MA=15
end


NREMinclMA_binary_vector = boutscore_vector==4 | boutscore_vector==15;
[NREMinclMA_onset, NREMinclMA_offset] = binary_to_OnOff(NREMinclMA_binary_vector);
NREMinclMA_duration = NREMinclMA_offset-NREMinclMA_onset;

total_duration = length(FullHypno); % second

NREM_percentage = sum(sws_duration)./total_duration *100;
Wake_percentage = sum(wake_woMA_duration)./total_duration*100;
REM_percentage = sum(REM_duration)./total_duration*100;
NREMinclMA_percentage = sum(NREMinclMA_duration)./total_duration*100;
MA_percentage = sum(MA_duration)./total_duration*100;

NREM_bouts = length(sws_duration)./(total_duration/3600); % #/per hour
Wake_bouts = length(wake_woMA_duration)./(total_duration/3600);
REM_bouts = length(REM_duration)./(total_duration/3600);
NREMinclMA_bouts = length(NREMinclMA_duration)./(total_duration/3600);
MA_bouts = length(MA_duration)./(total_duration/3600);

NREM_bout_duration = mean(sws_duration)/60; %min
Wake_bout_duration = mean(wake_woMA_duration)/60; %min
REM_bout_duration = mean(REM_duration)/60; %min
NREMinclMA_bout_duration = mean(NREMinclMA_duration)/60; %min
MA_bout_duration = mean(MA_duration); % second

MA_NREM = length(MA_duration)*60./sum(sws_duration);

Sleep_Composition = [ ...
    NREM_percentage, ...
    Wake_percentage, ...
    REM_percentage, ...
    NREMinclMA_percentage, ...
    MA_percentage, ...
    NREM_bouts, ...
    Wake_bouts, ...
    REM_bouts, ...
    NREMinclMA_bouts, ...
    MA_bouts, ...
    NREM_bout_duration, ...
    Wake_bout_duration, ...
    REM_bout_duration, ...
    NREMinclMA_bout_duration, ...
    MA_bout_duration, ...
    MA_NREM ...
];
