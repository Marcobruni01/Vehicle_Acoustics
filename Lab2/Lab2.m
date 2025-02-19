clc
clear
close all
%%

plot_og = true;
plot_filtered = true;
plot_TF = true;
plot_TF2 = true;
plot_avg = true;

%% SBN

%legge e importa i dati sperimentali per ogni microfono 
%per far funzionare il codice dovete inserire il percoso delle tabelle nel vostro pc

mic1 = readtable('students_data/SBN/Cal_Mic1_94dB_1kHz.txt');
mic2 = readtable('students_data/SBN/Cal_Mic2_94dB_1kHz.txt');
mic3 = readtable('students_data/SBN/Cal_Mic3_94dB_1kHz.txt');
mic4 = readtable('students_data/SBN/Cal_Mic4_94dB_1kHz.txt');
manichinoL = readtable('students_data/SBN/Cal_MicML_945dB_1kHz.txt');
manichinoR = readtable('students_data/SBN/Cal_MicMR_945dB_1kHz.txt');

%mics = [mic1, mic2, mic3, mic4, manichinoL, manichinoR];

% calibrazione
% variabili
Pref = 2*10^(-5); %Pa
Pref_dB = [94 94.5]; %dB

Prms = Pref .* 10.^(Pref_dB./20);

% root mean squared di V per confrontarlo con Prms:

% estraggo i valori dei vari V
V(1,:) = table2array(mic1(:,"Var1"));
V(2,:) = table2array(mic2(:,"Var2"));
V(3,:) = table2array(mic3(:,"Var5"));
V(4,:) = table2array(mic4(:,"Var6"));
V(5,:) = table2array(manichinoL(:,"Var3"));
V(6,:) = table2array(manichinoR(:,"Var4"));

% root mean square
V_rms = rms(V,2);

K = zeros(size(V_rms));

% calcoliamo la sensitività
K(1:4) = V_rms(1:4)/Prms(1);
K(5:6) = V_rms(5:6)/Prms(2);

clear mic1 mic2 mic3 mic4 manichinoL manichinoR

%%

empty = readtable('students_data/SBN/Test_Imp_EmptyCabin.txt');
full = readtable('students_data/SBN/Test_Imp_Cabin&passengers.txt');

% Parametri del filtro PASSA-ALTO
fs = 12500; % Frequenza di campionamento
order = 4; % Ordine del filtro
cutoff_freq = 10; % Frequenza di taglio in Hz

mic(1,:) = empty.Var1;
mic(2,:) = empty.Var2;
mic(3,:) = empty.Var5;
mic(4,:) = empty.Var6;
mic(5,:) = empty.Var3;
mic(6,:) = empty.Var4;

acc(1,:) = empty.Var7;   % Accelerometro 1 (Colonna 7)
acc(2,:) = empty.Var8;   % Accelerometro 2 (Colonna 8)
acc(3,:) = empty.Var9;   % Accelerometro 3 (Colonna 9)
acc(4,:) = empty.Var10;  % Accelerometro 4 (Colonna 10)
acc(5,:) = empty.Var11;  % Accelerometro 5 (Colonna 11)
acc(6,:) = empty.Var12;  % Accelerometro 6 (Colonna 12)

% Vibrometro (Laser)
vib = empty.Var14.'; % Vibrometro (Colonna 14)

% Segnale del martello (Hamm)
hamm = empty.Var15.';  % Segnale del martello (Colonna 15)

mic_full(1,:) = full.Var1;
mic_full(2,:) = full.Var2;
mic_full(3,:) = full.Var5;
mic_full(4,:) = full.Var6;
mic_full(5,:) = full.Var3;
mic_full(6,:) = full.Var4;

acc_full(1,:) = full.Var7;   % Accelerometro 1 (Colonna 7)
acc_full(2,:) = full.Var8;   % Accelerometro 2 (Colonna 8)
acc_full(3,:) = full.Var9;   % Accelerometro 3 (Colonna 9)
acc_full(4,:) = full.Var10;  % Accelerometro 4 (Colonna 10)
acc_full(5,:) = full.Var11;  % Accelerometro 5 (Colonna 11)
acc_full(6,:) = full.Var12;  % Accelerometro 6 (Colonna 12)

% Vibrometro (Laser)
vib_full = full.Var14.'; % Vibrometro (Colonna 14)

% Segnale del martello (Hamm)
hamm_full = full.Var15.';  % Segnale del martello (Colonna 15)

hamm_full(326217:326466) = zeros(length(hamm_full(326217:326466)), 1);

clear empty full

%%
t_axis = (0:1:length(hamm)-1)./fs;
t_axis_full = (0:1:length(hamm_full)-1)./fs;

if plot_og
    figure
    sgtitle("Empty Cabin")
    subplot(5,1,1)
    plot(t_axis, mic(1:4,:));
    title("Microphone signals")
    
    subplot(5,1,2)
    plot(t_axis, mic(5:6,:));
    title("Torso + Head signals")
    
    subplot(5,1,3)
    plot(t_axis, acc);
    title("Accelerometer signals")
    
    subplot(5,1,4)
    plot(t_axis, vib);
    title("Vibrometer signals")
    
    subplot(5,1,5)
    plot(t_axis, hamm);
    title("Impact hammer input")
    
    % full cabin

    figure
    sgtitle("Empty Cabin")
    subplot(5,1,1)
    plot(t_axis_full, mic_full(1:4,:));
    title("Microphone signals")
    
    subplot(5,1,2)
    plot(t_axis_full, mic_full(5:6,:));
    title("Torso + Head signals")
    
    subplot(5,1,3)
    plot(t_axis_full, acc_full);
    title("Accelerometer signals")
    
    subplot(5,1,4)
    plot(t_axis_full, vib_full);
    title("Vibrometer signals")
    
    subplot(5,1,5)
    plot(t_axis_full, hamm_full);
    title("Impact hammer input")
end
   
%%

% Design del filtro passa-alto
[b, a] = butter(order, cutoff_freq / (fs / 2), 'high');

% Filtra i segnali per la cabina vuota
filtered_hamm = filtfilt(b, a, hamm); % Segnale del martello filtrato
filtered_vib = filtfilt(b, a, vib); % Segnale del vibrometro filtrato

filtered_mic = zeros(size(mic));
filtered_acc = zeros(size(acc));

for i = 1:height(mic)
    filtered_mic(i,:) = filtfilt(b, a, mic(i,:)); % Microfono i filtrato
    filtered_acc(i,:) = filtfilt(b, a, acc(i,:)); % Accellerometro i filtrato
end

% Filtra i segnali per la cabina piena
filtered_hamm_full = filtfilt(b, a, hamm_full); % Segnale del martello filtrato
filtered_vib_full = filtfilt(b, a, vib_full); % Segnale del vibrometro filtrato

filtered_mic_full = zeros(size(mic_full));
filtered_acc_full = zeros(size(acc_full));

for i = 1:height(mic)
    filtered_mic_full(i,:) = filtfilt(b, a, mic_full(i,:)); % Microfono i filtrato
    filtered_acc_full(i,:) = filtfilt(b, a, acc_full(i,:)); % Accellerometro i filtrato
end

%%
if plot_filtered
    figure;
    sgtitle("Empty cabin")
    subplot(5,1,1);
    plot(t_axis, filtered_hamm);
    title('Hammer Filtered Signal');
    subplot(5,1,2);
    plot(t_axis, filtered_mic(1:4,:));
    title('Microphone Filtered Signals');
    subplot(5,1,3);
    plot(t_axis, filtered_mic(5:6,:));
    title('Head+Torso Filtered Signals');
    subplot(5,1,4);
    plot(t_axis, filtered_acc);
    title('Accelerometer Filtered Signals');
    subplot(5,1,5);
    plot(t_axis, filtered_vib);
    title('Vibrometer Filtered Signal');

    figure;
    sgtitle("Full cabin")
    subplot(5,1,1);
    plot(t_axis_full, filtered_hamm_full);
    title('Hammer Filtered Signal');
    subplot(5,1,2);
    plot(t_axis_full, filtered_mic_full(1:4,:));
    title('Microphone Filtered Signals');
    subplot(5,1,3);
    plot(t_axis_full, filtered_mic_full(5:6,:));
    title('Head+Torso Signals');
    subplot(5,1,4);
    plot(t_axis_full, filtered_acc_full);
    title('Accelerometer Filtered Signals');
    subplot(5,1,5);
    plot(t_axis_full, filtered_vib_full);
    title('Vibrometer Filtered Signal');
end

%% compute Transfer Functions unsing TFestimate and different windows
windnames = ["Rectangular", "Triangular", "Hamming", "Blackman"];

sz = 4096;

winds(1,:) = ones(sz,1);
% winds(2,:) = bartlett(sz);
% winds(3,:) = hamming(sz);
% winds(4,:) = blackman(sz);

% empty cabin

[mic_TFs, f_axis] = calculateInpulseTFs(filtered_hamm,filtered_mic./K,winds,length(winds)/2,fs);
[acc_TFs, ~] = calculateInpulseTFs(filtered_hamm,filtered_acc,winds,length(winds)/2,fs);
[vib_TFs, ~] = calculateInpulseTFs(filtered_hamm,filtered_vib,winds,length(winds)/2,fs);

mic_TFs = mic_TFs(:,:,:,1:end-1);
acc_TFs = acc_TFs(:,:,:,1:end-1);
vib_TFs = vib_TFs(:,:,:,1:end-1);

f_axis = f_axis(1:end-1);

mic_TFs = mic_TFs./max(abs(mic_TFs),[],4);
acc_TFs = acc_TFs./max(abs(acc_TFs),[],4);
vib_TFs = vib_TFs./max(abs(vib_TFs),[],4);
% full cabin

[mic_TFs_full, f_axis_full] = calculateInpulseTFs(filtered_hamm_full,filtered_mic_full./K,winds,length(winds)/2,fs);
[acc_TFs_full, ~] = calculateInpulseTFs(filtered_hamm_full,filtered_acc_full,winds,length(winds)/2,fs);
[vib_TFs_full, ~] = calculateInpulseTFs(filtered_hamm_full,filtered_vib_full,winds,length(winds)/2,fs);

f_axis_full = f_axis_full(1:end-1);

mic_TFs_full = mic_TFs_full(:,:,:,1:end-1);
acc_TFs_full = acc_TFs_full(:,:,:,1:end-1);
vib_TFs_full = vib_TFs_full(:,:,:,1:end-1);

mic_TFs_full = mic_TFs_full./max(abs(mic_TFs_full),[],4);
acc_TFs_full = acc_TFs_full./max(abs(acc_TFs_full),[],4);
vib_TFs_full = vib_TFs_full./max(abs(vib_TFs_full),[],4);

%% compute Transfer functions as Fourier(Out)/Fourier(In)

% empty cabin
[hamm_FTs, fft_axis] = computeFTs(filtered_hamm,filtered_hamm,sz);
[mic_FTs, ~] = computeFTs(filtered_hamm,filtered_mic./K,sz);
[acc_FTs, ~] = computeFTs(filtered_hamm,filtered_acc,sz);
[vib_FTs, ~] = computeFTs(filtered_hamm,filtered_vib,sz);

mic_TFs_2 = mic_FTs./hamm_FTs; 
acc_TFs_2 = acc_FTs./hamm_FTs; 
vib_TFs_2 = vib_FTs./hamm_FTs;

mic_TFs_2 = mic_TFs_2./max(abs(mic_TFs_2),[],3);
acc_TFs_2 = acc_TFs_2./max(abs(acc_TFs_2),[],3);
vib_TFs_2 = vib_TFs_2./max(abs(vib_TFs_2),[],3);

% full cabin
[hamm_FTs_full, fft_axis_full] = computeFTs(filtered_hamm_full,filtered_hamm_full,sz);
[mic_FTs_full, ~] = computeFTs(filtered_hamm_full,filtered_mic_full./K,sz);
[acc_FTs_full, ~] = computeFTs(filtered_hamm_full,filtered_acc_full,sz);
[vib_FTs_full, ~] = computeFTs(filtered_hamm_full,filtered_vib_full,sz);

mic_TFs_full_2 = mic_FTs_full./hamm_FTs_full; 
acc_TFs_full_2 = acc_FTs_full./hamm_FTs_full; 
vib_TFs_full_2 = vib_FTs_full./hamm_FTs_full; 

mic_TFs_full_2 = mic_TFs_full_2./max(abs(mic_TFs_full_2),[],3);
acc_TFs_full_2 = acc_TFs_full_2./max(abs(acc_TFs_full_2),[],3);
vib_TFs_full_2 = vib_TFs_full_2./max(abs(vib_TFs_full_2),[],3);

%% compute the coherence
for i = 1:height(winds)
    vib_coherence(:,:,i) = reshape(mscohere(squeeze(permute(vib_TFs(:,:,i,:), [4 2 3 1])), squeeze(permute(vib_TFs_2, [3 1 2]))),512,13,1);
    mic_coherence(:,:,i,:) = reshape(mscohere(squeeze(permute(mic_TFs(:,:,i,:), [4 2 3 1])), squeeze(permute(mic_TFs_2, [3 2 1]))),512,13,1,6);
    acc_coherence(:,:,i,:) = reshape(mscohere(squeeze(permute(acc_TFs(:,:,i,:), [4 2 3 1])), squeeze(permute(acc_TFs_2, [3 2 1]))),512,13,1,6);

    vib_coherence_full(:,:,i) = reshape(mscohere(squeeze(permute(vib_TFs_full(:,:,i,:), [4 2 3 1])), squeeze(permute(vib_TFs_full_2, [3 1 2]))),512,14,1);
    mic_coherence_full(:,:,i,:) = reshape(mscohere(squeeze(permute(mic_TFs_full(:,:,i,:), [4 2 3 1])), squeeze(permute(mic_TFs_full_2, [3 2 1]))),512,14,1,6);
    acc_coherence_full(:,:,i,:) = reshape(mscohere(squeeze(permute(acc_TFs_full(:,:,i,:), [4 2 3 1])), squeeze(permute(acc_TFs_full_2, [3 2 1]))),512,14,1,6);
end

%%
if plot_TF

    figure
    sgtitle("Transfer Function for hit number 8 (empty cabin)")
    for j = 1:height(mic_TFs)
        subplot(6,2,2*(j-1)+1);
        plot(f_axis, db(abs(squeeze(mic_TFs(j,8,:,:)))))
        hold on
        plot(fft_axis, db(abs(squeeze(mic_TFs_2(j,8,:)))), ':')
        xlim([0 500])
        title("Magnitude for Microphone: " + j)
        subplot(6,2,2*(j-1)+2);
        plot(f_axis, angle(squeeze(mic_TFs(j,8,:,:))))
        hold on
        plot(fft_axis, angle(squeeze(mic_TFs_2(j,8,:))), ':')
        xlim([0 500])
        title("Phase for Microphone: " + j)
    end
    legend(["tfestimate", "FT division"])

    figure
    sgtitle("Coherence of the 2 transfer functions for microphones (empty cabin)")
    for j = 1:height(mic_TFs)
        subplot(3,2,j);
        plot(f_axis(1:4:end), abs(squeeze(mic_coherence(:,8,:,j))))
        xlim([0 500])
        title("Microphone: " + j)
    end
    
    figure
    sgtitle("Transfer Function for hit number 8 (empty cabin)")
    for j = 1:height(acc_TFs)
        subplot(6,2,2*(j-1)+1);
        plot(f_axis, db(abs(squeeze(acc_TFs(j,8,:,:)))))
        hold on
        plot(fft_axis, db(abs(squeeze(acc_TFs_2(j,8,:)))), ':')
        xlim([0 500])
        title("Magnitude for Accelerometer: " + j)
        subplot(6,2,2*(j-1)+2);
        plot(f_axis, angle(squeeze(acc_TFs(j,8,:,:))))
        hold on
        plot(fft_axis, angle(squeeze(acc_TFs_2(j,8,:))), ':')
        xlim([0 500])
        title("Phase for Accelerometer: " + j)
    end
    legend(["tfestimate", "FT division"])

    figure
    sgtitle("Coherence of the 2 transfer functions for accelerometers (empty cabin)")
    for j = 1:height(mic_TFs)
        subplot(3,2,j);
        plot(f_axis(1:4:end), abs(squeeze(acc_coherence(:,8,:,j))))
        xlim([0 500])
        title("Accelerometer: " + j)
    end

    figure
    sgtitle("Vibrometer Transfer Function for hit number 8 (empty cabin)")
    subplot(2,1,1)
    title("Magnitude")
    plot(f_axis, db(abs(squeeze(vib_TFs(:,8,:,:)))))
    hold on
    plot(fft_axis, db(abs(squeeze(vib_TFs_2(:,8,:)))), ':')
    xlim([0 500])
    subplot(2,1,2)
    title("Phase")
    plot(f_axis, angle(squeeze(vib_TFs(:,8,:,:))))
    hold on
    plot(fft_axis, angle(squeeze(vib_TFs_2(:,8,:))), ':')
    xlim([0 500])
    legend(["tfestimate", "FT division"])

    figure
    sgtitle("Coherence of the 2 transfer functions for vibrometer (empty cabin)")
    plot(f_axis(1:4:end), abs(squeeze(vib_coherence(:,8,:))))
    xlim([0 500])

    % full cabin

    figure
    sgtitle("Transfer Function for hit number 9 (full cabin)")
    for j = 1:height(mic_TFs)
        subplot(6,2,2*(j-1)+1);
        plot(f_axis, db(abs(squeeze(mic_TFs_full(j,9,:,:)))))
        hold on
        plot(fft_axis, db(abs(squeeze(mic_TFs_full_2(j,9,:)))), ':')
        xlim([0 500])
        title("Magnitude for Microphone: " + j)
        subplot(6,2,2*(j-1)+2);
        plot(f_axis, angle(squeeze(mic_TFs_full(j,9,:,:))))
        hold on
        plot(fft_axis, angle(squeeze(mic_TFs_full_2(j,9,:))), ':')
        xlim([0 500])
        title("Phase for Microphone: " + j)
    end
    legend(["tfestimate", "FT division"])

    figure
    sgtitle("Coherence of the 2 transfer functions for microphones (full cabin)")
    for j = 1:height(mic_TFs)
        subplot(3,2,j);
        plot(f_axis(1:512), abs(squeeze(mic_coherence_full(:,9,:,j))))
        xlim([0 500])
        title("Microphone: " + j)
    end
    
    figure
    sgtitle("Transfer Function for hit number 9 (full cabin)")
    for j = 1:height(acc_TFs)
        subplot(6,2,2*(j-1)+1);
        plot(f_axis, db(abs(squeeze(acc_TFs_full(j,9,:,:)))))
        hold on
        plot(fft_axis, db(abs(squeeze(acc_TFs_full_2(j,9,:)))), ':')
        xlim([0 500])
        title("Magnitude for Accelerometer: " + j)
        subplot(6,2,2*(j-1)+2);
        plot(f_axis, angle(squeeze(acc_TFs_full(j,9,:,:))))
        hold on
        plot(fft_axis, angle(squeeze(acc_TFs_full_2(j,9,:))), ':')
        xlim([0 500])
        title("Phase for Accelerometer: " + j)
    end
    legend(["tfestimate", "FT division"])

    figure
    sgtitle("Coherence of the 2 transfer functions for accelerometers (full cabin)")
    for j = 1:height(mic_TFs)
        subplot(3,2,j);
        plot(f_axis(1:512), abs(squeeze(acc_coherence_full(:,9,:,j))))
        xlim([0 500])
        title("Accelerometer: " + j)
    end

    figure
    sgtitle("Vibrometer Transfer Function for hit number 9 (full cabin)")
    subplot(2,1,1)
    title("Magnitude")
    plot(f_axis, db(abs(squeeze(vib_TFs_full(:,9,:,:)))))
    hold on
    plot(fft_axis, db(abs(squeeze(vib_TFs_full_2(:,9,:)))), ':')
    xlim([0 500])
    subplot(2,1,2)
    title("Phase")
    plot(f_axis, angle(squeeze(vib_TFs_full(:,9,:,:))))
    hold on
    plot(fft_axis, angle(squeeze(vib_TFs_full_2(:,9,:))), ':')
    xlim([0 500])
    legend(["tfestimate", "FT division"])

    figure
    sgtitle("Coherence of the 2 transfer functions for vibrometer (full cabin)")
    plot(f_axis(1:512), abs(squeeze(vib_coherence_full(:,9,:))))
    xlim([0 500])

end

%% calculate the average for the TFs over the inpulses (using the TF estimate obtained ones (can be changed))
avg_mic_TFs = squeeze(sum(mic_TFs.^2,2)./width(mic_TFs)).^(1/2);
avg_acc_TFs = squeeze(sum(acc_TFs.^2,2)./width(acc_TFs)).^(1/2);
avg_vib_TFs = squeeze(sum(vib_TFs.^2,2)./width(vib_TFs)).^(1/2);
avg_mic_TFs_full = squeeze(sum(mic_TFs_full.^2,2)./width(mic_TFs_full)).^(1/2);
avg_acc_TFs_full = squeeze(sum(acc_TFs_full.^2,2)./width(acc_TFs_full)).^(1/2);
avg_vib_TFs_full = squeeze(sum(vib_TFs_full.^2,2)./width(vib_TFs_full)).^(1/2);

%%
if plot_avg
    figure
    sgtitle("Average Transfer Function for Microphones (empty cabin)")
    for j = 1:height(avg_mic_TFs)
        subplot(3,2,j);
        plot(f_axis, db(abs(squeeze(mic_TFs(j,:,:,:)))),'Color', '#666666')
        hold on 
        plot(f_axis, db(abs(squeeze(avg_mic_TFs(j,:,:)))), 'LineWidth', 1.5,'Color', 'green')
        xlim([0 500])
        title("Microphone: " + j)
        %legend(windnames)
    end
    
    figure
    sgtitle("Average Transfer Function for Accelerometers (empty cabin)")
    for j = 1:height(avg_acc_TFs)
        subplot(3,2,j);
        plot(f_axis, db(abs(squeeze(acc_TFs(j,:,:,:)))),'Color', '#666666')
        hold on 
        plot(f_axis, db(abs(squeeze(avg_acc_TFs(j,:,:)))), 'LineWidth', 1.5,'Color', 'green')
        xlim([0 500])
        title("Accelerometer: " + j)
        %legend(windnames)
    end
    
    figure
    sgtitle("Average Transfer Function for Vibrometer (empty cabin)")
    plot(f_axis, db(abs(squeeze(vib_TFs(:,:,:,:)))),'Color', '#666666')
    hold on 
    plot(f_axis, db(abs(squeeze(avg_vib_TFs(:,:,:)))), 'LineWidth', 1.5,'Color', 'green')
    xlim([0 1500])
    %legend(windnames)

    % full cabin

    figure
    sgtitle("Average Transfer Function for Microphones (full cabin)")
    for j = 1:height(avg_mic_TFs)
        subplot(3,2,j);
        plot(f_axis, db(abs(squeeze(mic_TFs_full(j,:,:,:)))),'Color', '#666666')
        hold on 
        plot(f_axis, db(abs(squeeze(avg_mic_TFs_full(j,:,:)))), 'LineWidth', 1.5,'Color', 'green')
        xlim([0 500])
        title("Microphone: " + j)
        %legend(windnames)
    end
    
    figure
    sgtitle("Average Transfer Function for Accelerometers (full  cabin)")
    for j = 1:height(avg_acc_TFs)
        subplot(3,2,j);
        plot(f_axis, db(abs(squeeze(acc_TFs_full(j,:,:,:)))),'Color', '#666666')
        hold on 
        plot(f_axis, db(abs(squeeze(avg_acc_TFs_full(j,:,:)))), 'LineWidth', 1.5,'Color', 'green')
        xlim([0 500])
        title("Accelerometer: " + j)
        %legend(windnames)
    end
    
    figure
    sgtitle("Average Transfer Function for Vibrometer (fulls cabin)")
    plot(f_axis, db(abs(squeeze(vib_TFs_full(:,:,:,:)))),'Color', '#666666')
    hold on 
    plot(f_axis, db(abs(squeeze(avg_vib_TFs_full(:,:,:)))), 'LineWidth', 1.5,'Color', 'green')
    xlim([0 1500])
    %legend(windnames)
end
%%
f_min = f_axis(2);
f_max = f_axis(end);

k = 0;
f_center = [];
while true
    f_current = f_min * (2^(k/3)); % Calculate center frequency for band index k
    if f_current > f_max
        break; % Stop if the current center frequency exceeds the maximum frequency
    end
    f_center = [f_center, f_current]; % Append the center frequency to the list
    k = k + 1; % Increment band index
end

mic_third = zeros(length(f_center),height(avg_mic_TFs));
acc_third = zeros(length(f_center),height(avg_acc_TFs));
vib_third = zeros(length(f_center),width(avg_vib_TFs));

mic_third_full = zeros(length(f_center),height(avg_mic_TFs_full));
acc_third_full = zeros(length(f_center),height(avg_acc_TFs_full));
vib_third_full = zeros(length(f_center),width(avg_vib_TFs_full));

for i = 1:length(f_center)
    % Compute lower and upper band edges
    f_lower = f_center(i) / (2^(1/6));
    f_upper = f_center(i) * (2^(1/6));
    
    % Find indices of frequencies within the current band
    idx = (f_axis >= f_lower) & (f_axis <= f_upper);
    idx_full = (f_axis_full >= f_lower) & (f_axis_full <= f_upper);
    
    % Calculate the average energy (or magnitude) in the current band
    if any(idx) % Ensure there are data points in this bands
        app1 = avg_mic_TFs(:,idx);
        app2 = avg_acc_TFs(:,idx);
        app3 = avg_vib_TFs(idx,:);
        mic_third(i,:) = sqrt(mean(app1.^2,2));
        acc_third(i,:) = sqrt(mean(app2.^2,2));
        vib_third(i,:) = sqrt(mean(app3.^2,1));
    else
    mic_third(i,:) = zeros(1,height(avg_mic_TFs));
    acc_third(i,:) = zeros(1,height(avg_acc_TFs));
    vib_third(i,:) = zeros(1,width(avg_vib_TFs));
    end

    if any(idx_full) % Ensure there are data points in this band
        app1 = avg_mic_TFs_full(:,idx_full);
        app2 = avg_acc_TFs_full(:,idx_full);
        app3 = avg_vib_TFs_full(idx_full,:);
        mic_third_full(i,:) = sqrt(mean(app1.^2,2));
        acc_third_full(i,:) = sqrt(mean(app2.^2,2));
        vib_third_full(i,:) = sqrt(mean(app3.^2,1));
    else
    mic_third_full(i,:) = zeros(1,height(avg_mic_TFs_full));
    acc_third_full(i,:) = zeros(1,height(avg_acc_TFs_full));
    vib_third_full(i,:) = zeros(1,width(avg_vib_TFs_full));
    end
end

%%

figure;
bar(1:length(f_center),abs([mic_third(:,1) mic_third_full(:,1)]));
xlabel('Frequency (Hz)');
ylabel('Band Magnitude (dB)');
title('One-Third Octave Band Analysis');
grid on;
xticks(1:length(f_center)); % Set x-ticks at integer indices
xticklabels(arrayfun(@(x) sprintf('%.1f Hz', x), f_center, 'UniformOutput', false));
xlim([8.5, 23.5]);
xtickangle(45);

%%

ramp = readtable('students_data/SBN/Test_Ramp.txt');

mic_ramp(1,:) = ramp.Var1;
mic_ramp(2,:) = ramp.Var2;
mic_ramp(3,:) = ramp.Var5;
mic_ramp(4,:) = ramp.Var6;
mic_ramp(5,:) = ramp.Var3;
mic_ramp(6,:) = ramp.Var4;

acc_ramp(1,:) = ramp.Var7;   
acc_ramp(2,:) = ramp.Var8;   
acc_ramp(3,:) = ramp.Var9;  
acc_ramp(4,:) = ramp.Var10;  
acc_ramp(5,:) = ramp.Var11;  
acc_ramp(6,:) = ramp.Var12;

clear ramp



figure
sgtitle("Microphones spectrogram")
for i = 1:height(mic_ramp)
    [Spec_mic, f_mic, t_mic] = spectrogram(mic_ramp(i,:), 1024,512,2048,fs);
    subplot(3,2,i)
    imagesc(t_mic, f_mic, db(abs(Spec_mic))); % Time and frequency axes
    axis xy
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title("Microphone " + i);
    set(gca, 'YDir', 'normal');
    ylim([0 500])
    xlim([7 20])
    colorbar;
end

figure
sgtitle("Accelerometer spectrogram")
for i = 1:height(acc_ramp)
    [Spec_acc, f_acc, t_acc] = spectrogram(acc_ramp(i,:), 1024,512,2048,fs);
    subplot(3,2,i)
    imagesc(t_acc, f_acc, db(abs(Spec_acc))); % Time and frequency axes
    axis xy
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title("Accelerometer " + i);
    set(gca, 'YDir', 'normal');
    ylim([0 500])
    xlim([7 20])
    colorbar;
end

%% ABN

test1 = readtable('students_data/ABN/Test1.txt');
test2 = readtable('students_data/ABN/Test2.txt');
test3 = readtable('students_data/ABN/Test3.txt');
test4 = readtable('students_data/ABN/Test4.txt');

test1_mics(1,:) = test1.Var1;
test1_mics(2,:) = test1.Var2;
test1_mics(3,:) = test1.Var5;
test1_mics(4,:) = test1.Var6;
test1_mics(5,:) = test1.Var3;
test1_mics(6,:) = test1.Var4;

test2_mics(1,:) = test2.Var1;
test2_mics(2,:) = test2.Var2;
test2_mics(3,:) = test2.Var5;
test2_mics(4,:) = test2.Var6;
test2_mics(5,:) = test2.Var3;
test2_mics(6,:) = test2.Var4;

test3_mics(1,:) = test3.Var1;
test3_mics(2,:) = test3.Var2;
test3_mics(3,:) = test3.Var5;
test3_mics(4,:) = test3.Var6;
test3_mics(5,:) = test3.Var3;
test3_mics(6,:) = test3.Var4;

test4_mics(1,:) = test4.Var1;
test4_mics(2,:) = test4.Var2;
test4_mics(3,:) = test4.Var5;
test4_mics(4,:) = test4.Var6;
test4_mics(5,:) = test4.Var3;
test4_mics(6,:) = test4.Var4;

Qsource(1,:) = test1.Var13;
Qsource(2,:) = test2.Var13;
Qsource(3,:) = test3.Var13;
Qsource(4,:) = test4.Var13;

clear test1 test2 test3 test4

[b, a] = butter(order, cutoff_freq / (fs / 2), 'high');

[test1_FTs, f_axis] = ffg(filtfilt(b,a,test1_mics.'), sz, 1/fs);
[test2_FTs, ~] = ffg(filtfilt(b,a,test2_mics.'), sz, 1/fs);
[test3_FTs, ~] = ffg(filtfilt(b,a,test3_mics.'), sz, 1/fs);
[test4_FTs, ~] = ffg(filtfilt(b,a,test4_mics.'), sz, 1/fs);
[Q_FTs, ~] = ffg(filtfilt(b,a,Qsource.'), sz, 1/fs);

test_FTs(1,:,:) = (test1_FTs./max(test1_FTs,[],1)).';
test_FTs(2,:,:) = (test2_FTs./max(test2_FTs,[],1)).';
test_FTs(3,:,:) = (test3_FTs./max(test3_FTs,[],1)).';
test_FTs(4,:,:) = (test4_FTs./max(test4_FTs,[],1)).';

clear test1_FTs test2_FTs test3_FTs test4_FTs

Q_FTs = (Q_FTs./max(Q_FTs,[],1)).';

test_TFs = zeros(size(test_FTs));

for i = 1:height(test_FTs)
    test_TFs(i,:,:) = squeeze(test_FTs(i,:,:))./Q_FTs(i,:);
    test_TFs(i,:,:) = test_TFs(i,:,:)./max(test_TFs(i,:,:),[],3);
end

%%
for i = 1:size(test_TFs,1)
    for j = 1:size(test_TFs,2)
        coherence(i,j,:) = mscohere(squeeze(test_FTs(i,j,:)),Q_FTs(i,:));
    end
end

%%

for i = 1:height(test_TFs)
    figure
    sgtitle("Transfer Functions for Test " + i)
    for j = 1:width(test_TFs)
        if(mod(j,2) == 0)
            amp_idx = 2*(j-2) + 2;
            ph_idx = 2*(j-2) + 4;
        else
            amp_idx = 2*(j-1) + 1;
            ph_idx = 2*(j-1) + 3;
        end
        subplot(6,2,amp_idx)
        plot(f_axis,db(abs(squeeze(test_TFs(i,j,:)))))
        title("TF for microphone " + j)
        xlim([0 5000])

        subplot(6,2,ph_idx)
        plot(f_axis,angle(squeeze(test_TFs(i,j,:))))
        xlim([0 5000])
    end

    figure
    sgtitle("Coherence with Qsource for Test " + i)
    for j = 1:width(test_TFs)
        subplot(3,2,j)
        plot(f_axis(1:4:end),squeeze(coherence(i,j,:)))
        title("coherence for microphone " + j)
        xlim([0 5000])

    end
end

%%
f_min = f_axis(2);
f_max = f_axis(end);

k = 0;
f_center = [];
while true
    f_current = f_min * (2^(k/3)); % Calculate center frequency for band index k
    if f_current > f_max
        break; % Stop if the current center frequency exceeds the maximum frequency
    end
    f_center = [f_center, f_current]; % Append the center frequency to the list
    k = k + 1; % Increment band index
end

third = zeros([size(test_TFs,[1 2]),length(f_center)]);

for i = 1:length(f_center)
    % Compute lower and upper band edges
    f_lower = f_center(i) / (2^(1/6));
    f_upper = f_center(i) * (2^(1/6));
    
    % Find indices of frequencies within the current band
    idx = (f_axis >= f_lower) & (f_axis <= f_upper);
    
    % Calculate the average energy (or magnitude) in the current band
    if any(idx) % Ensure there are data points in this bands
        app = test_TFs(:,:,idx);
        third(:,:,i) = sqrt(mean(app.^2,3));
    else
    third(:,:,i) = zeros([size(test_TFs,[1 2]),1]);
    end
end

%%

figure;
bar(1:length(f_center),abs([squeeze(third(3,1,:)) squeeze(third(3,5,:)) squeeze(third(3,6,:))]));
xlabel('Frequency (Hz)');
ylabel('Band Magnitude (dB)');
title('One-Third Octave Band Analysis');
grid on;
xticks(1:length(f_center)); % Set x-ticks at integer indices
xticklabels(arrayfun(@(x) sprintf('%.1f Hz', x), f_center, 'UniformOutput', false));
xlim([5.5, 33.5]);
xtickangle(45);
legend("Microphone 1", "Head and Torso Left", "Head and Tprsp Right")

%%

figure;
bar(1:length(f_center),abs([squeeze(third(3,1,:)) squeeze(third(4,1,:))]));
xlabel('Frequency (Hz)');
ylabel('Band Magnitude (dB)');
title('One-Third Octave Band Analysis');
grid on;
xticks(1:length(f_center)); % Set x-ticks at integer indices
xticklabels(arrayfun(@(x) sprintf('%.1f Hz', x), f_center, 'UniformOutput', false));
xlim([5.5, 33.5]);
xtickangle(45);
legend("Test 3", "Test 4")