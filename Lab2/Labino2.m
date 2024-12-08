clc
clear
close all
%%

plot_og = false;
plot_filtered = true;
plot_TF = false;
plot_TF2 = false;
plot_avg = true;

%%

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

%%

empty = readtable('students_data/SBN/Test_Imp_EmptyCabin.txt');
full = readtable('students_data/SBN/Test_Imp_Cabin&passengers.txt');

%%

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

clear empty full

%%
t_axis = (0:1:length(hamm)-1)./fs;
t_axis_full = (0:1:length(hamm_full)-1)./fs;

if plot_og
    figure 
    subplot(5,1,1)
    title("Microphone signals")
    plot(t_axis, mic(1:4,:));
    
    subplot(5,1,2)
    title("Torso + Head signals")
    plot(t_axis, mic(5:6,:));
    
    subplot(5,1,3)
    title("Accelerometer signals")
    plot(t_axis, acc);
    
    subplot(5,1,4)
    title("Vibrometer signals")
    plot(t_axis, vib);
    
    subplot(5,1,5)
    title("Impact hammer input")
    plot(t_axis, hamm);
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
    title('Head+Torso Signals');
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

winds(1,:) = ones(4096,1);
winds(2,:) = bartlett(4096);
winds(3,:) = hamming(4096);
winds(4,:) = blackman(4096);

% empty cabin

[mic_TFs, f_axis] = calculateInpulseTFs(filtered_hamm,filtered_mic./K,winds,length(winds)/2,fs);
[acc_TFs, ~] = calculateInpulseTFs(filtered_hamm,filtered_acc,winds,length(winds)/2,fs);
[vib_TFs, ~] = calculateInpulseTFs(filtered_hamm,filtered_vib,winds,length(winds)/2,fs);

% full cabin

[mic_TFs_full, f_axis_full] = calculateInpulseTFs(filtered_hamm_full,filtered_mic_full./K,winds,length(winds)/2,fs);
[acc_TFs_full, ~] = calculateInpulseTFs(filtered_hamm_full,filtered_acc_full,winds,length(winds)/2,fs);
[vib_TFs_full, ~] = calculateInpulseTFs(filtered_hamm_full,filtered_vib_fulls,winds,length(winds)/2,fs);

%% compute Transfer functions as Fourier(Out)/Fourier(In)
sz = 4096;

% empty cabin
hamm_FTs = computeFTs(filtered_hamm,filtered_hamm,sz);
mic_FTs = computeFTs(filtered_hamm,filtered_mic,sz);
acc_FTs = computeFTs(filtered_hamm,filtered_acc,sz);
vib_FTs = computeFTs(filtered_hamm,filtered_vib,sz);

fft_axis = (0:1:length(hamm_FTs)-1).*(fs/sz);

mic_TFs_2 = mic_FTs./hamm_FTs; 
acc_TFs_2 = acc_FTs./hamm_FTs; 
vib_TFs_2 = vib_FTs./hamm_FTs; 

% full cabin
hamm_FTs_full = computeFTs(filtered_hamm_full,filtered_hamm_full,sz);
mic_FTs_full = computeFTs(filtered_hamm_full,filtered_mic_full,sz);
acc_FTs_full = computeFTs(filtered_hamm_full,filtered_acc_full,sz);
vib_FTs_full = computeFTs(filtered_hamm_full,filtered_vib_full,sz);

fft_axis_full = (0:1:length(hamm_FTs_full)-1).*(fs/sz);

mic_TFs_full_2 = mic_FTs_full./hamm_FTs_full; 
acc_TFs_full_2 = acc_FTs_full./hamm_FTs_full; 
vib_TFs_full_2 = vib_FTs_full./hamm_FTs_full; 

%% compute the coherence
mic_coherence = zeros(512,13,4,6);
acc_coherence = zeros(512,13,4,6);
for i = 1:height(winds)
    %vib_coherence(:,:,i) = mscohere(squeeze(permute(vib_TFs(:,:,i,:), [4 2 3 1])), squeeze(permute(vib_TFs_2, [3 1 2])));
    vib_coherence(:,:,i) = mscohere(squeeze(permute(vib_TFs(:,:,i,:), [4 2 3 1])), squeeze(permute(vib_TFs_2, [3 1 2])));
    mic_coherence(:,:,i,:) = reshape(mscohere(squeeze(permute(mic_TFs(:,:,i,:), [4 2 3 1])), squeeze(permute(mic_TFs_2, [3 2 1]))),512,13,1,6);
    acc_coherence(:,:,i,:) = reshape(mscohere(squeeze(permute(acc_TFs(:,:,i,:), [4 2 3 1])), squeeze(permute(acc_TFs_2, [3 2 1]))),512,13,1,6);

    %vib_coherence_full(:,:,i) = mscohere(squeeze(permute(vib_TFs_full(:,:,i,:), [4 2 3 1])), squeeze(permute(vib_TFs_full_2, [3 1 2])));
    %mic_coherence_full(:,:,i) = mscohere(squeeze(permute(mic_TFs_full(:,:,i,:), [4 2 3 1])), squeeze(permute(mic_TFs_full_2, [3 2 1])));
    %acc_coherence_full(:,:,i) = mscohere(squeeze(permute(acc_TFs_full(:,:,i,:), [4 2 3 1])), squeeze(permute(acc_TFs_full_2, [3 2 1])));    
end

%%
if plot_TF
    for i = 1:width(mic_TFs)
        figure
        sgtitle("Transfer Function for hit number " + i)
        for j = 1:height(mic_TFs)
            subplot(3,2,j);
            plot(f_axis, abs(squeeze(mic_TFs(j,i,:,:))))
            xlim([0 1500])
            title("Microphone: " + j)
            legend(windnames)
        end
    end
    
    for i = 1:width(acc_TFs)
        figure
        sgtitle("Transfer Function for hit number " + i)
        for j = 1:height(acc_TFs)
            subplot(3,2,j);
            plot(f_axis, abs(squeeze(acc_TFs(j,i,:,:))))
            xlim([0 1500])
            title("Accelerometer: " + j)
            legend(windnames)
        end
    end
    
    for i = 1:width(vib_TFs)
        figure
        title("Vibrocimeter Transfer Function for hit number " + i)
        plot(f_axis, abs(squeeze(vib_TFs(:,i,:,:))))
        xlim([0 1500])
        legend(windnames)
    end
end

%% calculate the average for the TFs over the inpulses (using the TF estimate obtained ones (can be changed))
avg_mic_TFs = squeeze(sum(mic_TFs,2)./width(mic_TFs));
avg_acc_TFs = squeeze(sum(acc_TFs,2)./width(acc_TFs));
avg_vib_TFs = squeeze(sum(vib_TFs,2)./width(vib_TFs));

if plot_avg
    figure
    sgtitle("Average Transfer Function for Microphones")
    for j = 1:height(avg_mic_TFs)
        subplot(3,2,j);
        plot(f_axis, db(abs(squeeze(avg_mic_TFs(j,:,:)))))
        xlim([0 1500])
        title("Microphone: " + j)
        legend(windnames)
    end
    
    figure
    sgtitle("Average Transfer Function for Accelerometers")
    for j = 1:height(avg_acc_TFs)
        subplot(3,2,j);
        plot(f_axis, db(abs(squeeze(avg_acc_TFs(j,:,:)))))
        xlim([0 1500])
        title("Accelerometer: " + j)
        legend(windnames)
    end
    
    figure
    sgtitle("Average Transfer Function for Vibrometer")
    plot(f_axis, db(abs(squeeze(avg_vib_TFs(:,:)))))
    xlim([0 1500])
    legend(windnames)
end

%%

ramp = readtable('students_data/SBN/Test_Ramp.txt');

mic_ramp(1,:) = ramp.Var1;
mic_ramp(2,:) = ramp.Var2;
mic_ramp(3,:) = ramp.Var5;
mic_ramp(4,:) = ramp.Var6;
mic_ramp(5,:) = ramp.Var3;
mic_ramp(6,:) = ramp.Var4;

acc_ramp(1,:) = ramp.Var7;   % Accelerometro 1 (Colonna 7)
acc_ramp(2,:) = ramp.Var8;   % Accelerometro 2 (Colonna 8)
acc_ramp(3,:) = ramp.Var9;   % Accelerometro 3 (Colonna 9)
acc_ramp(4,:) = ramp.Var10;  % Accelerometro 4 (Colonna 10)
acc_ramp(5,:) = ramp.Var11;  % Accelerometro 5 (Colonna 11)
acc_ramp(6,:) = ramp.Var12;  % Accelerometro 6 (Colonna 12)

clear ramp

figure
sgtitle("Accelerometers")
for i = 1:height(acc_ramp)
    subplot(3,2,i)
    spectrogram(acc_ramp(i,:), 'yaxis');
end

figure
sgtitle("Mics")
for i = 1:height(acc_ramp)
    subplot(3,2,i)
    spectrogram(mic_ramp(i,:), 'yaxis');
end