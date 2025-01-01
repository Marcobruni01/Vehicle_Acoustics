clc 
clear 
close all
%%
ext_base = "students_data/exterior_noise";
int_base = "students_data/interior_noise";
ped_d = 65.5; %[m] distance between the 2 pedals
Fs = 16384;
train_len = 174;
T_spl = 0.125;
num_samps = ceil(Fs*T_spl);

extdirs = strsplit(convertCharsToStrings(ls(ext_base)));
extdirs(end) = [];

p0 = 1;

%%

for i = 1:length(extdirs)
    str_app = ext_base + "/" + extdirs(i);
    extcont = strsplit(convertCharsToStrings(ls(str_app)));
    extcont(end) = [];
    for j = 1:length(extcont)
        str_app = ext_base + "/" + extdirs(i) + "/" + extcont(j);
        var = load(str_app);
        Mic_ext(i,j).name = var.variabile.nome;
        Mic_ext(i,j).data = var.variabile.Dati;
        Mic_ext(i,j).fsamp = var.variabile.fsamp;
        Mic_ext(i,j).xlab = var.variabile.xlab;
        Mic_ext(i,j).ylab = var.variabile.ylab;
        weight_filt = weightingFilter('A-weighting', 'SampleRate', Mic_ext(i,j).fsamp);
        Mic_ext(i,j).data = weight_filt(Mic_ext(i,j).data);
        if(j < 5)
            for k = 1:length(Mic_ext(i,j).data)/num_samps
               if(k == length(Mic_ext(i,j).data)/num_samps)
                   cut = k*num_samps - length(Mic_ext(i,j).data);
               else
                   cut = 0;
                end
                Mic_ext(i,j).SPL(k) = (sum(Mic_ext(i,j).data((k-1)*num_samps + 1:(k*num_samps - cut)).^2)./num_samps).^(1/2);
            end
            % Mic_ext(i,j).SPL = Mic_ext(i,j).SPL./max(Mic_ext(i,j).SPL);
            app = db(Mic_ext(i,j).SPL);
            [val, idx] = max(app);
            Mic_ext(i,j).t1_idx = idx;
            Mic_ext(i,j).t2_idx = idx;
            while(app(Mic_ext(i,j).t1_idx) > val - 20)
                Mic_ext(i,j).t1_idx = Mic_ext(i,j).t1_idx - 1;
            end
            while(app(Mic_ext(i,j).t2_idx) > val - 20)
                Mic_ext(i,j).t2_idx = Mic_ext(i,j).t2_idx + 1;
            end
            Mic_ext(i,j).t1 = Mic_ext(i,j).t1_idx/Fs;
            Mic_ext(i,j).t2 = Mic_ext(i,j).t2_idx/Fs;
        end
     end
end

clear var idx val t1_idx t2_idx i j k

%%

for i = 1:length(Mic_ext)
    [~, idx1] = findpeaks(Mic_ext(i,5).data,'MinPeakHeight', 0.5);
    [~, idx2] = findpeaks(Mic_ext(i,6).data,'MinPeakHeight', 0.5);
    for j = 1:width(Mic_ext)
        Mic_ext(i,j).Tp(1) = (idx1(end) - idx1(1))/Fs;
        Mic_ext(i,j).Tp(2) = (idx2(end) - idx2(1))/Fs;
        Mic_ext(i,j).avg_vel(1) = train_len/Mic_ext(i,j).Tp(1);
        Mic_ext(i,j).avg_vel(2) = train_len/Mic_ext(i,j).Tp(2);
        Mic_ext(i,j).avg_vel(3) = ped_d*Fs/(idx2(1) - idx1(1));
        Mic_ext(i,j).Tp(3) = train_len/Mic_ext(i,j).avg_vel(3);
    end
end

clear idx1 idx2 i

%%

for i = 1:length(Mic_ext)
    for j = 1:(width(Mic_ext)-2)
        % Mic_ext(i,j).TEL = 10*log10(sum((Mic_ext(i,j).data(Mic_ext(i,j).t1_idx:Mic_ext(i,j).t2_idx).^2)./p0^2)/Mic_ext(i,j).Tp);
        SPL = [];
        for k = 1:length(Mic_ext(i,j).SPL)
            SPL = [SPL, ones(1,num_samps).*Mic_ext(i,j).SPL(k)];
        end
        Mic_ext(i,j).TEL = 10*log10(sum((Mic_ext(i,j).data(Mic_ext(i,j).t1_idx:Mic_ext(i,j).t2_idx).^2)./SPL(Mic_ext(i,j).t1_idx:Mic_ext(i,j).t2_idx).^2)./Mic_ext(i,j).Tp(1));
    end
end

clear SPL
%%
for i = 1:length(Mic_ext)
    app1(:,i) = Mic_ext(i,:).Tp;
    app2(:,i) = Mic_ext(i,:).avg_vel;
end

figure()
subplot(2,1,1)
bar(1:length(Mic_ext),app1)
grid on
legend("Tp from pad 1", "Tp from pad 1", "Tp from avg speed")
subplot(2,1,2)
bar(1:length(Mic_ext),app2)
legend("avg vel from pad 1", "avg vel from pad 1", "avg vel from distance of 2 pads")
grid on

clear app1 app2 i
%%

figure
for i = 1:4
    subplot(2,2,i)
    hold on
    grid on
    for j = 1:height(Mic_ext)
        scatter(Mic_ext(j,i).avg_vel(1), Mic_ext(j,i).TEL)

    end
end

%%
[Spec_near, f_near, t_near] = spectrogram(Mic_ext(1,1).data, 512,256,8192,Fs);
[Spec_far, f_far, t_far] = spectrogram(Mic_ext(1,3).data, 512,256,8192,Fs);
Spec_far = Spec_far./max(abs(Spec_far),[],1);
Spec_near = Spec_near./max(abs(Spec_near),[],1);
figure
sgtitle("Spectrogram comparison")
subplot(1,2,1)
imagesc(f_near, t_near, db(abs(Spec_near).')); % Time and frequency axes
axis xy
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title("Microphone near");
set(gca, 'YDir', 'normal');
%ylim([7 14])
colorbar;
sgtitle("Spectrogram comparison")
subplot(1,2,2)
imagesc(f_far, t_far, db(abs(Spec_far).')); % Time and frequency axes
axis xy
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title("Microphone far");
set(gca, 'YDir', 'normal');
%ylim([8 13])
colorbar;