clc
clear
close all
%%
fs = 12500; % Frequenza di campionamento
order = 4; % Ordine del filtro
cutoff_freq = 10; % Frequenza di taglio in Hz
sz = 4096;

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

test_FTs(1,:,:) = test1_FTs.';
test_FTs(2,:,:) = test2_FTs.';
test_FTs(3,:,:) = test3_FTs.';
test_FTs(4,:,:) = test4_FTs.';

clear test1_FTs test2_FTs test3_FTs test4_FTs

Q_FTs = Q_FTs.';

test_TFs = zeros(size(test_FTs));

for i = 1:height(test_FTs)
    test_TFs(i,:,:) = squeeze(test_FTs(i,:,:))./Q_FTs(i,:);
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
end
