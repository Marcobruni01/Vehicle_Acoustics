function FTs = computeFTs(hammer,in, size)

    [~,idx] = findpeaks(hammer, 'MinPeakHeight', 0.02);
    FTs = zeros(height(in),length(idx),size/2+1);
    for j = 1:length(idx)
        if j == length(idx)
            start = idx(j);
            finish = length(in);
        else
            start = idx(j);
            finish = idx(j+1);
        end
        imp_resp = in(:,start:finish);
        app = fft(imp_resp,size,2);
        FTs(:,j,:) = app(:,size/2:size);
    end

end