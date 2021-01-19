function [mu, P, w, group_idx] = kMeans(k, data)
%kMeans performs a naive N-D k-means implementation

%num data
Ndim = size(data,1);
Ndata = size(data,2);

%group each point belongs to
group_idx = zeros(Ndata,1);
group_idx_new = zeros(Ndata,1);

%get random initialization
mu_idx = randi([1 Ndata], k,1);

%means
mu = data(:, mu_idx);

%distance
distmat = zeros(Ndata,k);

%loop
loopmaxidx = 100;
looplogic = true;
loopidx = 1;
while looplogic
    
    %calculate distance to each mean of each data point
    for ii = 1:k
        for jj = 1:Ndata
            distmat(jj,ii) = norm(data(:,jj) - mu(:,ii));
        end
    end
    
    %assign group index based on min distance
    for ii = 1:Ndata
        [~, minidx] = min(distmat(ii,:));
        group_idx_new(ii) = minidx;
    end
    
    %calculate the number of changed points
    N_change = Ndata - sum(group_idx == group_idx_new);
    
    %update groups
    group_idx = group_idx_new;
    
    %update means
    for ii = 1:k
        mu(:,ii) = mean(data(:,group_idx == ii),2);
    end
    
    %loop managment
    loopidx = loopidx + 1;
    if(loopidx >= loopmaxidx)
        disp("Max k-Mean idx!")
        looplogic = false;
    end
    
    if(N_change == 0)
        looplogic = false;
    end
end

%calculate the variance
P = zeros(Ndim,Ndim,k);
for ii = 1:k
    P(:,:,ii) = cov(data(:,group_idx == ii)');
end

%calculate the weights
w = zeros(k,1);
for ii = 1:k
    w(ii) = sum(group_idx == ii)/Ndata;
end

end

