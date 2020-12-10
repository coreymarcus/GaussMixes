function [mu, P, w] = kMeans(k, data)
%kMeans performs a naive 1-D k-means implementation

%num data
N = length(data);

%group each point belongs to
group_idx = zeros(N,1);
group_idx_new = zeros(N,1);

%get random initialization
mu_idx = randi([1 N], k,1);

%means
mu = data(mu_idx);

%distance
distmat = zeros(N,k);

%loop
loopmaxidx = 100;
looplogic = true;
loopidx = 1;
while looplogic
    
    %calculate distance to each mean of each data point
    for ii = 1:k
        distmat(:,ii) = abs(data - mu(ii));
    end
    
    %assign group index based on min distance
    for ii = 1:N
        [~, minidx] = min(distmat(ii,:));
        group_idx_new(ii) = minidx;
    end
    
    %calculate the number of changed points
    N_change = N - sum(group_idx == group_idx_new);
    
    %update groups
    group_idx = group_idx_new;
    
    %update means
    for ii = 1:k
        mu(ii) = mean(data(group_idx == ii));
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
P = zeros(k,1);
for ii = 1:k
    P(ii) = var(data(group_idx == ii));
end

%calculate the weights
w = zeros(k,1);
for ii = 1:k
    w(ii) = sum(group_idx == ii)/N;
end

end

