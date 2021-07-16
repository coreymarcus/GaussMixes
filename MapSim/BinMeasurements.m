function fitidx = BinMeasurements(meas,gauss_list,R, bindist)
%BinMeasurements determines which element each measurement belongs to

%locals
Ndraw = size(meas,2);
Ngauss = length(gauss_list);
x_meas = meas(1,:);
y_meas = meas(2,:);

%initialize output
sigdist = zeros(Ndraw,Ngauss);
fitidx = zeros(Ndraw,1);
for jj = 1:Ndraw
    for kk = 1:Ngauss
        sigdist(jj,kk) = gauss_list{kk}.GaussSig(x_meas(jj),y_meas(jj),R);
    end
    
    %determine which one is the best fit
    [minval, minidx] = min(sigdist(jj,:));
    if(minval < bindist)
        fitidx(jj) = minidx(1);
    end
end

end

