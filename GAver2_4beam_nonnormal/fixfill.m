function threshold = fixfill(IN, nbins, fill)

%finding threshold that enforces the proper fill ratio
%see Notebook 1, p. 34 for details

N = numel(IN);
nsample = fill*N;
maxIN = max(real(reshape(IN,N,1)));
minIN = min(real(reshape(IN,N,1)));
range = maxIN - minIN;
%binsize = INmax / (nbins-1)

%h = histogram(real_part(IN),nbins=nbins,binsize=maxIN/nbins,min=0);

% IN must be a vector, NOT a 3D array

%h = fliplr(hist(real(reshape(IN,N,1)),minIN:range/nbins:maxIN));
%cum = cumsum(h)/N;
%cum = cumsum(fliplr(hist(real(IN),minIN:range/nbins:maxIN)))/N;
%DF = .25; %original held constant
%df = DF;
%last = 0;
%go = 1;
%l = 1;
%while go == 1
%    i = find( (cum > (fill-df)) .* (cum < (fill+df)) );
%    if numel(i) > 1
%        l = l + (last == 1);
%        df = df - DF/2^l;
%        last = -1;
%    elseif numel(i) < 1
%        l = l + (last == -1);
%        df = df + DF/2^l;
%        last = 1;
%    else
%        go = 0;
%        %disp('stop');
%    end
%end

h = hist(real(reshape(IN,N,1)),minIN:range/nbins:maxIN);
if numel(h) ~= nbins+1
    disp('nbins ~= 1');
end
if (minIN == maxIN) 
        disp('maxIN == minIN');
end
sum_hi = 0;
i = nbins+1; %searching the threshold from highest to lowest
go = 1;
while go == 1
    if i > numel(h) || i < 1
        a = 1;
    end
    sum_hi = sum_hi + h(i);
    if (sum_hi >= nsample)
       go = 0;
       sum_lo = sum_hi - h(i);
       if (nsample - sum_lo) < (sum_hi - nsample)
           i = i + 1;
       end
    else
           i = i - 1;
    end
end

%threshwq = float(i)/nbins*INmax %**GLOBAL CHANGE threshold as a real number
threshold = i/nbins*range+minIN; %in dimensionless units of Intensity

%threshold = maxIN - i/nbins*range;