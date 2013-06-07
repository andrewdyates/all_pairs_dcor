clear;
close all;

M = 10000; n = 500;
data = rand(M, n);

% tmpData = 1:500;
% data = [tmpData; tmpData*2; sin(tmpData/180*pi); cos(tmpData/180*pi); rand(1, 500); rand(1, 500)];
% data = [rand(1, 500); rand(1, 500)];

[M, n] = size(data);

meanD = zeros(M, n);

for i = 1 : n
    tmp = data - repmat(data(:, i), 1, n);
    meanD(:, i) = mean(abs(tmp), 2);
end;

meanTotal = mean(meanD, 2);

dCov = zeros(M, M);
tic
for i = 1 : n
    if mod(i, 10) == 0
        i
        toc
    end;
    A = abs(data - repmat(data(:, i), 1, n)) - repmat(meanD(:, i), 1, n)...
        - meanD + repmat(meanTotal, 1, n);
    dCov = dCov + A*A'/n/n;
end;

dCov = sqrt(dCov);

dVar = diag(dCov);

dCorr = dCov;
for i = 1 : M
    dCorr(i, :) = dCorr(i, :) ./ sqrt(dVar');
end;

for i = 1 : M
    dCorr(:, i) = dCorr(:, i) ./ sqrt(dVar);
end;
