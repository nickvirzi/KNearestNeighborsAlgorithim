clear all; clc;

t = readtable('C:\Users\virzina\Documents\data.txt', 'ReadVariableNames', false);
t.Properties.VariableNames = {'P', 'N'};

A = t{:,:};

Z1(1:199, 1) = 0;
Z2(1:200, 1) = 1;

B = [Z1; Z2];

Anew = [3.4699 6.6262];
A(2,:) = [];
B(2,:) = [];

correct0 = 0;
correct1 = 0;
incorrect0 = 0;
incorrect1 = 0;
MYaccuracy = [];


for k = 1:398
    [NearestToA,NearestToB,DistanceFromPoint] = FindKNeasrestNeighbor(A,B,Anew,k)
end

RDistances = [];
ARPoints = [];
BRPoints = [];

for R = 1:10
    for I = 1:398
        for J = 1:2
            if DistanceFromPoint(I) <= R
                RDistances(end + 1) = DistanceFromPoint(I);
                ARPoints(I,J) = NearestToA(I,J);
                BRPoints(I) = NearestToB(I);
            end
        end
    end
    for T = 1:length(BRPoints)
        if (BRPoints(T) == 0)
            correct0 = correct0 + 1;
        else
            incorrect1 = incorrect1 + 1;
        end
    end
    
    MYaccuracy(R) = correct0 / (correct0 + incorrect1);
end

kplot = [1:R];

plot(kplot,MYaccuracy);
title('Balanced Accuracy vs r');
xlabel('r');
ylabel('Balanced Accuracy');

function [NearestToA,NearestToB,DistanceFromPoint] = FindKNeasrestNeighbor(A,B,Anew,k)
[C,~,B] = unique(B,'stable');
%Get Distance%
A = repmat(Anew,size(A,1),1)-A;
DistanceFromPoint = sqrt(sum(A.^2,2));
[~,I] = sort(DistanceFromPoint);
NearestToA = A(I(1:k),:);
% Check the number of output arguments
if nargout > 1
    NearestToB = C(B(I(1:k)));
end
if nargout > 2
    DistanceFromPoint = DistanceFromPoint(I(1:k));
end
end