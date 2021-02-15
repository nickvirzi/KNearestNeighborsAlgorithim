clear all; clc;

t = readtable('C:\Users\virzina\Documents\data.txt', 'ReadVariableNames', false);
t.Properties.VariableNames = {'P', 'N'};

ATrain = t{:,:};

Z1(1:199, 1) = 0;
Z2(1:200, 1) = 1;

BTrain = [Z1; Z2];

ATest = [];
BTest = zeros(79,1);

for I = 1:79
    for J = 1:2
        if I > 0 && I < 82
            ATest(I,J) = ATrain(I + 160,J);
        end
    end
    BTest(I) = BTrain(I + 160);
end

for K = 161:239
    ATrain(K,:) = [];
    BTrain(K,:) = [];
end

m = 399;

[Weights,Bias,Iter] = PerecptronTrn(ATrain,BTrain);
Iterations=Iter

Errors=PerecptronTst(ATest,BTest,Weights,Bias);
disp(['Test_Errors=' num2str(Errors) '     Test Data Size= ' num2str(m-320)])

l=BTrain==0;
hold on
plot(ATrain(l,1),ATrain(l,2),'k.' );
plot(ATrain(~l,1),ATrain(~l,2),'g.');
[l,p]=size(ATrain);
plot([0,10],[15,0],'r-')
axis([0 15 0 20]), axis square, grid on
drawnow
title('Perceptron Decision Boundary');
function Errors=PerecptronTst(A,B,Weights,Bias)
%==========================================
% Testing phase
%==========================================
tic
[l,p]=size(A);
Errors=0;
for i=1:l
    AA=A(i,:);
    EB=AA*Weights+Bias;
    if EB>=0.5
        EB=1;
    else
        EB=0;
    end
    if B(i)~=EB
        Errors=Errors+1;
    end
end
toc
end

function [Weights,Bias,Iter]=PerecptronTrn(A,B)
tic
[l,p]=size(A);
Weights= [0; -1]; 
Bias=0;          
InitError=1;        
Iter=0;      
LearnRate=.001;        
MaxNorm=max(sqrt(sum(A)));
while InitError==1 
    InitError=0;
    TrainError=0; 
    for i=1:l  
        AA=A(i,:);
        EstY=AA*Weights+Bias;
        if EstY>=0.5
            EstY=1;
        else
            EstY=0;
        end
        if B(i)~=EstY
            Error=B(i)-EstY;
            Weights=Weights'+(Error*LearnRate)*A(i,:);
            TrainError=TrainError+1 ; 
            Weights=Weights';
        end
    end
    ee=TrainError;   
    if ee>0  
        InitError=1;
    end
    Iter=Iter+1; % stop after 10000 iterations
    if Iter==10000
        InitError=0;
        Iter=0;
    end
end
disp(['Training_Errors=' num2str(TrainError) '     Training data Size=' num2str(l)])
toc
end