% Matlab source code for Kmeans 

close all 
clear all 
clc 

Limit = 20; 

X = [10*randn(400,2) ; 10 * randn(400,2)] ; 
plot(X(:,1), X(:,2), 'k.')

k =1; 
for i = 1:length(X(:,1))
    if (sqrt(X(i,1)^2 + X(i,2)^2 )) > Limit 
        X(i,1) = 0 ; 
        X(i,2) = 0 ;
    else 
        Y(k,1) =  X(i,1) ;
        Y(k,2) =  X(i,2) ; 
        k = k+1; 
    end 
end 

figure
plot(Y(:,1), Y(:,2), 'k.')