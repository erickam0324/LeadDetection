%% orthogonal 

A = [-0.6779, -0.7352]; 
B = [0.7352, -0.6779];

x = [A(1);B(1)];
y = [A(2);B(2)];

figure, plot(x,y,'color','k','LineWidth',2)
hold on
normal = [A(1),A(2)] - null(A-B)';
plot([A(1),normal(1)],[A(2),normal(2)],'color','r','LineWidth',2)