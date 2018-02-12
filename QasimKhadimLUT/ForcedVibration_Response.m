clear all
close all

syms m k c real

assume      (m>0) % mass             'm'
assumeAlso  (k>0) % spring constant  'k'
assumeAlso  (c>0) % damping constant 'c'

%% Solving the system equation
%%
syms x(t) F F0 A B w f K X phi

x= A*cos(w*t)+B*sin(w*t);                   %stady state solution of forced vibration equation
xp = diff (x,t);                            %velocity
xpp= diff (xp,t);                           %acceleration
Eq=m*xpp+c*xp+k*x==F0*sin(w*t);             %equation of motion

%% comparing both sides of equation for sin(w*t)and cos(w*t)
%%
Eq_1= collect (Eq, [cos(w*t) sin(w*t)]); 

Eq_2= A*(k- m*w^2)+ B*c*w==0;               %coefficients of cos(w*t)
Eq_3= - (A*c*w) + B*(k - m*w^2)==F0/k;      %coefficients of sin(w*t)

%% Subsituting the variables k,m,c and w with r and zeta
%%
syms r zeta real

Eq_4= subs(Eq_2,{k-m*w^2,c*w},{1-r^2,2*zeta*r});
Eq_5= subs(Eq_3,{k-m*w^2,c*w},{1-r^2,2*zeta*r});


%% Solving the Eq_4 and Eq_5 for A and B
%%
a= [1-r^2 2*zeta*r;2*zeta*r 1-r^2];         %matrix a
b=[0;F0/k];                                 %matrix b
sol=a^-1*b;                                 %X=A^-1
A=sol(1);
B=sol(2);

%%  Particular solution is 
%%
xp= A*cos(w*t)+B*sin(w*t);                  %particular or steady state solution, same as above

%% By comparing with %xp=X*sin(w*t+phi)=X*(sin(w*t)*cos(phi)+cos(w*t)sin(phi)
%%

X*cos(phi)==A;
X*sin(phi)==B;

%%  Calculating the amplification factor
%% by using X= (sqrt(A^2+B^2)), amplitude of vibration is obtained. By dividing and simplifying the maximum transmitted force with the amplitude of force, amplification factor can be written as:


X= (sqrt(A^2+B^2));         
M =1./sqrt((1-r.^2).^2+(2*zeta*r).^2);  %amplification factor


%% plotting the results
%%
r=[0 2.5];

M=fplot(subs(M,zeta,[0,0.1,0.2,0.3,0.4,0.5,1]),r,'k','LineWidth',1.5)          %plotting the function


ylabel({'Amplification factor ($X\frac{k}{F_{0}}$)'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Times New Roman')
ylim([0 5])
xlabel({'Frequency ratio($r=\frac{f}{f{n}}$)'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Times New Roman')
title('Forced vibrations for an amplitude','FontUnits','points','FontWeight','normal','FontSize',12,'FontName','Times New Roman')
annotation('textarrow',r,M,'String',['\varsigma =0',sprintf('\n'),'No damping','\varsigma =0.1','\varsigma =0.2','\varsigma =0.3','\varsigma =0.4','\varsigma =0.5','\varsigma =0',sprintf('\n'),'Critical damping'],'FontWeight','normal','FontSize',12,'FontName','Times New Roman')
hold off

print('Z:\Qasim_Cmim\Assignment#03\pic01','-dmeta')






