clc
clear all

mu1= 2; mu2=2; mu3=2; mu4=2; mu5=2; mu6=2; 
%syms mu1 mu2 mu3 mu4 mu5 mu6
A1 = [-mu1];
alpha1 = [1];

A2 = [-mu2, 0.5*mu2, 0;
    0, - mu4, 0.5*mu4;
    0, mu5, -mu5];
alpha2= [1, 0, 0];


A3= [- mu2, mu2, 0;
    0, -mu4, 0.2*mu4;
    0, 0.7*mu5, -mu5];
alpha3= [1, 0, 0];


A4 = [-mu5, 0.7*mu5;
    0.2*mu4, -mu4];
alpha4 = [1, 0];

A5=[-mu3, 0.7*mu3, 0, 0;
    0, -mu2, mu2, 0;
    0, 0, -mu4, 0.2*mu4;
    0, 0, 0.7*mu5, -mu5] ;
alpha5= [1,0,0,0];

%calculate PH Matrix of maximum of times for D1->S1 and D1->S2
[B12, gamma12] = MaxPH(A1, A2, alpha1, alpha2);
[B34, gamma34]= MaxPH(A3, A4, alpha3, alpha4);
[B345, gamma345] = MaxPH(B34, A5, gamma34, alpha5);


tau12 = -gamma12*inv(B12)*ones(length(B12), 1);


tau345 = -gamma345*inv(B345)*ones(length(B345), 1);

tau = 0.5*tau12+ 0.5*tau345 + 1/mu6


