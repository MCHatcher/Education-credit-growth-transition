%Kitaura model, Economics Letters 2012
%Simulation of the transition dynamics after a relaxation of credit
%Case of binding constraint is denoted 'bind'
%Code written by M. Hatcher (Aug 2021, m.c.hatcher@soton.ac.uk)

clear; clc;
close all;

%Parameter values
rho = 0.3;
A = 1;
thetta = 1;
alfa = 0.3900005;
eta = 0.25;
psi = 0.205;
psi_new = 0.248;
T_reform = 400;
T = T_reform+300;

coef_k = rho*alfa*(1-alfa)*A*(1-eta)/( (1+rho)*(alfa + eta*(1-alfa)) );
coef_k_bind = rho*alfa*(1-alfa)*A*(1-psi)/( (1+rho)*(alfa + psi*(1-alfa)) );

coef_h = thetta*eta^eta*((1-alfa)/alfa)^eta;
coef_h_bind = thetta*psi^eta*((1-alfa)/alfa)^eta;

if psi >= eta
    coef_k_bind = coef_k;
    coef_h_bind = coef_h;
    psi = eta;
end

%Initial values
k0 = 500;
h0 = 500;
y0 = A*k0^alfa*h0^(1-alfa);

k(1) = coef_k*k0^alfa*h0^(1-alfa);
k_bind(1) = coef_k_bind*k0^alfa*h0^(1-alfa);

h(1) = coef_h*coef_k^eta*k0^(alfa*eta)*h0^(1-alfa*eta);
h_bind(1) = coef_h_bind*coef_k_bind^eta*k0^(alfa*eta)*h0^(1-alfa*eta);

x(1) = k(1)/h(1);
x_bind(1) = k_bind(1)/h_bind(1);

y(1) = A*k(1)^alfa*h(1)^(1-alfa);
y_bind(1) = A*k_bind(1)^alfa*h_bind(1)^(1-alfa);

growth(1) = y(1)/y0;
growth_bind(1) = y_bind(1)/y0;

for t=2:T

if t ~= T_reform 
    
   if t>T_reform
        psi = psi_new;
        coef_k = rho*alfa*(1-alfa)*A*(1-eta)/( (1+rho)*(alfa + eta*(1-alfa)) );
        coef_k_bind = rho*alfa*(1-alfa)*A*(1-psi)/( (1+rho)*(alfa + psi*(1-alfa)) );
        coef_h = thetta*eta^eta*((1-alfa)/alfa)^eta;
        coef_h_bind = thetta*psi^eta*((1-alfa)/alfa)^eta;
        
        if psi >= eta
            coef_k_bind = coef_k;
            coef_h_bind = coef_h;
            psi = eta;
        end   
    end
    
k(t) = coef_k*k(t-1)^alfa*h(t-1)^(1-alfa);
k_bind(t) = coef_k_bind*k_bind(t-1)^alfa*h_bind(t-1)^(1-alfa);
e(t-1) = eta*(1-alfa)/alfa*k(t);
e_bind(t-1) = psi*(1-alfa)/alfa*k_bind(t);
%e(t-1) = eta*(1-alfa)*y(t)/R(t);
%e_bind(t-1) = psi*(1-alfa)*y_bind(t)/R_bind(t);

h(t) = thetta*e(t-1)^eta*h(t-1)^(1-eta);
h_bind(t) = thetta*e_bind(t-1)^eta*h_bind(t-1)^(1-eta);
%h(t) = coef_h*coef_k^eta*k(t-1)^(alfa*eta)*h(t-1)^(1-alfa*eta);
%h_bind(t) = coef_h_bind*coef_k_bind^eta*k_bind(t-1)^(alfa*eta)*h_bind(t-1)^(1-alfa*eta);

x(t) = k(t)/h(t);
x_bind(t) = k_bind(t)/h_bind(t);

y(t) = A*k(t)^alfa*h(t)^(1-alfa);
y_bind(t) = A*k_bind(t)^alfa*h_bind(t)^(1-alfa);
s(t) = (rho/(1+rho))*(1-alfa)*(1-eta)*y(t);
s_bind(t) = (rho/(1+rho))*(1-alfa)*(1-psi)*y_bind(t);

w(t) = (1-alfa)*A*(k(t)/h(t))^alfa;
R(t) = alfa*A*(k(t)/h(t))^(alfa-1);
w_bind(t) = (1-alfa)*A*(k_bind(t)/h_bind(t))^alfa;
R_bind(t) = alfa*A*(k_bind(t)/h_bind(t))^(alfa-1);

k_next(t) = alfa/(alfa + eta*(1-alfa))*s(t);
k_next_bind(t) = alfa/(alfa + psi*(1-alfa))*s_bind(t);
e_next(t) = eta*(1-alfa)/alfa*k_next(t);
e_next_bind(t) = psi*(1-alfa)/alfa*k_next_bind(t);
h_next(t) = thetta*e_next(t)^eta*h(t)^(1-eta);
h_next_bind(t) = thetta*e_next_bind(t)^eta*h_bind(t)^(1-eta);

R_next(t) = alfa*A*(k_next(t)/h_next(t))^(alfa-1);
R_next_bind(t) = alfa*A*(k_next_bind(t)/h_next_bind(t))^(alfa-1);

c(t) = s(t)/rho;
c_bind(t) = s_bind(t)/rho;
d(t) = R(t)*s(t-1);
d_bind(t) = R_bind(t)*s_bind(t-1);

d_next(t) = R_next(t)*s(t); 
d_next_bind(t) = R_next_bind(t)*s_bind(t); 

U(t) = log(c(t)) + rho*log(d_next(t));
U_bind(t) = log(c_bind(t)) + rho*log(d_next_bind(t));
U_old(t) = log(c(t-1)) + rho*log(d(t)); 
U_old_bind(t) = log(c_bind(t-1)) + rho*log(d_bind(t));

growth(t) = y(t)/y(t-1);
growth_bind(t) = y_bind(t)/y_bind(t-1);

end

if t==T_reform
    k(t) = coef_k*k(t-1)^alfa*h(t-1)^(1-alfa);
    k_bind(t) = coef_k_bind*k_bind(t-1)^alfa*h_bind(t-1)^(1-alfa);
    e(t-1) = eta*(1-alfa)/alfa*k(t);
    e_bind(t-1) = psi*(1-alfa)/alfa*k_bind(t);
    %e(t-1) = eta*(1-alfa)*y(t)/R(t);
    %e_bind(t-1) = psi*(1-alfa)*y_bind(t)/R_bind(t);

    h(t) = thetta*e(t-1)^eta*h(t-1)^(1-eta);
    h_bind(t) = thetta*e_bind(t-1)^eta*h_bind(t-1)^(1-eta);
    %h(t) = coef_h*coef_k^eta*k(t-1)^(alfa*eta)*h(t-1)^(1-alfa*eta);
    %h_bind(t) = coef_h_bind*coef_k_bind^eta*k_bind(t-1)^(alfa*eta)*h_bind(t-1)^(1-alfa*eta);

    x(t) = k(t)/h(t);
    x_bind(t) = k_bind(t)/h_bind(t);

    y(t) = A*k(t)^alfa*h(t)^(1-alfa);
    y_bind(t) = A*k_bind(t)^alfa*h_bind(t)^(1-alfa);
    s(t) = (rho/(1+rho))*(1-alfa)*(1-eta)*y(t);
    s_bind(t) = (rho/(1+rho))*(1-alfa)*(1-psi)*y_bind(t);

    w(t) = (1-alfa)*A*(k(t)/h(t))^alfa;
    R(t) = alfa*A*(k(t)/h(t))^(alfa-1);
    w_bind(t) = (1-alfa)*A*(k_bind(t)/h_bind(t))^alfa;
    R_bind(t) = alfa*A*(k_bind(t)/h_bind(t))^(alfa-1);
    
    k_next(t) = alfa/(alfa + eta*(1-alfa))*s(t);
    k_next_bind(t) = alfa/(alfa + psi_new*(1-alfa))*s_bind(t);
    e_next(t) = eta*(1-alfa)/alfa*k_next(t);
    e_next_bind(t) = psi_new*(1-alfa)/alfa*k_next_bind(t);
    h_next(t) = thetta*e_next(t)^eta*h(t)^(1-eta);
    h_next_bind(t) = thetta*e_next_bind(t)^eta*h_bind(t)^(1-eta);

    R_next(t) = alfa*A*(k_next(t)/h_next(t))^(alfa-1);
    R_next_bind(t) = alfa*A*(k_next_bind(t)/h_next_bind(t))^(alfa-1);

    c(t) = s(t)/rho;
    c_bind(t) = s_bind(t)/rho;
    d(t) = R(t)*s(t-1);
    d_bind(t) = R_bind(t)*s_bind(t-1);

    d_next(t) = R_next(t)*s(t); 
    d_next_bind(t) = R_next_bind(t)*s_bind(t); 

    U(t) = log(c(t)) + rho*log(d_next(t));
    U_bind(t) = log(c_bind(t)) + rho*log(d_next_bind(t));
    U_old(t) = log(c(t-1)) + rho*log(d(t)); 
    U_old_bind(t) = log(c_bind(t-1)) + rho*log(d_bind(t)); 
    
    growth(t) = y(t)/y(t-1);
    growth_bind(t) = y_bind(t)/y_bind(t-1);
    
end
    
    
end

for t=T_reform-4:T
    %EV calcs 
    U_orig(t) = U(T_reform-1) + (1+rho)*(t+1-T_reform)*log(growth(T_reform-1));
    U_orig_bind(t) = U_bind(T_reform-1) + (1+rho)*(t+1-T_reform)*log(growth_bind(T_reform-1));
    U_orig_old(t) = U_old(T_reform-1) + (1+rho)*(t+1-T_reform)*log(growth(T_reform-1));
    U_orig_old_bind(t) = U_old_bind(T_reform-1) + (1+rho)*(t+1-T_reform)*log(growth_bind(T_reform-1));
    
    Lambda(t) = 100*( exp((U(t)-U_orig(t))/(1+rho)) -1);
    Lambda_bind(t) = 100*( exp((U_bind(t)-U_orig_bind(t))/(1+rho)) -1);
    Lambda_old(t) = 100*( exp((U_old(t)-U_orig_old(t))/(1+rho)) -1);
    Lambda_old_bind(t) = 100*( exp((U_old_bind(t)-U_orig_old_bind(t))/(1+rho)) -1);
    Generations(t) = t+1-T_reform; 
    
    g_bind(t) = growth_bind(t)/growth_bind(T_reform-1);
    ratio(t) = x_bind(t)/x_bind(T_reform-1);
    ratio1(t) = h_bind(t)/(h_bind(T_reform-1)*(growth_bind(T_reform-1))^(t+1-T_reform));
    %For plots
    zero(t) = 0;
    one(t) = 1;
    
end

EV = [Lambda_old_bind(T_reform), Lambda_bind(T_reform:T_reform+24)];
Gen = Generations(T_reform-1:T_reform+24); 
Zero = zero(T_reform-1:T_reform+24);

hold on, set(0,'DefaultLineLineWidth',0.75)
figure(1)
plot(Gen(1:11), one(T_reform-1:T_reform+9),'--k'), hold on, plot(Gen(1:11), g_bind(T_reform-1:T_reform+9),'r'), 
hold on, plot(Gen(1:11), ratio(T_reform-1:T_reform+9),'b'), hold on, plot(Gen(1:11), ratio1(T_reform-1:T_reform+9),'b'), 
title('Growth'),  xlabel('Time')
figure(2)
hold on, plot(Gen,Zero,'--k'), hold on, plot(Gen,EV,'b'), title('Welfare'), xlabel('Generations')
figure(3)
subplot(1,2,1), hold on, plot(Gen,Zero,'--k'), hold on, plot(Gen,EV,'k'), title('Equivalent variations'), xlabel('Generations'), ylabel('%')
subplot(1,2,2), plot(Gen(1:11), one(T_reform-1:T_reform+9),'--k'), hold on, plot(Gen(1:11), ratio1(T_reform-1:T_reform+9),'k'),
hold on, plot(Gen(1:11), g_bind(T_reform-1:T_reform+9),'k'), hold on, plot(Gen(1:11), ratio(T_reform-1:T_reform+9),'k'), 
title('Growth and capital'), xlabel('Time')




