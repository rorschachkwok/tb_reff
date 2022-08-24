
close all; clear all; clc;


%% step1： Define all symbolic variables and parameters involved in the differential equation of SEIpInR model ：
syms S_1 E_1 I_p1 I_n1 R_1;
syms S_2 E_2 I_p2 I_n2 R_2;
syms beta_11 beta_12 beta_22 beta_21 kappa m q omega_1 omega_2 lambda theta br dr;
syms mu_1 mu_2 alpha tau gamma_1 gamma_2 f_1 f_2;


%% step2：Define infected variables：
variables = [E_1 I_p1 I_n1 E_2 I_p2 I_n2];


%% step3：
dim = numel(variables); % dimension, number of infected variables
F = sym(zeros(dim,1)); % initialize the vector of newly infections
V = sym(zeros(dim,1)); % initialize the vector of transitions 


%% step4：Define vector F and V 
F(1) = beta_11*S_1*(I_p1+kappa*I_n1) + beta_21*S_1*(I_p2+kappa*I_n2);
F(4) = beta_22*S_2*(I_p2+kappa*I_n2) + beta_12*S_2*(I_p1+kappa*I_n1);
V(1) = m*theta*E_1 + dr*E_1 + q*omega_1*E_1 - lambda*mu_1*I_p1 + (1-q)*omega_2*E_1 - lambda*mu_2*I_n1 - alpha*tau*R_1;
V(2) = - q*omega_1*E_1 + lambda*mu_1*I_p1 + (1-lambda)*gamma_1*I_p1 + (dr+f_1)*I_p1;
V(3) = - (1-q)*omega_2*E_1 + lambda*mu_2*I_n1 + (1-lambda)*gamma_2*I_n1 + (dr+f_2)*I_n1;
V(4) = m*theta*E_2 + dr*E_2 + q*omega_1*E_2 - lambda*mu_1*I_p2 + (1-q)*omega_2*E_2 - lambda*mu_2*I_n2 - alpha*tau*R_2;
V(5) = - q*omega_1*E_2 + lambda*mu_2*I_p2 + (1-lambda)*gamma_2*I_p2 + (dr+f_2)*I_p2;
V(6) = - (1-q)*omega_2*E_2 + lambda*mu_2*I_n2 + (1-lambda)*gamma_2*I_n2 + (dr+f_2)*I_n2



%% step5：construct the next generation matrix： 
% compute the jacobi matrix
JF = sym(zeros(dim));
JV = sym(zeros(dim));
for i = 1:dim
    for j = 1:dim
        JF(i,j) = diff(F(i),variables(j));
        JV(i,j) = diff(V(i),variables(j));
    end
end


%% step6：construct the next-generation matrix Mat = F*V^(-1)：
invJV = inv(JV)
Mat = JF*invJV




%% step7：generate the equations of four interactive Reff between and within groups
idx1 = 1:3;
idx2 = 4:6;
M11 = Mat(idx1,idx1);
M12 = Mat(idx1,idx2);
M21 = Mat(idx2,idx1);
M22 = Mat(idx2,idx2);

R11 = eig(M11)
R22 = eig(M22)
R12 = eig(M21)
R21 = eig(M12)



%% step8：import parameters
kappa = 0.2;
m = 0.25;
q = 0.05;
omega_1 = 0.667;
omega_2 = 0.667;
lambda = 0.05 ;
theta = 0.01;
mu_1 = 0.005;
mu_2 = 0.005;
alpha = 0.006;
tau = 0.002;
gamma_1 = 0.071;
gamma_2 = 0.071;
f_1 = 0.01847;
f_2 = 0.01847;


B = readmatrix("tb_beta.xlsx",Sheet="jilin"); % Sheet can be changed
b = B(i, 2:end); 

beta_11 = b(1);
beta_12 = b(2);
beta_22 = b(3);
beta_21 = b(4);
D = readmatrix("br_dr.xlsx",Sheet="jilin");
br = D(i, 2);
dr = D(i, 3);
N = readmatrix("s1_s2.xlsx",Sheet="jilin");
S_1 = N(i, 3);
S_2 = N(i, 2);

Reff = zeros(size(B,1), 4); 
for i = 1:size(B,1) 
    b = B(i, 2:end); 
    % year = B(i,1);
    beta_11 = b(1);
    beta_12 = b(2);
    beta_22 = b(3);
    beta_21 = b(4);
    br = D(i, 2);
    dr = D(i, 3);
    S_1 = N(i, 3);
    S_2 = N(i, 2);
    Reff(i,1) = eval(R11(end)); 
    Reff(i,2) = eval(R12(end));
    Reff(i,3) = eval(R22(end));
    Reff(i,4) = eval(R21(end));
end 
 
%% step9：
% R0 = eval(R11(end))


