%% Script to Remake Under 60 yrs and Above 60 yrs COVID Figure
%% Define Variables
% Ensure these variables are defined before running the script
% h, S_1, S_2, I_1, I_2, R_1, R_2, D_1, D_2, N_1, N_2, mu

% Set N1 and N2
Npobla=[2456344 + 2866927+ 3570635, 2307647 ];
N_1 = Npobla(1); 
N_2 = Npobla(2);
 
mu = 0.0051; % based on delta from previous script

% Infected_cases
CfDe00a19=[1	1	1	2	7	7	0	2	10	3	6	3	8	8	3	8	3	5	4	3	5	2	3	6	1	4	9	5	3	5	3]';
CfDe20a39=[4	7	10	11	22	27	27	27	43	19	16	21	24	29	11	20	20	31	30	12	12	21	16	9	6	23	23	28	22	25	13]';
CfDe40a59=[7	10	13	11	16	16	32	28	24	22	8	16	21	16	6	15	14	32	28	2	12	22	8	10	5	16	24	19	14	20	12]';
CfDe60a99=[1	6	7	1	4	4	0	8	16	3	5	6	3	10	4	5	7	6	21	0	6	20	8	2	13	4	4	7	9	9	5]';
CfNinhos= [1	1	1	2	7	7	0	2	10	3	6	3	8	3	3	8	3	5	4	3	5	2	3	6	1	0	9	0	3	4	3]';

CFcov19_grupo_edades = [CfDe00a19+CfDe20a39+CfDe40a59,CfDe60a99];

% Set I_1 and I_2
I_1 = CFcov19_grupo_edades(:, 1);
I_2 = CFcov19_grupo_edades(:, 2); 

% Recovered_cases
AlDe00a19=[	0	0	2	0	5	5	1	0	2	1	5	5	1	2	5	5	3	9	3	7	5	2	4	6	12	2	3	3	1	3	1]';
AlDe20a39=[1	2	4	5	3	3	4	3	1	7	12	7	7	4	19	23	46	49	4	22	19	22	10	7	28	16	9	25	18	23	5]';
AlDe40a59=[2	1	3	4	6	7	6	4	3	4	13	6	5	10	8	26	41	54	5	21	21	11	17	7	14	12	12	17	24	16	4]';
AlDe60a99=[0	1	0	1	2	1	1	2	2	4	7	9	2	0	0	11	6	11	1	4	10	8	5	3	3	8	0	7	8	5	3]';
Almenor18=[	0	0	2	0	5	4	1	0	2	0	4	5	1	1	4	4	1	5	3	7	3	2	3	5	10	2	3	3	0	3	0]';
Alta19=   [0	0	0	0	0	1	0	0	0	1	1	0	0	1	1	1	2	4	0	0	2	0	1	1	2	0	0	0	1	0	1]';

ALcov19_grupo_edades = [AlDe00a19+AlDe20a39+AlDe40a59,AlDe60a99];

% Set R_1 and R_2
R_1 = ALcov19_grupo_edades(:, 1);
R_2 = ALcov19_grupo_edades(:, 2);

% Death_cases
FaDe00a19=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 0 0 0 0 0 0 0 0]';
FaDe20a39=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 0 0 0 0 0 0 0 0]';
FaDe40a59=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0 0 0 0 0 1 0 0 0]';
FaDe60a99=[0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	2	0 0 1 0 0 1 0 0 1]';

FAcov19_grupo_edades = [FaDe00a19+FaDe20a39+FaDe40a59,FaDe60a99];

% Set D_1 and D_2
D_1 = FAcov19_grupo_edades(:, 1);
D_2 = FAcov19_grupo_edades(:, 2);

% Need Susceptible Population numbers for each age group. 
% Set S_1 and S_2
S_1 = N_1 - I_1 - R_1 - D_1;
S_2 = N_2 - I_2 - R_2 - D_2;


%% Estimate Parameters
h = ones(length(S_1) - 1, 1); % Create an array for h with uniform time steps

%% S1 Population
num_points = length(S_1) - 1;
S1_matrixA = zeros(num_points, 2);
S1_matrixB = zeros(num_points, 1);
for i = 1:num_points
    S1_matrixA(i, :) = [-h(i) * S_1(i) * I_1(i) / N_1, -h(i) * S_1(i) * I_2(i) / N_1];
    S1_matrixB(i) = S_1(i+1) - S_1(i) + h(i) * mu * S_1(i) - h(i) * mu * N_1;
end
A = S1_matrixA;
b = S1_matrixB;
xhat = lsqnonneg(A, b); % Use lsqnonneg to ensure non-negative constraints
beta_11 = xhat(1);
beta_12 = xhat(2);

%% S2 Population
num_points = length(S_2) - 1;
S2_matrixA = zeros(num_points, 2);
S2_matrixB = zeros(num_points, 1);
for i = 1:num_points
    S2_matrixA(i, :) = [-h(i) * S_2(i) * I_2(i) / N_2, -h(i) * S_2(i) * I_1(i) / N_2];
    S2_matrixB(i) = S_2(i+1) - S_2(i) + h(i) * mu * S_2(i) - h(i) * mu * N_2;
end
A = S2_matrixA;
b = S2_matrixB;
xhat = lsqnonneg(A, b); % Use lsqnonneg to ensure non-negative constraints
beta_22 = xhat(1);
beta_21 = xhat(2);

%% I1 Population
num_points = length(I_1) - 1;
I1_matrixA = zeros(num_points, 3);
I1_matrixB = zeros(num_points, 1);
for i = 1:num_points
    I1_matrixA(i, :) = [h(i) * S_1(i) * I_1(i) / N_1, h(i) * S_1(i) * I_2(i) / N_1, -h(i) * I_1(i)];
    I1_matrixB(i) = I_1(i+1) - I_1(i) + h(i) * mu * I_1(i);
end
A = I1_matrixA;
b = I1_matrixB;
xhat = lsqnonneg(A, b); % Use lsqnonneg to ensure non-negative constraints
beta_11 = xhat(1);
beta_12 = xhat(2);
alpha1_plus_kappa1 = xhat(3);

%% I2 Population
num_points = length(I_2) - 1;
I2_matrixA = zeros(num_points, 3);
I2_matrixB = zeros(num_points, 1);
for i = 1:num_points
    I2_matrixA(i, :) = [h(i) * S_2(i) * I_2(i) / N_2, h(i) * S_2(i) * I_1(i) / N_2, -h(i) * I_2(i)];
    I2_matrixB(i) = I_2(i+1) - I_2(i) + h(i) * mu * I_2(i);
end
A = I2_matrixA;
b = I2_matrixB;
xhat = lsqnonneg(A, b); % Use lsqnonneg to ensure non-negative constraints
beta_22 = xhat(1);
beta_21 = xhat(2);
alpha2_plus_kappa2 = xhat(3);

%% R1 Population
num_points = length(R_1) - 1;
R1_matrixA = zeros(num_points, 1);
R1_matrixB = zeros(num_points, 1);
for i = 1:num_points
    R1_matrixA(i, :) = h(i) * I_1(i);
    R1_matrixB(i) = R_1(i+1) - R_1(i) + h(i) * mu * R_1(i);
end
A = R1_matrixA;
b = R1_matrixB;
xhat = lsqnonneg(A, b); % Use lsqnonneg to ensure non-negative constraints
alpha_1 = xhat(1);

%% R2 Population
num_points = length(R_2) - 1;
R2_matrixA = zeros(num_points, 1);
R2_matrixB = zeros(num_points, 1);
for i = 1:num_points
    R2_matrixA(i, :) = h(i) * I_2(i);
    R2_matrixB(i) = R_2(i+1) - R_2(i) + h(i) * mu * R_2(i);
end
A = R2_matrixA;
b = R2_matrixB;
xhat = lsqnonneg(A, b); % Use lsqnonneg to ensure non-negative constraints
alpha_2 = xhat(1);

%% D1 Population
num_points = length(D_1) - 1;
D1_matrixA = zeros(num_points, 1);
D1_matrixB = zeros(num_points, 1);
for i = 1:num_points
    D1_matrixA(i, :) = h(i) * I_1(i);
    D1_matrixB(i) = D_1(i+1) - D_1(i);
end
A = D1_matrixA;
b = D1_matrixB;
xhat = lsqnonneg(A, b); % Use lsqnonneg to ensure non-negative constraints
kappa_1 = xhat(1);

%% D2 Population
num_points = length(D_2) - 1;
D2_matrixA = zeros(num_points, 1);
D2_matrixB = zeros(num_points, 1);
for i = 1:num_points
    D2_matrixA(i, :) = h(i) * I_2(i);
    D2_matrixB(i) = D_2(i+1) - D_2(i);
end
A = D2_matrixA;
b = D2_matrixB;
xhat = lsqnonneg(A, b); % Use lsqnonneg to ensure non-negative constraints
kappa_2 = xhat(1);
