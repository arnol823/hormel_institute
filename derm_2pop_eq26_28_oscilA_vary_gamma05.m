%% % vary gamma value to be lower than 1/8 to simulate cancer for the oscil
% epidermis population model
% derm_2pop_eq26_28_oscilA.m
% first half = normal gamma / carrying capacity, second = increasing gamma,
% no overlay



% Initial values
N0_initial = 51;  
N1_initial = 0;  

% Initial population vector
N_initial = [N0_initial; N1_initial];

% Parameters
alpha_minus_alpha2 = 0.0481;    
k_0 = (59/213)/50;   % averaged 0 and 0.0026
gamma_initial = 1/8; % initial gamma value
tau = 20; % time constant for exponential decay
change_time = 50; % time at which gamma starts to decrease
alpha2 = 0.2332; % previous: 0.2121

% alpha oscillating function
alpha_0 = 1; % Set the average growth rate
alpha_1 = 0.05; % Set the amplitude of oscillation
omega = (2 * pi) / 12; % Set the angular frequency of oscillation

% Oscillating alpha function
alpha = @(t) alpha_0 + alpha_1 * sin(omega * t);

% Time span
t_span = [0 100]; % Time from 0 to 100

% Define a range of final gamma values to vary
final_gamma_values = [1/30, 1/60, 1/90];

% Plotting
figure;
hold on;

% Loop through the different final gamma values
for gamma_final = final_gamma_values
    % Gamma function with exponential decay
    gamma_func = @(t) (t < change_time) * gamma_initial + ...
                      (t >= change_time) * (gamma_final + (gamma_initial - gamma_final) * exp(-(t - change_time) / tau));

    % Solve the differential equations
    [t, N] = ode45(@(t, N) populations(t, N, alpha, alpha_minus_alpha2, k_0, gamma_func, alpha2), t_span, N_initial);

    % Unpack the solution
    N0 = N(:, 1);
    N1 = N(:, 2);

    % Calculate the total population
    N_total = N0 + N1;

    % Plot the results
    % plot(t, N0, 'LineWidth', 2, 'DisplayName', ['N_0 (gamma_{final} = ' num2str(gamma_final) ')']);
    plot(t, N0, 'LineWidth', 2, 'DisplayName', ['N_0 (gamma_{final} = ' num2str(gamma_final) ')'], 'Color', 'b');
    plot(t, N1, '--', 'LineWidth', 2, 'DisplayName', ['N_1 (gamma_{final} = ' num2str(gamma_final) ')']);
    plot(t, N_total, ':', 'LineWidth', 2, 'DisplayName', ['Total Population (gamma_{final} = ' num2str(gamma_final) ')']);
end

hold off;
% title('Population Dynamics with Varying Gamma');
xlabel('Time', 'FontSize', 24);
ylabel('Cell Populations', 'FontSize', 24);
% legend('show');
grid on;


%% Function definition
function dNdt = populations(t, N, alpha, alpha_minus_alpha2, k_0, gamma_func, alpha2)
    % Unpack the population variables
    N0 = N(1);
    N1 = N(2);
    
    % Calculate the current value of alpha(t)
    current_alpha = alpha(t);
    
    % Calculate the current value of gamma(t)
    gamma = gamma_func(t);
    
    % Define k_1 in terms of the current value of gamma
    k_1 = gamma / 10; % Example: k_1 is 1/10th of the current gamma
    
    % Define the differential equations
    dN0dt = (current_alpha - alpha2) * N0 - k_0 * N0^2;
    dN1dt = alpha2 * N0 + k_0 * N0^2 - gamma * N1 - k_1 * N1^2;
    
    % Return the derivatives as a column vector
    dNdt = [dN0dt; dN1dt];
end
