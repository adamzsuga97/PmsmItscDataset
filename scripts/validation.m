% Always run both Simulink models before running this script

t = i_algebraic.Time;

% Figure 1: Comparison of i_u currents
figure(1);
plot(t,i_algebraic.Data(1:end,1),'LineWidth',2);
hold on;
plot(t,i_inverse.Data(1:end,1),'LineWidth',2)
title('Comparison of i_u'); xlabel('t[s]'); ylabel('i[A]');
legend('algebraic','inverse','FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'LineWidth', 1.5); % Increase the font size and axis line width
hold off;

% Figure 2: Comparison of i_v currents
figure(2);
plot(t,i_algebraic.Data(1:end,2),'LineWidth',2);
hold on;
plot(t,i_inverse.Data(1:end,2),'LineWidth',2)
title('Comparison of i_v'); xlabel('t[s]'); ylabel('i[A]');
legend('algebraic','inverse');
hold off;

% Figure 3: Comparison of i_w currents
figure(3);
plot(t,i_algebraic.Data(1:end,3),'LineWidth',2);
hold on;
plot(t,i_inverse.Data(1:end,3),'LineWidth',2)
title('Comparison of i_w'); xlabel('t[s]'); ylabel('i[A]');
legend('algebraic','inverse','FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'LineWidth', 1.5); % Increase the font size and axis line width
hold off;

% Figure 4: Differences between algebraic and inverse currents
figure(4);
subplot(1,3,1);
plot(t,i_algebraic.Data(1:end,1) - i_inverse.Data(1:end,1),'LineWidth',2);
title('Difference of i_u'); xlabel('t[s]'); ylabel('\Delta i[A]');

subplot(1,3,2);
plot(t,i_algebraic.Data(1:end,2) - i_inverse.Data(1:end,2),'LineWidth',2);
title('Difference of i_b'); xlabel('t[s]'); ylabel('\Delta i[A]');

subplot(1,3,3);
plot(t,i_algebraic.Data(1:end,3) - i_inverse.Data(1:end,3),'LineWidth',2);
title('Difference of i_w'); xlabel('t[s]'); ylabel('\Delta i[A]');

% Figure 5: Fault current comparison
figure(5);
plot(t, dqf_currents_algebraic.Data(1:end,3),'LineWidth', 2);
hold on;
plot(t, dqf_currents_inverse.Data(1:end,3),'LineWidth', 2);
title('Fault current comparison'); xlabel('t[s]'); ylabel('i[A]');
legend('algebraic','inverse','FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'LineWidth', 1.5); % Increase the font size and axis line width
hold off;

% Figure 6: Fault current difference
figure(6);
plot(t, dqf_currents_inverse.Data(1:end,3) - dqf_currents_algebraic.Data(1:end,3));
title('Fault current difference'); xlabel('t[s]'); ylabel('\Delta i[A]');

% Figure 7: Torque comparison
figure(7);
plot(t, m_algebraic.Data(1:end,1),'LineWidth',2);
hold on;
plot(t, m_inverse.Data(1:end,1),'LineWidth',2);
title('Torque comparison'); xlabel('t[s]'); ylabel('\Delta M[Nm]');
legend('algebraic','inverse','FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'LineWidth', 1.5); % Increase the font size and axis line width
hold off;

% Figure 8: Torque difference
figure(8);
plot(t, m_algebraic.Data(1:end,1) - m_inverse.Data(1:end,1),'LineWidth',2);
title('Torque difference'); xlabel('t[s]'); ylabel('\Delta M[Nm]');

% Figure 9: Three-phase currents of the inverse map model
figure(9);
plot(t, i_inverse.Data(1:end,1:end), 'LineWidth', 2);
title('Three-phase currents of the inverse map model'); xlabel('t[s]'); ylabel('i[A]');
legend('i_u','i_v','i_w','FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'LineWidth', 1.5); % Increase the font size and axis line width

% Assuming t is the time vector and dqf_currents_algebraic.Data is the data matrix

% Extract the signal for algebraic model
signal_algebraic = dqf_currents_algebraic.Data(1:end,1);

% Extract the signal for inverse model
signal_inverse = dqf_currents_inverse.Data(1:end,1);

% Compute the sampling frequency
Fs = 1 / mean(diff(t)); % Sampling frequency

% Compute the FFT for algebraic model
L = length(signal_algebraic); % Length of signal
Y_algebraic = fft(signal_algebraic); % Compute the FFT

% Compute the two-sided spectrum P2 for algebraic model
P2_algebraic = abs(Y_algebraic/L);

% Compute the single-sided spectrum P1 for algebraic model
P1_algebraic = P2_algebraic(1:L/2+1);
P1_algebraic(2:end-1) = 2*P1_algebraic(2:end-1);

% Compute the FFT for inverse model
Y_inverse = fft(signal_inverse); % Compute the FFT

% Compute the two-sided spectrum P2 for inverse model
P2_inverse = abs(Y_inverse/L);

% Compute the single-sided spectrum P1 for inverse model
P1_inverse = P2_inverse(1:L/2+1);
P1_inverse(2:end-1) = 2*P1_inverse(2:end-1);

% Define the frequency domain f
f = Fs*(0:(L/2))/L;

% Plot the frequency spectrum
figure(10);
plot(f, P1_algebraic, 'LineWidth', 2);
hold on;
plot(f, P1_inverse, 'LineWidth', 2);
title('Single-Sided Amplitude Spectrum of i_d current');
xlabel('f (Hz)');
ylabel('|P1(f)|');
legend('Algebraic Model', 'Inverse Model');
hold off;
