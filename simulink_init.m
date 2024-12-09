clear
clc

addpath('dataset')
addpath('scripts')

load dataset/faulted_PMSM_data.mat



ModelParam.Udc=500;      %[V]            %DC link voltage
ModelParam.UNom = 215;   %[V]            %Motor nominal voltage
ModelParam.INom = 250;   %[A]            %Motor nominal current
ModelParam.IQmax = ModelParam.INom *(2/3);
ModelParam.FNom = 245;   %[Hz]           %Motor nominal electrical frequency


% Common parameters for both model
ModelParam.Ld  = 0.00190051776107054;       % [H]
ModelParam.Lq  = 0.00567347930411143;       % [H]
ModelParam.Rs  =  0.0523;                  % [Ohm]  
ModelParam.u = 1/72;
ModelParam.Rf = 0.01;
ModelParam.Psi = 0.169954396290924;        % [Vs]       %flux induced by magnet
ModelParam.p   = 4;               %[-]        %number of pole pairs

% Model parameters for algebraic contstraint model
ModelParam.MachineModel.PsidMap = psid_map;
ModelParam.MachineModel.PsiqMap = psiq_map;
ModelParam.MachineModel.PsifMap = psif_map;
ModelParam.MachineModel.IdGrid  = id_grid;
ModelParam.MachineModel.IqGrid  = iq_grid;
ModelParam.MachineModel.IfGrid  = if_grid;
ModelParam.MachineModel.ThetaGrid = theta_grid;

% Model parameters for inverse flux map model
ModelParam.MachineModel.IdMap   = id_map;
ModelParam.MachineModel.IqMap   = iq_map;
ModelParam.MachineModel.IfMap   = if_map;
ModelParam.MachineModel.PsidGrid= psid_grid;
ModelParam.MachineModel.PsiqGrid= psiq_grid;
ModelParam.MachineModel.PsifGrid= psif_grid;

% Initial values of the integrators
ModelParam.MachineModel.PsidInit= psid_map(6,6,6,1);
ModelParam.MachineModel.PsiqInit= psiq_map(6,6,6,1);
ModelParam.MachineModel.PsifInit= psif_map(6,6,6,1);

% Torque map
ModelParam.MachineModel.TorqueMap = torque_map;

% Mechanical data
ModelParam.Theta=0.0334;    %[kg*m2]        %inertia
ModelParam.D=0.0;         %[Nm/(rad/s)]   %damping coefficient of viscous friction
ModelParam.Mload=100;     %[Nm]           %load torq
% Tsl=0.0001;    %[Nm]           %load sample time(commented out because it is contionus model)

%% Controller parameters

ControllerParam.Kpid = 300*ModelParam.Rs;
ControllerParam.Ipid = ModelParam.Ld / ModelParam.Rs;

ControllerParam.Kpiq = ControllerParam.Kpid;
ControllerParam.Ipiq = ModelParam.Lq / ModelParam.Rs;

% create grids for proper visualization
[psid_mesh,psiq_mesh,psif_mesh] = ndgrid(psid_grid,psiq_grid,psif_grid);
[id_mesh,iq_mesh,if_mesh] = ndgrid(id_grid,iq_grid,if_grid);

if_idx = 5;
% Plot the results for visualization
figure(5);
subplot(2,3,1)
surf(psid_mesh(:,:,if_idx),psiq_mesh(:,:,if_idx),squeeze(id_map(:,:,if_idx,1)))
title('I_d(\psi_d,\psi_q,\psi_f) map', 'FontSize', 16)
xlabel('\psi_d [Vs]', 'FontSize', 14)
ylabel('\psi_q [Vs]', 'FontSize', 14)
zlabel('I_d [A]', 'FontSize', 14)

subplot(2,3,2)
surf(psid_mesh(:,:,if_idx),psiq_mesh(:,:,if_idx),iq_map(:,:,if_idx,1))
title('I_q(\psi_d,\psi_q,\psi_f) map', 'FontSize', 16)
xlabel('\psi_d [Vs]', 'FontSize', 14)
ylabel('\psi_q [Vs]', 'FontSize', 14)
zlabel('I_q [A]', 'FontSize', 14)

subplot(2,3,3)
surf(psid_mesh(:,:,if_idx),psiq_mesh(:,:,if_idx),if_map(:,:,if_idx,1))
title('I_f(\psi_d,\psi_q,\psi_f) map', 'FontSize', 16)
xlabel('\psi_d [Vs]', 'FontSize', 14)
ylabel('\psi_q [Vs]', 'FontSize', 14)
zlabel('I_f [A]', 'FontSize', 14)

subplot(2,3,4)
surf(squeeze(iq_mesh(if_idx,:,:)),squeeze(if_mesh(if_idx,:,:)),squeeze(psid_map(if_idx,:,:,1)));
title('\psi_d map at id = -50 A', 'FontSize', 16)
xlabel('i_q [A]', 'FontSize', 14)
ylabel('i_f [A]', 'FontSize', 14)
zlabel('[Wb]', 'FontSize', 14)

subplot(2,3,5)
surf(squeeze(iq_mesh(if_idx,:,:)),squeeze(if_mesh(if_idx,:,:)),squeeze(psiq_map(if_idx,:,:,1)));
title('\Psi_Q map at id = -50 A', 'FontSize', 16)
xlabel('i_q [A]', 'FontSize', 14)
ylabel('i_f [A]', 'FontSize', 14)
zlabel('[Wb]', 'FontSize', 14)

subplot(2,3,6)
surf(squeeze(iq_mesh(if_idx,:,:)),squeeze(if_mesh(if_idx,:,:)),squeeze(psif_map(if_idx,:,:,1)));
title('\Psi_F map at id = -50 A', 'FontSize', 16)
xlabel('i_q [A]', 'FontSize', 14)
ylabel('i_f [A]', 'FontSize', 14)
zlabel('[Wb]', 'FontSize', 14)

clearvars -except ModelParam && ControllerParam