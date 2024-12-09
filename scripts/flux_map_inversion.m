%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PMSM ITSC model current table fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% psid,psiq,psif tables extracted from FEM simulations%
% Tables are recorded in -250 to 250 A in case of id,iq and -1250 to 1250 A
% in case of if. This scripts maps the inverse table of the flux tables:
% psid(id,iq,if), psiq(id,iq,if) , psif(id,iq,if) ==>
% id(psid,psiq,psif), iq(psid,psiq,psif), if0(psid,psiq,psif)


%% Load flux tables and extract data
clc 
disp("Initialize flux table inversion to generate current tables...")

load dataset/faulted_PMSM_data.mat


psid = permute(psid_map,[1 2 3 4]);
psiq = permute(psiq_map,[1 2 3 4]);
psif = permute(psif_map,[1 2 3 4]);

% Get table dimension
d_dim       = size(id_grid,2);
q_dim       = size(iq_grid,2);
f_dim       = size(if_grid,2);
theta_dim   = size(theta_grid,2);


%% Load data into 3D matrices


[id_mesh,iq_mesh,if_mesh] = ndgrid(id_grid,iq_grid,if_grid);

min_i = min(id_mesh(:));
max_i = max(iq_mesh(:));


min_if = min(if_mesh(:));
max_if = max(if_mesh(:));

extrap_ratio = 15/11; % to have 15 element matrices

id_grid_new = linspace(min_i*extrap_ratio,max_i*extrap_ratio,round(d_dim*extrap_ratio));
iq_grid_new = id_grid_new';
if_grid_new = linspace(min_if*extrap_ratio,max_if*extrap_ratio,round(d_dim*extrap_ratio));

[id_mesh_new,iq_mesh_new,if_mesh_new] = ndgrid(id_grid_new,iq_grid_new,if_grid_new);

psid_map_extrap = zeros(round(d_dim*extrap_ratio),round(q_dim*extrap_ratio),round(f_dim*extrap_ratio),theta_dim);
psiq_map_extrap = zeros(round(d_dim*extrap_ratio),round(q_dim*extrap_ratio),round(f_dim*extrap_ratio),theta_dim);
psif_map_extrap = zeros(round(d_dim*extrap_ratio),round(q_dim*extrap_ratio),round(f_dim*extrap_ratio),theta_dim);

%% GOverwrite dimensions for the extrapolated tables
d_dim = round(d_dim*extrap_ratio);
q_dim = round(q_dim*extrap_ratio);
f_dim = round(f_dim*extrap_ratio);

%Extrapolation for better inversion at the surfaces edges
for i = 1:theta_dim

    psid_act = psid(:,:,:,i);
    psiq_act = psiq(:,:,:,i);
    psif_act = psif(:,:,:,i);

    PSID = griddedInterpolant(id_mesh,iq_mesh,if_mesh,psid_act,'linear','linear');
    PSIQ = griddedInterpolant(id_mesh,iq_mesh,if_mesh,psiq_act,'linear','linear');
    PSIF = griddedInterpolant(id_mesh,iq_mesh,if_mesh,psif_act,'linear','linear');

    psid_map_extrap(:,:,:,i) = PSID(id_mesh_new,iq_mesh_new,if_mesh_new);
    psiq_map_extrap(:,:,:,i) = PSIQ(id_mesh_new,iq_mesh_new,if_mesh_new);
    psif_map_extrap(:,:,:,i) = PSIF(id_mesh_new,iq_mesh_new,if_mesh_new);
end


%% Visualize extrapolated current space
points_Id = reshape(id_mesh_new, [], 1);
points_Iq = reshape(iq_mesh_new, [], 1);
points_If = reshape(if_mesh_new, [], 1);

figure(1);
scatter3(points_Id, points_Iq, points_If);
xlabel('I_d[A]');
ylabel('I_q[A]');
zlabel('I_f[A]]');
title('The current domain of fluxes');
grid on;

%% Visualize extrapolated flux maps
figure(2);
subplot(2, 3, 1);
surf(squeeze(iq_mesh_new(4,:,:)),squeeze(if_mesh_new(4,:,:)),squeeze(psid_map_extrap(4,:,:,1))); 
title('Normalized flux map of original f_{\Psi_d} at I_d = -50 A'); xlabel('i_{q,norm}'); ylabel('i_{f,norm}'); zlabel('\Psi_{d,norm}');

subplot(2, 3, 2);
surf(squeeze(iq_mesh_new(4,:,:)),squeeze(if_mesh_new(4,:,:)),squeeze(psiq_map_extrap(4,:,:,1))); 
title('Normalized flux map of original f_{\Psi_q} at I_d = -50 A'); xlabel('i_{q,norm}'); ylabel('i_{f,norm}'); zlabel('\Psi_{q,norm}');

subplot(2, 3, 3);
surf(squeeze(iq_mesh_new(4,:,:)),squeeze(if_mesh_new(4,:,:)),squeeze(psif_map_extrap(4,:,:,1))); 
title('Normalized flux map of original f_{\Psi_f} at I_d = -50 A'); xlabel('i_{q,norm}'); ylabel('i_{f,norm}'); zlabel('\Psi_{f,norm}');


subplot(2, 3, 4);
surf(id_mesh_new(:,:,6),iq_mesh_new(:,:,6),psid_map_extrap(:,:,6,1)); 
title('Normalized flux map of original f_{\Psi_d} at I_f = 0 A'); xlabel('i_{d,norm}'); ylabel('i_{q,norm}'); zlabel('\Psi_{d,norm}');

subplot(2, 3, 5);
surf(id_mesh_new(:,:,6),iq_mesh_new(:,:,6),psiq_map_extrap(:,:,6,1));
title('Normalized flux map of original f_{\Psi_q} at I_f = 0 A'); xlabel('i_{d,norm}'); ylabel('i_{q,norm}'); zlabel('\Psi_{q,norm}');

subplot(2, 3, 6);
surf(id_mesh_new(:,:,6),iq_mesh_new(:,:,6),psif_map_extrap(:,:,6,1));
title('Normalized flux map of original f_{\Psi_f} at I_f = 0 A'); xlabel('i_{d,norm}'); ylabel('i_{q,norm}'); zlabel('\Psi_{f,norm}');

%% Invertibility Check of Flux Maps
disp("Checking invertibility of the flux mappings...");

% Preallocate matrices for jacobian determinant calculation
jacobian_det = zeros(d_dim, q_dim, f_dim,theta_dim);

dpsid_did = zeros(d_dim,q_dim,f_dim,theta_dim);
dpsid_diq = zeros(d_dim,q_dim,f_dim,theta_dim);
dpsid_dif = zeros(d_dim,q_dim,f_dim,theta_dim);

dpsiq_did = zeros(d_dim,q_dim,f_dim,theta_dim);
dpsiq_diq = zeros(d_dim,q_dim,f_dim,theta_dim);
dpsiq_dif = zeros(d_dim,q_dim,f_dim,theta_dim);

dpsif_did = zeros(d_dim,q_dim,f_dim,theta_dim);
dpsif_diq = zeros(d_dim,q_dim,f_dim,theta_dim);
dpsif_dif = zeros(d_dim,q_dim,f_dim,theta_dim);

% Loop through each point in the grid and calculate the Jacobian determinant
for m = 1:1
    for d = 1:d_dim
        for q = 1:q_dim
            
            for f = 1:f_dim
                [dpsid_diq(:,:,:,m), dpsid_did(:,:,:,m), dpsid_dif(:,:,:,m)] = gradient(permute(psid_map_extrap(:,:,:,m),[2 1 3 4]),id_grid_new, iq_grid_new, if_grid_new);
                [dpsiq_diq(:,:,:,m), dpsiq_did(:,:,:,m), dpsiq_dif(:,:,:,m)] = gradient(permute(psiq_map_extrap(:,:,:,m),[2 1 3 4]),id_grid_new, iq_grid_new, if_grid_new);
                [dpsif_diq(:,:,:,m), dpsif_did(:,:,:,m), dpsif_dif(:,:,:,m)] = gradient(permute(psif_map_extrap(:,:,:,m),[2 1 3 4]),id_grid_new, iq_grid_new, if_grid_new);

                % Form the Jacobian matrix at each point
                J = [dpsid_did(d, q, f, m), dpsid_diq(d, q, f, m), dpsid_dif(d, q, f, m);
                     dpsiq_did(d, q, f, m), dpsiq_diq(d, q, f, m), dpsiq_dif(d, q, f, m);
                     dpsif_did(d, q, f, m), dpsif_diq(d, q, f, m), dpsif_dif(d, q, f, m)];
                
                % Calculate the determinant of the Jacobian
                jacobian_det(d, q, f, m) = det(J);
            end
            
        end
    end
end

figure(3);
surf(iq_grid_new,if_grid_new,squeeze(jacobian_det(7,:,:,1)));
title('jacobian determinant with normalized values for iq,if dependency');xlabel('i_{d,norm}'); ylabel('i_{q,norm}'); zlabel('det(J)');

%% Invert flux tables with fmincon
% Lower and upper bounds
lb = [min(id_mesh_new(:))*1, min(iq_mesh_new(:))*1, min(if_mesh_new(:))*1];
ub = [max(id_mesh_new(:))*1, max(iq_mesh_new(:))*1, max(if_mesh_new(:))*1];

% Constraints: id, iq, if should be within the bounds
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = @(xyz) deal([], []);

% Solve the optimization problem using fmincon with interior point method
options = optimoptions('fmincon', 'Algorithm','interior-point','OptimalityTolerance',1e-8,'Display', 'off');

start_value = [0,0,0];

% Preallocate memory for current maps
id_map = zeros(d_dim,q_dim,f_dim,theta_dim);
iq_map = zeros(d_dim,q_dim,f_dim,theta_dim);
if_map = zeros(d_dim,q_dim,f_dim,theta_dim);

nr_of_iters = d_dim*q_dim*theta_dim;
iterations = nr_of_iters;
progress = 0;

disp("Start inversion process...")
fprintf('\rProgress: 0.00%', progress);

% Find the minimum and maximum values for psid, psiq, and psif
psid_min = min(psid_map_extrap(:));
psiq_min = min(psiq_map_extrap(:));
psif_min = min(psif_map_extrap(:));
psid_max = max(psid_map_extrap(:));
psiq_max = max(psiq_map_extrap(:));
psif_max = max(psif_map_extrap(:));

% Create grids for psid, psiq, and psif using linspace
psid_grid = linspace(psid_min,psid_max,d_dim);
psiq_grid = linspace(psiq_min,psiq_max,q_dim);
psif_grid = linspace(psif_min,psif_max,f_dim);

psid_act = zeros(d_dim,q_dim,f_dim,theta_dim);
psiq_act = zeros(d_dim,q_dim,f_dim,theta_dim);
psif_act = zeros(d_dim,q_dim,f_dim,theta_dim);

% Generate 3D grids for psid, psiq, and psif
[psid_mesh,psiq_mesh,psif_mesh] = ndgrid(psid_grid,psiq_grid,psif_grid);

tic
for m = 1:theta_dim
    for i = 1:d_dim
        for j = 1:q_dim
            tic
            parfor k = 1:f_dim
                % extract actual value
                psid_act = psid_map_extrap(:,:,:,m);
                psiq_act = psiq_map_extrap(:,:,:,m);
                psif_act = psif_map_extrap(:,:,:,m);

                psid_act_mid = (abs(max(psid_act(:))) + abs(min(psid_act(:))))/2;
                psiq_act_mid = (abs(max(psiq_act(:))) + abs(min(psiq_act(:))))/2;
                psif_act_mid = (abs(max(psif_act(:))) + abs(min(psif_act(:))))/2;

                FD = griddedInterpolant(id_mesh_new,iq_mesh_new,if_mesh_new, psid_act,'makima','makima');
                FQ = griddedInterpolant(id_mesh_new,iq_mesh_new,if_mesh_new, psiq_act,'makima','makima');
                FF = griddedInterpolant(id_mesh_new,iq_mesh_new,if_mesh_new, psif_act,'makima','makima');

                FID = scatteredInterpolant(psid_act(:),psiq_act(:),psif_act(:), id_mesh_new(:),'linear','linear');
                FIQ = scatteredInterpolant(psid_act(:),psiq_act(:),psif_act(:), iq_mesh_new(:),'linear','linear');
                FIF = scatteredInterpolant(psid_act(:),psiq_act(:),psif_act(:), if_mesh_new(:),'linear','linear');

                start_value_id = FID(psid_mesh(i,j,k) , psiq_mesh(i,j,k) , psif_mesh(i,j,k));
                start_value_iq = FIQ(psid_mesh(i,j,k) , psiq_mesh(i,j,k) , psif_mesh(i,j,k));
                start_value_if = FIF(psid_mesh(i,j,k) , psiq_mesh(i,j,k) , psif_mesh(i,j,k));

                % % 
                psi_j            = @(xyz)           sqrt(((FD(xyz(1), xyz(2), xyz(3)) - psid_mesh(i,j,k))/psid_act_mid)^2 +  ...
                                                         ((FQ(xyz(1), xyz(2), xyz(3)) - psiq_mesh(i,j,k))/psiq_act_mid)^2 +  ...
                                                         ((FF(xyz(1), xyz(2), xyz(3)) - psif_mesh(i,j,k))/psif_act_mid)^2); 

                solution = fmincon(psi_j, [start_value_id,start_value_iq,start_value_if], A, b, Aeq, beq, lb, ub, nonlcon, options);
              
                id_map(i,j,k,m)      = solution(1);
                iq_map(i,j,k,m)      = solution(2);
                if_map(i,j,k,m)      = solution(3);
   
            end
            toc
            iterations = iterations - 1;
            progress = (1 - iterations/nr_of_iters) * 100;
            fprintf('\rProgress:  %.2f%% and %.2f iterations in d-direction', progress,iterations);
        end
    end
end
toc


% Define the filename for the .mat file
filename = 'output_data.mat';

% Save the variables to the .mat file
save(filename, 'psid_grid', 'psiq_grid', 'psif_grid', 'theta_grid', 'id_map', 'iq_map', 'if_map');

disp('Variables saved to output_data.mat');

%% Validation of current map
% Preallocate memory for checking the flux maps
psid_check = zeros(d_dim,q_dim,f_dim,theta_dim);
psiq_check = zeros(d_dim,q_dim,f_dim,theta_dim);
psif_check = zeros(d_dim,q_dim,f_dim,theta_dim);

id_check = zeros(d_dim,q_dim,f_dim,theta_dim);
iq_check = zeros(d_dim,q_dim,f_dim,theta_dim);
if_check = zeros(d_dim,q_dim,f_dim,theta_dim);    

mse_id = zeros(theta_dim,1);
mse_iq = zeros(theta_dim,1);
mse_if = zeros(theta_dim,1);

error_id = zeros(d_dim,q_dim,f_dim,theta_dim);
error_iq = zeros(d_dim,q_dim,f_dim,theta_dim);
error_if = zeros(d_dim,q_dim,f_dim,theta_dim);

[psid_mesh,psiq_mesh,psif_mesh] = ndgrid(psid_grid,psiq_grid,psif_grid);

iterations = d_dim*q_dim*1
% Loop through each angle (theta)
for m = 1:1 % Currently only processes the first angle. You might want to change this to i = 1:theta_dim
    for i = 1:d_dim
        for j = 1:q_dim
            for k = 1:f_dim
                id_act = id_map(:,:,:,m);
                iq_act = iq_map(:,:,:,m);
                if_act = if_map(:,:,:,m);

                % Extract the extended flux maps for the current angle
                psid_check_act = psid_map_extrap(:,:,:,m);
                psiq_check_act = psiq_map_extrap(:,:,:,m);
                psif_check_act = psif_map_extrap(:,:,:,m);

                % Create scattered interpolants for Id, Iq, and If as functions of psid, psiq, and psif
                FD = griddedInterpolant(psid_mesh,psiq_mesh,psif_mesh,id_act,'linear','linear');
                FQ = griddedInterpolant(psid_mesh,psiq_mesh,psif_mesh,iq_act,'linear','linear');
                FF = griddedInterpolant(psid_mesh,psiq_mesh,psif_mesh,if_act,'linear','linear');

                % Interpolate the current values onto the flux grids
                id_check(i,j,k,m) = FD(psid_check_act(i,j,k),psiq_check_act(i,j,k),psif_check_act(i,j,k));
                iq_check(i,j,k,m) = FQ(psid_check_act(i,j,k),psiq_check_act(i,j,k),psif_check_act(i,j,k));
                if_check(i,j,k,m) = FF(psid_check_act(i,j,k),psiq_check_act(i,j,k),psif_check_act(i,j,k));

                % Calculate relative error in percentage
                error_id(i,j,k,m) = abs((id_check(i,j,k,m) - id_mesh_new(i,j,k)));
                error_iq(i,j,k,m) = abs((iq_check(i,j,k,m) - iq_mesh_new(i,j,k)));
                error_if(i,j,k,m) = abs((if_check(i,j,k,m) - if_mesh_new(i,j,k)));
            end
            iterations = iterations-1;
            fprintf('\rProgress:  %.2f iterations in d-direction',iterations);
        end
    end
end

for i = 1:1
    % Calculate Mean Squared Error for each flux component

    D_square = error_id(:,:,:,i).^2;
    mse_d = mean(D_square(:));
    rmse_d = sqrt(mse_d);


    D_square = error_iq.^2;
    mse_q = mean(D_square(:));
    rmse_q = sqrt(mse_q);

    D_square = error_if.^2;
    mse_f = mean(D_square(:));
    rmse_f = sqrt(mse_f);


end

if_idx = 8;
% Plot the results for visualization
figure(5);
subplot(2,3,1)
surf(psid_mesh(:,:,if_idx),psiq_mesh(:,:,if_idx),squeeze(id_map(:,:,if_idx,1)))
title('I_d(\Psi_d,\Psi_q,\Psi_f) map')
xlabel('\psi_d [Vs]')
ylabel('\psi_q [Vs]')
zlabel('I_d [A]')

subplot(2,3,2)
surf(psid_mesh(:,:,if_idx),psiq_mesh(:,:,if_idx),iq_map(:,:,if_idx,1))
title('I_q(\Psi_d,\Psi_q,\Psi_f) map')
xlabel('\psi_d [Vs]')
ylabel('\psi_q [Vs]')
zlabel('I_q [A]')

subplot(2,3,3)
surf(psid_mesh(:,:,if_idx),psiq_mesh(:,:,if_idx),if_map(:,:,if_idx,1))
title('I_f(\Psi_d,\Psi_q,\Psi_f) map')
xlabel('\psi_d [Vs]')
ylabel('\psi_q [Vs]')
zlabel('I_f [A]')

subplot(2,3,4)
surf(squeeze(iq_mesh_new(if_idx,:,:)),squeeze(if_mesh_new(if_idx,:,:)),squeeze(error_id(if_idx,:,:,1)));
title('i_d error')
xlabel('i_q [Wb]')
ylabel('i_f [Wb]')
zlabel('[A]')

subplot(2,3,5)
surf(squeeze(iq_mesh_new(if_idx,:,:)),squeeze(if_mesh_new(if_idx,:,:)),squeeze(error_iq(if_idx,:,:,1)));
title('i_q error')
xlabel('i_q [Wb]')
ylabel('i_f [Wb]')
zlabel('[A]')


subplot(2,3,6)
surf(squeeze(iq_mesh_new(if_idx,:,:)),squeeze(if_mesh_new(if_idx,:,:)),squeeze(error_if(if_idx,:,:,1)));
title('i_f error')
xlabel('i_q [Wb]')
ylabel('i_f [Wb]')
zlabel('[A]')

