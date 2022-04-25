clear; clc;

addpath('datas')
addpath('D:\MATLAB\functions')

%%
mesh = myMESH('MDPS_bearing_hexa.inp');
surf_mesh = myMESH('MDPS_bearing_surf.inp');

material.E = 210000; % [Mpa]
material.nu = 0.3; % Poisson
material.rho = 7.85e-09; % [ton/mm^3];

model = myFEM();
model.set_mesh(mesh);
model.set_material(material);
model.compute_MK('saved_MK_bearing.mat');

%%
% 1-56,	   mount 1 0,101,-15
% 57-112,  mount 2 -87.469,-50.5,-15
% 113-168, mount 3 87.468,-50.5,15
% 169-264, bearing 1 0,0,-5(위부터   -> 1~6
% 265-564, bearing 2 0,0,-57.5      -> 7~12
% 565-984, bearing 3 0,0,-152.5     -> 13~18
% 985, sensing 1 -> 19~21
% 986, sensing 2 -> 22~24
% 987, sensing 3 -> 25~27
% 988, 예비 (바닥 중심)
% 172-8997, surf
% 8998-22840, res

%%
    slave_node_idx{1} = 1:56;
    slave_node_idx{2} = 57:112;
    slave_node_idx{3} = 113:168;
    slave_node_idx{4} = 169:264;
    slave_node_idx{5} = 265:564;
    slave_node_idx{6} = 565:984;
    master_node_xyz = [0,101,-15;
                       -87.469,-50.5,-15;
                        87.468,-50.5,15;
                        0,0,-5;
                        0,0,-57.5;
                        0,0,-152.5];
model.compute_RBE2(slave_node_idx, master_node_xyz);
% model.load_RBE2('saved_Trbe_simple.mat');


%%
%     full_dof = model.num_dof;
%     BC_dof = 1:18;
%     org_dof = 1:full_dof;
%     org_dof(BC_dof) = [];
% model.compute_BC(BC_dof);
% [eig_vec_rbe,eig_val_rbe] = model.compute_eig(10);

% U_rbe = zeros(full_dof,10);
% U_rbe(org_dof,:) = eig_vec_rbe;
% 
% anime = myPLOT(surf_mesh);
% anime.fig_init([-120 120 -120 120 -250 20])
% view([-150,30])
% 
% mode_num = 10;
% disp = model.T_rbe(1:size(surf_mesh.node,1)*3,:) * U_rbe(:,mode_num);
% for i = 1:1000
% tic
% 
%     anime.fig_next(disp * sin(i * 1/100 * 2*pi),1)
%     title(sprintf('Frame# = %d, FPS = %.1f',i,1/toc))
%     
% %     pause(0.2-toc)
% end

%% Crag - Bampton
ROM_cb = myROM(model);
ROM_cb.CB(1:48,100); % 1-36 RBE / 37-48 Ref

model_cb = myFEM();
model_cb.set_MK(ROM_cb.M_reduction, ROM_cb.K_reduction);

    full_dof = model_cb.num_dof;
    BC_dof = 1:18;
    org_dof = 1:full_dof;
    org_dof(BC_dof) = [];

model_cb.compute_BC(BC_dof);

%% Time integration

NM = myNewmark(model_cb,1/100);
NM.compute_init();

iter_end = 3000;
f_vec = zeros(model_cb.num_dof,1);

tic
for i = 1:1000
%     f_vec(1) = 300 * (i/1000) * cos(i * 1/300 * 2*pi);
%     f_vec(2) = 300 * (i/1000) * sin(i * 1/300 * 2*pi);
    f_vec(7) = 800 * (i/1000) * cos(i * 1/300 * 2*pi);
    f_vec(8) = 800 * (i/1000) * sin(i * 1/300 * 2*pi);
    f_vec(13) = 1000 * (i/1000) * cos(i * 1/300 * 2*pi);
    f_vec(14) = 1000 * (i/1000) * sin(i * 1/300 * 2*pi);
    NM.compute_next(f_vec);
    
 if(rem(i,1000) == 0) 
 fprintf('- Compute Time Integration %.1f %% (%.1f sec)\n',i * 100 / iter_end, toc) 
 end
         
end

for i = 1001:iter_end
%     f_vec(1) = 300 * cos(i * 1/300 * 2*pi);
%     f_vec(2) = 300 * sin(i * 1/300 * 2*pi);
    f_vec(7) = 800 * cos(i * 1/300 * 2*pi);
    f_vec(8) = 800 * sin(i * 1/300 * 2*pi);
    f_vec(13) = 1000 * cos(i * 1/300 * 2*pi);
    f_vec(14) = 1000 * sin(i * 1/300 * 2*pi);
    NM.compute_next(f_vec);
    
 if(rem(i,1000) == 0) 
 fprintf('- Compute Time Integration %.1f %% (%.1f sec)\n',i * 100 / iter_end, toc) 
 end
         
end

%% Virtual Sensing

NM_inv = myNewmark(model_cb,1/100);
m_dof = [7,8,13,14];%19:27;
f_dof = [7,8,13,14];%[1,2,7,8,13,14];
NM_inv.compute_inv_init(m_dof,f_dof,0);

result = NM.data_log;
tic
for i = 1:iter_end

    NM_inv.compute_inv(result.U_log(m_dof,i));
    
 if(rem(i,1000) == 0) 
 fprintf('- Compute Time Integration %.1f %% (%.1f sec)\n',i * 100 / iter_end, toc) 
 end
         
end

result_2 = NM_inv.data_log;

% subplot(2,1,1)
% plot(result.F_log(1,:))
% title('Forward')
% subplot(2,1,2)
% plot(result_2.F_log(1,:))
% title('Inverse')

%%
fig = figure();
fig.Color = [1 1 1];

title_array = ["Bearing 1 / X", "Bearing 1 / Y", "Bearing 2 / X", "Bearing 2 / Y"];
for i = 1:numel(f_dof)
subplot(numel(f_dof)/2,2,i)
plot(result.t_log, result.F_log(f_dof(i),:),'r','LineWidth',1.5)
hold on
plot(result_2.t_log, result_2.Ff_log(i,:),'--k','LineWidth',1.5)
title(title_array(i))
xlabel("Force [N]",'FontSize',8); ylabel("Time [Sec]",'FontSize',8);
legend("Ref.", "Estimated")
end 


%%
fig = figure();
fig.Color = [1 1 1];

for i = 1:numel(m_dof)
subplot(numel(m_dof)/3,3,i)
plot(result.t_log, result.U_log(m_dof(i),:),'r','LineWidth',2)
hold on
plot(result_2.t_log, result_2.Um_log(i,:),'--k','LineWidth',2)
end 


%%
result = NM.data_log;
U_nonbc = zeros(full_dof,size(result.U_log,2));
U_nonbc(org_dof,:) = result.U_log;

i = 1; clear disp
for t = 1:10:size(U_nonbc,2)
    
    disp(:,i) = model.T_rbe(1:size(surf_mesh.node,1)*3,:) * ROM_cb.T_mat * U_nonbc(:,t);
    i = i+1
end

%% Plot Init.
anime = myPLOT(surf_mesh);
anime.fig_init([-120 120 -120 120 -250 20]); hold on; view([-150,30])
anime.fig_now.Position = [2000 100 700 600];
anime.fig_now.Color = [1 1 1];

%% Plot Loop
force_arrow =...
    quiver3(master_node_xyz(4:6,1),master_node_xyz(4:6,2),master_node_xyz(4:6,3),...
     result.F_log(1:6:18,1),result.F_log(2:6:18,1),result.F_log(3:6:18,1),...
     0.3,'r','LineWidth',2,'MaxHeadSize',2);
 
for i = 1:size(disp,2)
    tic

    delete(force_arrow)

    force_arrow =...
        quiver3(master_node_xyz(4:6,1),master_node_xyz(4:6,2),master_node_xyz(4:6,3),...
         result.F_log(1:6:18,10*i-9),result.F_log(2:6:18,10*i-9),result.F_log(3:6:18,10*i-9),...
         0.3,'r','LineWidth',2,'MaxHeadSize',2);

    anime.fig_next(disp(:,i),1000)
%     title(sprintf('Frame# = %d, FPS = %.1f',i,1/toc))
    title(sprintf('Frame # = %d',i))
end