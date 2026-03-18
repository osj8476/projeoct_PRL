%기본 파라미터
function params = parameters()
    params.m = 43;
    params.Ib=diag([0.41,2.1,2.1]);
    params.g=[0;0;-9.81];
    params.myu = 0.6;
    params.tau_max = 250;
    params.theta_weight = 1;
    params.z_weight = 50;
    params.yaw_rate_weight = 1;
    params.v_weight = 1;
    params.alpha = 1*10^-6;
    params.f_min = 0; % -> 솔버에서 해당 제약조건 f_min =10 일때 오류코드 -2반환, 0으로 강제 (이유는?)
    params.f_max = 666;
    params.l1 = 0.23; % thigh
    params.l2 = 0.23; % shank
    params.m_l1 = 0.6; 
    params.m_l2 = 0.4; % 전체질량 = 43kg, 다리 4개의 질량 = 4.3kg(10%) , 다리 하나의 질량 = 약 1kg 
    params.width_body = 0.25;
    params.length_body = 0.6;
    params.q_init = [0,0,0,0;0.8,0.8,0.8,0.8;-1.6,-1.6,-1.6,-1.6]; % hip roll, hip pitch , knee pitch 
end
params = parameters();

function cmd = command()
    cmd.v_x = 1;
    cmd.v_y = 1;
    cmd.yaw_rate = 0;
    cmd.yaw = atan2(cmd.v_y,cmd.v_x);
end
cmd = command();

% dcm of psi
function r = f_r(angle) 
    r = [cos(angle),-sin(angle),0;sin(angle),cos(angle),0;0,0,1];
end 

% skew matrix (외적곱 -> 행렬곱 연산)
function s = skew(r)
%skew matrix 외적 곱 -> 행렬곱으로 연산
s = [0,-r(3),r(2);r(3),0,-r(1);-r(2),r(1),0];
end

% A_C, B_C
function [A_c,B_c]=acbc(params,x,r_all)
  
    psi = x(3);
    rz=f_r(psi);
    rzT=rz';
    

    I_hat = rz*params.Ib*rz';
    I_hat_inv=inv(I_hat);

    %A_c
    zc= zeros(3);
    i3 = eye(3);
    A_c12 = [zc,zc,rzT,zc;zc,zc,zc,i3;zc,zc,zc,zc;zc,zc,zc,zc]; 
    A_c = blkdiag(A_c12,0);
    A_c(10:12, 13) = params.g; % -> v_dot = g_vector + 1/m(f_1+f_2+f_3+f_4)

    %B_c
    B_c=zeros(13,12);

    for i = 1:4
    ri_skew = skew(r_all(:,i));
    
    B_c(7:9,1+3*(i-1):3*i)=I_hat_inv*ri_skew;
    B_c(10:12,1+3*(i-1):3*i)=eye(3)/ params.m;
    end

end

%p_local : hip hinge <-> foot tip / r_body : center of mass (COM) <-> foot tip
function [p_local,r_body] = q2p(q,params)
    l1 = params.l1;
    l2 = params.l2;
    w =  params.width_body;
    l = params.length_body;
    
    %offset of FL,FR,HL,HR
    offset =[l/2,l/2,-l/2,-l/2;w/2,-w/2,w/2,-w/2;0,0,0,0];
    
    p_local = zeros(3,4);
    for i=1:4
        q1=q(1,i); q2=q(2,i); q3=q(3,i);
        x = l1*sin(q2) + l2*sin(q2+q3);
        y = l1*cos(q2)*sin(q1) + l2*cos(q2+q3)*sin(q1);
        z = -l1*cos(q2)*cos(q1) - l2*cos(q2+q3)*cos(q1);
        p_local(:,i)=[x;y;z];
    end    

    r_body= p_local + offset;
end

% r_all(world frame) 행렬
function r_all = update_r_i(q_init,x_future, x_current,params)
    [~ , r_body] = q2p(q_init, params); 
    
    p_com_current = x_current(4:6);    % 현재 COM 위치 [x; y; z]
    psi = x_current(3);        % 현재 Yaw 각
    Rz = f_r(psi);     % 현재 회전 행렬 (3x3)
    
    p_foot_world_fixed = p_com_current + (Rz * r_body);
    p_com_future = x_future(4:6); % 예측된 COM 위치 using mpc
    r_all = p_foot_world_fixed - p_com_future;
end

%discretize 이산화
function [A_d, B_d] = discretize_sys(A_c, B_c, Ts)
    
    sys_c = ss(A_c, B_c, eye(size(A_c)), zeros(size(A_c,1), size(B_c,2)));
    sys_d = c2d(sys_c, Ts);
    
    A_d = sys_d.A;
    B_d = sys_d.B;
end

function [A_qp, B_qp] = qp_matrices(A_d, B_d, k_h)
    
    n_row = 13;
    n_col = 12;
    A_qp = zeros(n_row*k_h,n_row);
    B_qp = zeros(n_row*k_h,n_col*k_h);
    
    for i=1:k_h
        A_qp(n_row*(i-1)+1:n_row*i,:) = A_d^(i);
        
        for j=1:i
            if j==i
                B_qp(n_row*(i-1)+1:n_row*i,n_col*(j-1)+1:n_col*j)=B_d;
            else        
                B_qp(n_row*(i-1)+1:n_row*i,n_col*(j-1)+1:n_col*j)=A_d^(i-j)*B_d;
            end
        end
    end
end

% Reference Trajectory
function x_ref = refer(horizon,x_current,cmd,pos_int,Ts)

n_row= 13;

x_ref = zeros(n_row*horizon,1);

% 현재값
pos_current = x_current(4:6);
yaw_heading = cmd.yaw;


% 목표값
v_des = [cmd.v_x;cmd.v_y;0];
yaw_rate_des = cmd.yaw_rate;
z_des = pos_int(3);

    for k=1:horizon
        kt = k*Ts;
        yaw_ref = yaw_heading;
        
        pos_ref = pos_current + v_des*kt;
        pos_ref(3)=z_des;

        x_ref_k=zeros(n_row,1);
        x_ref_k(3) = yaw_ref;
        x_ref_k(4:6) = pos_ref;
        x_ref_k(9) = yaw_rate_des;
        x_ref_k(10:12) = v_des;
        x_ref_k(13) = 1;

        x_ref(n_row*(k-1)+1:n_row*k) =  x_ref_k;
    end
end

function [C,D,L,K] = matrix_of_weights(params,horizon)
    n_row = 13;
    n_col = 12;
    myu=params.myu;

    Q = diag([ params.theta_weight, params.theta_weight, params.theta_weight,0,0,params.z_weight,0,0, ...
params.yaw_rate_weight,params.v_weight,params.v_weight,params.v_weight,0]);
    R = eye(n_col)*params.alpha;
    
    C = zeros(6*4*horizon,3*4*horizon);
    D = zeros(6*4*horizon,1);
    L = zeros(n_row*horizon,n_row*horizon);
    K = zeros(n_col*horizon,n_col*horizon);
   
    c_i = [1,0,-myu;-1,0,-myu;0,1,-myu;0,-1,-myu;0,0,1;0,0,-1];
    d_i = [0;0;0;0;params.f_max;-params.f_min];
    for i=1:horizon

        L(n_row*(i-1)+1:n_row*i,n_row*(i-1)+1:n_row*i) = Q;
        K(n_col*(i-1)+1:n_col*i,n_col*(i-1)+1:n_col*i) = R;
        
        for n=1:4
            row_ik = 24*(i-1) + 6*(n-1) +1;  
            col_ik = 12*(i-1) + 3*(n-1) +1;
            
            C(row_ik:row_ik+5,col_ik:col_ik+2)=c_i;
            D(row_ik:row_ik+5,1)=d_i;
        end
    end
end

%%
%parameter
Ts=0.05;%sampling time Ts
horizon = 10; % horizon k
T_sim =5; % simulation time
n_sim = round(T_sim/Ts);%No of step 


% 상태 x 초기값 정의
euler_angle_int = [0;0;cmd.yaw];
pos_int = [0;0;0.32]; % l1,l2,q2,q3 초기값으로 계산 z_int : 32cm
omega_int = [0;0;0];
vel_int = [0;0;0];
gravity_const = 1;

x=[euler_angle_int;pos_int;omega_int;vel_int;gravity_const]; % 13x1 x matrix 

% simulation loop

x_history = zeros(13,n_sim);
u_history = zeros(12,n_sim);

for t=1:n_sim

    x_ref = refer(horizon,x,cmd,pos_int,Ts);
    
    x_future = x; % 근사화
    
    r_all = update_r_i(params.q_init,x_future, x,params); % 고민해야할 부분->x_future = x로 가정 (짧은 duration). 
    
    [A_c,B_c]=acbc(params,x,r_all);
    [A_d, B_d] = discretize_sys(A_c, B_c, Ts);
    [A_qp, B_qp] = qp_matrices(A_d, B_d,horizon);
    [C,D,L,K] = matrix_of_weights(params,horizon);
    
    H=2*(B_qp'*L*B_qp+K);
    H = (H + H') / 2;
    g=2*B_qp'*L*(A_qp*x-x_ref);
    
    %등식제약조건
    %A_eq *U = B_eq
    %특정 지면반력 ->0으로 강제 (swing)
    A_eq = zeros(12*horizon,12*horizon);
    B_eq = zeros(12*horizon,1);

    for k=1:horizon
        index_sw = t + k ; 
        standard_sw = mod(floor(index_sw/4),2);

        if standard_sw == 0
          gait_matrix = [1,0,0,1];
        else
          gait_matrix = [0,1,1,0];
        end
        
        for j=1:4
            if gait_matrix(j)==0
                row_g = 12*(k-1) + 3*(j-1) +1;
                col_g = 12*(k-1) + 3*(j-1) +1;
                A_eq(row_g:row_g+2,col_g:col_g+2) = eye(3);
            end
        end
    end

    empty_rows = all(A_eq == 0, 2);
    A_eq(empty_rows, :) = [];
    B_eq(empty_rows, :) = [];



    options = optimoptions('quadprog', 'Algorithm', 'interior-point-convex','Display', 'off');
    
    
    [u_star, fval, exitflag] = quadprog(H, g, C, D, A_eq, B_eq, [], [], [], options);
    
    if exitflag > 0 && ~isempty(u_star)
        u_current = u_star(1:12);
    else
        % 솔버 실패 시 경고 메시지 출력 및 안전한 기본값(0) 대입
        warning('QP 솔버가 해를 찾지 못했습니다. Exitflag: %d', exitflag);
        u_current = zeros(12, 1); 
    end
    
   

    x = A_d*x+B_d*u_current; % x업데이트 x(i+1)

    u_history(:,t)= u_current ;
    x_history(:,t)= x ;

end

%% 시뮬레이션 결과 시각화

t_vec = (0:n_sim-1) * Ts;
figure('Name', 'MPC Quadruped Control Results', 'Color', 'w', 'Position', [100, 100, 1000, 800]);

% --- 1. 3D Trajectory ---
subplot(3, 2, 1);
plot3(x_history(4,:), x_history(5,:), x_history(6,:), 'b', 'LineWidth', 2);
grid on; hold on;
plot3(x_history(4,1), x_history(5,1), x_history(6,1), 'ro', 'MarkerFaceColor', 'r');
plot3(x_history(4,end), x_history(5,end), x_history(6,end), 'gx', 'MarkerSize', 10, 'LineWidth', 2); 
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('3D World Trajectory');
axis equal;

% --- 2. Vertical Forces (Fz)  ---
subplot(3, 2, 2);

plot(t_vec, u_history(3:3:12, :), 'LineWidth', 1.2);
grid on;
title('Vertical Ground Reaction Forces (F_z)');
legend('FL', 'FR', 'HL', 'HR');
xlabel('Time [s]'); ylabel('Force [N]');
ylim([0, params.f_max * 1.1]);

% --- 3. Body Height (Z) ---
subplot(3, 2, 3);
plot(t_vec, x_history(6,:), 'r', 'LineWidth', 1.5);
line([t_vec(1) t_vec(end)], [0.32 0.32], 'Color', 'k', 'LineStyle', '--'); 
grid on;
xlabel('Time [s]'); ylabel('Z [m]');
title('Body Height (Z)');
legend('Actual', 'Target');
ylim([0, 0.32 * 1.5]);

% --- 4. Yaw Angle ---
subplot(3, 2, 4);
plot(t_vec, rad2deg(x_history(3,:)), 'k', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]'); ylabel('Yaw [deg]');
title('Body Yaw Angle');

% --- 5. Forward/Side Velocity (Vx, Vy) ---
subplot(3, 2, 5);
plot(t_vec, x_history(10,:), 'b', 'LineWidth', 1.5,'DisplayName', 'Actual Vx'); hold on;
plot(t_vec, x_history(11,:), 'g', 'LineWidth', 1.5,'DisplayName', 'Actual Vy');
grid on;
hold on;
plot(t_vec, repmat(cmd.v_x, 1, n_sim), 'r--', 'DisplayName', 'Desired Vx');
xlabel('Time [s]'); ylabel('Vel [m/s]');
legend;
title('Body Linear Velocity Tracking Performance');
ylim([0, 1.2]);

sgtitle('MIT Cheetah 3 Convex MPC Simulation Analysis', 'FontSize', 14, 'FontWeight', 'bold');



%% --- 3D Animation Implementation (Updated: 3D Box Body) ---
fig = figure('Name', 'Quadruped MPC - 3D Box Body', 'Color', 'w', 'Position', [100, 100, 800, 600]);
view(3); grid on; axis equal;
hold on;

% 1. 시간 표시 
h_time = text(0.05, 0.95, 'Time: 0.00 s', 'Units', 'normalized', ...
              'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');

% 2. 로봇 몸체 (직육면체)
wb = params.width_body;
lb = params.length_body;
hb = 0.15; % 직육면체의 높이 

% 직육면체의 8개 꼭짓점 정의 (로컬 좌표계, COM 기준)
body_vertices = [-lb/2, -wb/2, -hb/2; 
                  lb/2, -wb/2, -hb/2;  
                  lb/2,  wb/2, -hb/2;  
                 -lb/2,  wb/2, -hb/2;  
                 -lb/2, -wb/2,  hb/2;  
                  lb/2, -wb/2,  hb/2;  
                  lb/2,  wb/2,  hb/2;  
                 -lb/2,  wb/2,  hb/2]; 

% 직육면체의 6개 면 정의
body_faces = [1 2 3 4;  % 밑면
              5 6 7 8;  % 윗면
              1 2 6 5;  % 우면
              2 3 7 6;  % 전면
              3 4 8 7;  % 좌면
              4 1 5 8]; % 후면

% 면별 색상 정의 (기본 회색, 전면만 주황색)
face_colors = repmat([0.5 0.5 0.5], 6, 1); % 모든 면 회색 초기화
face_colors(4, :) = [1 0.5 0];           % 전면만 주황색으로 설정

% patch 함수로 3D 직육면체 생성
h_body_3d = patch('Vertices', body_vertices, 'Faces', body_faces, ...
                  'FaceVertexCData', face_colors, ... 
                  'FaceColor', 'flat', ...           
                  'FaceAlpha', 0.8, ...              
                  'EdgeColor', 'k', 'LineWidth', 1);

h_com = plot3(0,0,0, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8); % COM 

% 3. 다리 및 지면 반력 시각화 
leg_colors = {'r', 'g', 'b', 'm'}; 
h_legs = cell(1, 4);
h_force = cell(1, 4);
for i = 1:4
    h_legs{i} = plot3([0 0], [0 0], [0 0], leg_colors{i}, 'LineWidth', 3);
    h_force{i} = quiver3(0,0,0,0,0,0, 0, 'Color', leg_colors{i}, 'MaxHeadSize', 0.4, 'LineWidth', 1.5);
end

% 4. 격자판(Ground Grid) 시각화
[gx, gy] = meshgrid(-10:0.5:10, -10:0.5:10);
gz = zeros(size(gx));
plot3(gx, gy, gz, 'Color', [0.8 0.8 0.8], 'LineWidth', 0.1, 'HandleVisibility', 'off');

% --- 속도 설정 ---
slow_factor = 3; % 3배 느리게
frame_delay = Ts * slow_factor; 

% --- Animation Loop ---
for t = 1:1:n_sim
    tic;
    
    current_time = (t-1) * Ts;
    set(h_time, 'String', sprintf('Time: %.2f s', current_time));
    
    % 현재 상태 및 회전 행렬
    psi = x_history(3, t);
    pos = x_history(4:6, t);
    Rz = f_r(psi);
    
    % 핵심: 3D 직육면체 꼭짓점 좌표 업데이트 (회전 및 평행이동 반영)
    % 로컬 꼭짓점 좌표를 월드 좌표계로 변환
    body_v_world = (Rz * body_vertices')' + pos';
    set(h_body_3d, 'Vertices', body_v_world); 
    
    set(h_com, 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3));
    
    % 다리 및 힘 업데이트
    [~, r_body] = q2p(params.q_init, params);
    for i = 1:4
        hip_world = pos + Rz * (r_body(:,i));
        fz = u_history(3*i, t); 
        scale = 0.0022; 
        
        if fz > 1.0 
            foot_pos = [hip_world(1:2); 0]; 
            set(h_legs{i}, 'XData', [hip_world(1), foot_pos(1)], ...
                           'YData', [hip_world(2), foot_pos(2)], ...
                           'ZData', [hip_world(3), foot_pos(3)], 'Visible', 'on');
            set(h_force{i}, 'XData', foot_pos(1), 'YData', foot_pos(2), 'ZData', foot_pos(3), ...
                            'UData', 0, 'VData', 0, 'WData', fz * scale, 'Visible', 'on');
        else
            set(h_legs{i}, 'Visible', 'off'); 
            set(h_force{i}, 'Visible', 'off');
        end
    end
    
    % 로봇 중심 시점
    view_range = 1.2; 
    axis([pos(1)-view_range, pos(1)+view_range, ...
          pos(2)-view_range, pos(2)+view_range, ...
          0, 1]); 
    
    drawnow;
    
    loop_time = toc;
    pause(max(0.01, frame_delay - loop_time)); 
end


%foot jacobian matrix ,j 
% function jacobian_i = f_jacobian_i(parmas,q_i)
% 
%     l1 = parmas.l1;
%     l2 = parmas.l2;
% 
%     q1=q_i(1); q2=q_i(2); q3=q_i(3);
%     s1 = sin(q1); c1 = cos(q1);
%     s2 = sin(q2); c2 = cos(q2);
%     s23 = sin(q2+q3); c23 = cos(q2+q3);
% 
%     jacobian_i = zeros(3,3);
% 
%     jacobian_i(1,1) = 0;
%     jacobian_i(1,2) = l1*c2 + l2*c23;
%     jacobian_i(1,3) = l2*c23;
% 
%     % 행 2: dy/dq1, dy/dq2, dy/dq3
%     jacobian_i(2,1) = l1*c2*c1 + l2*c23*c1;
%     jacobian_i(2,2) = -l1*s2*s1 - l2*s23*s1;
%     jacobian_i(2,3) = -l2*s23*s1;
% 
%     % 행 3: dz/dq1, dz/dq2, dz/dq3
%     jacobian_i(3,1) = l1*c2*s1 + l2*c23*s1;
%     jacobian_i(3,2) = l1*s2*c1 + l2*s23*c1;
%     jacobian_i(3,3) = l2*s23*c1;
% end

%lamda matrix 
% function lamda_i=f_lamda_i(params,q_i)
%     l1 = params.l1; l2 = params.l2;
%     m_l1 = params.m_l1; m_l2 = params.m_l2; 
% 
%     I_knee = m_l2*(l2/2)^2;
%     I_hip_pitch = m_l1*(l1/2)^2 + m_l2*(l2/2)^2;
%     I_hip_roll = I_hip_pitch;
% 
%     M = diag([I_hip_pitch,I_hip_roll,I_knee]);
% 
%     j_i = f_jacobian_i(params,q_i);
%     lamda_inv = j_i*(M\j_i');
%     lamda_i=inv(lamda_inv+eye(3)*1e-6);
% end