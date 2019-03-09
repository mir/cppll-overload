% parameters
omega_free = 0;
T_ref = 10^-3;
R = 1000;
C = 10^-6;
K_vco = 1000;
I_p = 10^-3;

% recalculated values
tau_2N = R*C/T_ref
K_N = I_p*R*K_vco*T_ref
F_N = 1/(2*pi)*sqrt(K_N/tau_2N)
dzeta = sqrt(K_N*tau_2N)/2

% initial data
v_1 = 0;
tau_1 = -0.1*T_ref;
if (tau_1 < -T_ref)
    errordlg('Impossible tau_1. Can not be lower than -T_ref.');
end

% number of steps (tau_k) to simulate
max_step = 10000;

% initialize PFD output with initial data
pfd_output = zeros((max_step-1)*4,2);
pfd_output(1,:) = [0 0];
if (tau_1 >= 0)
    pfd_output(2,:) = [0 I_p];
    pfd_output(3,:) = [tau_1 I_p];
    pfd_output(4,:) = [tau_1 0];
    t_k_middle = tau_1;
    initial_vco_phase = 1 - (...
           (omega_free + K_vco*v_1 + K_vco*I_p*R - K_vco*tau_1*I_p/C)*tau_1 ...
           + K_vco*I_p/2/C*tau_1^2 ...
        );        
    % overload
    if (-K_vco*I_p/C*tau_1 + K_vco*v_1 < 0)
        if (-K_vco*I_p/C*tau_1 + K_vco*v_1 +K_vco*I_p*R < 0)
           if (omega_free + K_vco*v_1 + K_vco*I_p*R < 0)
                errordlg('Impossible initial condition. v_1 is too small.');
           end
           initial_vco_phase = 1 - (...
                (omega_free + K_vco*v_1 + K_vco*I_p*R)^2/2/(K_vco*I_p/C) ...
            );
        else
           initial_vco_phase = 1 - (...
                tau_1^2*K_vco*I_p/2/C+ ...
                tau_1*(K_vco*I_p*R + K_vco*v_1 + omega_free - tau_1*K_vco*I_p/C)...
            );
        end
    end
    if (initial_vco_phase < 0)
        errordlg('Impossible initial condition. v_1 or tau_1 is too big.');
    end
    initial_ref_phase = 0;
else
    pfd_output(2,:) = [0 -I_p];
    pfd_output(3,:) = [-tau_1 -I_p];
    pfd_output(4,:) = [-tau_1 0];
    t_k_middle = -tau_1;
    initial_vco_phase = 0;
    initial_ref_phase = 1 + tau_1/T_ref;
end
initial_filter_state = v_1 - I_p*tau_1/C;
index = 4;


tau_v = zeros(max_step,2);
tau_v(1,:) = [tau_1 v_1];
tau_k = tau_1;
v_k = v_1;
for step = 2:(max_step - 1)  
    [tau_k1,v_k1,tau_k_zero] = righthand(tau_k,v_k ,...
                                K_vco, T_ref, I_p, C, R, omega_free);

    %check for VCO overload
    if ((tau_k > 0 ...
            && (v_k+omega_free/K_vco - I_p/C*tau_k) < 0)...
        ||...
        (tau_k < 0 ...
            && v_k+omega_free/K_vco - I_p*R < 0))       
        [tau_k_o,v_k_o,tau_k_zero] = righthand_overload(tau_k,v_k ,...
                                             tau_k1,v_k1,...
                                K_vco, T_ref, I_p, C, R, omega_free);
        tau_k1 = tau_k_o;
        v_k1 = v_k_o;
    end
                            
    tau_v(step,:) = [tau_k1 v_k1];
    t_k1 = t_k_middle + tau_k_zero;
    index = index + 1;
    pfd_output(index,:) = [t_k1 0];
    
    if (tau_k1 ~= 0)
        index = index + 1;
        pfd_output(index,:) = [t_k1 I_p*sign(tau_k1)];
        
        t_k1_middle = t_k1 + abs(tau_k1);
        index = index + 1;
        pfd_output(index,:) = [t_k1_middle I_p*sign(tau_k1)];
        index = index + 1;
        pfd_output(index,:) = [t_k1_middle 0];
    end
    
    t_k = t_k1;
    t_k_middle = t_k1_middle;
    tau_k = tau_k1;
    v_k = v_k1;
end
[tau_k1,v_k1,tau_k_zero] = righthand(tau_k,v_k ,...
                                K_vco, T_ref, I_p, C, R, omega_free);
tau_v(max_step,:) = [tau_k1,v_k1];

% truncate trailing zeros
last_non_zero = find(pfd_output(:,1),1,'last');
pfd_output = pfd_output(1:last_non_zero,:);

plot(pfd_output(:,1),pfd_output(:,2));
ylim([-1.1*I_p 1.1*I_p]);
% paemel_simulation;
sim('pfd_overload_simulink');
                            