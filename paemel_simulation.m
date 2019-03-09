% check that omega_free is zero
if (omega_free > 0)
    errordlg('Please set omega_free to zero');
end

% initialize PFD output with initial data
paemel_output = zeros((max_step-1)*4,2);
paemel_output(1,:) = [0 0];
if (tau_1 >= 0)
    paemel_output(2,:) = [0 I_p];
    paemel_output(3,:) = [tau_1 I_p];
    paemel_output(4,:) = [tau_1 0];
    t_k_middle = tau_1;
    initial_vco_phase = 1 - (...
           (K_vco*v_1)*tau_1 ...
           + K_vco*I_p/2/C*tau_1^2 ...
           - K_vco*I_p/C*tau_1^2 ...
        );        
    initial_ref_phase = 0;
else
    paemel_output(2,:) = [0 -I_p];
    paemel_output(3,:) = [-tau_1 -I_p];
    paemel_output(4,:) = [-tau_1 0];
    t_k_middle = -tau_1;
    initial_vco_phase = 0;
    initial_ref_phase = 1 + tau_1/T_ref;
end
initial_filter_state = v_1 - I_p*tau_1/C;
index = 4;


tau_v_paemel = zeros(max_step,2);
tau_v_paemel(1,:) = [tau_1 v_1];
tau_k = tau_1;
v_k = v_1;
for step = 2:(max_step - 1)  
    [tau_k1,v_k1,tau_k_zero] = paemel_righthand(tau_k,v_k ,...
                                K_vco, T_ref, I_p, C, R);

    % Compute case 6
    if (tau_k1 < -T_ref)
        if (step <= 2)
            errordlg('ERROR in Paemel algo case 6), no v(-1)');
            paemel_output = zeros((max_step-1)*4,2);
            break;
        end        
        v_n_1 = tau_v_paemel(step-2,2);         
        t_sum = 0;
        while (t_sum < abs(tau_k))
            t_n = (v_n_1 - I_p*R - sqrt((v_n_1 - I_p*R)^2 - 2*I_p/C/K_vco)...
              )/(I_p/C);
            t_sum = t_sum + t_n;
            v_n_prev = v_n_1;
            v_n_1 = v_n_1 - I_p/C*t_n;            
        end 
        t_sum = t_sum - t_n;
        t_a = -tau_k - t_sum;
        t_b = (1/K_vco - t_a*(v_n_prev - I_p*R) + I_p/C*t_a^2/2)/v_k;
        tau_k1 = t_b - T_ref;
        
        v_k1 = v_k + tau_k1*I_p/C;
        tau_k_zero = t_b;
    end                            
                            
    tau_v_paemel(step,:) = [tau_k1 v_k1];
    t_k1 = t_k_middle + tau_k_zero;
    index = index + 1;
    paemel_output(index,:) = [t_k1 0];
    
    if (tau_k1 ~= 0)
        index = index + 1;
        paemel_output(index,:) = [t_k1 I_p*sign(tau_k1)];
        
        t_k1_middle = t_k1 + abs(tau_k1);
        index = index + 1;
        paemel_output(index,:) = [t_k1_middle I_p*sign(tau_k1)];
        index = index + 1;
        paemel_output(index,:) = [t_k1_middle 0];
    end
    
    % case 5
    if (tau_k1 >= T_ref)
        tau_k1 = rem(tau_k1,T_ref);
    end
    
    t_k = t_k1;
    t_k_middle = t_k1_middle;
    tau_k = tau_k1;
    v_k = v_k1;
end
[tau_k1,v_k1,tau_k_zero] = righthand(tau_k,v_k ,...
                                K_vco, T_ref, I_p, C, R, omega_free);
tau_v_paemel(max_step,:) = [tau_k1,v_k1];

% truncate trailing zeros
last_non_zero = find(paemel_output(:,1),1,'last');
if (last_non_zero > 1)
    paemel_output = paemel_output(1:last_non_zero,:);
end
% plot(pfd_output(:,1),pfd_output(:,2));
% ylim([-1.1*I_p 1.1*I_p]);
                            