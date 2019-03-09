function [tau_k1,v_k1,tau_k_zero] = paemel_righthand( ...
    tau_k,v_k,...
    K_vco, T_ref, I_p, C, R)
%righthandside Paemel 
if(tau_k >= 0)
    % case 1
    a = I_p/(2*C);
    b = v_k + I_p*R;
    c = (T_ref - tau_k)*v_k-1/K_vco;
    tau_k1 = (-b + sqrt(b^2 - 4*a*c))/(2*a);
    tau_k_zero = T_ref - tau_k;
    if (tau_k1 < 0) 
        % case 3 Paemel
        tau_k1 = 1/(K_vco*v_k) - T_ref + tau_k;
        tau_k_zero = 1/(K_vco*v_k);
    end
else
    % case 2 Paemel
    tau_k1 = (...
        1/K_vco  - I_p*R*tau_k - I_p*tau_k^2/(2*C)...
        )/v_k - T_ref + tau_k;
    
    lk = -tau_k;
    S_lk = (K_vco*v_k-I_p*R*K_vco)*lk+(K_vco*I_p*lk^2)/(2*C);
    tau_k_zero = (1-rem(S_lk,1))/(K_vco*v_k);
    
    if( tau_k1 > 0)
        % case 4 Paemel
        tau_k1 = (...
            -I_p*R -v_k + sqrt((I_p*R + v_k)^2 + (2*I_p/C)*v_k*tau_k1) ...
            )/(I_p/C);
            
        tau_k_zero = T_ref;
    end
end
v_k1 = v_k + tau_k1*I_p/C;
end