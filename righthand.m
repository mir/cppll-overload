function [tau_k1,v_k1,tau_k_zero] = righthand( ...
    tau_k,v_k,...
    K_vco, T_ref, I_p, C, R, omega_free)
%righthandside Corrected 
if(tau_k >= 0)
    c = (T_ref - rem(tau_k,T_ref))*(omega_free+K_vco*v_k)-1;
    if (c <= 0)
%         tau(k+1) > 0, case 1) 
        a = K_vco*I_p/(2*C);
        b = omega_free + K_vco*v_k + K_vco*I_p*R;
        tau_k1 = (-b + sqrt(b^2 - 4*a*c))/(2*a);
        tau_k_zero = T_ref - rem(tau_k,T_ref);
    else
%         tau(k+1) < 0, case 2)
        tau_k1 = 1/(omega_free + K_vco*v_k) - T_ref + rem(tau_k,T_ref);
        tau_k_zero = 1/(omega_free + K_vco*v_k);
    end
else
    lk = -tau_k;
    S_lk = (K_vco*v_k-I_p*R*K_vco+omega_free)*lk+(K_vco*I_p*lk^2)/(2*C);
    S_la = rem(S_lk,1);
    S_lb = 1-S_la;
    lb = S_lb/(K_vco*v_k + omega_free);
    if lb <= T_ref
        %                 tau(k+1) < 0
        l_k1 = T_ref - lb;
        tau_k1 = -l_k1;
        tau_k_zero = lb;
    else
        %                 tau(k+1) >= 0
        S_Tref = T_ref*(K_vco*v_k +omega_free);
        c = S_la + S_Tref -1;
        b = omega_free +K_vco*v_k +K_vco*I_p*R;
        a = K_vco*I_p/(2*C);
        tau_k1 = (-b+sqrt(b^2-4*a*c))/(2*a);
        tau_k_zero = T_ref;
    end
end
v_k1 = v_k + tau_k1*I_p/C;
end