function [tau_k_o,v_k_o,tau_k_zero] = righthand_overload( ...
    tau_k,v_k,...
    tau_k_1,v_k_1,...
    K_vco, T_ref, I_p, C, R, omega_free)
%righthandside overload 
root2 = @(a,b,c) ((-b + sqrt(b^2 - 4*a*c))/2/a);

omega_vco = omega_free + K_vco*v_k;
a = K_vco*I_p/2/C;
b = omega_vco + K_vco*I_p*R;

if(tau_k < 0)
    l_x = min([-C/I_p*(v_k + omega_free/K_vco - I_p*R) -tau_k]);
    S = K_vco*(tau_k + l_x)^2*I_p/2/C;
    Sla = rem(S,1);
    if (omega_vco > 0)        
        if (tau_k_1 < 0)
            % O1            
            l_b = (1-Sla)/omega_vco;
            tau_k_o = -(T_ref - l_b);
            
            tau_k_zero = l_b;
        else
            % O2
            Sref = T_ref*omega_vco;
            d = Sla + Sref - 1;
            tau_k_o = root2(a,b,d);
            
            tau_k_zero = T_ref;
        end
    else
        if (v_k + omega_free/K_vco + I_p*R < 0)
            % O3
            l_b_0 = C/I_p*(-omega_free/K_vco - I_p*R - v_k);
            l_b_plus = sqrt((1 - Sla)*2*C/K_vco/I_p);
            tau_k_o = l_b_0 + l_b_plus;
            
            tau_k_zero = T_ref;
        else
            % O4            
            d_0 = Sla - 1;            
            tau_k_o = root2(a,b,d_0);
            
            tau_k_zero = T_ref;
        end
    end
else
    if (omega_vco <= 0)
        if (tau_k_1 > 0)
            % O5            
            tau_k_o = root2(a,b,-1);
            
            tau_k_zero = T_ref - rem(tau_k,T_ref);
        end
    else
        if (tau_k_1 > 0)
            % O6 equal to case 1) no overload
            c = (T_ref - rem(tau_k,T_ref))*(omega_free+K_vco*v_k)-1;
            tau_k_o = root2(a,b,c);            
            
            tau_k_zero = T_ref - rem(tau_k,T_ref);
        else
            % O7 equal to case 2) no overload
            tau_k_o = 1/(omega_free + K_vco*v_k) - T_ref + rem(tau_k,T_ref);
        
            tau_k_zero = T_ref - rem(tau_k,T_ref) + tau_k_o;
        end
    end
end
v_k_o = v_k + tau_k_o*I_p/C;
end