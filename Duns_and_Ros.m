%Inputs
Qg = 0.1992  ;                                      % Gas flow rate             Unit = ft^3/sec 
QL = 0.1290  ;                                      % Liquid flow rate          Unit = ft^3/sec
miu_L = 18   ;                                      % Liquid viscosity          Unit = cp
rho_L = 56.6   ;                                    % Liquid density            Unit = lb/ft^3
sigma = 30   ;                                      % Interfacial tension       Unit = dyne/cm
d = 0.249    ;                                      % Tubing ID                 Unit = ft
eps = 0.0006   ;                                    % Relative roughness        eps = epsilon/d
L = 10000   ;                                       % True Length               Unit = ft
P_tf = 250  ;                                       % Top flow pressure         Unit = psi
T_tf = 100  ;                                       % Top flow temperature      Unit = F
T_wf = 200  ;                                       % Well flow temperature     Unit = F
MW = 21     ;                                       % Molecular Weight          Unit = lb/lb mole
T_c = 407   ;                                       % Critical Temoerature      Unit = R
P_c = 607   ;                                       % Critical Pressure         Unit = psi
no_i = 10000   ;                                    % Number of element

A = pi * (d^2) / 4   ;                              % Tubing area               Unit = ft^2
V_sG = Qg / A        ;                         % Superfacial liquid velocity    Unit = ft/sec
V_sL = QL / A        ;                         % Superfacial gas velocity       Unit = ft/sec
V_m = V_sL + V_sG    ;                              % mixture velocity          Unit = ft/sec
N_Lv = 1.938 * V_sL * (rho_L / sigma)^(1/4) ;       % Liquid Velocity number                     
N_gv = 1.938 * V_sG * (rho_L / sigma)^(1/4) ;       % Gas Velocity number
N_d = 120.872 * d * (rho_L / sigma)^(1/4)   ;       % Pipe diameter number
N_L = 0.15726 * miu_L * (1 / (rho_L * sigma))^(1/4) ;  % Liquid viscosity number
Re_L = (1488 * rho_L * V_sL * d) / miu_L      ;        % Liquid Reynolds number
L_s = 50 + 36 * N_Lv ;
L_m = 75 + 84 * (N_Lv)^(3/4) ;
L_1 = L1 (N_d)         ;
L_2 = L2 (N_d)         ;

dL = L / no_i  ;
P_list = [P_tf]  ;                                  % Pressure list
L_list = [0]    ;                                   % Length list
l = 0           ;                                   % initial length
P_1 = P_tf      ;                                   % initial pressure

for i = 1 : no_i
    T1 = g(l)  ;                                    % Temperature of up of element
    l = l + dL ;
    L_list = [L_list , l]  ;
    T2 = g(l)   ;                                   % Temperature of down of element
    T_avg = (T1 + T2) / 2  ;
    
    Delta_P2 = 0         ;
    Delta_P = 5   ;                                % Unit = psi
    
    while abs(Delta_P - Delta_P2) > 0.001
        Delta_P2 = Delta_P   ;
        P_2 = P_1 + Delta_P  ;
        P_avg = (P_1 + P_2) / 2   ;
        dP_dl = flow_regime(N_gv,P_avg,T_avg)   ;
        Delta_P = dP_dl * dL    ;
    end 
    P_1 = P_1 + Delta_P  ;
    P_list = [P_list , P_1]  ;
end

plot(P_list,L_list)
xlabel('Pressure')
ylabel('Measured depth')
set(gca,'XAxisLocation','top','YAxisLocation','left', 'YDir','reverse')

function [dp_dL] = flow_regime (Ngv,P,T)
Qg = .1992  ;                                       % Gas flow rate             Unit = ft^3/sec 
QL = .1290  ;                                       % Liquid flow rate          Unit = ft^3/sec
miu_L = 18  ;                                       % Liquid viscosity          Unit = cp
MW = 21  ;
rho_g = SRK (P,T)  ;
miu_g = miu(T , MW , rho_g)   ;
rho_L = 56.6   ;                                    % Liquid density            Unit = lb/ft^3
sigma = 30  ;                                       % Interfacial tension       Unit = dyne/cm
d = .249    ;                                       % Tubing ID                 Unit = ft
eps = 0.0006   ;                                    % Relative roughness        eps = epsilon/d

A = pi * (d^2) / 4   ;                              % Tubing area               Unit = ft^2
V_sG = Qg / A        ;                         % Superfacial liquid velocity    Unit = ft/sec
V_sL = QL / A        ;                         % Superfacial gas velocity       Unit = ft/sec
V_m = V_sL + V_sG    ;                              % mixture velocity          Unit = ft/sec
N_Lv = 1.938 * V_sL * (rho_L / sigma)^(1/4) ;       % Liquid Velocity number                     
N_gv = 1.938 * V_sG * (rho_L / sigma)^(1/4) ;       % Gas Velocity number
N_d = 120.872 * d * (rho_L / sigma)^(1/4)   ;       % Pipe diameter number
N_L = 0.15726 * miu_L * (1 / (rho_L * sigma))^(1/4) ;  % Liquid viscosity number
Re_L = (1488 * rho_L * V_sL * d) / miu_L      ;        % Liquid Reynolds number
L_s = 50 + 36 * N_Lv ;
L_m = 75 + 84 * (N_Lv)^(3/4) ;
L_1 = L1 (N_d)         ;
L_2 = L2 (N_d)         ;

%Bubble
if (0 <= Ngv) && (Ngv <= L_1 + L_2 * N_Lv)
   
    %elevation
    F_1 = F1 (N_L)         ;
    F_2 = F2 (N_L)         ;
    F_3 = F3 (N_L)         ;
    F_4 = F4 (N_L)         ;
    S = F_1 + F_2 * N_Lv + (F_3 - (F_4) / N_d) * (Ngv / (1 + N_Lv))^2  ;
    V_s = (S / 1.938) * (sigma / rho_L)^(1/4)   ;
    H_L = (V_s - V_m + sqrt((V_m - V_s)^2 + 4 * V_s * V_sL)) / (2 * V_s)  ;
    H_g = 1 - H_L   ;
    dp_ele = (rho_L * H_L) + (rho_g * H_g)   ;
    
    %friction
    f_1 = Moody (Re_L,eps)   ;
    f_2 = f2 ((f_1 * V_sg * N_d^(2/3)) / V_sL)   ;
    f_3 = 1 + f_1 * sqrt(V_sG / (50 * V_sL))   ;
    f_tp = (f_1 * f_2) / f_3   ;
    dp_friction = (f_tp * rho_L * V_sL * V_m) / (2 * 32.2 * d)   ;
    
    %acc negligible
    
    dp_dL = (dp_ele + dp_friction) / 144  ;
end

%Slug
if (L_1 + L_2 * N_Lv < Ngv) &&  (Ngv <= L_s)
    
    %elevation
    F_5 = F5 (N_L)         ;
    F_6 = F6 (N_L)         ;
    F_7 = F7 (N_L)         ;
    S = (1 + F_5) * ((Ngv^(0.982) + 0.029 * N_d + F_6) / (1 + F_7 * N_Lv)^2)   ;
    V_s = (S / 1.938) * (sigma / rho_L)^(1/4)   ;
    H_L = (V_s - V_m + sqrt((V_m - V_s)^2 + 4 * V_s * V_sL)) / (2 * V_s)  ;
    H_g = 1 - H_L   ;
    dp_ele = (rho_L * H_L) + (rho_g * H_g)   ;
    
    %friction
    f_1 = Moody (Re_L,eps)   ;
    f_2 = f2 ((f_1 * V_sG * N_d^(2/3)) / V_sL)   ;
    f_3 = 1 + f_1 * sqrt(V_sG / (50 * V_sL))   ;
    f_tp = (f_1 * f_2) / f_3   ;
    dp_friction = (f_tp * rho_L * V_sL * V_m) / (2 * 32.2 * d)   ;
    
    %acc negligible
    
    dp_dL = (dp_ele + dp_friction) / 144  ;
end

%Annular
if (L_m < Ngv)
    
    %elevation
    Lambda_L = QL / (QL + Qg)  ;
    Lambda_g = 1 - Lambda_L    ;
    dp_ele = (rho_L *  Lambda_L) + (rho_g * Lambda_g)   ;
    rho_n = dp_ele  ;
    
    %friction
    e = eps * d ;                                          % e = epsilon
    N_we = 452.84 * (rho_g * (V_sG)^2 * e) / sigma   ;
    N_miu = 6.23*10^(-3) * (((miu_L)^2) / (rho_L * sigma * e))  ;
    iter = 0 ;
    d_avg = d ;
    V_sG_avg = V_sG   ;
    while abs(e - iter) > 10^(-6)
        iter = e ;
        w = N_we * N_miu  ;
        if w <= 0.005
            e_d_avg = (0.0749 * sigma) / (452.84 * rho_g * (V_sG_avg)^2 * d_avg)  ;
        end
    
        if w > 0.005
            e_d_avg = (0.3713 * sigma * w^(0.302)) / (452.84 * rho_g * (V_sG_avg)^2 * d_avg)  ;
        end
    
        e = e_d_avg * d_avg  ;
        d_avg = d - e ; 
        V_sG_avg = V_sG * (d^2) / (d_avg)^2    ;
        N_we = 452.84 * (rho_g * (V_sG_avg)^2 * e) / sigma   ;
        N_miu = 6.23*10^(-3) * (((miu_L)^2) / (rho_L * sigma * e))  ;
    end

    Re_g = (1488 * rho_g * V_sG_avg * d_avg) / miu_g  ;
    f_tp = Moody (Re_g,eps)         ;
    dp_friction = (f_tp * rho_g * (V_sG_avg)^2) / (2 * 32.2 * d_avg) ;
    
    %acc
    E_K = (V_m * V_sG * rho_n) / (32.2 * P) ;
    
    dp_dL = ((dp_ele + dp_friction) / (1 - E_K)) / 144  ;
end

%Transition
if (L_s < Ngv) && (Ngv <= L_m)
    A = (L_m - N_gV) / (L_m - L_s)  ;
    B = 1 - A   ;
    
    %slug_section
    %elevation
    S = (1 + F_5) * ((Ngv^(0.982) + 0.029 * N_d + F_6) / (1 + F_7 * N_Lv)^2)   ;
    V_s = (S / 1.938) * (sigma / rho_L)^(1/4)   ;
    H_L = (V_s - V_m + sqrt((V_m - V_s)^2 + 4 * V_s * V_sL)) / (2 * V_s)  ;
    H_g = 1 - H_L   ;
    dp_ele = (rho_L * H_L) + (rho_g * H_g)   ;
    %friction
    f_3 = 1 + f_1 * sqrt(V_sG / (50 * V_sL))   ;
    f_tp = (f_1 * f_2) / f_3   ;
    dp_friction = (f_tp * rho_L * V_sL * V_m) / (2 * 32.2 * d)   ;
    %acc negligible

    dp_slug = dp_ele + dp_friction  ;
    
    %Annular_section
    %elevation
    Lambda_L = QL / (QL + Qg)  ;
    Lambda_g = 1 - Lambda_L    ;
    dp_ele = (rho_L *  Lambda_L) + (rho_g * Lambda_g)   ;
    rho_n = dp_ele  ;
    %friction
    rho_g = (rho_g * N_gv) / L_m  ;
    e = eps * d ;                                          % e = epsilon
    N_we = 452.84 * (rho_g * (V_sG)^2 * e) / sigma   ;
    N_miu = 6.23*10^(-3) * (((miu_L)^2) / (rho_L * sigma * e))  ;
    iter = 0 ;
    d_avg = d ;
    V_sG_avg = V_sG   ;
    while abs(e - iter) > 10^(-6)
        iter = e ;
        w = N_we * N_miu  ;
        if w <= 0.005
            e_d_avg = (0.0749 * sigma) / (452.84 * rho_g * (V_sG_avg)^2 * d_avg)  ;
        end
    
        if w > 0.005
            e_d_avg = (0.3713 * sigma * w^(0.302)) / (452.84 * rho_g * (V_sG_avg)^2 * d_avg)  ;
        end
    
        e = e_d_avg * d_avg  ;
        d_avg = d - e ; 
        V_sG_avg = V_sG * (d^2) / (d_avg)^2    ;
        N_we = 452.84 * (rho_g * (V_sG_avg)^2 * e) / sigma   ;
        N_miu = 6.23*10^(-3) * (((miu_L)^2) / (rho_L * sigma * e))  ;
    end

    Re_g = (1488 * rho_g * V_sG_avg * d_avg) / miu_g  ;
    f_tp = Moody (Re_g,eps)         ;
    dp_friction = (f_tp * rho_g * (V_sG_avg)^2) / (2 * 32.2 * d_avg) ;
    %acc
    E_K = (V_m * V_sG * rho_n) / (32.2 * P) ;
    
    dp_annular = (dp_ele + dp_friction) / (1 - E_K)  ;
    
    dp_dL = ((A * dp_slug) + (B * dp_annular)) / 144  ;
    
end
end


function [L__1] = L1 (x) 
% L_1 vs N_d
if (10 <= x) && (x <= 21) 
    L__1 = -0.0004*x^2 + 0.0146*x + 1.8968 ;
end

if (21 < x) && (x <= 81)
    L__1 = -6*10^(-9)*x^5 + 10^(-6)*x^4 - 2*10^(-5)*x^3 - 0.0033*x^2 + 0.1503*x + 0.3329 ;
end

if (81 < x) && (x <= 340)
    L__1 = 1 ;
end
end

function [L__2] = L2 (x)
% L_2 vs N_d
if (17 <= x) && ( x <= 72)
    L__2 = -2*10^(-6)*x^3 + 10^(-4)*x^2 + 0.0158*x + 0.1457 ;
end

if (72 < x) && (x <= 340)
    L__2 = -5*10^(-15)*x^6 + 7*10^(-12)*x^5  -4*10^(-9)*x^4 +10^(-6)*x^3 - 0.0002*x^2 + 0.012*x +0.7897 ;
end
end

function [F__1] = F1 (x)
% F_1 vs N_L
if (0.005 <= x) && (x <= 0.01)
    F__1 = 497.07*x^2 - 5.1602*x + 1.1649 ;
end

if (0.01 < x) && (x <= 0.08)
    F__1 = -5*10^(6)*x^5 + 10^(6)*x^4 - 141303*x^3 + 6668.7*x^2 - 109.37*x + 1.7507 ;
end

if (0.08 < x) && (x <= 1.63)
    F__1 = 0.6228*x^4 - 2.4804*x^3 + 3.7667*x^2 - 3.1892*x^ + 2.4517 ;
end
end

function [F__2] = F2 (x)
% F_2 vs N_L
if (0.002 <= x) && (x <= 0.006)
    F__2 = 7594.9*x^2 - 43.771*x + 0.9297 ;
end

if (0.006 < x) && (x <= .06)
    F__2 = -6*10^(7)*x^5 + 10^(7)*x^4 -681222*x^3 + 18984*x^2 - 135.87*x + 1.2035 ;
end

if (0.06 < x) && (x <= 0.6)
    F__2 = 21.826*x^3 - 21.273*x^2 + 6.9208*x + 3.3198 ;
end
end

function [F__3] = F3 (x)
% F_3 vs N_L
if (0.0025 <= x) && (x <= 0.0163)
    F__3 = 2*10^(6)*x^4 - 63124*x^3 + 657.22*x^2 - 2.6559*x + 0.2448 ;
end

if (0.0163 < x) && (x <= 0.146)
    F__3 = -332074*x^5 + 144063*x^4 - 23100*x^3 + 1607.4*x^2 - 36.145*x + 0.05047 ;
end

if (0.146 < x) && (x <= 1.75)
    F__3 = 0.1751*x^2 - 0.5066*x + 1.11 ;
end
end

function [F__4] = F4 (x)
% F_4 vs N_L
if (0.002 < x) && (x <= 0.00716)
    F__4 = 2*10^(8)*x^3 - 4*10^(6)*x^2 + 29924*x - 65.363 ;
end

if (0.00716 < x) && (x <= 0.0487)
    F__4 = -3*10^(7)*x^4 + 4*10^(6)*x^3 - 194107*x^2 + 5235.2*x -11.875 ;
end

if (0.0487 < x) && (x <= 1)
    F__4 = 106.03*x^5 - 329.3*x^4 + 385.45*x^3 - 209.26*x^2 + 49.336*x + 53.821 ;
end
end

function [F__5] = F5 (x)
% F_5 vs N_L
if (0.002 <= x) && (x <= 0.0457)
    F__5 = 78011*x^4 - 9305.9*x^3 + 403.08*x^2 - 8.5677*x + 0.1965 ;
end

if (0.0457 < x) && (x <= 0.118)
    F__5 = -15119*x^4 + 5208*x^3 - 635.02*x^2 + 31.067*x - 0.4274 ;
end

if (0.118 < x) && (x <= 2)
    F__5 = 0.137*x^6 - 0.8169*x^5 + 1.9053*x^4 - 2.1833*x^3 + 1.237*x^2 - 0.2515*x + 0.0376 ;
end
end

function [F__6] = F6 (x)
% F_6 vs N_L
if (0.002 <= x) && (x <= 0.0755)
    F__6 = 4*10^(8)*x^6 - 10^(8)*x^5 + 10^(7)*x^4 - 694041*x^3 + 21115*x^2 - 292.41*x + 1.3814 ;
end

if (0.0755 < x) && (x <= 0.224)
    F__6 = -6143.7*x^4 + 4465.8*x^3 - 1196.1*x^2 + 137.09*x - 3.6098 ;
end
if (0.224 < x) && (x <= 2)
    F__6 = 0.4959*x^6 - 3.4669*x^5 + 9.7738*x^4 - 14.201*x^3 + 11.202*x^2 - 4.456*x + 2.3787 ;
end
end

function [F__7] = F7 (x)
% F_7 vs N_L
if (0.002 <= x) && (x <= 0.03)
    F__7 = 86088*x^4 - 7145.5*x^3 +271.95*x^2 - 6.6153*x + 0.1081 ;
end

if (0.03 < x) && (x <= 0.139)
    F__7 = 420.44*x^4 - 166.29*x^3 + 24.416*x^2 - 1.6323*x + 0.0619 ;
end

if (0.139 < x) && (x <= 2.6)
    F__7 = 0.0009*x^4 - 0.0006*x^3 + 0.0137*x^2 - 0.0135*x + 0.0183 ;
end
end

function [f__1] = Moody (Re,eps)
% Moody diagram - friction factor vs Re
% Re = Reynolds number                eps = epsilon / d
if Re <= 2100
    % Laminar flow
    f__1 = 64 / Re ;
end

if Re > 2100
   % Turbulent flow
   if eps == 0
       %smooth
       f__1 = 0.0056 + 0.5*Re^(-0.32) ;
   end
   if eps ~= 0
       f__1 = 1 / (1.14 - 2*log10(eps + 21.25 / Re^.9))^2 ;
   end
end     
end

function [f__2] = f2 (x)
% f_2 vs (f_1 * V_sg * N_d^(2/3)) / V_sL
if (0.001 <= x) && (x <= 0.0124)
    f__2 = 9*10^(9)*x^5 -3*10^(8)*x^4 + 4*10^(6)*x^3 - 18578*x^2 + 32.466*x + .09824 ;
end

if (0.0124 < x) && (x <= 0.574)
    f__2 = 3.4224*x^3 - 4.2965*x^2 + 1.8128*x + 0.9695 ;
end

if (0.574 < x) && (x <= 7.84)
    f__2 = -0.0009*x^4 + 0.015*x^3 - 0.0733*x^2 - 0.0362*x + 1.2893 ;
end

if (7.84 < x) && (x <= 55.18)
    f__2 = -10^(-8)*x^5 + 2*10^(-6)*x^4 - 10^(-4)*x^3 + 0.0041*x^2 - 0.0732*x + 0.8468 ;
end

if (55.18 < x) && (x <= 1000)
    f__2 = -2*10^(-15)*x^5 + 5*10^(-12)*x^4 - 6*10^(-9)*x^3 + 3*10^(-6)*x^2 - 0.0007*x + 0.2745 ;
end
end

function [T1] = g (l)
L = 10000 ;
T_tf = 100  ;                                    
T_wf = 200  ; 
T1 = ((T_wf - T_tf) / L) * l + T_tf ;
end

function [rho] = SRK (P,T)
T_c = 407    ;                                       % Critical Temoerature      Unit = R
P_c = 607   ;
MW = 21  ;
R = 10.732  ;
T = T + 459.67  ;
P_r = P / P_c   ;
T_r = T / T_c   ;
z = (0.1 * log10(T_r) + 0.73)^P_r + 0.1 * P_r   ;
rho = P * MW / (z * R*T)   ;
end

function [miu__g] = miu (T , MW , rho_g)
T = T + 459.67      ;                         % F to R
rho_g = rho_g * 453.59 / 30.48^3      ;        % lb / ft^3 to g/cc
k = ((9.379 + 0.01607 * MW) * T^(1.5)) / (209.2 + 19.26 * MW + T)   ;
x = 3.448 + ( 986.4 ./ T ) + 0.01009 * MW   ;
y = 2.447 - 0.2224 * x   ;
miu__g = 10^(-4) * k * exp(x * (rho_g)^y)  ;
end