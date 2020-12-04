%% CLEANING MEMORY AND CLOSING WINDOWS
clc; clear; close all;

%% PURPOSE: The code aims to compute the maximum displacement and stresses 
% values and the energy norm given by the First Order Shear Deformation
% AUTHOR: Renan Miranda Portela
% SUBJECT: 

%% MATERIAL INPUT
material = 'orthotropic'; % isotropic or orthotropic

a = 1; b = 1;       % Plate length and width
h = 1/10;           % Plate thickness
q = 1e-3;           % Load intensity
fprintf('a/h = %4.1f \n', a/h)

Ks = 5/6; %shear correction factor
fprintf('shear correction factor = %4.2f \n', Ks)

switch material
    case 'isotropic'
        E = 10920; E2 = E;
        nu = 0.25;
        Q11 = E/(1-nu^2);
        Q22 = E2/(1-nu^2);
        Q12 = E*nu/(1-nu^2);
        Q66 = (Q11-Q12)/2;
        Q55 = Q66;
        Q44 = Q66;  
        A11 = E*h/(1-nu^2);
        A12 = nu*A11;
        A22 = A11;
        A66 = (1-nu)/2*A11;
        A44 = Ks*A66;
        A55 = Ks*A66;
        D11 = E*h^3/(12*(1-nu^2));
        D12 = nu*D11;
        D22 = D11; D16 = 0; D26 = 0;
        D66 = (1-nu)/2*D11;
    case 'orthotropic'
        E1 = 25e3;
        E2 = 1e3;
        G23 = 0.2e3;
        G13 = 0.5e3;
        G12 = 0.5e3;
        nu12 = 0.25; nu21 = nu12*E2/E1;
        
       Q11 = E1/(1-nu12*nu21);
       Q22 = E2/(1-nu12*nu21);
       Q12 = E1*nu21/(1-nu12*nu21);
       Q66 = G12;
       Q55 = G23;
       Q44 = G13;        

        A11 = E1*h/(1 - nu12*nu21);
        A22 = E2*h/(1 - nu12*nu21);
        A12 = E1*h*nu21/(1 - nu12*nu21);
        A66 = G12*h;
        A44 = Ks*G23*h;
        A55 = Ks*G13*h;
        
        D11 = E1*h^3/(12*(1 - nu12*nu21));
        D22 = E2*h^3/(12*(1 - nu12*nu21));
        D12 = nu21*E1*h^3/(12*(1 - nu12*nu21));
        D16 = 0; D26 = 0;
        D66 = G12*h^3/12;
        
            case 'Laminate (0/90/0)'
        E1 = 25e3;
        E2 = 1e3;
        G23 = 200;
        G13 = 500;
        G12 = 500;
        nu12 = 0.25; nu21 = nu12*E2/E1;
        
        geo = [1 -h/2 -h/2 + h/3 0;
               2 -h/2 + h/3 -h/2 + 2*h/3 pi/2;
               3 -h/2 + 2*h/3 h/2 0];
           
       A11 = 0; A12 = 0; A22 = 0; A44 = 0; A55 = 0; A66 = 0;
       D11 = 0; D12 = 0; D22 = 0; D44 = 0; D55 = 0; D66 = 0;
       F11 = 0; F12 = 0; F22 = 0; F44 = 0; F55 = 0; F66 = 0;
       H11 = 0; H12 = 0; H22 = 0; H44 = 0; H55 = 0; H66 = 0;
           
       for ii = 1 : size(geo,1)
           Q(1,1) = E1/(1-nu12*nu21);
           Q(2,2) = E2/(1-nu12*nu21);
           Q(1,2) = E1*nu21/(1-nu12*nu21); Q(2,1) = Q(1,2);
           Q(3,3) = G12;
           Q(4,4) = G23;
           Q(5,5) = G13; Q(6,6) = 1;
           
           zbot = geo(ii,2);
           ztop = geo(ii,3);
           theta = geo(ii,4);
           
           m = cos(theta);
           n = sin(theta);
           
      T = [m^2, n^2,   0,   0,   0,   2*m*n;
          n^2, m^2,   0,   0,   0,  -2*m*n;
            0,   0,   1,   0,   0,       0;
            0,   0,   0,   m,  -n,       0;
            0,   0,   0,   n,   m,       0;
         -m*n, m*n,   0,   0,   0, m^2-n^2];  
     
           Q_bar = T\Q*T; 
           Q11 = Q_bar(1,1);
           Q12 = Q_bar(1,2);
           Q22 = Q_bar(2,2);
           Q66 = Q_bar(3,3);
           Q44 = Q_bar(4,4);
           Q55 = Q_bar(5,5);
           
           A11 = A11 + Q11*(ztop - zbot);
           A22 = A22 + Q22*(ztop - zbot);
           A12 = A12 + Q12*(ztop - zbot);
           A66 = A66 + Q66*(ztop - zbot);
           A44 = A44 + Ks*Q44*(ztop - zbot);
           A55 = A55 + Ks*Q55*(ztop - zbot);
           
           D11 = D11 + Q11*(ztop^3 - zbot^3)/3;
           D12 = D12 + Q12*(ztop^3 - zbot^3)/3;
           D22 = D22 + Q22*(ztop^3 - zbot^3)/3;
           D44 = D44 + Q44*(ztop^3 - zbot^3)/3;
           D55 = D55 + Q55*(ztop^3 - zbot^3)/3;
           D66 = D66 + Q66*(ztop^3 - zbot^3)/3;
           
           F11 = F11 + Q11*(ztop^5 - zbot^5)/5;
           F12 = F12 + Q12*(ztop^5 - zbot^5)/5;
           F22 = F22 + Q22*(ztop^5 - zbot^5)/5;
           F44 = F44 + Q44*(ztop^5 - zbot^5)/5;
           F55 = F55 + Q55*(ztop^5 - zbot^5)/5;
           F66 = F66 + Q66*(ztop^5 - zbot^5)/5;
           
           H11 = H11 + Q11*(ztop^7 - zbot^7)/7;
           H12 = H12 + Q12*(ztop^7 - zbot^7)/7;
           H22 = H22 + Q22*(ztop^7 - zbot^7)/7;
           H44 = H44 + Q44*(ztop^7 - zbot^7)/7;
           H55 = H55 + Q55*(ztop^7 - zbot^7)/7;
           H66 = H66 + Q66*(ztop^7 - zbot^7)/7;
       end
       case 'Laminate (0/90/90/0)'
        E1 = 25e3;
        E2 = 1e3;
        G23 = 200;
        G13 = 500;
        G12 = 500;
        nu12 = 0.25; nu21 = nu12*E2/E1;
        
        geo = [1 -h/2 -h/4 0;
               2 -h/4 0 pi/2;
               3  0   h/4 pi/2;
               4  h/4 h/2 0];
           
       A11 = 0; A12 = 0; A22 = 0; A44 = 0; A55 = 0; A66 = 0;
       D11 = 0; D12 = 0; D22 = 0; D44 = 0; D55 = 0; D66 = 0;
       F11 = 0; F12 = 0; F22 = 0; F44 = 0; F55 = 0; F66 = 0;
       H11 = 0; H12 = 0; H22 = 0; H44 = 0; H55 = 0; H66 = 0;
           
       for ii = 1 : size(geo,1)
           Q(1,1) = E1/(1-nu12*nu21);
           Q(2,2) = E2/(1-nu12*nu21);
           Q(1,2) = E1*nu21/(1-nu12*nu21); Q(2,1) = Q(1,2);
           Q(3,3) = G12;
           Q(4,4) = G23;
           Q(5,5) = G13; Q(6,6) = 1;
           
           zbot = geo(ii,2);
           ztop = geo(ii,3);
           theta = geo(ii,4);
           
           m = cos(theta);
           n = sin(theta);
           
      T = [m^2, n^2,   0,   0,   0,   2*m*n;
          n^2, m^2,   0,   0,   0,  -2*m*n;
            0,   0,   1,   0,   0,       0;
            0,   0,   0,   m,  -n,       0;
            0,   0,   0,   n,   m,       0;
         -m*n, m*n,   0,   0,   0, m^2-n^2];  
     
           Q_bar = T\Q*T; 
           Q11 = Q_bar(1,1);
           Q12 = Q_bar(1,2);
           Q22 = Q_bar(2,2);
           Q66 = Q_bar(3,3);
           Q44 = Q_bar(4,4);
           Q55 = Q_bar(5,5);
           
           A11 = A11 + Q11*(ztop - zbot);
           A22 = A22 + Q22*(ztop - zbot);
           A12 = A12 + Q12*(ztop - zbot);
           A66 = A66 + Q66*(ztop - zbot);
           A44 = A44 + Ks*Q44*(ztop - zbot);
           A55 = A55 + Ks*Q55*(ztop - zbot);
           
           D11 = D11 + Q11*(ztop^3 - zbot^3)/3;
           D12 = D12 + Q12*(ztop^3 - zbot^3)/3;
           D22 = D22 + Q22*(ztop^3 - zbot^3)/3;
           D44 = D44 + Q44*(ztop^3 - zbot^3)/3;
           D55 = D55 + Q55*(ztop^3 - zbot^3)/3;
           D66 = D66 + Q66*(ztop^3 - zbot^3)/3;
           
           F11 = F11 + Q11*(ztop^5 - zbot^5)/5;
           F12 = F12 + Q12*(ztop^5 - zbot^5)/5;
           F22 = F22 + Q22*(ztop^5 - zbot^5)/5;
           F44 = F44 + Q44*(ztop^5 - zbot^5)/5;
           F55 = F55 + Q55*(ztop^5 - zbot^5)/5;
           F66 = F66 + Q66*(ztop^5 - zbot^5)/5;
           
           H11 = H11 + Q11*(ztop^7 - zbot^7)/7;
           H12 = H12 + Q12*(ztop^7 - zbot^7)/7;
           H22 = H22 + Q22*(ztop^7 - zbot^7)/7;
           H44 = H44 + Q44*(ztop^7 - zbot^7)/7;
           H55 = H55 + Q55*(ztop^7 - zbot^7)/7;
           H66 = H66 + Q66*(ztop^7 - zbot^7)/7;
       end
end

%% ENERGY NORM 

syms x y m n Wn O_x O_y

w = Wn*sin(m*pi*x/a)*sin(n*pi*y/b);
phi_x = O_x*cos(m*pi*x/a)*sin(n*pi*y/b);
phi_y = O_y*sin(m*pi*x/a)*cos(n*pi*y/b);

dwdx = diff(w,x);
dwdy = diff(w,y);

dphi_xdx = diff(phi_x,x);
dphi_xdy = diff(phi_x,y);

dphi_ydx = diff(phi_y,x);
dphi_ydy = diff(phi_y,y);

df1 = A55*dwdx*dwdx - A55*phi_x*dwdx + A44*dwdy*dwdy - A44*phi_y*dwdy;
df2 = -D11*dphi_xdx*dphi_xdx - D12*dphi_ydy*dphi_xdx - D16*(dphi_xdy*dphi_xdx + dphi_ydx*dphi_xdx) + ...
      -D16*dphi_xdx*dphi_xdy - D26*dphi_ydy*dphi_xdy - D66*(dphi_xdy*dphi_xdy + dphi_ydx*dphi_xdy) + ...
      A55*dwdx*phi_x - A55*phi_x*phi_x;
df3 = -D16*dphi_xdx*dphi_ydx - D26*dphi_ydy*dphi_ydx - D66*(dphi_xdy*dphi_ydx + dphi_ydx*dphi_ydx) + ...
      -D12*dphi_xdx*dphi_ydy - D22*dphi_ydy*dphi_ydy - D26*(dphi_xdy*dphi_ydy + dphi_ydx*dphi_ydy) + ...
      A44*dwdy*phi_y - A44*phi_y*phi_y;  
df = df1 + df2 + df3;

Ea = 0;

%% SINUSOIDAL DISTRIBUTED LOAD

APX = 1;
odd = 1:2:2*APX;
w = 0;
strxx = 0; stryy = 0; strxy = 0; strxz = 0; stryz = 0;

for ii = 1 : APX
    for jj = 1 : APX
        K = zeros(3);
        R = zeros(3,1);
        mi = odd(ii);
        ni = odd(jj);
        alfa = mi*pi/a;
        beta = ni*pi/b;
        R(1,1) = q;
        K(1,1) = alfa^2*A55 + beta^2*A44;
        K(1,2) = alfa*A55;
        K(1,3) = beta*A44;
        K(2,1) = K(1,2);
        K(2,2) = D11*alfa^2 + D66*beta^2 + A55;
        K(2,3) = (D12 + D66)*alfa*beta;
        K(3,1) = K(1,3);
        K(3,2) = K(2,3);
        K(3,3) = D66*alfa^2 + D22*beta^2 + A44;
        U = K\R;
        F = subs(df, [Wn O_x O_y m n], [U(1) U(2) U(3) mi ni]);
        Ea = Ea + double(int(int(F, x, [0 1]), y, [0 1]));
        fprintf('Ea = %d \n', Ea)
        w = w + U(1)*sin(mi*pi/2)*sin(ni*pi/2);    
        strxx = strxx - h/2*(Q11*alfa*U(2) + Q12*beta*U(3));
        stryy = stryy - h/2*(Q12*alfa*U(2) + Q22*beta*U(3));
        strxy = strxy - h/2*(Q66*(beta*U(2) + alfa*U(3)));
        stryz = stryz + Q55*(U(3) + beta*U(1));
        strxz = strxz + Q44*(U(2) + alfa*U(1));
    end
end

fprintf('Wmn (Sinusoidal) = %4.4f \n', w*D11/q/a*100)

%% STRESSES

fprintf('S_xx = %4.4f \n', strxx*h^2/a^2/q)
fprintf('S_yy = %4.4f \n', stryy*h^2/a^2/q)
fprintf('S_xy = %4.4f \n', strxy*h^2/a^2/q)
fprintf('S_xz = %4.4f \n', strxz*h/a/q)
fprintf('S_yz = %4.4f \n', stryz*h/a/q)

%% UNIFORMLY DISTRIBUTED LOAD

APX = 19;
odd = 1:2:2*APX;
w = 0; strxx = 0; stryy = 0; strxy = 0; strxz = 0; stryz = 0;

for ii = 1 : APX
    for jj = 1 : APX
        K = zeros(3);
        R = zeros(3,1);
        
        mi = odd(ii);
        ni = odd(jj);
        alfa = mi*pi/a;
        beta = ni*pi/b;
        
        R(1,1) = 16*q/pi^2/mi/ni;
        
        K(1,1) = alfa^2*A55 + beta^2*A44;
        K(1,2) = alfa*A55;
        K(1,3) = beta*A44;
        
        K(2,1) = K(1,2);
        K(2,2) = D11*alfa^2 + D66*beta^2 + A55;
        K(2,3) = (D12 + D66)*alfa*beta;
        
        K(3,1) = K(1,3);
        K(3,2) = K(2,3);
        K(3,3) = D66*alfa^2 + D22*beta^2 + A44;
        
        U = K\R;
        w = w + U(1)*sin(mi*pi/2)*sin(ni*pi/2);
        
        strxx = strxx - h/2*(Q11*alfa*U(2) + Q12*beta*U(3))*sin(mi*pi/2)*sin(ni*pi/2);
        stryy = stryy - h/2*(Q12*alfa*U(2) + Q22*beta*U(3))*sin(mi*pi/2)*sin(ni*pi/2);
        strxy = strxy - h/2*(Q66*(beta*U(2) + alfa*U(3)))*cos(mi*pi)*cos(ni*pi);
        stryz = stryz + Q55*(U(3) + beta*U(1))*sin(mi*pi/2)*cos(ni*pi);
        strxz = strxz + Q55*(U(2) + alfa*U(1))*cos(mi*pi)*sin(ni*pi/2);        
    end
end

fprintf('Wmn (Uniform) = %4.4f \n', w*h^3*E2/q/a^4)
%% STRESSES

fprintf('S_xx = %4.4f \n', strxx*h^2/a^2/q)
fprintf('S_yy = %4.4f \n', stryy*h^2/a^2/q)
fprintf('S_xy = %4.4f \n', strxy*h^2/a^2/q)
fprintf('S_xz = %4.4f \n', strxz*h/a/q)
fprintf('S_yz = %4.4f \n', stryz*h/a/q)