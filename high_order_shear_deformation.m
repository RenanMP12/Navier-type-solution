%% CLEANING MEMORY AND CLOSING WINDOWS
clear; clc;

%% PURPOSE: The code aims to compute the maximum displacement and stresses 
% values and the energy norm given by the High Order Shear Deformation
% AUTHOR: Renan Miranda Portela

%% MATERIAL INPUT
material = 'isotropic'; % isotropic or orthotropic

a = 1; b = 1; q = 1e-3; h = 1/10;  
fprintf('a/h = %4.1f \n', a/h)
c1 = 4/3/h^2; c2 = 3*c1;

switch material
    case 'isotropic'
        E1 = 10920; E2 = E1;
        nu = 0.3;
        
        geo = [1 -h/2 h/2 0];
           
        A = 0; B = 0; D = 0; E = 0; F = 0; H = 0; 
           
       for ii = 1 : size(geo,1)
           C11 = E1/(1-nu^2);
           C12 = nu*E1/(1-nu^2);
           Q = [C11 C12 0 0 0 C12;
                C12 C11 0 0 0 C12;
                0  0 (C11-C12)/2 0 0 0;
                0  0 0 (C11-C12)/2 0 0;
                0  0 0 0 (C11-C12)/2 0;
                C12 C12 0 0 0 C11]; 
           
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
            
           A = A + Q_bar*(ztop - zbot); 
           A11 = A(1,1); A12 = A(1,2); A16 = A(1,3); A22 = A(2,2); 
           A26 = A(2,3); A44 = A(4,4); A55 = A(5,5); A66 = A(3,3);
           
           B = B + Q_bar*(ztop^2 - zbot^2)/2;
           B11 = B(1,1); B12 = B(1,2); B16 = B(1,3); B22 = B(2,2); 
           B26 = B(2,3); B44 = B(4,4); B55 = B(5,5); B66 = B(3,3);
           
           D = D + Q_bar*(ztop^3 - zbot^3)/3;
           D11 = D(1,1); D12 = D(1,2); D16 = D(1,3); D22 = D(2,2); 
           D26 = D(2,3); D44 = D(4,4); D55 = D(5,5); D66 = D(3,3); D45 = D(4,5);          
           
           E = E + Q_bar*(ztop^4 - zbot^4)/4;
           E11 = E(1,1); E12 = E(1,2); E16 = E(1,3); E22 = E(2,2); 
           E26 = E(2,3); E44 = E(4,4); E55 = E(5,5); E66 = E(3,3);           
           
           F = F + Q_bar*(ztop^5 - zbot^5)/5;
           F11 = F(1,1); F12 = F(1,2); F16 = F(1,3); F22 = F(2,2); 
           F26 = F(2,3); F44 = F(4,4); F55 = F(5,5); F66 = F(3,3);           
           
           H = H + Q_bar*(ztop^7 - zbot^7)/7;
           H11 = H(1,1); H12 = H(1,2); H16 = H(1,3); H22 = H(2,2); 
           H26 = H(2,3); H44 = H(4,4); H55 = H(5,5); H66 = H(3,3);
       end
        
    case 'orthotropic'
        E1 = 25e3;
        E2 = 1e3;
        G23 = 0.2e3;
        G13 = 0.5e3;
        G12 = 0.5e3;
        nu12 = 0.25; nu21 = nu12*E2/E1;
        
        geo = [1 -h/2 h/2 0];
           
        A = 0; B = 0; D = 0; E = 0; F = 0; H = 0; 
           
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
            
           A = A + Q_bar*(ztop - zbot); 
           A11 = A(1,1); A12 = A(1,2); A16 = A(1,3); A22 = A(2,2); 
           A26 = A(2,3); A44 = A(4,4); A55 = A(5,5); A66 = A(3,3);
           
           B = B + Q_bar*(ztop^2 - zbot^2)/2;
           B11 = B(1,1); B12 = B(1,2); B16 = B(1,3); B22 = B(2,2); 
           B26 = B(2,3); B44 = B(4,4); B55 = B(5,5); B66 = B(3,3);
           
           D = D + Q_bar*(ztop^3 - zbot^3)/3;
           D11 = D(1,1); D12 = D(1,2); D16 = D(1,3); D22 = D(2,2); 
           D26 = D(2,3); D44 = D(4,4); D55 = D(5,5); D66 = D(3,3); D45 = D(4,5);          
           
           E = E + Q_bar*(ztop^4 - zbot^4)/4;
           E11 = E(1,1); E12 = E(1,2); E16 = E(1,3); E22 = E(2,2); 
           E26 = E(2,3); E44 = E(4,4); E55 = E(5,5); E66 = E(3,3);           
           
           F = F + Q_bar*(ztop^5 - zbot^5)/5;
           F11 = F(1,1); F12 = F(1,2); F16 = F(1,3); F22 = F(2,2); 
           F26 = F(2,3); F44 = F(4,4); F55 = F(5,5); F66 = F(3,3);           
           
           H = H + Q_bar*(ztop^7 - zbot^7)/7;
           H11 = H(1,1); H12 = H(1,2); H16 = H(1,3); H22 = H(2,2); 
           H26 = H(2,3); H44 = H(4,4); H55 = H(5,5); H66 = H(3,3);
       end
        
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
           
        A = 0; B = 0; D = 0; E = 0; F = 0; H = 0; 
           
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
            
           A = A + Q_bar*(ztop - zbot); 
           A11 = A(1,1); A12 = A(1,2); A16 = A(1,3); A22 = A(2,2); 
           A26 = A(2,3); A44 = A(4,4); A55 = A(5,5); A66 = A(3,3);
           
           B = B + Q_bar*(ztop^2 - zbot^2)/2;
           B11 = B(1,1); B12 = B(1,2); B16 = B(1,3); B22 = B(2,2); 
           B26 = B(2,3); B44 = B(4,4); B55 = B(5,5); B66 = B(3,3);
           
           D = D + Q_bar*(ztop^3 - zbot^3)/3;
           D11 = D(1,1); D12 = D(1,2); D16 = D(1,3); D22 = D(2,2); 
           D26 = D(2,3); D44 = D(4,4); D55 = D(5,5); D66 = D(3,3);           
           
           E = E + Q_bar*(ztop^4 - zbot^4)/4;
           E11 = E(1,1); E12 = E(1,2); E16 = E(1,3); E22 = E(2,2); 
           E26 = E(2,3); E44 = E(4,4); E55 = E(5,5); E66 = E(3,3);           
           
           F = F + Q_bar*(ztop^5 - zbot^5)/5;
           F11 = F(1,1); F12 = F(1,2); F16 = F(1,3); F22 = F(2,2); 
           F26 = F(2,3); F44 = F(4,4); F55 = F(5,5); F66 = F(3,3);           
           
           H = H + Q_bar*(ztop^7 - zbot^7)/7;
           H11 = H(1,1); H12 = H(1,2); H16 = H(1,3); H22 = H(2,2); 
           H26 = H(2,3); H44 = H(4,4); H55 = H(5,5); H66 = H(3,3);
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
           
         A = 0; B = 0; D = 0; E = 0; F = 0; H = 0; 
           
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
            
           A = A + Q_bar*(ztop - zbot); 
           A11 = A(1,1); A12 = A(1,2); A16 = A(1,3); A22 = A(2,2); 
           A26 = A(2,3); A44 = A(4,4); A55 = A(5,5); A66 = A(3,3);
           
           B = B + Q_bar*(ztop^2 - zbot^2)/2;
           B11 = B(1,1); B12 = B(1,2); B16 = B(1,3); B22 = B(2,2); 
           B26 = B(2,3); B44 = B(4,4); B55 = B(5,5); B66 = B(3,3);
           
           D = D + Q_bar*(ztop^3 - zbot^3)/3;
           D11 = D(1,1); D12 = D(1,2); D16 = D(1,3); D22 = D(2,2); 
           D26 = D(2,3); D44 = D(4,4); D55 = D(5,5); D66 = D(3,3);           
           
           E = E + Q_bar*(ztop^4 - zbot^4)/4;
           E11 = E(1,1); E12 = E(1,2); E16 = E(1,3); E22 = E(2,2); 
           E26 = E(2,3); E44 = E(4,4); E55 = E(5,5); E66 = E(3,3);           
           
           F = F + Q_bar*(ztop^5 - zbot^5)/5;
           F11 = F(1,1); F12 = F(1,2); F16 = F(1,3); F22 = F(2,2); 
           F26 = F(2,3); F44 = F(4,4); F55 = F(5,5); F66 = F(3,3);           
           
           H = H + Q_bar*(ztop^7 - zbot^7)/7;
           H11 = H(1,1); H12 = H(1,2); H16 = H(1,3); H22 = H(2,2); 
           H26 = H(2,3); H44 = H(4,4); H55 = H(5,5); H66 = H(3,3);
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

df1 = c2*(D45*phi_y*dwdx + D45*dwdy*dwdx + D55*phi_x*dwdx + D55*dwdx*dwdx) + ...
      c2*(D45*phi_y*dwdx + D45*dwdy*dwdx + D55*phi_x*dwdx + D55*dwdx*dwdx +...
      c2*(F45*phi_y*dwdx + F45*dwdy*dwdx + F55*phi_x*dwdx + F55*dwdx*dwdx)) +...
      c2*(D44*phi_y*dwdy + D44*dwdy*dwdy + D45*phi_x*dwdy + D45*dwdx*dwdy) + ...
      c2*(D44*phi_y*dwdy + D44*dwdy*dwdy + D45*phi_x*dwdy + D45*dwdx*dwdy +...
      c2*(F44*phi_y*dwdy + F44*dwdy*dwdy + F45*phi_x*dwdy + F45*dwdx*dwdy)) +...

%% SINUSOIDAL

APX = 1;
odd = 1:2:2*APX;
w = 0; strxx = 0; stryy = 0; strxy = 0; strxz = 0; stryz = 0;

for ii = 1 : APX
    for jj = 1 : APX
        S = zeros(3);
        R = zeros(3,1);
        
        m = odd(ii);
        n = odd(jj);
        alfa = m*pi/a;
        beta = n*pi/b;
        
        R(1,1) = q;
        
        S(1,1) = alfa^2*A55 + beta^2*A44 - 8/h^2*( alfa^2*D55 + beta^2*D44 ) + (4/h^2)^2*(alfa^2*F55 + beta^2*F44) + (4/3/h^2)^2*( alfa^4*H11 + 2*(H12 + 2*H66)*alfa^2*beta^2 + beta^4*H22 );
        S(1,2) = alfa*A55 - 8/h^2*alfa*D55 + (4/h^2)^2*alfa*F55 - 4/3/h^2*(alfa^3*F11 + alfa*beta^2*(F12+2*F66)) + (4/3/h^2)^2*( alfa^3*H11 + alfa*beta^2*(H12 +2*H66) ); S(2,1) = S(1,2);
        S(1,3) = beta*A44 - 8/h^2*beta*D44 + (4/h^2)^2*beta*F44 - 4/3/h^2*(alfa^2*beta*(F12+2*F66)+beta^3*F22) + (4/3/h^2)^2*(alfa^2*beta*(H12+2*H66)+beta^3*H22); S(3,1) = S(1,3);
        S(2,2) = A55 + alfa^2*D11 + beta^2*D66 - 8/h^2*D55 + (4/h^2)^2*F55 - 8/3/h^2*(alfa^2*F11 + beta^2*F66) + (4/3/h^2)^2*(alfa^2*H11 + beta^2*H66);
        S(2,3) = alfa*beta*(D12 + D66 - 8/3/h^2*(F12+F66) + (4/3/h^2)^2*(H12 + H66)); S(3,2) = S(2,3);
        S(3,3) = A44 + alfa^2*D66 + beta^2*D22 - 8/h^2*D44 + (4/h^2)^2*F44 - 8/3/h^2*(alfa^2*F66 + beta^2*F22) + (4/3/h^2)^2*(beta^2*H22 + alfa^2*H66);
        
        U = S\R;
        w = w + U(1)*sin(m*pi/2)*sin(n*pi/2);
        Z = h/2; 
        strxx = strxx + Q_bar(1,1)*( -Z*alfa*U(2) + Z^3*c1*(alfa*U(2) + alfa^2*U(1)))*sin(m*pi/2)*sin(n*pi/2) + Q_bar(1,2)*(-Z*beta*U(3) + Z^3*c1*(beta*U(3)+alfa^2*U(1)))*sin(m*pi/2)*sin(n*pi/2);
%         Z = h/4; 
        strxy = strxy - Q_bar(3,3)*(Z*(beta*U(2) + alfa*U(3)) - Z^3*c1*(beta*U(2) + alfa*U(3) + 2*alfa*beta*U(1)))*cos(m*pi)*cos(n*pi);
        Z = 0;
        stryz = stryz - Q_bar(5,5)*(1 - c2*Z^2)*(U(3) + beta*U(1))*sin(m*pi/2)*cos(n*pi);
        strxz = strxz - Q_bar(4,4)*(U(2) + alfa*U(1))*cos(m*pi)*sin(m*pi/2);
    end
end

fprintf('Nondimensional displacement (Sinusoidal) = %4.4f\n', w*h^3*E2/q/a^4*100)
fprintf('Nondimensional stress xx = %4.4f\n', strxx*h^2/q/a^2)
fprintf('Nondimensional stress xy = %4.4f\n', strxy*h^2/q/a^2)
fprintf('Nondimensional stress yz = %4.4f\n', stryz*h/q/a)
fprintf('Nondimensional stress xz = %4.4f\n', strxz*h/q/a)

%% UNIFORM

APX = 29;
odd = 1:2:2*APX;
w = 0; strxx = 0; stryy = 0; strxy = 0; strxz = 0; stryz = 0;

for ii = 1 : APX
    for jj = 1 : APX
        S = zeros(3);
        R = zeros(3,1);
        
        m = odd(ii);
        n = odd(jj);
        alfa = m*pi/a;
        beta = n*pi/b;
        
        R(1,1) = 16*q/pi^2/m/n;
        
        S(1,1) = alfa^2*A55 + beta^2*A44 - 8/h^2*( alfa^2*D55 + beta^2*D44 ) + (4/h^2)^2*(alfa^2*F55 + beta^2*F44) + (4/3/h^2)^2*( alfa^4*H11 + 2*(H12 + 2*H66)*alfa^2*beta^2 + beta^4*H22 );
        S(1,2) = alfa*A55 - 8/h^2*alfa*D55 + (4/h^2)^2*alfa*F55 - 4/3/h^2*(alfa^3*F11 + alfa*beta^2*(F12+2*F66)) + (4/3/h^2)^2*( alfa^3*H11 + alfa*beta^2*(H12 +2*H66) ); S(2,1) = S(1,2);
        S(1,3) = beta*A44 - 8/h^2*beta*D44 + (4/h^2)^2*beta*F44 - 4/3/h^2*(alfa^2*beta*(F12+2*F66)+beta^3*F22) + (4/3/h^2)^2*(alfa^2*beta*(H12+2*H66)+beta^3*H22); S(3,1) = S(1,3);
        S(2,2) = A55 + alfa^2*D11 + beta^2*D66 - 8/h^2*D55 + (4/h^2)^2*F55 - 8/3/h^2*(alfa^2*F11 + beta^2*F66) + (4/3/h^2)^2*(alfa^2*H11 + beta^2*H66);
        S(2,3) = alfa*beta*(D12 + D66 - 8/3/h^2*(F12+F66) + (4/3/h^2)^2*(H12 + H66)); S(3,2) = S(2,3);
        S(3,3) = A44 + alfa^2*D66 + beta^2*D22 - 8/h^2*D44 + (4/h^2)^2*F44 - 8/3/h^2*(alfa^2*F66 + beta^2*F22) + (4/3/h^2)^2*(beta^2*H22 + alfa^2*H66);
        
        U = S\R;
        w = w + U(1)*sin(m*pi/2)*sin(n*pi/2);
        Z = h/2; 
        strxx = strxx + Q_bar(1,1)*( -Z*alfa*U(2) + Z^3*c1*(alfa*U(2) + alfa^2*U(1)))*sin(m*pi/2)*sin(n*pi/2) + Q_bar(1,2)*(-Z*beta*U(3) + Z^3*c1*(beta*U(3)+alfa^2*U(1)))*sin(m*pi/2)*sin(n*pi/2);
        strxy = strxy - Q_bar(3,3)*(Z*(beta*U(2) + alfa*U(3)) - Z^3*c1*(beta*U(2) + alfa*U(3) + 2*alfa*beta*U(1)))*cos(m*pi)*cos(n*pi);
        Z = 0;
        stryz = stryz - Q_bar(5,5)*(1 - c2*Z^2)*(U(3) + beta*U(1))*sin(m*pi/2)*cos(n*pi);
        strxz = strxz - Q_bar(4,4)*(U(2) + alfa*U(1))*cos(m*pi)*sin(m*pi/2);
    end
end

fprintf('Nondimensional displacement (Uniform N = %d) = %4.4f\n', APX, w*h^3*E2/q/a^4*100)
fprintf('Nondimensional stress xx = %4.4f\n', strxx*h^2/q/a^2)
fprintf('Nondimensional stress xy = %4.4f\n', strxy*h^2/q/a^2)
fprintf('Nondimensional stress yz = %4.4f\n', stryz*h/q/a)
fprintf('Nondimensional stress xz = %4.4f\n', strxz*h/q/a)