% % ------------------------------------------------------------------------------------------
% % ------------------------------------------rayonnement 3D------------------------------------------------
% % ------------------------------------------------------------------------------------------
Phi_s=0;
 Theta_s=pi/6; % changer Theta et phi (steering)


Phi=0;
Theta=pi/6;  % changer Theta et phi (steering)
 a = 0; b =d_x; c =0; d =d_y;
m = 8; n = 12;
dx = (b - a)/m; dy = (d - c)/n;
i = 1 : m; j = 1 : n;
u = a + (i - 1/2)*dx;
v = c + (j - 1/2)*dy;
[u,v] = meshgrid(u,v);
intgralmax=sum( sum(abs(calcul_Ee_J(u,v,A,Xs,alpha,beta,d_x,d_y,MB,NB,Ny,y1,y2,x1,x2,L,k,epsilon_r,omega,mu_0,epsilon_0,w)).*exp(1i.*k.*((u.*sin(Theta).*cos(Phi)+v.*sin(Theta).*sin(Phi))-(u.*sin(Theta_s).*cos(Phi_s)+v.*sin(Theta_s).*sin(Phi_s)))) ))*dx*dy;
radiationmax=intgralmax;
theta=-pi :0.1 : pi;
phi=0 : 0.2:2*pi;
[Theta, Phi] = meshgrid(theta, phi); % standard spherical co-ordinates
size_of_phi = size(Phi);
for i = 1: size_of_phi(1)
    for j = 1: size_of_phi(2)
 Su(i,j)=sin(Theta(i,j))*cos(Phi(i,j));
 Sv(i,j)=sin(Theta(i,j))*sin(Phi(i,j));
        if ((Su(i,j)^2) +(Sv(i,j)^2) <1)
        Integral(i,j)=sum(sum(abs(calcul_Ee_J(u,v,A,Xs,alpha,beta,d_x,d_y,MB,NB,Ny,y1,y2,x1,x2,L,k,epsilon_r,omega,mu_0,epsilon_0,w)).*exp(1i.*k.*((u.*sin(Theta(i,j)).*cos(Phi(i,j))+v.*sin(Theta(i,j)).*sin(Phi(i,j)))-(u.*sin(Theta_s).*cos(Phi_s)+v.*sin(Theta_s).*sin(Phi_s)))) ))*dx*dy;
        E_panel(i,j) = Integral(i,j)./radiationmax;
        W(i,j)=20.*log10(abs( E_panel(i,j)));
       else
         W(i,j)=20.*log10(abs(eps/10));
       end
    end   
end

% --------------------------Directivity (avec 3D reperesentation)-------------------------------------------
dth=theta(2)-theta(1) ;       
dph=phi(2)-phi(1);

Prad=sum(sum(abs(E_panel).^2.*sin(Theta)*dth*dph));
Directivity=4*pi*max(max(abs(E_panel).^2))/Prad;
Directivity_db=10*log10(Directivity)
Directivity_linear_value=((Directivity));
% ---------------------------------------------------------------------
% -------------Representation 3D--------------------------------------------------------
% ---------------------------------------------------------------------
  Su=sin(Theta).*cos(Phi);
Sv=sin(Theta).*sin(Phi);
dBmin =60;
 W(W<-dBmin) = NaN;
 figure(1);
h = surfc(Su,Sv, W);
set(h,'edgecolor','none')
xlabel('U');ylabel('V');colorbar
zlim([-dBmin dBmin])
figure(2)
surfc(Su,Sv,W)
shading flat
grid
axis image
xlabel('V')
ylabel('U');
axis([-1 1 -1 1]);colorbar
% % ------------------------------------------------------------------------------------------
% % ------------------------------------------rayonnement 3D------------------------------------------------
% % ------------------------------------------------------------------------------------------
Phi_s=0;
 Theta_s=pi/6; % changer Theta et phi (steering)


Phi=0;
Theta=pi/6;  % changer Theta et phi (steering)
 a = 0; b =d_x; c =0; d =d_y;
m = 8; n = 12;
dx = (b - a)/m; dy = (d - c)/n;
i = 1 : m; j = 1 : n;
u = a + (i - 1/2)*dx;
v = c + (j - 1/2)*dy;
[u,v] = meshgrid(u,v);
intgralmax=sum( sum(abs(calcul_Ee_J(u,v,A,Xs,alpha,beta,d_x,d_y,MB,NB,Ny,y1,y2,x1,x2,L,k,epsilon_r,omega,mu_0,epsilon_0,w)).*exp(1i.*k.*((u.*sin(Theta).*cos(Phi)+v.*sin(Theta).*sin(Phi))-(u.*sin(Theta_s).*cos(Phi_s)+v.*sin(Theta_s).*sin(Phi_s)))) ))*dx*dy;
radiationmax=intgralmax;
theta=-pi :0.1 : pi;
phi=0 : 0.2:2*pi;
[Theta, Phi] = meshgrid(theta, phi); % standard spherical co-ordinates
size_of_phi = size(Phi);
for i = 1: size_of_phi(1)
    for j = 1: size_of_phi(2)
 Su(i,j)=sin(Theta(i,j))*cos(Phi(i,j));
 Sv(i,j)=sin(Theta(i,j))*sin(Phi(i,j));
        if ((Su(i,j)^2) +(Sv(i,j)^2) <1)
        Integral(i,j)=sum(sum(abs(calcul_Ee_J(u,v,A,Xs,alpha,beta,d_x,d_y,MB,NB,Ny,y1,y2,x1,x2,L,k,epsilon_r,omega,mu_0,epsilon_0,w)).*exp(1i.*k.*((u.*sin(Theta(i,j)).*cos(Phi(i,j))+v.*sin(Theta(i,j)).*sin(Phi(i,j)))-(u.*sin(Theta_s).*cos(Phi_s)+v.*sin(Theta_s).*sin(Phi_s)))) ))*dx*dy;
        E_panel(i,j) = Integral(i,j)./radiationmax;
        W(i,j)=20.*log10(abs( E_panel(i,j)));
       else
         W(i,j)=20.*log10(abs(eps/10));
       end
    end   
end

% --------------------------Directivity (avec 3D reperesentation)-------------------------------------------
dth=theta(2)-theta(1) ;       
dph=phi(2)-phi(1);

Prad=sum(sum(abs(E_panel).^2.*sin(Theta)*dth*dph));
Directivity=4*pi*max(max(abs(E_panel).^2))/Prad;
Directivity_db=10*log10(Directivity)
Directivity_linear_value=((Directivity));
% ---------------------------------------------------------------------
% -------------Representation 3D--------------------------------------------------------
% ---------------------------------------------------------------------
  Su=sin(Theta).*cos(Phi);
Sv=sin(Theta).*sin(Phi);
dBmin =60;
 W(W<-dBmin) = NaN;
 figure(1);
h = surfc(Su,Sv, W);
set(h,'edgecolor','none')
xlabel('U');ylabel('V');colorbar
zlim([-dBmin dBmin])
figure(2)
surfc(Su,Sv,W)
shading flat
grid
axis image
xlabel('V')
ylabel('U');
axis([-1 1 -1 1]);colorbar
