%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% last update 6th December 2020, lne %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -> This program computes the Schrodinger-Poisson equation in heterostructures with
% the Dirichlet boundary conditions. This means that the doping level on the left 
% and right side of the structure are the boundary condition since they are setting
% the Fermi level.
% -> The quantum structure is sandwitch between contacts with spacers. The Schrodinger 
% solver is working ONLY in this domain.
% -> In order to keep the code fast but still usefull, the mass is kept constant
% all over the structure. It means that meff should be set at the value of the
% well. Obviously, the non-parabolicity of the bands is also not considered in 
% the Schrodinger solver and the density of states.
% -> Schottky contact can be simulated by setting the doping of the contact
% at zero and the bandgap energy of the contact material.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the code doesn t converge:
% -> increase the amount of loops, Nloops
% -> increase the damping, tau0
% -> increase the resolution dz
% -> increase the resolution dE
% -> increase the temperature (T=0K is very bad while T=10K is already much better)
% -> decrease the doping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h    = 6.62606896E-34;              %% Planck constant [J.s]
hbar = h/(2*pi);
e    = 1.602176487E-19;             %% electron charge [C]
m0   = 9.10938188E-31;              %% electron mass [kg]
Epsi0= 8.854187817620E-12;          %% Vaccum dielectric constant [F/m]
kB   = 1.3806488E-23;               %% Boltzmann's constant [J/K]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloops = 30;                 % number of loops
tau0   = 30;                  % damping coefficiant for convergency
n      = 3;                   % number of quantum state solutions
ScF    = 0.05;                % scaling factor to plot the wave function [Without Dimension]
dz     = 0.5e-9;              % resolution of the grid [m]
T      = 300;                 % Temperature [Kelvin], react on the Fermi function and the population only

plot_density     = 1;         % Activate the plot 0 or 1
plot_convergence = 1;         % Activate the plot 0 or 1
plot_field       = 0;         % Activate the plot 0 or 1
plot_Vbending    = 0;         % Activate the plot 0 or 1
plot_charges     = 0;         % Activate the plot 0 or 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input_file;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% NOTHING TO CHANGE ANYMORE !!! %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Grabbing the parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zt   = M(:,2)*1e-9;         % conversion of the length from Angstrom to meter
Dopt = M(:,3)*1e18*1e6;     % n doping conversion from cm-3 to m-3
CBOt = M(:,1);              % Conduction Band Offset [eV]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Discretisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here, I descretize the grid z, the potential V0 and the values that are needed

z=0; V0=CBOt(1); Dop=Dopt(1);

for i=1:length(zt)
    t=zt(i);
    zv= (z(end)+dz): dz : (z(end)+dz)+t;
    z=[z zv];
    V0  = [ V0     ones(size(zv)) * CBOt(i)  ];
    Dop = [ Dop    ones(size(zv)) * Dopt(i)  ];
end

V0=V0-min(V0);   % Shift the band, but not utlimatly necessary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% indexes in which the Schrodinger equations are solved %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idxL = find(abs(z-zt(1)-zt(2))<dz/2,1);
idxR = find(abs(z-z(end)+zt(end-1)+zt(end))<dz/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Boundary conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EgL = 1.4;                              % bandgap of the layer on the left [eV]
EfL = Fermi3D_n(EgL,T,meff,Dopt(1));    % Fermi level on the left

EgR = 1.4;                              % bandgap of the layer on the right [eV]
EfR = Fermi3D_n(EgR,T,meff,Dopt(end));  % Fermi level on the right

EfL=EfL+V0(1);  EfR=EfR+V0(end);
Ef=z*0+EfL;
Fbi=(EfL-EfR)/(z(end)-z(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Electron Energy grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take care of the grid! It is optimized to use the minimum amount of points.
% If you do something special, you might have to increase the range and the
% resolution

Emin = min( min(EfL,EfR) , min(V0) ) - 0.2;
Emax = max(EfL,EfR) + 0.2;
En = linspace( Emin , Emax, 300 );
[ZZ,EEn]=meshgrid(z,En);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Meshgrid of densities matrices %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ro3D_const = (1/(2*pi^2)) * ( (2*e*meff*m0/(hbar^2)).^(3/2) ) + z*0;
ro3D_const(idxL:idxR) = ro3D_const(idxL:idxR)*0;

[ro3D_const_M] = meshgrid(ro3D_const,En); % put the vector Mass_n in a matrix En-long
ro3D = ro3D_const_M .* sqrt(  EEn );

ro2Dcst = e*meff*m0/(pi*(hbar)^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fermi level construction %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Ef_M,EEn]=meshgrid(Ef,En);
FEc = 1./(1+exp((EEn-Ef_M)/(kB*T/e))) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Starting of the Poisson s loop %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vs=z*0; Vsold=Vs; ntot=0; nloop=1;
ErrVec=1; sumVtotVec=1;

if Dopt==0
    Nloops=2;
    plot_density     = 0;
    plot_convergence = 0;
    plot_field       = 0;
    plot_Vbending    = 0;
    plot_charges     = 0;
end

while nloop<Nloops
    
    nloop
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    Vbending=Vs;
    Vbitot=V0+Vbending;
    tau = tau0*(1 + 2^((nloop - Nloops*0.8 )/10)); % tau will increase at each loop
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% matrix density calcul %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [Vbitot_M]=meshgrid(Vbitot,En);
   
    %%%%%%%%%%%%%%%%%%%%%%%%% 3D electrons density %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    idx  = ( EEn - Vbitot_M ) > 0;
    ro3D = ro3D_const_M .* sqrt( EEn -Vbitot_M).*idx;
           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Schrodinger solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    [Ec,psic] = Schroed1D_FEM_f(z(idxL:idxR),Vbitot(idxL:idxR),meff,n);
    PSI=zeros(length(z),n);
    PSI(idxL:idxR,:)=psic;
    
    for ii=1:length(Ec)
      PSI_M(:,:,ii)=meshgrid(PSI(:,ii),En);   
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%% 2D electrons density %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ro2D = zeros(length(En),length(z),length(Ec)) + ro2Dcst;
    for ii=1:length(Ec)
       ro2D( En<Ec(ii),:,ii) = 0;
       ro2D(:,:,ii) = ro2D(:,:,ii) .* abs(PSI_M(:,:,ii).^2);
    end
    
    ro2D = sum(ro2D,3);

    %%%%%%%%%%%%%%%%%%%%%%%%%% sum of all the densities %%%%%%%%%%%%%%%%%%%%%%%%
    
    ro = (ro3D + ro2D) .* FEc;       %% bulk + well
    Ntot    = trapz(En,ro);
    NtotX   = Ntot-Dop;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% double integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ntot = ntot + (NtotX-ntot)/tau;        % It add slowly the total number of electrons in order to converge
    
    F  = e*cumtrapz(z,ntot)./(Epsi0*Epsi); % Electrical Field
    MF = trapz(z,F)/(z(end)-z(1));         % MF is the mean(F) function on a nonlinear grid z
    F  = F - MF - Fbi ;
              
    Vsold=Vs;                              % storing the old value
    Vs  = -cumtrapz(z,F);                  % integal on a nonlinear grid z
           
    %%%%%%%%%%%%%%%%%%%%%%% Convergence analysis/plot %%%%%%%%%%%%%%%%%%%%%%%%%%

    Err = abs(  1 - sumVtotVec(end)/sum(Vs)  );
    sumVtotVec(nloop) = sum(Vs);
    ErrVec = [ErrVec Err];

    nloop=nloop+1;
    
    if Err<1e-10
       break 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Scaling and shifting the wavefunctions %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(Ec)
   PSIc(:,i)=abs(psic(:,i)).^2/max(abs(psic(:,i)).^2)*ScF + Ec(i); % normalisation for the plotting
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[10 100 1000 700],'color','w');
subplot(1,1,1,'fontsize',15)
hold on;grid on;box on;
col=colormap(jet);

if plot_density==1
    grid off
    pcolor(ZZ*1e9,EEn,ro*1e-6)
    set(gca,'color',col(1,:))
    shading flat
    hcb=colorbar;
    title(hcb,'\fontsize{8}cm-3')
    
    plot(z*1e9,V0,'b--','linewidth',1)
    plot(z*1e9,Vbitot,'w','linewidth',2)
    
elseif plot_density==0
    plot(z*1e9,V0,'b--','linewidth',1)
    plot(z*1e9,Vbitot,'b','linewidth',2)
end

for i=1:length(Ec)
    plot(z(idxL:idxR)*1e9,PSIc(:,i),'color','r','linewidth',1)
end

plot(z*1e9,Ef,'g-','linewidth',1)
text(z(1)*1e9,Ef(1)-0.01,'\color{green}Fermi')
text(z(end)*1e9*0.95,Ef(end)+0.01,'\color{green}Fermi')

xlabel('z (nm)')
ylabel('Energy (eV)')
xlim([z(1) z(end)]*1e9)
ylim([min(Vbitot)-0.1 max(Vbitot)+0.1])
title(strcat('\fontsize{12}T=',num2str(T),'K; meff=',num2str(meff),'; Epsilon=',num2str(Epsi),'; dz=',num2str(dz*1e9),'nm; dE=',num2str( (En(2)-En(1))*1e3,'%.1f'),'meV'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_convergence==1

    figure
    semilogy(1:nloop,ErrVec,'bo-')
    hold on; grid on;box on;
    xlabel('Cycles')
    ylabel('Convergence (norm. units)')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_field==1
    
    figure
    hold on;grid on;box on;
    [AX,H1,H2]=plotyy(z*1e9,F*1e-2*1e-3,z*1e9,Dop*1e-18*1e-6);
        
    set(H1,'color','r')
    set(H2,'color','b')
    
    xlabel('z (nm)')
    ylabel(AX(1),'E- field (kV/cm)','color','red')
    ylabel(AX(2),'Doping (1e18 cm-3)','color','blue')
    
    set(AX(1),'ycolor','red')
    set(AX(2),'ycolor','blue')
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_Vbending==1

    figure
    hold on; grid on;box on;

    [AX,H1,H2]=plotyy(z*1e9,Vbending,z*1e9,(ntot+Dop)*1e-18*1e-6);

    set(H1,'color','r')
    set(H2,'color','b')
    
    xlabel('z (nm)')
    ylabel(AX(1),'Vbending (eV)','color','red')
    ylabel(AX(2),'ntot (1e18 cm-3)','color','blue')
    
    set(AX(1),'ycolor','red')
    set(AX(2),'ycolor','blue')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_charges==1

    figure
    hold on; grid on;box on;
    
    plot(z*1e9,NtotX*1e-18*1e-6,'r')
    plot(z*1e9,ntot*1e-18*1e-6,'b--')
        
    xlabel('z (nm)')
    ylabel('Charges density (1e18 cm-3)')
    legend('NtotX','ntot')
    title('Simulations good when both curves are on each other')
    
    % From that graph, one can see if the convergence parameters are well 
    % set. In order to make the convergence, the charges are added 
    % progressively and not all are injected from the begining. As a 
    % results, it can be that the simulation finishes but not all the 
    % charges were injected.
    % NtotX is the full charge density while ntot is the one used during
    % the simulation. If both do not match, either nloop have to be
    % increased or tau0 has to be decreased.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%