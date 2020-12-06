function[Ef]=Fermi3D_n(Eg,T,meff,N3D)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h    = 6.62606896E-34;              %% Planck constant [J.s]
hbar = h/(2*pi);
e    = 1.602176487E-19;             %% electron charge [C]
m0   = 9.10938188E-31;              %% electron mass [kg]
kB   = 1.3806488E-23;               %% Boltzmann's constant [J/K]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if N3D==0
   Ef=-Eg/2;
   return
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Fermi level at T=0K %%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ef0 = ( hbar^2 / (2*e*meff*m0) )  *  (3*pi^2* N3D )^(2/3);
if T==0
    Ef=Ef0;
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%% 3D electrons density %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here, I try to optimize the meshing
if Ef0<0.3
    Emax=0.3;
else
    Emax=Ef0+0.3;
end
En1=linspace(0,Emax,1e4);
En2=linspace(Emax,3,50);
En=sort([En1 En2]);

ro3Dn = (1/(2*pi^2))*( (2*e*meff*m0/(hbar^2))^(3/2) ) * sqrt(En) ;

%%%%%%%%%%%%%%%%%%%%% Fermi level at any temperature %%%%%%%%%%%%%%%%%%%%%%

Ef=Ef0;
FEc= 1./(1+exp((En-Ef)/(kB*T/e))) ;
ro3DEfn=ro3Dn.*FEc;
NtotX=trapz(En,ro3DEfn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the Fermi level is obiously at T=0K the max possible for n doped
% so now, it scans down with big step (ddE) 
% to find the Fermi level at T different from zero

ddE=0.005;   % step of the scan in eV to pass below the Fermi level
 
while ( ((NtotX - N3D)) > 0)
    
    Ef = Ef - ddE;    
    FEc= 1./(1+exp((En-Ef)/(kB*T/e))) ;

    ro3DEfn = ro3Dn.*FEc;
    NtotX = trapz(En,ro3DEfn);
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now, it will try to get as close as posible to the real Ef with an
% error of Epsilon by dichotomy

Ef1 = Ef;
Ef2 = Ef1+ddE;
Epsilon = 1e-10;
            
while ( (abs((NtotX - N3D)/N3D)) > Epsilon)  % find the Fermi level at any temperature
   
    if ( NtotX > (N3D))
        Ef  = Ef - abs(Ef1-Ef2)/2 ;  
        Ef1 = Ef ;
        FEc = 1./(1+exp((En-Ef)/(kB*T/e))) ;
    else
        Ef  = Ef + abs(Ef1-Ef2)/2 ;
        Ef2 = Ef ;
        FEc = 1./(1+exp((En-Ef)/(kB*T/e))) ;
    end

    ro3DEfn = ro3Dn.*FEc;
    NtotX = trapz(En,ro3DEfn);
    
end

if Ef<-Eg/2
   Ef=-Eg/2;
   return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure
% hold on;grid on;box on;
% 
% plot(En,ro3Dn*1e-6,'b')
% plot([1 1]*Ef,[0 max(ro3DEfn)]*1e-6,'g')
% plot([1 1]*Ef0,[0 max(ro3DEfn)]*1e-6,'g--')
% plot(En,ro3DEfn*1e-6,'ro-')
% 
% xlabel('Energy (eV)')
% ylabel('States density (cm-3.eV-1)')
% title(strcat('NtotX=',num2str(NtotX*1e-6,'%.2e'),'cm-3; Ef=',num2str(Ef),'eV'))

end
