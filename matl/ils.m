% resolution probleme non lineaire de diffusion
% schema implicite linearise
clear;
N=51; L=1; gamma=10;
dx=L/(N-1);X=[0:dx:L];
% parametre
K0=0.01; sigma=1; beta=300;
% nds internes
I=2:N-1;
% second membre
delta=0.2; Q=(X<delta)*beta;
% loi K
K=inline('t.^2'); 
% C.I.
Un=ones(1,N); dUn=zeros(1,N); Fn=zeros(1,N);
% iterations
nitmax=1000; epsilon=1.0e-5;
for it=1:nitmax
    % calcul du dt
    Umax=max(Un);
    dt=gamma*2/(4*sigma*Umax^3+4*K0*K(Umax)/dx^2);
    % calcul du 2nd membre
    Fn(I)=Un(I)+dt*(Q(I)+sigma);
    % C.L.
    Fn(1)=Un(1)+dt*(Q(1)+sigma);
    Fn(N)=1;
    % calcul de la matrice 3Diag 
    Kn12=0.5*(K(Un(1:N-1))+K(Un(2:N))); 
    A=[1+dt*(K0/dx^2*(2*Kn12(1))+sigma*Un(1)^3),1+dt*(K0/dx^2*(Kn12(I-1)+Kn12(I))+sigma*Un(I).^3),1];
    B=[-dt*K0/dx^2*(2*Kn12(1)), -dt*K0/dx^2*(Kn12(I)),0];
    C=[0,-dt*K0/dx^2*(Kn12(I-1)),0]; 
    J=[C;A;B];
    % resolution
    Un1=tridiag(J,Fn);
    % erreur
    Err=norm((Un1-Un)/dt)/sqrt(N);
    % fin
    Un=Un1; 
    if (Err<epsilon)  break; end;
  
end;