% resolution probleme non lineaire de diffusion
% schema Newton Ralphson
clear;
N=51; L=1;
dx=L/(N-1); X=[0:dx:L];
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
    % calcul du 2nd membre -Fi
    % nds internes
    Kn12=0.5*(K(Un(1:N-1))+K(Un(2:N))); 
    Fn(I)=K0/dx^2*(Kn12(I-1).*(Un(I-1)-Un(I))+Kn12(I).*(Un(I+1)-Un(I)))-sigma*(Un(I).^4-1.0)+Q(I);
    % C.L en 0
    Fn(1)=K0/dx^2*(Kn12(1)*(Un(2)-Un(1))+Kn12(1)*(Un(2)-Un(1)))-sigma*(Un(1).^4-1.0)+Q(1);
    Fn(N)=0;
    % calcul de la matrice Jacobienne
    A=[K0/dx^2*(2*Kn12(1))+4*sigma*Un(1)^3,K0/dx^2*(Kn12(I-1)+Kn12(I))+4*sigma*Un(I).^3,1];
    B=[-K0/dx^2*(2*Kn12(1)), -K0/dx^2*(Kn12(I)),0];
    C=[0,-K0/dx^2*(Kn12(I-1)),0]; 
    J=[C;A;B];
    % resolution
    dUn=tridiag(J,Fn);
    % erreur
    Err=norm(Fn)/sqrt(N);
    % fin
    Un=Un+dUn;  
    if (Err<epsilon) break; end;
end;