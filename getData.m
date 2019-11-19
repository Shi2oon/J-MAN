function [X, Y, A, E, S, W, H, a, SF] = getData (Maps)
% function that get the displacement defromation graident as a matrix
% A(position, matrix), strain, X, Y, Map width (w), map hight (h) and crack
% distance in x (a)
X        = Maps.X(1,:)*1e3;    	Y = Maps.Y(:,1)*1e3; % in mm
A(:,:,1) = Maps.A11;      	A(:,:,2) = Maps.A12;            A(:,:,3) = Maps.A13;
A(:,:,4) = Maps.A21;      	A(:,:,5) = Maps.A22;            A(:,:,6) = Maps.A23;
A(:,:,7) = Maps.A31;      	A(:,:,8) = Maps.A32;            A(:,:,9) = Maps.A33;

E(:,:,1) = Maps.E11;       	E(:,:,2) = Maps.E22;            E(:,:,3) = Maps.E12;    
S(:,:,1) = Maps.S11;       	S(:,:,2) = Maps.S22;            S(:,:,3) = Maps.S12; 
W = max(X);               	H = max(Y)/2;           

imagesc(X,Y,log10(Maps.GND));        hold on
set(gca,'Ydir','normal');       axis equal;     axis tight;   
set(gca,'CLim',[14 15.5]);      c = colorbar;                  	
title({'\rho_G_N_D_s';''});
set(gcf,'position',[500,100,1200,800]);  colormap(jet)
xlabel('x[\mum]');          	ylabel('y[\mum]');
  
[a,~] = ginput(1);
close all;


SF = Maps.Stiffness; % in Pa