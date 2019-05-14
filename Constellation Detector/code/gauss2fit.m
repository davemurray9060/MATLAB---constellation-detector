function [K,X0,Y0,Sx,Sy,Phi,chi2,zmodel,error]=gauss2fit(x,y,z,K,X0,Y0,Sx,Sy,Phi,tol,graphout,mask)

% [K,X0,Y0,Smj,Smi,Phi,chi2,err]=
%			gauss2fit(x,y,z,K,X0,Y0,Smj,Smi,Phi,tol,graphout,mask)
%
% Fits a 2-D Gaussian to the (x,y,z) dataset using the Levenberg-Marquart
% algorithm. The parameters for the fit are the height of the Gaussian (K), 
% its center (X0,Y0), the standard deviations along the major and minor axis 
% (Smj,Smi), and the position angle of the major axis (Phi). Other inputs are 
% the tolerance of the fit (tol), defined as the percentual change in the 
% sum of the square of the residuals (chi2). The fiting process is shown 
% graphically using contour-plots when the input graphout is nonzero.
% The mathematical expression for the 2-D gaussian used here is
%
% model = K exp(-1/2 ((+(x-X0) cos(Phi) + (y-Y0) sin(Phi))^2/Sx^2 + ...
%                     (-(x-X0) sin(Phi) + (y-Y0) cos(Phi))^2/Sy^2))
%
% An initial guess has to be supplied for starting the algorithm.
% When the initial guess is not good enough and the algorithm fails to
% converge to the specified tolerance a flag is returned (error==1).

if (nargin==11),
	mask=ones(1,6);
end
if (any(mask==0)),
	disp('Attempting constrained gaussian fit:');
	if (mask(1)==0), disp('	- Fixed peak amplitude'); end
	if (mask(2)==0), disp('	- Fixed center X coordinate'); end
	if (mask(3)==0), disp('	- Fixed center Y coordinate'); end
	if (mask(4)==0), disp('	- Fixed size in X'); end
	if (mask(5)==0), disp('	- Fixed size in Y'); end
	if (mask(6)==0), disp('	- Fixed orientation'); end
end
PA=Phi*pi/180;
L=0.001;
peak=max(max(z));
levs=[1:5]*peak/5;
%% Inverse squared weighting matrix...
%W=abs(z)/sum(sum(abs(z)));
W=ones(size(z)); W=W/sum(sum(W));

A=+(x-X0).*cos(PA)+(y-Y0).*sin(PA);
B=-(x-X0).*sin(PA)+(y-Y0).*cos(PA);
zmodel=K*exp(-0.5*((A/Sx).^2+(B/Sy).^2));
delta=(z-zmodel);
deltaw=delta.*W;
chi2=sum(sum(W.*delta.^2));

DK= zmodel/K;
DSx=zmodel.*A.^2/Sx^3;
DSy=zmodel.*B.^2/Sy^3;
DX0=zmodel.*(A*cos(PA)/Sx^2-B*sin(PA)/Sy^2);
DY0=zmodel.*(A*sin(PA)/Sx^2+B*cos(PA)/Sy^2);
DPA=zmodel.*(A.*B)*(1/Sx^2+1/Sy^2);

betaK= sum(sum(deltaw.*DK))*mask(1);
betaX0=sum(sum(deltaw.*DX0))*mask(2);
betaY0=sum(sum(deltaw.*DY0))*mask(3);
betaSx=sum(sum(deltaw.*DSx))*mask(4);
betaSy=sum(sum(deltaw.*DSy))*mask(5);
betaPA=sum(sum(deltaw.*DPA))*mask(6);

beta=[betaK;betaSx;betaSy;betaX0;betaY0;betaPA];

hold off
if (graphout>-1)
  contour(x,y,z,levs,'-g')
  axis('square')
  hold on
end
if (graphout>0),
	contour(x,y,zmodel,levs,'-r')
	contour(x,y,zmodel,[0,K/2],'--r')
	plot(X0,Y0,'xr')
	drawnow
end

converge=0;
error=0;
while (~converge),

	alpha=[sum(sum(W.*DK.*DK))*(1+L),sum(sum(W.*DK.*DSx)),sum(sum(W.*DK.*DSy)),...
	       sum(sum(W.*DK.*DX0)),sum(sum(W.*DK.*DY0)),sum(sum(W.*DK.*DPA))];
	alpha=[alpha;alpha(1,2),sum(sum(W.*DSx.*DSx))*(1+L),sum(sum(W.*DSx.*DSy)),...
	       sum(sum(W.*DSx.*DX0)),sum(sum(W.*DSx.*DY0)),sum(sum(W.*DSx.*DPA))];
	alpha=[alpha;alpha(1,3),alpha(2,3),sum(sum(W.*DSy.*DSy))*(1+L),...
	       sum(sum(W.*DSy.*DX0)),sum(sum(W.*DSy.*DY0)),sum(sum(W.*DSy.*DPA))];
	alpha=[alpha;alpha(1,4),alpha(2,4),alpha(3,4),...
	       sum(sum(W.*DX0.*DX0))*(1+L),sum(sum(W.*DX0.*DY0)),sum(sum(W.*DX0.*DPA))];
	alpha=[alpha;alpha(1,5),alpha(2,5),alpha(3,5),...
	       alpha(4,5),sum(sum(W.*DY0.*DY0))*(1+L),sum(sum(W.*DY0.*DPA))];
	alpha=[alpha;alpha(1,6),alpha(2,6),alpha(3,6),...
	       alpha(4,6),alpha(5,6),sum(sum(W.*DPA.*DPA))*(1+L)];
	
	if (rcond(alpha)<1e-15),
		disp('Ill conditioned data')
		error=2;
		return
	end

	da=inv(alpha)*beta;
	newK =K +da(1);
	newSx=abs(Sx+da(2));
	newSy=abs(Sy+da(3));
	newX0=X0+da(4);
	newY0=Y0+da(5);
	newPA=PA+da(6);

	newA=+(x-newX0).*cos(newPA)+(y-newY0).*sin(newPA);
	newB=-(x-newX0).*sin(newPA)+(y-newY0).*cos(newPA);
	newzmodel=newK*exp(-0.5*((newA/newSx).^2+(newB/newSy).^2));
	newdelta=(z-newzmodel);
	newdeltaw=newdelta.*W;
	newchi2=sum(sum(W.*newdelta.^2));

	if (newchi2>=chi2),
		L=L*10;
		converge=(L>100000);
		if (converge),
			error=1;
			disp('Exceeding L limit');
		end
	else
		converge=((chi2-newchi2)/newchi2);
		converge=(converge<tol);

		L=L/10;
		K=newK;
		Sx=newSx;
		Sy=newSy;
		X0=newX0;
		Y0=newY0;
		PA=newPA;
		A=newA;
		B=newB;
		zmodel=newzmodel;
		chi2=newchi2;
		deltaw=newdeltaw;

		DK= zmodel/K;
		DSx=zmodel.*A.^2/Sx^3;
		DSy=zmodel.*B.^2/Sy^3;
		DX0=zmodel.*(A*cos(PA)/Sx^2-B*sin(PA)/Sy^2);
		DY0=zmodel.*(A*sin(PA)/Sx^2+B*cos(PA)/Sy^2);
		DPA=zmodel.*(A.*B)*(1/Sx^2+1/Sy^2);

		betaK= sum(sum(deltaw.*DK))*mask(1);
		betaX0=sum(sum(deltaw.*DX0))*mask(2);
		betaY0=sum(sum(deltaw.*DY0))*mask(3);
		betaSx=sum(sum(deltaw.*DSx))*mask(4);
		betaSy=sum(sum(deltaw.*DSy))*mask(5);
		betaPA=sum(sum(deltaw.*DPA))*mask(6);

		beta=[betaK;betaSx;betaSy;betaX0;betaY0;betaPA];
	end

	if (graphout>0),
		hold off
		contour(x,y,z,levs,'-g')
		axis('square')
		hold on
		contour(x,y,z,[0,peak/2],':g')
		contour(x,y,zmodel,levs,'-r')
		contour(x,y,zmodel,[0,K/2],'--r')
		plot(X0,Y0,'xr')
		drawnow
	end

end

Phi=PA*180/pi;
if (Sx<Sy),
	Phi=90+Phi;
	S=Sy;
	Sy=Sx;
	Sx=S;
end
if (Phi<0),
	Phi=Phi+180;
end

% "Normalizing" the Chi^2...

[a,b]=size(x);
chi2=chi2/(a*b);

if (graphout>-1)
clf
pcolor(x,y,z-zmodel); shading flat; hold on
contour(x,y,z,[0,peak/2],':g')
contour(x,y,zmodel,levs,'-r')
contour(x,y,zmodel,[0,K/2],'--r')
plot(X0,Y0,'xr')
a=axis;
xlen=a(2)-a(1);
ylen=a(4)-a(3);
fc=sqrt(log(256));
text(a(1)+0.05*xlen,a(4)-0.05*ylen,['K=' num2str(K)]);
text(a(1)+0.05*xlen,a(4)-0.10*ylen,['FWHMmj=' num2str(fc*Sx)]);
text(a(1)+0.05*xlen,a(4)-0.15*ylen,['FWHMmi=' num2str(fc*Sy)]);
text(a(1)+0.05*xlen,a(4)-0.20*ylen,['X0=' num2str(X0)]);
text(a(1)+0.05*xlen,a(4)-0.25*ylen,['Y0=' num2str(Y0)]);
text(a(1)+0.05*xlen,a(4)-0.30*ylen,['Phi=' num2str(Phi)]);
text(a(1)+0.05*xlen,a(4)-0.35*ylen,['Chi2=' num2str(chi2)]);
drawnow;
end
if ((error==1)&(graphout>-1)),
	text(a(1)+0.05*xlen,a(4)-0.40*ylen,'Bad fit: L limit');
end

return

