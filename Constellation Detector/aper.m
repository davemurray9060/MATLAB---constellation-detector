 function [flx,err]=aper(im,col,row,rad,ir,or,Kccd,saturation)

% APER(im,col,row,rad,ir,or,Kccd) Do aperture photometry of image "im"
% for a star centered at the "row,col" coordinates, using aperture of
% radius "rad" and sky annulus of inner/outer radii "ir/or" with CCD
% gain of Kccd ADU/electron. Optionally, an 8th parameter can be passed
% with the saturation value for the CCD.

if (nargin==7), saturation=Inf; end

[a,b]=size(im);
[xx,yy]=meshgrid(1:b,1:a);

ixsrc=(xx-col).^2+(yy-row).^2<=rad^2;
ixsky=((xx-col).^2+(yy-row).^2<=or^2)&((xx-col).^2+(yy-row).^2>=ir^2);

sky=median(im(ixsky));                        % sky value
pix=im(ixsrc)-sky;                            % source without sky
sig=sqrt(im(ixsrc)/Kccd);                     % photon noise per pixel
ssig=std(im(ixsky))/sqrt(length(ixsky))/Kccd; % sky noise in average
flx=sum(pix)/Kccd;                            % flux
err=sqrt(sum(sig).^2+ssig^2);                 % total error
if (max(im(ixsrc))>saturation), err=NaN; end

fw=or;
ix=find((xx>=(col-2*fw))&(xx<=(col+2*fw))&(yy>=(row-2*fw))&(yy<=(row+2*fw)));
aa=length(find((xx(1,:)>=(col-2*fw))&(xx(1,:)<=(col+2*fw))));
bb=length(find((yy(:,1)>=(row-2*fw))&(yy(:,1)<=(row+2*fw))));
px=reshape(xx(ix),bb,aa);
py=reshape(yy(ix),bb,aa);
pz=reshape(im(ix),bb,aa);
clf;
imagesc(px(1,:),py(:,1),pz);
if (~isempty(im(ixsrc))),
  caxis([sky max(im(ixsrc))]);
end
axis image
p=(0:360)*pi/180;
xc=cos(p);
yc=sin(p);
hold on
plot(col+rad*xc,row+rad*yc,'w');
plot(col+ir*xc,row+ir*yc,'r');
plot(col+or*xc,row+or*yc,'r');
if (isnan(err)),
    ht=text(col,row,'SATURATED');
    set(ht,'horizontalalign','center','color','w','fontweig','bold');
end
hold off



