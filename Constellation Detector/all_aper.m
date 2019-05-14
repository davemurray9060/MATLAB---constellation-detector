function [flux,ferr]=all_aper(fout,fpos,im)

% ALL_APER(fout,fpos,im) Do aperture photometry of image "im" using a star 
% position file produced by FINDSTARS

st=load(fpos); % File with star positions

Kccd=16;    % Gain of the CCD
inrad=10;   % Inner radius for sky annulus
outrad=13;  % Outer radius for sky annulus
saturation=20000; % Maximum counts in ADU before image is saturated
photrad=median(st(:,4)); % Aperture radius for photometry 

flux=zeros(length(st),1);
ferr=flux;
fp=fopen(fout,'w');
fprintf(fp,'%%  X     \t   Y    \t   Flux \t  Error  \t Flag \n');
for i=1:length(st),
    [flx,err]=aper(im,st(i,1),st(i,2),photrad,inrad,outrad,Kccd,saturation);
    flux(i)=flx;
    ferr(i)=err;
    mess=sprintf('Photometry good (1) or bad (0) for source %d?',i);
    r=input(mess);
    if (isnan(err)), r=0; end
    fprintf(fp,'%8.2f \t %8.2f \t %8.2f \t %8.2f \t %d\n',st(i,1),st(i,2),flx,err,r);
end
fclose(fp);

    
