function st=findstars(im,mv,fw)

% FINDSTARS(fname,im,mv,fw) Finds sources in image "im" down to
% a peak level of "mv" counts. The parameter "fw" controls the size
% of the patch blanked around each source. Set it equal to or larger than 
% the radius of the largest sources for best results.
%im = im2double(im);
%im = fitsread('m67vc.fits');
%im = imread('many_stars.jpg');
im = im2double(im);
[v,ixm]=max(size(im(:)));
[a,b]=size(im);
[xx,yy]=meshgrid(1:b,1:a);
xmv=xx(ixm);
ymv=yy(ixm);
imr=im;
i=0;
bgk=median(im(:));
while (v>=(mv+bgk)),
    i=i+1;
    ix=find((xx>=(xmv-2*fw))&(xx<=(xmv+2*fw))&(yy>=(ymv-2*fw))&(yy<=(ymv+2*fw)));
    aa=length(find((xx(1,:)>=(xmv-2*fw))&(xx(1,:)<=(xmv+2*fw))));
    bb=length(find((yy(:,1)>=(ymv-2*fw))&(yy(:,1)<=(ymv+2*fw))));
    px=reshape(xx(ix),bb,aa);
    py=reshape(yy(ix),bb,aa);
    pz=reshape(im(ix),bb,aa);
    bg=median(pz(:));
    err=1;
    nt=0;
    while (err&&(nt<10)),
        fwt=(1+(rand-0.5)*0.3);
        [K,X,Y,Sj,Si,P,chi2,zm,err]=gauss2fit(px,py,pz-bg,v-bg,xmv,ymv,fwt,fwt,0,1e-3,-1);
        nt=nt+1;
    end
    ixb=(xx-xmv).^2+(yy-ymv).^2<2*fw^2;
    st(i,:)=[X,Y,K,2.35*Sj,2.35*Si,P];
    imr(ixb)=bg;
    fprintf('Source %d: ampl=%0.1f, x=%0.1f, y=%0.1f\n',i,K,X,Y);
    [v,ixm]=max(imr(:));
    xmv=xx(ixm);
    ymv=yy(ixm);
%    imagesc(imr); axis image; colorbar; drawnow
end
fprintf(1,'Total %d sources\n',i);
fprintf(1,'Returning table of: \n');
fprintf(1,'center X, center Y, amplitude, major axis, minor axis, position angle\n');
%fprintf(1,'See also file %s\n',im);
save('findstars_data');
clf
imagesc(xx(1,:),yy(:,1),im);
hold on
caxis([min(im(:)),0.05*(max(im(:))-min(im(:)))+min(im(:))])
plot(st(:,1),st(:,2),'ow','markersize',12,'linewidth',1.5);
hold off



