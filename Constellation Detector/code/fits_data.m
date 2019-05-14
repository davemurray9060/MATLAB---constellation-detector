vc=rfits('myfile.fits');
imagesc(vc.data); 
axis image

%fitswrite('many_stars.jpg','many_stars_fits');