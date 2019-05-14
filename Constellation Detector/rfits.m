function im=rfits(file,bval,bval2)

% RFITS	Read in FITS
%	A simple FITS image reader.
%
%	RFITS('file') reads a file written in FITS format and returns
%	all the relevant information in a structure, including the header
%	values, the data, and coordinate arrays for each axis.
%
%	RFITS('file',bval) fills replaces the blank pixels in the image with
%	the value bval. If bval is not specified, it defaults to NaN.
%
%	RFITS('file','options'), where the second parameter is a string
%	rather than a number, allows you to control the workings of RFITS.
%	The currently implemented options are:
%	- head or header: return just the header information, do not read
%	        	  in the data (faster and less memory-intensive than
%			  reading the whole thing).
%	- nowcs: do not produce "World Coordinate System" matrices, which
%		 is very time consuming and only important if the image or
%		 cube are big enough for the type of projection to be 
%		 important.
%   - xten#: specify which extension number to read (right now just one
%            digit
%
%	RFITS('file','options',bval) specifies a blanking value.


rlen=2880;	% FITS record length
if (nargin==1),
	bval=NaN;
end
if (ischar(bval)),
    head=~isempty(findstr(bval,'head'));
    nowcs=~isempty(findstr(bval,'nowcs'));
    silent=~isempty(findstr(bval,'silent'));
    silent=isempty(findstr(bval,'verbose'));
    ix=~isempty(findstr(bval,'xten'));
    if (ix>0),
        xten=str2double(bval(ix+4));
    else
        xten=0;
    end
    if (nargin==1),
        bval=NaN;
    elseif (nargin>2),
        bval=bval2;
    end
else
    head=0;
	nowcs=0;
    xten=0;
    silent=1;
end
h=fopen(file,'r','b'); % Big endian IEEE format
if (h<0),
	error('File not found');
end
if (fscanf(h,'SIMPLE =%s')~='T'),
	error('File is not in FITS format');
end
xti=0;
st=0;
if (xten),
    nr=0;
    xt=[];
    im=rfits(file,'head'); % now read header for entire file even if soliciting xtension
    while (isempty(xt)&&(st==0)),
        nr=nr+1;
        st=fseek(h,nr*rlen,-1);
        xt=fscanf(h,'XTENSION=%s');
    end
    im.xtension=xt;
    xti=rlen*nr;
end
if (st<0),
    error(ferror(h));
end    
fseek(h,xti+80,-1);
bitpix=fscanf(h,'BITPIX =%d');
if isempty(bitpix),
	error('BITPIX keyword not found');
end
fseek(h,xti+2*80,-1);
naxis=fscanf(h,'NAXIS =%d');
if isempty(naxis),
	error('NAXIS keyword not found');
end
im.naxis=naxis;
im.bitpix=bitpix;
im.numpt=zeros(1,naxis);
im.crval=zeros(naxis,1);
im.crpix=zeros(naxis,1);
im.cdelt=ones(naxis,1);
im.crota=zeros(naxis,1);
im.ctype=cell(naxis,1);
im.cunit=cell(naxis,1);
im.bscale=1;
im.bzero=0;
im.bunit='';
im.blank=NaN;
npt=1;
nh=1;
nc=1;
ln=0;
for ax=1:naxis,
	ln=2+ax;
	fseek(h,xti+ln*80,-1);
	fmt=['NAXIS' num2str(ax) ' =%d'];
	im.numpt(ax)=fscanf(h,fmt,1);
	npt=npt*im.numpt(ax);
	if isempty(im.numpt(ax)),
		error(['NAXIS' num2str(ax) ' keyword not found']);
	end
end
flg=1;
while (flg),                                    % Header parsing (new, AB 11/21/2007)
	ln=ln+1;
	fseek(h,xti+ln*80,-1); 			% Reposition file at start of new line
    lnstr=fscanf(h,'%c',80); 	 	% Read line in
    field=lower(deblank(lnstr(1:8)));
	ixf=findstr(field,'-');
	field(ixf)='_';
	fend=findstr(lnstr,'/'); 		% Find scanning range (look for comment)
	if (isempty(fend)),
	   fend=80;
	else
	   fend=fend(end)-1;			% Sometimes there can be more than one '/'
	end
	vstr=lnstr(10:fend);
	ixq=findstr(vstr,''''); 		% Is it a string or a number?
	if (length(ixq)>1),
	   vstr=vstr((ixq(1)+1):(ixq(2)-1));
	   value=deblank(vstr);  		% "value" always contains the field value
       if (~isempty(value)),
           if (value(end)=='&'),        % Some strings are continued on the next line using this convention
               if (iscontin==0), cfield=field; end
               iscontin=1;
               value=value(1:(end-1));  % drop trailing '&'
           else
               iscontin=0;
           end
       end
	else
	   value=sscanf(vstr,'%f',1);
	end
	froot=field(1:(end-1));
    if (~isempty(strmatch(froot,{'crval','crpix','cdelt','crota'},'exact'))), % Special treatment for these fields
        ax=str2double(field(6));
        im.(field(1:5))(ax)=value;		
    elseif (~isempty(strmatch(froot,{'ctype','cunit'},'exact'))), % CTYPE,CUNIT are cell arrays
        ax=str2double(field(6));
        im.(field(1:5)){ax}=value;		
    elseif (strcmp(field,'history')),       % HISTORY also gets special treatment (can be repeated)
        im.history{nh}=deblank(lnstr(9:end));
        nh=nh+1;
    elseif (strcmp(field,'comment')),       % COMMENT is like history
        im.comment{nc}=deblank(lnstr(9:end));
        nc=nc+1;
	elseif (strcmp(field,'end')),		% END of course signals the end
        flg=0;
    elseif (strcmp(field,'continue')),  % CONTINUE statement, append to previous value
        im.(cfield)=[im.(cfield),value];
	elseif (~isempty(field)),		% otherwise just assign field
        im.(field)=value;
    end
end

if (isfield(im,'cd1_1')),
   if (~isfield(im,'cd1_2')), % force define cd1_2 for some funky headers
      	im.cd1_2=0;
   end
   if (~isfield(im,'cd2_1')),
        im.cd2_1=0;
   end
end

if (head),              			% Return now if only requesting header
  fclose(h);
  return; 
end

im.x=cell(1,naxis);
if (isfield(im,'cd1_1')),                % CD format for distortion coeffs.
    xv1=(1:im.numpt(1))-im.crpix(1);
    xv2=(1:im.numpt(2))-im.crpix(2);
    [xx,yy]=meshgrid(xv1,xv2);
    im.x(1)={im.cd1_1*xx+im.cd1_2*yy};
    im.x(2)={im.cd2_1*xx+im.cd2_2*yy};
    dcp=1;
elseif (isfield(im,'pc1_1')),            % PC format for distortion coeffs.
    xv1=(1:im.numpt(1))-im.crpix(1);
    xv2=(1:im.numpt(2))-im.crpix(2);
    [xx,yy]=meshgrid(xv1,xv2);
    im.x(1)={im.cdelt(1)*(im.pc1_1*xx+im.pc1_2*yy)};
    im.x(2)={im.cdelt(2)*(im.pc2_1*xx+im.pc2_2*yy)};
    dcp=1;
else                                    % usual axis definition
    for i=1:naxis,
        im.x(i)={((1:im.numpt(i))-im.crpix(i))*im.cdelt(i)};
    end
    dcp=0;
end

if ((~nowcs)&&(~isempty(im.ctype))),
ctype=char(im.ctype{1});
if (length(ctype)>=8),
  if (~silent), disp('Projection found, working on WCS...'); end
  if (dcp==0),
    x1=double(im.x{1});
    x2=double(im.x{2});  
% speed up calculations for big arrays
% AB 3/24/2004
    [yy,xx]=meshgrid(x2,x1);
  else
    xx=double(im.x{1});
    yy=double(im.x{2});
  end
  switch (ctype(6:8)),
    case 'SIN',
	phi=atan2(yy,xx)+pi/2;
	theta=acos(pi/180*sqrt(xx.^2+yy.^2));
	wcs=1;
    case 'TAN',
	phi=atan2(yy,xx)+pi/2;
	theta=atan(180./(pi*sqrt(xx.^2+yy.^2)));	
	wcs=1;
    otherwise,
    disp(['Unsupported projection ' ctype]);
	wcs=0;
  end
else
  wcs=0;
end
if (wcs),
% make up some space for big arrays
% AB 3/24/2004
    clear xx yy	
    dc=im.crval(2)*pi/180;
    ac=im.crval(1)*pi/180;
    sind=sin(theta).*sin(dc)-cos(theta).*cos(phi).*cos(dc);
    cosdsina=cos(theta).*sin(phi);
    cosdcosa=sin(theta).*cos(dc)+cos(theta).*cos(phi).*sin(dc);
    im.wcs1=180/pi*(ac+atan2(cosdsina,cosdcosa));
    im.wcs2=180/pi*atan2(sind,sqrt(cosdsina.^2+cosdcosa.^2));
end
end
for i=1:naxis,
    im.x(i)={im.crval(i)+im.x{i}};
end

cp=ftell(h);
np=rlen*ceil(cp/rlen);
fseek(h,np,-1);		% Position at the begining of the data records

if (bitpix==32),
	prec='int32';
elseif (bitpix==16),
	prec='int16';
elseif (bitpix==-32),
	prec='float32';
elseif (bitpix==-64),
	prec='float64';
else
	prec='uint8';
end

nval=ceil(npt);
if (~silent), fprintf(1,'Header read. Getting %d data points...',nval); end
if (length(im.numpt)>2),
    [v,cnt]=fread(h,nval,prec);
    im.data=reshape(v,im.numpt);
    clear v
elseif (isempty(im.numpt)),
    return
else
    [im.data,cnt]=fread(h,im.numpt,prec);
end
if (~silent), disp('Data read...'); end
if (cnt<nval),
	disp(['File too short at pixel ' num2str(cnt)]);
else
    if (~ischar(bval)), 
        if (~silent), disp('Finding blanks in cube...'); end
        ix=find(im.data==im.blank); 
        if (~isnan(im.blank)), ix=[ix;find(isnan(im.data))]; end
    end
    if ((im.bscale~=1)||(im.bzero~=0)),
        if (~silent), disp('Rescaling data according to header...'); end
        im.data=im.data*im.bscale+im.bzero;
    end
	if (~ischar(bval)),
        if (~silent), disp('Reassigning blanking values...'); end
        im.data(ix)=bval;
        im.blank=bval;
    end
%    im.data=reshape(v,im.numpt);
end

fclose(h);
return

