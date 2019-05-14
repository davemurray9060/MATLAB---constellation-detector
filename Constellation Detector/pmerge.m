function pmerge(fout,fname1,xy1,fname2,xy2)

% PMERGE(fout,fname1,xy1,fname2,xy2) Merge two photometry files "fname1" and
% "fname2" produced by ALL_APER. A common star needs to be identified.
% The coordinates of that star in image 1 are specified in "xy1".
% Similarly, its coordinates in image 2 are specified in "xy2". If both
% images have different exposures their integration times should be 
% specified in as the third value in "xy1" and "xy2".


f1=load(fname1);
f2=load(fname2);

xs=f2(:,1)+xy1(1)-xy2(1);
ys=f2(:,2)+xy1(2)-xy2(2);
if (length(xy1)>2),
    fs1=100/xy1(3); % scaling factors to rescale flux to a common
    fs2=100/xy2(3); % integration time of 100s
else
    fs1=1;
    fs2=1;
end
u=zeros(length(f2(:,1)),1);
k=0;
for i=1:length(f1(:,1)),
    if (f1(i,5)),           % if phot1 valid
        k=k+1;
        d=sqrt((f1(i,1)-xs).^2+(f1(i,2)-ys).^2);
        [dm,ix]=sort(d);    % find nearest match
        id=ix(1);
        if (dm(1)>=3),      % if farther than 3 pixels
            fm(k,:)=f1(i,:);% assign a new source
	    fm(k,[3,4])=fm(k,[3,4])*fs1; % rescale if necessary
            f1p(k)=i;       % keep track of where it comes from
            f2p(k)=0;
        else
            u(id)=1;        % mark as found in phot2
            if (f2(id,5)),  % if phot2 valid, average with phot1
                fm(k,1)=(f1(i,1)+xs(id))/2;
                fm(k,2)=(f1(i,2)+ys(id))/2;
                fm(k,3)=(f1(i,3)/(fs1*f1(i,4)^2)+f2(id,3)/(fs2*f2(id,4)^2))/(1/(fs1*f1(i,4))^2+1/(fs2*f2(id,4))^2);
                fm(k,4)=1/sqrt(1/(fs1*f1(i,4))^2+1/(fs2*f2(id,4))^2);
                fm(k,5)=1;
                f1p(k)=i;   % keep track of where it comes from
                f2p(k)=id;
            else            % if phot2 invalid, ignore
                fm(k,:)=f1(i,:);
		fm(k,[3,4])=fm(k,[3,4])*fs1;
                f1p(k)=i;
                f2p(k)=0;
            end
        end
    end
end
ix=find((u==0)&(f2(:,5)>0)); % list of sources only in phot2 and valid
for i=1:length(ix),
    fm(k+i,:)=f2(ix(i),:);
    fm(k+i,[3,4])=fm(k+i,[3,4])*fs2;
    f1p(k+i)=0;
    f2p(k+i)=ix(i);
end
fprintf(1,'Found %d duplicated sources valid in both files\n',sum(u));
fprintf(1,'Found %d valid sources unique to %s\n',k,fname1);
fprintf(1,'Found %d valid sources unique to %s\n',length(ix),fname2);
                
fo=fopen(fout,'w');
fprintf(fo,'%%  X     \t   Y    \t   Flux \t  Error  \t Flag \t File1 \t File2 \n');
for i=1:length(fm(:,1)),
    fprintf(fo,'%8.2f \t %8.2f \t %8.2f \t %8.2f \t %d \t %d \t %d \n',fm(i,1),fm(i,2),fm(i,3),fm(i,4),fm(i,5),f1p(i),f2p(i));
end
fclose(fo);
