function scolor(fout,fname1,xy1,fname2,xy2)

% SCOLOR(fout,fname1,xy1,fname2,xy2) Compute the stellar colors from
% two photometry files produced by "pmerge" or "all_aper"


f1=load(fname1);
f2=load(fname2);

xs=f2(:,1)+xy1(1)-xy2(1);
ys=f2(:,2)+xy1(2)-xy2(2);
u=zeros(length(f2(:,1)),1);
k=0;
for i=1:length(f1(:,1)),
    if (f1(i,5)),           % if phot1 valid
        d=sqrt((f1(i,1)-xs).^2+(f1(i,2)-ys).^2);
        [dm,ix]=sort(d);    % find nearest match
        id=ix(1);
        if (dm(1)<3),       % if closer than 3 pixels
            u(id)=1;        % mark as found in phot2
            if (f2(id,5)),  % if phot2 valid, compute color
                k=k+1;
                fm(k,1)=(f1(i,1)+xs(id))/2;   % X
                fm(k,2)=(f1(i,2)+ys(id))/2;   % Y
                fm(k,3)=-2.5*log10(f1(i,3));  % mag1
                fm(k,4)=-2.5*log10(f2(id,3)); % mag2
                fm(k,5)=fm(k,3)-fm(k,4);      % color
                f1p(k)=i;   % keep track of where it comes from
                f2p(k)=id;
            end
        end
    end
end
fprintf(1,'Found valid colors for %d sources\n',k);
                
fo=fopen(fout,'w');
fprintf(fo,'%%  X     \t   Y    \t    Mag1 \t   Mag2  \t    Color \t File1 \t File2 \n');
for i=1:k,
    fprintf(fo,'%8.2f \t %8.2f \t %8.2f \t %8.2f \t %8.2f \t %d \t %d \n',fm(i,1),fm(i,2),fm(i,3),fm(i,4),fm(i,5),f1p(i),f2p(i));
end
fclose(fo);

clf;
plot(fm(:,5),fm(:,3),'o');
xlabel('Mag1-Mag2');
ylabel('Mag1');
set(gca,'ydir','rev');
