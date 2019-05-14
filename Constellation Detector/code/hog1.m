function h=hog(I,list,nbins)
r = 4;
nI=zeros(ceil(size(I)/r)*r);      % alloc space for nI with zeros
nI(1:size(I,1),1:size(I,2)) = I;  % set nI(1:#rows,1:#col) = [I] pad-matrix
hx = [1; 0; -1];
hy = hx';
dx = imfilter(nI,hx);           % sliding window thru x
size(dx)
dy = imfilter(nI,hy);           % sliding window thru y
size(dy)
angle = mod(atan2(dy,dx),2*pi);   % direction of gradient
rho = ((dy.^2)+(dx.^2)).^.5+0.02;      % magnitude of difference angle

negative_angles = angle<=0; % negative values
angle(negative_angles) = mod(angle(negative_angles)+2*pi,2*pi); % [0,2*pi]

binsize = (2*pi)/(nbins);   % set size of each theta-bin for angles, theta--> binsize
theta = mod(angle+binsize/2,2*pi); % 128 x 171
tbin = ceil(theta/binsize);
theta(isnan(theta))=0;
rho(isnan(rho))=0;
%%%%%% L2-normalization
rho=rho/sqrt(norm(rho)^2);
for z=1:length(rho)
  if rho(z)>0.4
   rho(z)=0.4;
  end
end

rhos = zeros([size(tbin) nbins]);
[I,J] = meshgrid(1:size(tbin,1), 1:size(tbin,2));
K = sub2ind(size(tbin),I(:),J(:));
%size(rho)
rhos(sub2ind(size(rhos),I(:),J(:),tbin(K))) = rho(sub2ind(size(rho),I(:),J(:)));  % 128 x 171 x 16
%size(rhos)
[x,y] = meshgrid(-r:r,-r:r);    % x: 9x9: -4 to 4    y: 9x9: -4 to 4
xi = reshape(kron(list(:,1),ones(numel(x),1)),[size(x),size(list,1)]); % 9 x 9 x 100
yi = reshape(kron(list(:,2),ones(numel(y),1)),[size(y),size(list,1)]); % 9 x 9 x 100
xi = xi+repmat(x,[1,1,size(list,1)]); % 9 x 9 x 100
yi = yi+repmat(y,[1,1,size(list,1)]); % 9 x 9 x 100
xi = reshape(xi,[numel(x),size(list,1)]);     % (r+1)^2  x  maxheight/nbins
                                         % numel(x) = prod(size(x)).  if
                                         % size(x) = 2 x 3 x 4, then prod(size(x)) = 2*3*4 = 24  
yi = reshape(yi,[numel(x),size(list,1)]);     % (r+1)^2  x  maxheight/nbins
%size(xi)                                
hii = reshape(kron(1:nbins,ones(numel(xi),1)),[size(xi) nbins]);  % reshape to 81 x 100 x 16
%size(hii)
xii = repmat(xi,[1,1,nbins]);     % 81 x 100 x 16
size(xii);
yii = repmat(yi,[1,1,nbins]);     % 81 x 100 x 16
save('h_star_data');
rhoss = rhos(sub2ind(size(rhos),xii,yii,hii));    % 81 x 100 x 16
h = rhoss;
save('h_star_data');
end