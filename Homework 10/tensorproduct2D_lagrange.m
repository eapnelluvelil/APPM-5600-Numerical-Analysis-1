function z = tensorproduct2D_lagrange(xy, pointx, pointy, pointz)

x = xy(1,:);
y = xy(2,:);
nxy = size(xy,2);

% create the 1D interpolation matrices
Lx = lagrange1D(x,pointx);
Ly = lagrange1D(y,pointy);
nxi = length(pointx);
nyi = length(pointy);

L = ones(nxy,nxi*nyi);

for i = 1:nxi
    L(:,(i-1)*nyi+(1:nyi)) = (Lx(:,i)*ones(1,nyi)).*Ly;
end

z = L*pointz';
end