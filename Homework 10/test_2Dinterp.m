function test_2Dinterp

% f = @(x,y) (x+y)./(1+25*(x.^2+y.^2));
% f = @(x,y) sin(x+y);
f = @(x, y) (x+y)./(1 + 25*(x+y).^(2));

xeval = [-1:0.1:1];
nn = length(xeval);
yeval = xeval;
[Xeval,Yeval] = meshgrid(xeval,yeval);
xxeval = reshape(Xeval,1,nn^2);
yyeval = reshape(Yeval,1,nn^2);
xyeval = [xxeval;yyeval];

fxeval = f(xxeval,yyeval);


N = 30;
% equispaced nodes
% xtmp = linspace(-1,1,N);
% Chebychev nodes
xtmp = cos((2*(0:N-1)+1)*pi/(2*(N)));

ytmp = xtmp;
[X,Y] = meshgrid(xtmp,xtmp);
x = reshape(X,1,N^2);
y = reshape(Y,1,N^2);
fxy = f(x,y);

z = tensorproduct2D_lagrange(xyeval,xtmp,ytmp, fxy);

norm(fxeval-z')

figure(1)
Z = reshape(z,nn,nn);
surf(Xeval,Yeval,Z)
title('approximation')


figure(2)
surf(Xeval,Yeval, abs(Z-f(Xeval,Yeval)))
title('error')


% keyboard



return

function z = tensorproduct2D_lagrange(xy, pointx,pointy, pointz)

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



return



function L=lagrange1D(x,pointx)
%
%
n=size(pointx,2);
L2=ones(n,size(x,2));
   for i=1:n
      for j=1:n
         if (i~=j)
            L2(i,:)=L2(i,:).*(x-pointx(j))/(pointx(i)-pointx(j));
         end
      end
   end
   L = L2.';  

   
return

