function L=lagrange1D(x,pointx)
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
end