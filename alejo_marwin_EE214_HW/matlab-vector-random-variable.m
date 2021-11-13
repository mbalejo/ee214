N = 100000;
bw = 0.04;
xbins = [-1:bw:1];
ybins = [-1:bw:1];
iter = 100;
M = length(xbins);
Nsamples = zeros(M);
count = 0;

for ii = 1:iter
x = 2*rand(1,N)-1;
y = 2*rand(1,N)-1;

X = [];
Y = [];

for k=1:N
if x(k)^2+y(k)^2<1
X = [X x(k)];
Y = [Y y(k)];
end
end

count = count+length(X);

for m=1:length(xbins)
    for n=1:length(ybins)
        temp1=(abs(X - xbins(m))<bw/2);
        temp2=(abs(Y - ybins(n))Nbw/2);
        Nsamples(m,n)=Nsamples(m,n)+sum(temp1.*temp2);
    end
end

end

PDFest = Nsamples/(count*bw^2);
mesh(xbins, ybins, PDFest);
xlabel('x'); ylabel('y');
zlabel('Joint PDF');


