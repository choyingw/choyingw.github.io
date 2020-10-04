function y = linear_sg(x,gamma,lambda)

k1=0.2; k2=0; kh=0.3;
y1=zeros(size(x)); y2=y1; y3=y1; %id = 1:length(x);
id1=find(x<=gamma(1)); id3=find(x>=gamma(2));% id2=setdiff(1:length(x),[id1',id3']);  
id2=find(x>gamma(1) & x<gamma(2)); 
x(x>gamma(3))=gamma(3);
y1(id1)=x(id1)*(2-kh)/(-gamma(1))+2;
y2(id2)=x(id2)*(kh-k1)/(gamma(1)-gamma(2))+kh-gamma(1)*(kh-k1)/(gamma(1)-gamma(2));
y3(id3)=x(id3)*(k1-k2)/(gamma(2)-gamma(3))+k1-gamma(2)*(k1-k2)/(gamma(2)-gamma(3));
y=lambda*(y1+y2+y3);
%y(y==0)=0.1;
