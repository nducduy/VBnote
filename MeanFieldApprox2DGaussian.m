mu = [0,0]; %// data
sigma = [1 0.8; 0.8 1]; %// data
 x = linspace(-5,5,100); %// x axis
 y = linspace(-4,4,100); %// y axis
 
 [X Y] = meshgrid(x,y); %// all combinations of x, y


%mean-field approximation
m1=0.1;
m2=0.1;
Lambda11=sigma(1,1);
Lambda22=sigma(2,2);

n=50;
for i=1:n
    m1=mu(1)-sigma(1,2)/sigma(1,1)*(m2-mu(2));
    m2=mu(2)-sigma(2,1)/sigma(2,2)*(m1-mu(1));
end


q1=@(x)normpdf(x,m1,1/sqrt(sigma(1,1)));
q2=@(y)normpdf(y,m2,1/sqrt(sigma(2,2)));

q=@(x,y)q1(x)*q2(y);
Q=q1(X(:)).*q2(Y(:));
Q = reshape(Q,size(X));

% % define fxy
% f=@(x,y)exp(-0.5.*([x;y]-mu')'*inv(sigma)*([x;y]-mu'))./(2*pi*sqrt(det(sigma)));
% F=f(X(:),Y(:));
% F = reshape(F,size(X));
% figure
% 
% fcontour(q,'r');
% hold on
% fcontour(f,'b');
% hold off;


XY = cat(3,X,Y);
% subtract mu
XYmmu = bsxfun(@minus,XY,shiftdim(mu(:),-2));

isigXY = squeeze(sum(bsxfun(@times,shiftdim(inv(sigma),-2),XYmmu),3));
XYisXY = sum(isigXY .* XYmmu,3);

F = (1/(2*pi*sqrt(det(sigma)))) * exp(-0.5 * XYisXY);

figure
hold on

[~, h1] = contour(X,Y,Q, '-r');
h1_ = plot(NaN, '-r');
[~, h2] = contour(X,Y,F, '-b');
h2_ = plot(NaN, '-b');

L = legend([h1_ h2_], 'Mean-field approximation: q(z)', 'Exact Posterior: p(z)','Interpreter','latex','FontSize',36);
