%% System Formation 
clear
clc
syms omega phi kapa xl yl zl s x1 y1 z1 x2 y2 z2
adindan =load('adindan.txt');
adindan =adindan(:,2:4);
wgs84 =load('WGS84.txt');
wgs84 =wgs84(:,2:4);

Rx = [1 0 0; 0 cos(omega) -sin(omega);0 sin(omega) cos(omega)]
Ry = [cos(phi) 0 sin(phi); 0 1 0;-sin(phi) 0 cos(phi)]
Rz = [cos(kapa) -sin(kapa) 0;sin(kapa) cos(kapa) 0;0 0 1]
Rxyz(omega,phi,kapa) = Rz*Rx*Ry
dv = [omega phi kapa xl yl zl s]
T =[xl;yl;zl]
X1 = [x1;y1;z1]
X2 = [x2;y2;z2]
objf(x1,y1,z1,x2,y2,z2)= (X2 - s*Rxyz*X1+T).^2
fobjf =0;
for i=1:8
    X22 =adindan(i,:)';
    X11 =wgs84(i,:)';
    fobjf =fobjf+objf(X11(1),X11(2),X11(3),X22(1),X22(2),X22(3));
end
fobjf(omega,phi,kapa,xl,yl,zl,s) =sum(fobjf)
gv = gradient(fobjf,dv)

of =matlabFunction(fobjf);
gv =matlabFunction(gv);
H =hessian(fobjf,dv);
H =matlabFunction(H);
%% Initialization 
solold = [0,0,0,0,0,0,1]';
%ftol = 1E-10; 
%gtol =1E-10;
lamda=1;
iter=0;
 fold = of(solold(1),solold(2),solold(3),solold(4),solold(5),solold(6),solold(7));
%% Iteration 
while 1
   iter=iter+1;
   gvf = gv(solold(1),solold(2),solold(3),solold(4),solold(5),solold(6),solold(7));
   Hx =H(solold(1),solold(2),solold(3),solold(4),solold(5),solold(6),solold(7));
   %search direction
   sd = -inv(Hx +lamda*eye(size(dv,1)))*gvf;
   solnew = solold+sd;
   fnew =of(solnew(1),solnew(2),solnew(3),solnew(4),solnew(5),solnew(6),solnew(7));
   if fnew<fold
      lamda =lamda/2;
   else
       lamda =2*lamda;
   end
   if abs(sd(1))>1E-9|abs(sd(2))>1E-9|abs(sd(3))>1E-9 |abs(sd(4))>1E-6|abs(sd(5))>1E-6|abs(sd(6))>1E-6| abs(sd(7))>1E-12 
       fprintf('Iter : %.0f,F Loss:%.10f \n',iter,fnew)
   else
       disp('Final Solution' )
       disp(solnew);
       disp(fnew);
       break;
   end
solold=solnew;
fold =fnew;
end
%%