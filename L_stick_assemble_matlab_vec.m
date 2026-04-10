function [L]=L_stick_assemble_matlab_vec(nSticks,nNodes,NN,G,radius,npg,Wint,n_thread)
L=zeros(nSticks,nSticks);
%%
[gauss_P,gauss_W] = gauleg(npg);
%%
for hh=1:nSticks
    %self-inductancE
    PPh=NN(1:3,G(1:2,hh));
    [PPgh,ll_h]=Gauss_line_nvar(PPh,gauss_P,npg);
    aa=log(ll_h/radius(hh)+sqrt((ll_h/radius(hh))^2+1));
    bb=sqrt(1+(radius(hh)/ll_h)^2);
    cc=radius(hh)/ll_h;
    L(hh,hh)=4*ll_h*(aa-bb+cc+(Wint)*0.25);
    ut_h=(PPh(1:3,2)-PPh(1:3,1))/ll_h;
    %mutual-inductance
    for kk=(hh+1):nSticks
        PPk=NN(1:3,G(1:2,kk));
        ll_k=norm(PPk(1:3,2)-PPk(1:3,1));
        ut_k=(PPk(1:3,2)-PPk(1:3,1))/ll_k;
        integ =0;
        for ii=1:npg
            ri=fun_my_norm([PPgh(1,ii)-PPk(1,1);...
                            PPgh(2,ii)-PPk(2,1);...
                            PPgh(3,ii)-PPk(3,1)]);
            rf=fun_my_norm([PPgh(1,ii)-PPk(1,2);...
                            PPgh(2,ii)-PPk(2,2);...
                            PPgh(3,ii)-PPk(3,2)]);
            eps=ll_k./(ri+rf).';
            log_eps=log((1+eps)./(1-eps));
            integ=integ+gauss_W(ii)*log_eps.*fun_my_dot(ut_h,ut_k).';
        end 
        L(hh,kk)=ll_h*integ.';
	    L(kk,hh)=L(hh,kk).';
    end 
end 
L=1.0d-7*0.5*L;
end
%% ********************************************
function [C] = fun_my_dot(A,B)
C = A(1)*B(1,:)+A(2)*B(2,:)+A(3)*B(3,:);
end
%% ********************************************
function [B] = fun_my_norm(A)
B = sqrt(A(1,:).^2 +A(2,:).^2 + A(3,:).^2);
end
%% ************* GAUSS-LEGENDRE *********************
function [xabsc,weig] = gauleg(ngp)
%%
xabsc=zeros(ngp,1);
weig=zeros(ngp,1);
%%
EPS2=3.0e-15;
m=(ngp+1)/2;
%Roots are symmetric in the interval so only need to find half of them
for i=1:m
z=cos(pi*(i-0.25d0)/(ngp+0.5d0)); %starting approximation
z1=2;
%Newton's method
while max(abs(z-z1)) > EPS2
    p1 = 1;
    p2 = 0;
    %Loop up the recurrence relation to get the Legendre
    %polynomial evaluated at z
    for j=1:ngp
        p3=p2;
        p2=p1;
        p1=((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j;
    end
    %p1 is now the desired Legendre polynomial. We next compute pp,
    %its derivative, by a standard relation involving also p2, the
    %polynomial of one lower order.
    pp=ngp*(z*p1-p2)/(z*z-1.0d0);
    z1=z;
    z=z1-p1/pp;             %Newton's Method
end
xabsc(i)=-z;                    		% Roots will be bewteen -1.0 & 1.0
xabsc(ngp+1-i)=+z;            		% and symmetric about the origin
weig(i)=2.0d0/((1.0d0-z*z)*pp*pp);	% Compute the weight and its
weig(ngp+1-i)=weig(i);             % symmetric counterpart
end 
end
%% ************** Gauss_line_nvar ****************** 
function [PPg,ll]=Gauss_line_nvar(NN,gauss_P,np)
ll = norm(NN(1:3,2)-NN(1:3,1));
PPg(1,1:np) = NN(1,1)+0.5*(NN(1,2)-NN(1,1))*(1.0+gauss_P).';
PPg(2,1:np) = NN(2,1)+0.5*(NN(2,2)-NN(2,1))*(1.0+gauss_P).';
PPg(3,1:np) = NN(3,1)+0.5*(NN(3,2)-NN(3,1))*(1.0+gauss_P).';
end