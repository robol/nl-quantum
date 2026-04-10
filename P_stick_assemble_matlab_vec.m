function [P]= P_stick_assemble_matlab_vec(nSticks,nNodes,NN,G,radius,npg,ne_cap_max,Cap_Elem,n_thread)
%%
[gauss_P,gauss_W] = gauleg(npg);
%%
epsilon0=8.85418781762e-12;
P=zeros(nNodes,nNodes);
for hh=1:nNodes
    ne_cap_hh=Cap_Elem(1,hh);
    idE_hh=Cap_Elem(2:ne_cap_max+1,hh);
    P_self_jj_jj=0.0;
    P_self_jj_ii=0.0;
    ll_tot_hh=0.0;
    for jj=1:ne_cap_hh
        NN_edge=NN(1:3,G(1:2,idE_hh(jj)));
        NN_jj(1:3,1) = 0.5*(NN_edge(1:3,1)+NN_edge(1:3,2));
        NN_jj(1:3,2) = NN(1:3,hh);
        [alpha_jj]=fun_1_R_2stick_self(NN_jj,radius(idE_hh(jj)));
        P_self_jj_jj=P_self_jj_jj+alpha_jj;
        [PPg_jj,ll_jj]=Gauss_line_nvar(NN_jj,gauss_P,npg);
        ll_tot_hh=ll_tot_hh+ll_jj;
        for ii=(jj+1):ne_cap_hh
            NN_edge=NN(1:3,G(1:2,idE_hh(ii)));
            NN_ii(1:3,1) = 0.5*(NN_edge(1:3,1)+NN_edge(1:3,2));
            NN_ii(1:3,2) = NN(1:3,hh);
            ll_ii=norm(NN_ii(1:3,2)-NN_ii(1:3,1));
            integ_self=0.0;
            for ff=1:npg
                [log_eps]=fun_1_R_stick(NN_ii,PPg_jj(1:3,ff),ll_ii);
                integ_self=integ_self+gauss_W(ff)*log_eps;
            end
            P_self_jj_ii=P_self_jj_ii+0.5*ll_jj*integ_self;
        end
    end
	P(hh,hh)=(P_self_jj_jj+2*P_self_jj_ii)/(4*pi*epsilon0*ll_tot_hh^2);
end
for hh=1:nNodes
    ne_cap_hh=Cap_Elem(1,hh);
    idE_hh=Cap_Elem(2:ne_cap_max+1,hh);
    glob_P=0.0;
    for kk=hh+1:nNodes
        ne_cap_kk=Cap_Elem(1,kk);
        idE_kk=Cap_Elem(2:ne_cap_max+1,kk);
        ll_tot_hh=0.0;
        P_mutual_jj_ii = 0.0;
            for jj = 1:ne_cap_hh
                NN_edge=NN(1:3,G(1:2,idE_hh(jj)));
                NN_jj(1:3,1) = 0.5*(NN_edge(1:3,1)+NN_edge(1:3,2));
                NN_jj(1:3,2) = NN(1:3,hh);          
                [PPg_jj,ll_jj]=Gauss_line_nvar(NN_jj,gauss_P,npg); 
                ll_tot_hh=ll_tot_hh+ll_jj; 
                ll_tot_kk=0.0;		  
                for ii = 1:ne_cap_kk 
                    NN_edge=NN(1:3,G(1:2,idE_kk(ii)));
                    NN_ii(1:3,1) = 0.5*(NN_edge(1:3,1)+NN_edge(1:3,2));
                    NN_ii(1:3,2) = NN(1:3,kk);
                    ll_ii=norm(NN_ii(1:3,2)-NN_ii(1:3,1));
                    ll_tot_kk=ll_tot_kk+ll_ii; 			 
                    integ_mutual=0.0;
                        for ff=1:npg
                            [log_eps]=fun_1_R_stick(NN_ii,PPg_jj(1:3,ff),ll_ii);
                            integ_mutual=integ_mutual+gauss_W(ff)*log_eps;			 
                        end
                    P_mutual_jj_ii=P_mutual_jj_ii+0.5*ll_jj*integ_mutual;
                end 		  
            end 
        P(hh,kk)=P_mutual_jj_ii/(ll_tot_kk*ll_tot_hh*4*pi*epsilon0);
        P(kk,hh)=P(hh,kk);
    end 
end 
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
%% ************** fun_1_R_2stick_self ****************
function [alpha] = fun_1_R_2stick_self(NN,radius)
ll=norm(NN(1:3,2)-NN(1:3,1));
aa=log(ll/radius+sqrt((ll/radius)^2+1));
bb=sqrt(1+(radius/ll)^2);
cc=radius/ll;
alpha=2.0*ll*(aa-bb+cc);
end    
%% ************** fun_1_R_stick ***************
function [log_eps]=fun_1_R_stick(PP,PP0,ll)
ri=norm(PP0-PP(1:3,1));
rf=norm(PP0-PP(1:3,2));
eps=ll/(ri+rf);
log_eps=log((1+eps)/(1-eps));
end