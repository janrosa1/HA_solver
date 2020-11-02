% Huggett 1993 Replication
% Macro- ECON 516 Midterm
% Aditi Singh, Jan Rosa, Sudipta Ghosh 


tic
clear all

%% 1. Parameters

% endowment
eh=1; 
el=0.1;
pihh=0.925;
pihl=0.5;
trans_mat=[pihh pihl; 1-pihh 1-pihl]; % this is the transition matrix [phh phl;plh pll]

beta=0.99322;
sigma=1.5; % risk aversion para
sigma1=3;
% asset grid
%amin=-2; % can change here
amin_c=[-2 -4 -6 -8];
amax=8;
aobs=200; % how many grid points we have

% tolarence levels
ctol=1e-7;
qtol=2.5e-4;

%% Try function for different sigmas

[r_15, q_15] = HuggettSolve(amin_c,amax,aobs,eh,el,beta, sigma,ctol,qtol,trans_mat);
[r_3, q_3] = HuggettSolve(amin_c,amax,aobs,eh,el,beta, sigma1,ctol,qtol,trans_mat);

 
%% 2.For an initial guess on prices, get policy function I construct the following function
% I use endowgenous grid method here
function[r_out, q_out]= HuggettSolve(amin_c,amax,aobs,eh,el,beta, sigma,ctol,qtol,trans_mat)
% result matrices
q_out=zeros(1,4);
r_out=zeros(1,4);


for m=1:4
    amin=amin_c(m);
    agrid_1=makegrid(2,aobs,amin,amax)';
    agrid  =repmat(agrid_1,1,2);
    egrid  =repmat([eh el],aobs,1);
    
    % qguess
    qmin=beta+1e-4;
    qmax=1/0.9;  % initially 0.9
    q0=1/beta+0.002; % shouldnt have beta/q>1
    c0=agrid+egrid+q0*8; % initial guess for c

    k=0;
    q_iter=0;
    while q_iter==0
        c_iter=0;
        while c_iter==0
            k=k+1;
            RHS=beta/q0*uti(sigma,c0)*trans_mat; % this is ith rhs of ee
            ctil=(RHS.^(-1/sigma)); % this is the ctil(a',y)
            for j=1:2
                for i=1:aobs
                    if ctil(i,j)==0
                        ctil(i,j)=1e-4;
                    else
                        ctil(i,j)=ctil(i,j);
                    end
                end
            end
            astar=ctil+q0*agrid-egrid; % astar as function of a' and y

            alstar=astar(1,:); % the asset state induce borrowing constraint binding
            ahstar=astar(aobs,:); % highest value
            c1=zeros(aobs,2);
            for j=1:2
                for i=1:aobs
                    if agrid(i,j)>=alstar(j) && agrid(i,j)<ahstar(j)
                        %c1(i,j)=interp1q(astar(:,j),ctil(:,j),agrid(i,j)); % case borrowing constraint slack
                        % interploration
                        a=agrid(i,j);
                        l=min(find(astar(:,j)-a>0));% the one bigger than agrid

                        cl=ctil(l-1,j);
                        ch=ctil(l,j);
                        al=astar(l-1,j);
                        ah=astar(l,j);
                        c1(i,j)=cl+(ch-cl)*(a-al)/(ah-al);
                    elseif agrid(i,j)>=ahstar(j)
                        c1(i,j)=0;
                    else
                        c1(i,j)=agrid(i,j)+egrid(i,j)-q0*amin;% bc binding
                    end
                end
            end

            dif=abs(c1-c0);
            if max(max(dif))<=ctol
                c_iter=1;
            else
                c0=c1;
            end

        end

        % get policy a'

        ap=(agrid+egrid-c1)/q0;

        % iteration on stationary dist
        phi1=zeros(aobs,2);

        phi0=zeros(aobs,2);% this is the initial guess for dist over a and e
        for i=1:aobs
            phi0(i,:)=i/(2*aobs)*ones(1,2);
        end
        distol=1e-7;
        p_iter=0;
        while p_iter==0 
            phi=[0 0; phi0*trans_mat'];% this is phi(a',e') but a' not necessarily on grid
            for j=1:2
                for i=1:aobs
                    if max(ap(:,j))>=agrid(i,j)
                        l=min(find(ap(:,j)-agrid(i,j)>=0));
                    else
                        l=aobs;
                    end
                    if l<aobs && l>1
                        phih=phi(l+1,j);
                        phil=phi(l,j);
                        app=[amin amin;ap];
                        ah=app(l+1,j);
                        al=app(l,j);
                        phi1(i,j)=phil+(agrid(i,j)-al)/(ah-al)*(phih-phil);
                    elseif l>=aobs
                        phi1(i,j)=phi(aobs,j);
                    else
                        phi1(i,j)=phi(1,j);
                    end
                end
            end
            dif_p=max(max(abs(phi0-phi1)));
            if dif_p<=distol
                p_iter=1;
                phi_s=phi1;
            else
                phi0=phi1;
            end
        end

        % excess demand and suply
        p_s=diff([0 0; phi_s]);
        A=ap'*p_s;
        AD=A(1,1)+A(2,2); % this is aggregate demand, if >0 excess dem (over saving, should inc q), <0 sup
        if abs(AD)<qtol
            q_iter=1;
            q_e=q0;
        elseif AD>0
            qmin=q0;
            q0=0.5*(qmin+qmax);
        else
            qmax=q0;
            q0=0.5*(qmin+qmax);
        end


    end
    q_out(m)=q_e;
    r_out(m)=(1/q_e)^6;
end
toc
disp(['q1   q2   q3   q4'])
[q_out(1)  q_out(2)  q_out(3)  q_out(4)]
disp(['r1   r2   r3   r4'])
[r_out(1)-1  r_out(2)-1  r_out(3)-1  r_out(4)-1]
end

function u=uti(sig,c)
u=c.^(-sig);
end

function y=makegrid(k,n,min,max)

if k==1
    y=linspace(min,max,n);  % equidistant grid
elseif k==2
    y=linspace(0,log(1+max-min),n); % exponential grid, more points at low assets, makes sense if curvature (or action) is there
    y=exp(y)-1;
    y=y+min;
else
    y=linspace(0,log(log(1+max-min)+1),n);  % exponential grid, even more points at low assets, same reason
    y=exp(exp(y)-1)-1;
    y=y+min;
end
end