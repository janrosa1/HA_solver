%% 2.For an initial guess on prices, get policy function I construct the following function
% I use endowgenous grid method here
function[r_result, q_result]= HuggettSolveEGM(amin_c,amax,grid_len,eh,el,beta, sigma,c_tol,q_tol,trans_mat)
% result matrices
q_result=zeros(1,4);
r_result=zeros(1,4);


for m=1:4
    amin=amin_c(m);
    agrid_1=makegrid(2,grid_len,amin,amax)';
    agrid  =repmat(agrid_1,1,2);
    egrid  =repmat([eh el],grid_len,1);
    
    % qguess
    q_min=beta+1e-4;
    q_max=1/0.9;  % initially 0.9
    q0=1/beta+0.002; % shouldnt have beta/q>1
    c0=agrid+egrid+q0*8; % initial guess for c

    k=0;
    iter_q=0;
    while iter_q==0
        iter_cons=0;
        while iter_cons==0
            k=k+1;
            b=beta/q0*c0.^(-sigma)*trans_mat; % this is ith rhs of ee
            ctil=(b.^(-1/sigma)); % this is the ctil(a',y)
            for j=1:2
                for i=1:grid_len
                    ctil(i,j)= max(ctil(i,j),1e-4);
                end
            end
            astar=ctil+q0*agrid-egrid; % astar as function of a' and y

            alstar=astar(1,:); % the asset state induced binding borrowing constraint 
            ahstar=astar(grid_len,:); % highest value
            c_new =zeros(grid_len,2);
            for j=1:2
                for i=1:grid_len
                    if agrid(i,j)>=alstar(j) && agrid(i,j)<ahstar(j)
                        % interploration
                        a=agrid(i,j);
                        l=min(find(astar(:,j)-a>0)); % next greater value 

                        clow=ctil(l-1,j);
                        chigh=ctil(l,j);
                        alow=astar(l-1,j);
                        ahigh=astar(l,j);
                        c_new(i,j)=clow+(chigh-clow)*(a-alow)/(ahigh-alow);
                    elseif agrid(i,j)>=ahstar(j)
                        c_new(i,j)=0;
                    else
                        c_new(i,j)=agrid(i,j)+egrid(i,j)-q0*amin;
                    end
                end
            end

            delta=abs(c_new-c0);
            if max(max(delta))<=c_tol
                iter_cons=1;
            else
                c0=c_new;
            end

        end

        % get policy a'

        apolicy=(agrid+egrid-c_new)/q0;

        % iteration on stationary dist
        phi1=zeros(grid_len,2);

        phi0=zeros(grid_len,2);% this is the initial guess for dist over a and e
        for i=1:grid_len
            phi0(i,:)=i/(2*grid_len)*ones(1,2);
        end
        distol=1e-7;
        iter_pdf=0;
        while iter_pdf==0 
            phi=[0 0; phi0*trans_mat'];% this is phi(a',e') but a' not necessarily on grid
            for j=1:2
                for i=1:grid_len
                    if max(apolicy(:,j))>=agrid(i,j)
                        l=min(find(apolicy(:,j)-agrid(i,j)>=0));
                    else
                        l=grid_len;
                    end
                    if l<grid_len && l>1
                        phih=phi(l+1,j);
                        phil=phi(l,j);
                        policy_n=[amin amin;apolicy];
                        ahigh=policy_n(l+1,j);
                        alow=policy_n(l,j);
                        phi1(i,j)=phil+(agrid(i,j)-alow)/(ahigh-alow)*(phih-phil);
                    elseif l>=grid_len
                        phi1(i,j)=phi(grid_len,j);
                    else
                        phi1(i,j)=phi(1,j);
                    end
                end
            end
            delta_p=max(max(abs(phi0-phi1)));
            if delta_p<=distol
                iter_pdf=1;
                phi_s=phi1;
            else
                phi0=phi1;
            end
        end

        % excess demand and suply
        p_s=diff([0 0; phi_s]);
        A=apolicy'*p_s;
        AD=A(1,1)+A(2,2); % aggregate demand: if >0 excess demand (over saving:q should increase), <0 excess supply
        if abs(AD)<q_tol
            iter_q=1;
            q_equilm=q0;
        elseif AD>0
            q_min=q0;
            q0=0.5*(q_min+q_max);
        else
            q_max=q0;
            q0=0.5*(q_min+q_max);
        end


    end
    Dist = ( phi1(1:grid_len,1) )' + ( phi1(1:grid_len,2) )';
    q_result(m)=q_equilm;
    r_result(m)=((1/q_equilm)^6)-1;
if m==1   % Plot only for the borrowing constraint =-2  

% SOME ADDITIONAL RESULTS FOR GRAPHS
popSS=reshape(phi1',grid_len,2);                                           % density: col1 E ; col2 U
popCUM=cumsum(popSS);
indRich=find(popSS(:,1)>0,1,'first');                                  % identifies the richest agent
ub=indRich;
phi2 = cumsum(phi1);   
    
close all

 

% (GROSS) SAVINGS
figure(1)
plot(agrid_1,apolicy(:,1),'-',agrid_1,apolicy(:,2),'-',agrid_1,agrid_1,'--')
xlim([-2 8])
ylim([-2 8])
title('Optimal Decision Rule:(amin=-2) $\sigma$=',num2str(sigma),'Interpreter','latex')
xlabel('Current Period Asset')
ylabel('Next Period Asset')
legend('SAVINGS EMPLOYED','45 DEGREE LINE','SAVINGS UNEMPLOYED','Location','northwest')
saveas(figure(1),['optimal_policy',num2str(sigma),'.png'])

% WEALTH DISTRIBUTION

figure(2)
plot(agrid_1(1:grid_len),phi1(1:grid_len,1),'-',agrid_1,phi1(:,2),'-')
xlim([-2 0.97])
title('Stationary Distribution: $\sigma=$',num2str(sigma),'fontsize',8,'Interpreter','latex')
xlabel('Credit Stock')
ylabel('Frequency')
legend('CDF EMPLOYED','CDF UNEMPLOYED','Location','NorthWest')

saveas(figure(2),['Distribution',num2str(sigma),'.png'])
end

end
toc
disp(['q1   q2   q3   q4'])
[q_result(1)  q_result(2)  q_result(3)  q_result(4)]
disp(['r1   r2   r3   r4'])
[r_result(1)  r_result(2)  r_result(3)  r_result(4)]

end



