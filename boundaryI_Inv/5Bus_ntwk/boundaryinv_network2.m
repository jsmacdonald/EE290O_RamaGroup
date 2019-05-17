% NEW VERSION
% passing in knowns cell array now, deleted off knowns assignment for each
% bus type
function dxdt = boundaryinv_network2(t,x,c,params,Ts,invBus,knowns);
V=knowns{1}; del=knowns{2}; P=knowns{3}; Q=knowns{4};

dxdt=[]; % will fill this
% input c is case struct for power network from MATPOWER
N = size(c.bus,1);
L = size(c.branch,1);
%[Ybus,Yf,Yt]=makeYbus(c); % Yt almost the same as Y

% Known values, hard to generalize by pulling from c
%  Qmeas(4)=-0.5365;
%  Pmeas=[2.0988 -3.0021 0.2373 -3.9497 4.666]';
% vmagVec=[1 0.989 1 1 1]';
% thetaVec_deg=[3.273 -0.759 -0.492  0 4.112]'; % degrees
% thetaVec=thetaVec_deg*(pi/180); % convert to radians

   [M_k,Mbar_k,Y_k,Ybar_k,Y_lm,Ybar_lm,Y_ml,Ybar_ml]=computeYmats(c);
 %% PF equation matrix formulas
    h_V=@(x,M_k) trace(M_k*x*x'); % note Tr(Mk*x*x')=V^2, not V!
    h_del=@(x,Mbar_k) trace(Mbar_k*x*x'); % note Tr(Mk*x*x')=V^2, not V!
    h_P=@(x,Y_k) trace(Y_k*x*x');
    h_Q=@(x,Ybar_k) -trace(Ybar_k*x*x');
    h_Pline=@(x,Y_lm) trace(Y_lm*x*x'); % can take in l-->m or m-->l Y matrices
    h_Qline=@(x,Ybar_lm) -trace(Ybar_lm*x*x'); % add neg to make sign convention agree with runpf and vYv formulation
    %x=[real(V)' imag(V)']'

    
%% Collect knowns from case
% vmagVec=c.bus(:,8).*c.bus(:,10); % units=V, pu*base=actual
% thetaVec_deg=c.bus(:,9); % units=degrees
% thetaVec=thetaVec_deg*(pi/180); % convert to radians

%% build dxdt equations
    nInvState=17-4;
    ofs=nInvState; % offset state assignment to not overwrite inv states

    x_QVdroop=x(1:3); % internal states gi, no specific name
    Qcmd=x(4);
    Iqcmd=x(5); 
    Ipcmd=x(6);
    w=x(7);
    Pcmd=x(8);
    x_phys=x(9:10); % internal states gi
    Ipterm=x(11);
    Iqterm=x(12);
%     Vterm=x(13); 
%     Vterm_theta=x(14);
%     Pt=x(15);
%     Qt=x(16);
    Vref=x(13);
    
    Pt=x(ofs+2*N+1);
    Qt=x(ofs+2*N+2);
    Vterm=x(ofs+2*N+3);
    Vterm_theta=x(ofs+2*N+4);
    
    dxdt=[dxdt;...
      %% Now add inverter eqns
    QVdroop(x_QVdroop,Vterm,Vref,Qcmd,params); % 4 diff eq, g1/g2/g3/Qcmd
    current_control2(Qcmd,Pt,Vterm,Iqcmd,Ipcmd,Pcmd,params); % 2 diff eq, 3 alg, mixed ordering
    physConv(x_phys,Ipcmd,Iqcmd,Ipterm,Iqterm,params); % 2 diff eq, 2 alg, g1/g2/Ipterm/Iqterm
    0; % d(Vref)=0
];

%% add in PF eqns assoc with network
    for k = 1:N % creat 2*N bus equations, after this state vector is (ofs+2N)x1
            busType=c.bus(k,2);           
            switch busType % take different measurements at diff bus types
                case 1 % slack bus
                    % fill in unknowns
                    P(k)=x(ofs+2*k-1);
                    Q(k)=x(ofs+2*k); 
                    
                    % setup eqns with rest of unknowns
                    dxdt=[dxdt; -P(k)+h_P([V; del],Y_k{k})]; % create with possibilities of unknowns
                    dxdt=[dxdt; -Q(k)+h_Q([V; del],Ybar_k{k})]; % create with possibilities of unknowns

                case 2 % gen bus 
                    % fill in unknowns
                    del(k)=x(ofs+2*k-1);
                    Q(k)=x(ofs+2*k); 
                    % setup eqns with rest of unknowns
                    dxdt=[dxdt; -del(k)+h_del([V; del],Mbar_k{k})]; % create with possibilities of unknowns
                    dxdt=[dxdt; -Q(k)+h_Q([V; del],Ybar_k{k})]; % create with possibilities of unknowns

                case 3 % load bus
                    % fill in unknowns
                    V(k)=x(ofs+2*k-1);
                    del(k)=x(ofs+2*k);
                    % fill in knowns
                    if k==invBus
                        % need eqns for setting Qt,Pt,Vterm, Vterm_theta
                       eqn=[Vterm*Ipterm-P(k)-Pt;...
                           Vterm*Iqterm-Q(k)-Qt;...
                            V(k)-Vterm;...
                            del(k)-Vterm_theta]; % save for after this loop                           
                       %pause
                    end
                    % solve for V & del
                    dxdt=[dxdt; -V(k)+h_V([V; del],M_k{k})]; % create with possibilities of unknowns
                    dxdt=[dxdt; -del(k)+h_del([V; del],Mbar_k{k})]; % create with possibilities of unknowns
            end
    end
    dxdt=[dxdt; eqn]; % append to dxdt at end so can align with x asignments
end
