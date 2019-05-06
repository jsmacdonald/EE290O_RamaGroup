function f = current_control(x_Ictrl,Qcmd,Iqterm,Vterm,Iqcmd,Ipcmd,Vterm_theta,Pcmd,w,params)
% Function returns the sys in state space form, to be concatenated with
% other subsystems of in the inverter
% Inputs: [Qcmd, Qgen]
    % Qgen=Iqterm*Vterm, the whole inverter Qoutput
    % Vterm is a constant, but here treating as an input to the func
% Outputs: [Iqcmd, Ipcmd]

% Params that form the "boundary current" limits
% Iqmin=params.Iqmin % not used yet
% Iqmax=params.Iqmax % not used yet
% Ipmax=params.Ipmax % not used yet
Kvi=params.Kvi;
Kqi=params.Kqi;
ws=params.ws;
Tfrq=params.Tfrq;
Kw=params.Kw;
theta_inf=params.theta_inf;

% % See handwritten work for derivation of state space form from GE PV
% % inverter paper "Solar Photovoltaic (PV) Plant Models in PSLF"
% A=[0 0;...
%     Kvi 0];
% B=[KQi -KQi 0;...
%     0 0 0];
% C=[0 1];
% D=[0 0 1/Vterm];
% 
% % naming is needed for concatenation
% mySys=ss(A,B,C,D,'InputName',{'cur_Qcmd','cur_Qgen','cur_Pord'},'OutputName',{'Iqcmd','Ipcmd'});

g1=x_Ictrl(1);
kphi=60; % associated with 60Hz, see equations for derivation
f=[
    % Differential:
    %d(g1)/dt=
    Kqi*(Qcmd-Vterm*Iqterm); % Qgen=Vterm*Iqterm
    
    % Algebraic:
    %0=
    Pcmd/Vterm-Ipcmd; % Set Ipcmd
    
    % Differential:
    % d(Iqcmd)/dt=
    Kvi*(g1-Vterm);

    % Algebraic:
    % 0=
    kphi*(w-ws)-(Vterm_theta-theta_inf); % set w, w and ws are pu
        % change in torque (Milano eq 15.5 part 2), 377=2*pi*60 because w and ws are in pu   
        
    % 0=
%    (-1/Kw)*(w-ws)-Pcmd; % Pcmd=
    % d(Pcmd)/dt= 
    Kw*(w-ws); 
];
end

% referring to pg 4 of this paper for w eqns: https://arxiv.org/pdf/1206.5033.pdf
%  originally tried d(w)dt=g2, but this makes DAE second order