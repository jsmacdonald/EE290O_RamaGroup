% Create simple circuit, inverter, RLC load, gnd
function inverter_dxdt=bounaryinv_infBus(x,inverter_params);
%% Need to set these...
% Vmeas
% Vref
% Vslow
% Qgen
% Pord
% Vterm

% Load params
% foo

% Internal states of each block
x_QVdroop=x(1:3)
x_Ictrl=x(4:5)
x_phys=x(6:7)

% I/O exchanged between each block
Qcmd
Iqcmd
Ipcmd
Iterm

Pord=1; % Pord seems to be the raw power from PV panel, 
Vterm=1;
Vref=1;

inverter_dxdt=[
% DAEs
    QVdroop(x_QVdroop,Vmeas, Vref,Qcmd,params); % 4 diff eq
    current_control(x_Ictrl,Qcmd, Qgen, Vterm,Iqcmd, Ipcmd, inverter_params); % 2 diff eq, 1 alg
    physConv(x_phys,Ipcmd,Iqcmd,Vterm,Iterm,params) % 2 diff eq, 1 alg
    ];


%% Create combined inv model
% using series MATLAB func reference: https://www.mathworks.com/help/control/ref/series.html
sys1= QVdroop(Vmeas, Vref, Vslow, params)
sys2= current_control(Qcmd, Qgen, Pord, Vterm, params)
sys3= physConv(Ipcmd,Iqcmd,Vterm, params)

% double check I/O to make sure connecting right terminals
sys1.InputName
sys1.OutputName
sys2.InputName
sys2.OutputName
sys3.InputName
sys3.OutputName
sys12=series(sys1,sys2,[1],[1]) % connect 1st output to 1st input
combinedSys=series(sys12,sys3,[1 2],[1 2]) % connect the only 2 outputs to the only 2 inputs

% Once ready, should be able to solve combinedSys with ODE45:
% https://www.mathworks.com/matlabcentral/answers/146782-solve-state-space-equation-by-ode45