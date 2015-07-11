% Creates a data struct, as it is used in the function CreateDatFile. Each
% data struct create by this function is guranteed to converge to a
% solution in MUSCOD.
% In the basic version, values are chosen to replicate the results of my
% dissertation/ the ICRA paper on costfunctions.
%
% legJointType is either
%   'PEA' or 'SEA'
%
% Possible uses:
% CreateDatFile('Data/HOPPING_1D_PL_BASE',CreateDatStructMonoped1D('PEA'))
% CreateDatFile('Data/HOPPING_1D_SL_BASE',CreateDatStructMonoped1D('SEA'))
function datStruct = CreateDatStructRAMoneWalk(costFunction,v_avg,sigma,const_sel)
    % First, define all properties that only depend on the gait type:
    sd = zeros(25,1);
    sd(1) = 1;      % x
    sd(2) = v_avg;  % dx
    sd(3) = 0.4;    % y
    sd(4) = 0;      % dy
    sd(5) = 0;      % phi
    sd(6) = 0;      % dphi
    sd(7) = 0;      % alphaL
    sd(8) = 0;      % dalphaL
    sd(9) = 0;      % ualphaL
    sd(10)= 0;      % dualphaL
    sd(11)= pi/4;   % alphaR
    sd(12)= 0;      % dalphaR
    sd(13)= pi/4;   % ualphaR
    sd(14)= 0;      % dualphaR
    sd(15)= -pi/4;  % betaL
    sd(16)= 0;      % dbetaL
    sd(17)= -pi/4;  % ubetaL
    sd(18)= 0;      % dubetaL
    sd(19)= -pi/4;  % betaR
    sd(20)= 0;      % dbetaR
    sd(21)= -pi/4;  % ubetaR
    sd(22)= 0;      % dubetaR
    sd(23)= 0;      % posActWork
    sd(24)= 0;      % posElWork
    sd(25)= 0;      % totElLoss
    
    p = zeros(4,1);
    p(1) = costFunction;      % Cost function selection index (cost_fct_sel) [1; 2; 3]
    p(2) = v_avg;  % Average speed (v_avg) [m/s]
    p(3) = sigma;   % Knee spring regularization constant (sigma) [rad]
    p(4) = const_sel;      % Constraint selection index (const_sel) [binary flags]

    u  = zeros(4,1);
    h  = [0.5;0;0.5];  
    nshoot = [50;1;50];
    
    typeString = 'Walking';
    datStruct.libmodel = 'libAB_Walking';
    datStruct.libind = {'ind_rkf45'; 'ind_strans'; 'ind_rkf45'};
    datStruct.h_name = {'single stance right'; 'touchdown'; 'double stance'};
    datStruct.h_comment = {'Single stance until touchdown'; 'Touchdown'; 'Double stance until liftoff'};
    datStruct.nshoot = nshoot;
    datStruct.h      = h;
    datStruct.h_sca  = [1;0;1];
    datStruct.h_min  = [0.01;0.0;0.01];
    datStruct.h_max  = [2.0;0.0;2.0];
    datStruct.h_fix  = [0;1;0];
    
    % Parameterized properties:
    datStruct.header = {'/***********************************************';
                          '*';
                          ['*  Gait optimization for an articulating biped in a ',typeString];
                          '*  gait with series elastic legs';
                          '*';
                          '*';
                          ['*   0: ',num2str(p(1),'% 8.5f'),'  // cost_fct_sel  [1; 2; 3]          Cost function selector (1 = posActWorkCOT; 2 = posElWorkCOT; 3 = totElLossCOT)'];
                          ['*   1: ',num2str(p(2),'% 8.5f'),'  // v_avg         [m/s]              Enforced average speed'];
                          ['*   2: ',num2str(p(3),'% 8.5f'),'  // sigma         [rad]              Knee spring regularization constant'];
                          ['*   3: ',num2str(p(4),'% 8.5f'),'  // const_sel     [binary flags]     1: GroundClearance 2: TorqueLimits 4: SpeedLimits'];
                          '*';
                          '***********************************************/'};

    %% States:
    datStruct.xd_name = {'x'; 'dx'; 'y'; 'dy'; 'phi'; 'dphi'; 
                         'alphaL'; 'dalphaL'; 'ualphaL'; 'dualphaL'; 
                         'alphaR'; 'dalphaR'; 'ualphaR'; 'dualphaR'; 
                         'betaL'; 'dbetaL'; 'ubetaL'; 'dubetaL';
                         'betaR'; 'dbetaR'; 'ubetaR'; 'dubetaR';
                         'Act Work'; 'El Work'; 'El Loss'};
    datStruct.xd_comment = { 'x          [m]         Main body x position'; 
                             'dx         [m/s]       Main body x velocity'; 
                             'y          [m]         Main body y position'; 
                             'dy         [m/s]       Main body y velocity'; 
                             'phi        [rad]       Main body pitch'; 
                             'dphi       [rad/s]     Main body pitch velocity';
                             'alphaL     [rad]       Left hip angle';
                             'dalphaL    [rad/s]     Left hip angular velocity';
                             'ualphaL    [rad]       Left hip actuator angle';
                             'dualphaL   [rad/s]     Left hip actuator angular velocity';
                             'alphaR     [rad]       Right hip angle';
                             'dalphaR    [rad/s]     Right hip angular velocity';
                             'ualphaR    [rad]       Right hip actuator angle';
                             'dualphaR   [rad/s]     Right hip actuator angular velocity';
                             'betaL      [rad]       Left knee angle';
                             'dbetaL     [rad/s]     Left knee angular velocity';
                             'ubetaL     [rad]       Left knee actuator angle';
                             'dubetaL    [rad/s]     Left knee actuator angular velocity';
                             'betaR      [rad]       Right knee angle';
                             'dbetaR     [rad/s]     Right knee angular velocity';
                             'ubetaR     [rad]       Right knee actuator angle';
                             'dubetaR    [rad/s]     Right knee actuator angular velocity';
                             'posActWork [J]         Positive actuator work (must initially remain 0)'; 
                             'posElWork  [J]         Positive electrical work (must initially remain 0)'; 
                             'totElLoss  [J]         Total electrical losses (must initially remain 0)'};
                         
    datStruct.sd          = sd;
    datStruct.sd_sca      = ones(25,1);
    
    % States are unlimited apart from the body height (which is bounded by [0.6,0.2])
    % pitch ([-pi/4,pi/4]), hip joints ([-pi/2, pi/2]),  knee actuators 
    % ([-3*pi/4, pi/4]), and knee joints ([-3*pi/4, -0.05]).
    datStruct.sd_min          = -100*ones(25,1);
    datStruct.sd_min(3)       = 0.2; %y
    datStruct.sd_min(5)       = -pi/4; %phi
    datStruct.sd_min([7,11])  = -pi/2; %alphaL,alphaR
    datStruct.sd_min(15:2:21) = -3*pi/4; %betaL,ubetaL,betaR,ubetaR
    
    datStruct.sd_max          = 100*ones(25,1);
    datStruct.sd_max(3)       = 0.6; %y
    datStruct.sd_max(5)       = pi/4; %phi
    datStruct.sd_max([7,11])  = pi/2; %alphaL,alphaR
    datStruct.sd_max([15,19]) = -0.05; %betaL,betaR
    datStruct.sd_max([17,21]) = pi/4; %ubetaL,betaR
    
    % Some states fixed at the start of the simulation:
    datStruct.sd_fix_init = zeros(25,1);
    datStruct.sd_fix_init(1) = 1; % Fix initial x position
    datStruct.sd_fix_init(23:end) = 1; % Fix work terms
    
    % Everything can change thereafter:
    datStruct.sd_fix      = zeros(25,1);
    
    %% Parameters:
    datStruct.p_name = {'cost_fct_sel';
                        'v_avg';
                        'sigma';
                        'const_sel'};
        
    datStruct.p_comment = { 'cost_fct_sel  [0; 1; 2; 3]       Cost function selector (1 = posActWorkCOT; 2 = posElWorkCOT; 3 = totElLossCOT)';
                            'v_avg         [m/s]              Enforced average speed';
                            'sigma         [rad]              Knee spring regularization constant';
                            'const_sel     [binary flags]     1: GroundClearance 2: TorqueLimits 4: SpeedLimits'};
    datStruct.p      = p;
    datStruct.p_sca  = ones(4,1);  %scale factor
    datStruct.p_min  = -100*ones(4,1);
    datStruct.p_max  = 100*ones(4,1);
    % All Parameters are fixed by default
    datStruct.p_fix  = ones(4,1);

    %% Actuation
    datStruct.u_name = {'TalphaL'; 'TalphaR'; 'TbetaL'; 'TbetaR'};
    datStruct.u_comment = {'TalphaL     [N*m]   Left hip motor torque';
                           'TalphaR     [N*m]   Right hip motor torque';
                           'TbetaL      [N*m]   Left knee motor torque';
                           'TbetaR      [N*m]   Right knee motor torque'};
    datStruct.u_type = 2*ones(4,1); % piecewise linear
    datStruct.u      = u;
    datStruct.u_sca  = ones(4,1);
    
    % This limit should never go into effect, since the actuators are
    % limited within the dynamic functions. Yet, it helps stabilize the
    % optimization
    datStruct.u_min  = -20.0*ones(4,1);
    datStruct.u_max  = -datStruct.u_min;
    datStruct.u_fix  = zeros(4,1);

    %% Misc entries
    datStruct.rd_scaStart = 1.0;
    datStruct.rd_scaEnd   = 1.0;
    datStruct.rc_sca      = 1.0;

    datStruct.of_name = 'Energy';
    datStruct.of_sca  = 0.5;
    datStruct.of_min  = 0.0;
    datStruct.of_max  = 1.0;
    datStruct.nhist   = 100;

    datStruct.libhessian = 'hess_update';
    datStruct.libsolve   = 'solve_slse';
    datStruct.libcond    = 'cond_std';
    datStruct.libtchk    = 'tchk';
    datStruct.libmssqp   = 'mssqp_standard';
    datStruct.libeval    = 'eval_ind';
    datStruct.libqps     = 'qps_qpopt';
    datStruct.libplot    = 'plot_pgplot';

    datStruct.options_acc           = 1e-6;
    datStruct.options_ftol          = -1.0;
    datStruct.options_itol          = -1.0;
    datStruct.options_rfac          = 0.0;
    datStruct.options_levmar        = 0.0;
    datStruct.options_qp_featol     = 1.0e-8;
    datStruct.options_qp_relax      = 1.1;
    datStruct.options_nhtopy        = 0;
    datStruct.options_frstart       = 0;
    datStruct.options_frmax         = 0;
    datStruct.options_itmax         = 1000;
    datStruct.options_plevel_screen = 0;
    datStruct.options_plevel_file   = 1;
    datStruct.options_plevel_matlab = 0;
    datStruct.options_bflag         = -1;
    datStruct.options_qp_itmax      = 10000;
    datStruct.options_qp_expand     = 99999999;
    datStruct.options_sflag         = 0;
    datStruct.options_options_wflag = 0;
    datStruct.options_cflag         = 0;
    datStruct.options_output_ps     = 0;
    datStruct.options_output_gif    = 0;
end