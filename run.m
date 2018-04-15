function varargout = run(numSteps, choice, pauseLen, da, method)
% RUN PS2 EKF Feature-Based SLAM
%   RUN(ARG)
%   RUN(ARG, CHOICE, PAUSELEN)
%   RUN(ARG, CHOICE, PAUSELEN, DA)
%      ARG - is either the number of time steps, (e.g. 100 is a complete
%            circuit) or a data structure from a previous run.
%      CHOICE - is either 'sim' or 'vp' for simulator or Victoria Park
%               data set, respectively. Or 'prevp' preprocessed Victoria
%               Park data set.
%      PAUSELEN - set to `inf`, to manually pause, o/w # of seconds to wait
%                 (e.g., 0.3 is the default)
%      DA - data assocation, is one of either:
%           'known' - only available in simulator
%           'nn'    - incremental maximum likelihood nearest neighbor
%           'nndg'  - nn double gate on landmark creation
%                     (throws away ambiguous observations)
%           'jcbb'  - joint compatability branch and bound
%
%   DATA = RUN(ARG, CHOISE, PAUSELEN, DA)
%      DATA - is an optional output and contains the data array generated
%             and/or used during the simulation.
%
%   Note: more parameters can be controlled in the run.m file itself via
%   fields of the Param structure.

%   (c) 2009-2015
%   Ryan M. Eustice
%   University of Michigan
%   eustice@umich.edu

addpath('./slamsim');
addpath('./vicpark');
addpath('./seans_nn_functions')

if ~exist('pauseLen', 'var') || isempty(pauseLen)
    pauseLen = [];
end

clear global Param State Data;
global Param;
global State;
global Data;


% select which data association method to use in ekfupdate.m, choices are:
%   known - only available in simulator
%   nn    - incremental maximum likelhood nearest neighbor
%   nndg  - nn double gate on landmark creation (throws away ambiguous observations)
%   jcbb  - joint compatability branch and bound
if ~exist('da','var') || isempty(da)
    da = 'known';
end
Param.dataAssociation = da;


if ~exist('method','var') || isempty(method)
   method = 'EKF' 
end%% select which method to use
% EKF - EKF slam
% iSAM - iSAM slam
% EKFISAM = comare EKF and ISAM together on the sim env
Param.method = method

% select which update method to use in ekfupdate.m, choices are:
%   batch  - batch updates
%   seq    - sequential updates
% Param.updateMethod = 'batch';
Param.updateMethod = 'seq';

% size of bounding box for VP data set plotting
Param.bbox = 0; % bbox = 20 [m] speeds up graphics

% Structure of global State variable
%===================================================
State.Ekf.t     = 0;          % time
State.Ekf.mu    = zeros(3,1); % robot initial pose
State.Ekf.Sigma = zeros(3,3); % robot initial covariance
State.Ekf.iR    = 1:3;        % 3 vector containing robot indices
State.Ekf.iM    = [];         % 2*nL vector containing map indices
State.Ekf.iL    = {};         % nL cell array containing indices of landmark i
State.Ekf.sL    = [];         % nL vector containing signatures of landmarks
State.Ekf.nL    = 0;          % scalar number of landmarks
State.Ekf.rPoses = [];        % save all robot poses.
%===================================================

% Structure of global State variable for iSAM
%====================================================
State.iSAM.t = 0; % time step
State.iSAM.A = []; % we may not need to keep this matrix
State.iSAM.b = []; 
State.iSAM.R = []; % upper triangle matrix.
State.iSAM.d = [];
State.iSAM.rM = []; % robot poses indices in R matrix, 1x3nR
State.iSAM.nR = 0; % number of robot poses.
State.iSAM.nL = 0; % number of landmarks
State.iSAM.lM = []; % landmark indices in R matrix. 1x2nL landmark indices are different from signatures. The signature is the id of the landmark if given known data association. landmark indiex is given by the order of creating the landmarks.
State.iSAM.sL = []; % nL vector containning signatures of landmarks
State.iSAM.iL = {}; % nL cell array containing indices of landmark i
State.iSAM.rlin = []; % robot pose linearization points
State.iSAM.llin = []; % landmark linearization points
State.iSAM.Lu = []; % block diag lu matrix for robot
State.iSAM.Lz = []; % block diag lz matrix for landmark
State.iSAM.uList = []; % list of controls
State.iSAM.zList = []; % list of observations
State.iSAM.max_node = 0; % the maximum index of the node
State.iSAM.rNode = []; % the node index of the robot poses
State.iSAM.lNode = []; % the node index of the landmarks
%====================================================

switch lower(choice)
    case 'sim'
        Param.vp=false;
        Data = runsim(numSteps,pauseLen);
        if nargout > 0
            varargout{1} = Data;
        end
    case 'vp'
        Param.vp=true;
        runvp(numSteps,pauseLen);
    case 'prevp'
        Param.vp=true;
        runpreprocessedvp(numSteps, pauseLen);
    otherwise
        error('unrecognized selection "%s"', choice);
end
