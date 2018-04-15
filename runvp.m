function runvp(nSteps,pauseLen)

global Param;
global State;
global Data;

if ~exist('nSteps','var') || isempty(nSteps)
    nSteps = inf;
end

if ~exist('pauseLen','var')
    pauseLen = 0; % seconds
end

Data = load_vp_si();

% Initalize Params
%===================================================
% vehicle geometry
Param.a = 3.78; % [m]
Param.b = 0.50; % [m]
Param.L = 2.83; % [m]
Param.H = 0.76; % [m]
Param.initialStateMean = [Data.Gps.x(2), Data.Gps.y(2), 36*pi/180]';
% 2x2 process noise on control input
sigma.vc = 0.02; % [m/s]
sigma.alpha = 2*pi/180; % [rad]
Param.Qu = diag([sigma.vc, sigma.alpha].^2);

% 3x3 process noise on model error
sigma.x = 0.1; % [m]
sigma.y = 0.1; % [m]
sigma.phi = 0.5*pi/180; % [rad]
Param.Qf = diag([sigma.x, sigma.y, sigma.phi].^2);

% 2x2 observation noise
sigma.r = 0.05; % [m]
sigma.beta = 1*pi/180; % [rad]
Param.R = diag([sigma.r, sigma.beta].^2);
%===================================================

% Initialize State
%===================================================
State.Ekf.mu = [Data.Gps.x(2), Data.Gps.y(2), 36*pi/180]';
State.Ekf.Sigma = zeros(3);
State.Ekf.true_mu = [];
State.Ekf.estimated_mu = [];
global AAr;
AAr = [0:360]*pi/360;

% Save Video
makeVideo = false;
alg = 'task3';
% pauseLen = 0.1

if makeVideo
    try
        votype = 'avifile';
        video_name = [alg,'.avi']
        vo = avifile(video_name, 'fps', min(5, 1/pauseLen));
    catch
        votype = 'VideoWriter';
        video_name = [alg]
        vo = VideoWriter(video_name, 'MPEG-4');
        set(vo, 'FrameRate', min(5, 1/pauseLen));
        open(vo);
    end
end

figure(1); clf;
axis equal;
xlim([-120,20])
ylim([-120,20])
hold on
ci = 1; % control index
t = min(Data.Laser.time(1), Data.Control.time(1));

prediction_time_list = []; % first row: prediction time, second row: time
update_time_list = []; % first row: update time, second row: time
nL_list = []; % first row: number of landmarks, second row: time
%% iSAM intilization
iSAMInit();

for k=1:min(nSteps, length(Data.Laser.time))
    
    while (Data.Control.time(ci) < Data.Laser.time(k))
       % control available
       dt = Data.Control.time(ci) - t;
       t = Data.Control.time(ci);
       u = [Data.Control.ve(ci), Data.Control.alpha(ci)]';
       tic;
       switch lower(Param.method)
           case 'isam'
               iSAMOdometryUpdate_vp(u,dt);
           case 'ekf'
               ekfpredict_vp(u, dt);
           case 'sam'
       end
       pred_time = toc;
       prediction_time_list = [prediction_time_list, [pred_time;t]];
       ci = ci+1;
    end
    
    % observation available
    dt = Data.Laser.time(k) - t;
    t = Data.Laser.time(k);
    z = detectTreesI16(Data.Laser.ranges(k,:));
    
    % Get ground truth gps data
    index = find(Data.Gps.time <= t);
    if (~isempty(index))
        gps_index = index(end);
        State.Ekf.true_mu = [State.Ekf.true_mu, [Data.Gps.x(gps_index);Data.Gps.y(gps_index)]];
    end
    % Modify observation, as the beearing has pi/2 difference.
    tic
    z_pi = z;
    z_pi(2,:) = z_pi(2,:) - pi/2;
    for ind = 1:1:size(z_pi,2)
        z_pi(2,ind) = minimizedAngle(z_pi(2,ind));
    end
    switch lower(Param.method)
        case 'isam'
            iSAMMeasurementUpdate(z_pi);
            if mod(k,10) ==1
                iSAMRelinearization();
            end
        case 'ekf'
            ekfupdate(z_pi);
        case 'sam'
    end
    update_time = toc;
    update_time_list = [update_time_list, [update_time;t]];
    nL_list = [nL_list,[State.Ekf.nL;t]];
    % load estimated robot pose path
    State.Ekf.estimated_mu = [State.Ekf.estimated_mu, State.Ekf.mu(1:2)];
    
    
    switch lower(Param.method)
        case 'isam'
            %% plot robot traj
            plot(State.iSAM.rlin(1:3:end), State.iSAM.rlin(2:3:end),'r-');
            hold on
            %% plot landmark pose
            plot(State.iSAM.llin(1:2:end), State.iSAM.llin(2:2:end),'r*');
            xlim([-120,20])
            ylim([-120,20])
            hold on
        case 'ekf'
            %% Draw the robot uncertainty
            plotcov2d(State.Ekf.mu(1), State.Ekf.mu(2), State.Ekf.Sigma(1:2,1:2), 'red',false,[],[],3)
            %% Draw landmark uncertainty
            for ind = 1:1:State.Ekf.nL
                plotcov2d(State.Ekf.mu((3+2*ind-1)), State.Ekf.mu((3+2*ind)), State.Ekf.Sigma((3+2*ind-1):(3+2*ind),(3+2*ind-1):(3+2*ind)), 'm',false,[],[],3)
            end
            doGraphics(z);
            hold on;
        case 'sam'
    end
    plot(State.Ekf.true_mu(1,:), State.Ekf.true_mu(2,:),'k-') %% Plot the robot true path
    hold on
%     plot(State.Ekf.estimated_mu(1,:), State.Ekf.estimated_mu(2,:),'b-')%% Plot the estimated robot pose path
%     plotcov2d( State.Ekf.mu(1), State.Ekf.mu(2), State.Ekf.Sigma, 'blue', 0, 'blue', 0, 3); % Plot robot uncertainty

    hold off;
        drawnow;

    if pauseLen > 0
        pause(pauseLen);
    end
    
        %% Save video
    if makeVideo
        F = getframe(gcf);
        switch votype
          case 'avifile'
            vo = addframe(vo, F);
          case 'VideoWriter'
            writeVideo(vo, F);
          otherwise
            error('unrecognized votype');
        end
    end
end
if makeVideo
    fprintf('Writing video...');
    switch votype
      case 'avifile'
        vo = close(vo);
      case 'VideoWriter'
        close(vo);
      otherwise
        error('unrecognized votype');
    end
    fprintf('done\n');
end
disp('finish exp')
%% plot analysis figure
h1 = figure
hold on
plot(prediction_time_list(2,:), prediction_time_list(1,:),'b-','DisplayName','Prediction time')
plot(update_time_list(2,:),update_time_list(1,:),'r-','DisplayName','Update time')
legend('Location','NorthWest')
xlabel('time')
ylabel('CPU \Delta time')
saveas(h1,'task3_time','epsc')
h2 = figure
hold on
plot(nL_list(2,:), nL_list(1,:),'m-','DisplayName','Number of landmarks')
xlabel('time')
ylabel('Number of landmarks')
legend('Location','NorthWest')
saveas(h2,'task3_nL','epsc')

%==========================================================================
function doGraphics(z)
% Put whatever graphics you want here for visualization
%
% WARNING: this slows down your process time, so use sparingly when trying
% to crunch the whole data set!

global Param;
global State;
% plot the robot and 3-sigma covariance ellipsoid
plotbot(State.Ekf.mu(1), State.Ekf.mu(2), State.Ekf.mu(3), 'black', 1, 'blue', 1);
hold on
plotcov2d( State.Ekf.mu(1), State.Ekf.mu(2), State.Ekf.Sigma, 'blue', 0, 'blue', 0, 3);
% restrict view to a bounding box around the current pose
% BB=20;
% axis([[-BB,BB]+State.Ekf.mu(1), [-BB,BB]+State.Ekf.mu(2)]);
xlim([-110 20])
ylim([-110 20])

% project raw sensor detections in global frame using estimate pose
xr = State.Ekf.mu(1);
yr = State.Ekf.mu(2);
tr = State.Ekf.mu(3);
for k=1:size(z,2)
    r = z(1,k);
    b = z(2,k);
    xl = xr + r*cos(b+tr-pi/2);
    yl = yr + r*sin(b+tr-pi/2);
    plot([xr; xl], [yr; yl],'r',xl,yl,'r*');
end

% Plot landmark
for ind = 1:1:State.Ekf.nL
    plotcov2d(State.Ekf.mu((3+2*ind-1)), State.Ekf.mu((3+2*ind)), State.Ekf.Sigma((3+2*ind-1):(3+2*ind),(3+2*ind-1):(3+2*ind)), 'm',false,[],[],3)
end

hold off;

