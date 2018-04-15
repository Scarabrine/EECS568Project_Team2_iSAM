function [] = runpreprocessedvp(nSteps, pauseLen)
% run and test methods for on the preprocessed vp data set
global Param;
global State;
global Data;
if ~exist('pauseLen','var')
    pauseLen = 0.3; % seconds
end

%========================================
%Load preprocessed data
%========================================
% first column: 0-odometry, 1-landmark
% second and third columns: two node index, this row is the factor between these two nodes.
% For odometry, column 4:6 is the delta_x, delta_y, delta_theta, can be
% treated as control, H(u) = H(p_(t))^-1 H(p_(t+1))
% For landmark, column 4:5 is the landmark position in the robot current
% pose frame, z = R(r)^-1 (l-r)
load('./VP_Data.mat');
Data = load_vp_si();
%===========================
% Initialize States
%===========================
Param.Inf_sqrt_u = diag([100,500,500]); % Might have problem with this order
Param.Inf_sqrt_z = diag([1.581139, 1.581139]);
Param.R = inv(Param.Inf_sqrt_z*Param.Inf_sqrt_z) ;
Param.M = inv(Param.Inf_sqrt_u*Param.Inf_sqrt_u) *100;
Param.initialStateMean = [Data.Gps.x(2), Data.Gps.y(2), 36*pi/180]';
Param.initialStateSigma = eye(3)/10;
Param.initialStateSigma(3,3) = deg2rad(2);

% SAVE VIDEO
makeVideo = false;
alg = 'task2';
pauseLen = 0.1;
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
xlim([-120,120])
ylim([-120,120])
% hold on

switch lower(Param.method)
    case 'isam'
        %% Method for iSAM
        %% Initialize the model
        iSAMInit_pre();
        tic
        for t = 1:1:nSteps
            State.iSAM.t = t;
            %% Get new observation or control
            if VP_Data(t,1) == 0
               %% It is odometry data
                u = VP_Data(t, 4:6);
                nodes = VP_Data(t, 2:3);
                
                iSAMOdometryUpdate_pre(u, nodes);
            else
               %% It is landmark data
                z = VP_Data(t, 4:5);
                nodes = VP_Data(t, 2:3);
                iSAMMeasurementUpate_pre(z, nodes);
            end
            if mod(t,100) ==0
                iSAMRelinearization_pre();
        %         save('vp_all_pre_data_fast_given_rotation.mat','State');
            end

                %% plot robot traj
            plot(State.iSAM.rlin(1:3:end), State.iSAM.rlin(2:3:end),'r-');
            xlim([-200,100])
            ylim([-100,200])
            hold on
            %% plot landmark pose
            plot(State.iSAM.llin(1:2:end), State.iSAM.llin(2:2:end),'r*');
            xlim([-200,100])
            ylim([-100,200])
            hold off
            pause(0.01)
        end
        toc
        save('vp_all_pre_data.mat','State');
%         tic
    case 'ekf'
        tic
        State.iSAM.lNode = [];
        State.iSAM.max_node = 0;
        State.Ekf.mu = Param.initialStateMean;
        State.Ekf.Sigma = Param.initialStateSigma;
        PastPose = [];
        for t = 1:1:nSteps
            %% Get new observation or control
            if VP_Data(t,1) == 0
               %% It is odometry data
                u = VP_Data(t, 4:6);
                nodes = VP_Data(t, 2:3);
                ekf_pre_predict(u);
            else
               %% It is landmark data
                z = VP_Data(t, 4:5);
                nodes = VP_Data(t, 2:3);
                ekf_pre_update(z, nodes);
            end
            % plot robot traj
            PastPose = [PastPose,State.Ekf.mu(1:3)];
            clf;
            plot(PastPose(1,:), PastPose(2,:),'r-');
            xlim([-200,100])
            ylim([-100,200])
            hold on
            % plot landmark pose
            plot(State.Ekf.mu(4:2:end), State.Ekf.mu(5:2:end),'r*');
            xlim([-200,100])
            ylim([-100,200])
%            variance ellipse of last step
            plotcov2d(State.Ekf.mu(1), State.Ekf.mu(2), State.Ekf.Sigma(1:2,1:2), 'g', 1, 'g', 0.1, 3);
            for i=4:2:length(State.Ekf.mu)
                plotcov2d(State.Ekf.mu(i), State.Ekf.mu(i+1), State.Ekf.Sigma(i:i+1,i:i+1), 'y', 1, 'y', 0.1, 3);
            end
            pause(0.01)
        end
        toc
        save('vp_all_pre_data_ekf.mat','State','Param','PastPose')

    case 'sam'
        %% Initialize model
        iSAMInit_pre();
        tic
        for t = 1:1:nSteps
            State.iSAM.t = t
            if VP_Data(t,1) == 0
               %% It is odometry data
                u = VP_Data(t, 4:6);
                nodes = VP_Data(t, 2:3);
                SAMOdometryInit_pre(u, nodes);
            else
               %% It is landmark data
                z = VP_Data(t, 4:5);
                nodes = VP_Data(t, 2:3);
                SAMMeasurementInit_pre(z, nodes);
            end
        end
        disp('finish SAM init:')
        toc
        SAM_pre()
        disp('finish SAM:')
        toc
        plot(State.iSAM.rlin(1:3:end), State.iSAM.rlin(2:3:end),'r-');
        xlim([-200,100])
        ylim([-100,200])
        hold on
        %% plot landmark pose
        plot(State.iSAM.llin(1:2:end), State.iSAM.llin(2:2:end),'r*');
        xlim([-200,100])
        ylim([-100,200])
        hold off
        save('vp_all_pre_data_SAM_result.mat','State')
end
end

