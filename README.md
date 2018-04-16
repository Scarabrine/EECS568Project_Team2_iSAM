# EECS568Project_Team2_iSAM
EECS 568 Final project for team 2 on iSAM implementation in Victoria Park
Content:
A. iSAM and SAM
  1. run.m
     run is to start the SLAM. The input can define the process's dataset (choice), SLAM (method), total steps 
     (numSteps), data association (da) method and pause length (pauseLen). 
     RUN(ARG):
        RUN(ARG, CHOICE, PAUSELEN)
        RUN(ARG, CHOICE, PAUSELEN, DA)
        ARG - is either the number of time steps, (e.g. 100 is a complete
              circuit) or a data structure from a previous run.
        CHOICE - is either 'sim', 'vp' or 'vp_pre' for simulator, Victoria Park data set or the Victoria Park data set pre- 
                 processed by MIT respectively. Note*: 'vp' is not applicable for SAM and iSAM. 'vp_pre' is for SAM and iSAM.
        PAUSELEN - set to `inf`, to manually pause, o/w # of seconds to wait
                   (e.g., 0.3 is the default)
        DA - data assocation, is one of either:
             'known' - only available in simulator
             'nn'    - incremental maximum likelihood nearest neighbor
             'nndg'  - nn double gate on landmark creation
                       (throws away ambiguous observations)
             'jcbb'  - joint compatability branch and bound

        DATA = RUN(ARG, CHOISE, PAUSELEN, DA)
        DATA - is an optional output and contains the data array generated
               and/or used during the simulation.
     Recommended command: run(500,'prevp',[],'known','SAM'); -> Run SAM in prevp
                          run(2000,'prevp',[],'known','iSAM'); -> Run iSAM in prevp
  2. runsim.m
     runsim.m is called by run.m function which deals with the simluation's dataset when the 'Choice' is 'sim'. This function   
     has no output and input are total steps and pause length which are not supposed to be changed directly. For the 'isam' 
     case, iSAMOdometryUpdate_sim.m, iSAMMeasurementUpdate.m and iSAMRelinearization.m are called. For 'ekf' case,  
     ekfpredict_sim.m, ekfupdate.m are called. There is no implementation for 'sam' case in 'ekf'. Note*, to compare the 
     performance of EKF and iSAM in simulation dataset, there is the case 'ekfisam' which runs the ekf and isam at the same 
     time based on the same noise.

  3. runvp.m
     runvp.m is called by run.m function which deals with the victoria park's dataset when the 'Choice' is 'vp'. This function   
     has no output and input are total steps and pause length which are not supposed to be changed directly. 
     iSAMOdometryUpdate_vp.m, iSAMMeasurementUpdate.m and iSAMRelinearization.m are called for 'sim' case. ekfupdate.m, 
     ekfpredict_vp.m are called for 'ekf' case. For 'sam' case, there is no option for SAM in the normal Victoria Park 
     dataset.
     
  4. runpreprocessedvp.m is called by run.m function which deals with the preprocessed victoria park's dataset when the 
     'Choice' is 'prevp'. This function has no output and input are total steps and pause length which are not supposed to be 
     changed directly. In this function, the iSAMOdometryUpdate_pre.m, iSAMMeasurementUpate_pre.m and 
     iSAMRelinearization_pre.m are called for 'isam' case and ekf_pre_predict.m, ekf_pre_update.m are called for 'ekf' case.
     SAMOdometryInit_pre.m and SAMMeasurementInit_pre.m are called for 'sam' case
     
  1. Details of iSAM
  the iSAM mainly has three subfunctions: odemotry update and measurement update are to add new measurement and control input   
  to the iSAM. relinearization is to relinearize the motion and measurement models. The initilization is to get start with 
  iSAM. Each of them is written in the simulation dataset and and preprocessed vp dataset.

The JCBB implementation relies on four files:
1. da_jcbb.m
2. JCBB_R.m
3. joint_compability.m
4. individual_compatibility.m

Modify the JCBB.alpha in da_jcbb.m to define the significance level needed for assessing whether a measurement is associated with a new landmark. Currently the alpha level is set at .001. The Mahalanobis Distance threshold is chi2inv(1-JCBB.alpha,DOF), where DOF = (Number of measuremeents in current hypothesis) * 2. You can change the Mahalanobis distance threshold in joint_compatibility.m.

The function da_jcbb.m calls the recursive function JCBB_R.m, which implements the branch and bound joint-compability data association approach described in Jose Neira's paper.

joint_compatibility.m is a function that takes in a set of measurements and a hypothesis of the landmark ids associated with those measurements. It returns true if the assignment hypothesis is compatible--based on the alpha threshold.

individual_compatibility.m is a special case function that checks if a single measurement is individually compatible with a single hypothesis of the landmark associated with that measurement.
