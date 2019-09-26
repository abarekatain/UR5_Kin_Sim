%% Initialization
% q1 = q(1);
% q2 = q(2);
% q3 = q(3);
% q4 = q(4);
% q5 = q(5);
% q6 = q(6);

%% Create a connection To V-REP using remote API
vrep = remApi('remoteApi');
vrep.simxFinish(-1);
clientID=vrep.simxStart('127.0.0.1',19999,true,true,5000,5);    

if clientID > -1
    t = timer('ExecutionMode', 'FixedRate','Period', 0.05,'TimerFcn', {@tmrFunc}, 'StartDelay', 1);
    Q = {q1 q2 q3 q4 q5 q6 clientID vrep t};
    t.UserData = Q;
    start(t)
%     pause(0.9)
%     stop(t)
end


%% 
% if clientID > -1
%     disp('Connected!')
%     handle = zeros(6,1);
%     for i = 1:6
%         %Get Joint Handles
%         [returnCode,handle(i)]=vrep.simxGetObjectHandle(clientID,['UR5_joint',num2str(i)],vrep.simx_opmode_blocking);
%     end
% 
%     t = tout;
%     animate_step = 50;
%     for i = 1:animate_step:size(t)
%         %Calculate Joint Positions
%         %     jnt1 = [0; 0; 89.2];
%         %     jnt2 = [-425.0.*cos(q1).*cos(q2); -425.0.*cos(q2).*sin(q1); 425.0.*sin(q2) + 89.2];
%         %     jnt3 = [392.0.*cos(q1).*sin(q2).*sin(q3) - 392.0.*cos(q1).*cos(q2).*cos(q3) - 425.0.*cos(q1).*cos(q2); 392.0.*sin(q1).*sin(q2).*sin(q3) - 425.0.*cos(q2).*sin(q1) - 392.0.*cos(q2).*cos(q3).*sin(q1); 425.0.*sin(q2) + 392.0.*cos(q2).*sin(q3) + 392.0.*cos(q3).*sin(q2) + 89.2];
%         %     jnt4 = [392.0.*cos(q1).*sin(q2).*sin(q3) - 425.0.*cos(q1).*cos(q2) - 392.0.*cos(q1).*cos(q2).*cos(q3) - 109.3.*sin(q1); 109.3.*cos(q1) - 425.0.*cos(q2).*sin(q1) + 392.0.*sin(q1).*sin(q2).*sin(q3) - 392.0.*cos(q2).*cos(q3).*sin(q1); 425.0.*sin(q2) + 392.0.*cos(q2).*sin(q3) + 392.0.*cos(q3).*sin(q2) + 89.2];
%         %     jnt5 = [94.75.*sin(q4).*(cos(q1).*cos(q2).*cos(q3) - 1.0.*cos(q1).*sin(q2).*sin(q3)) - 425.0.*cos(q1).*cos(q2) - 109.3.*sin(q1) + 94.75.*cos(q4).*(cos(q1).*cos(q2).*sin(q3) + cos(q1).*cos(q3).*sin(q2)) - 392.0.*cos(q1).*cos(q2).*cos(q3) + 392.0.*cos(q1).*sin(q2).*sin(q3); 109.3.*cos(q1) - 425.0.*cos(q2).*sin(q1) - 94.75.*sin(q4).*(sin(q1).*sin(q2).*sin(q3) - 1.0.*cos(q2).*cos(q3).*sin(q1)) + 94.75.*cos(q4).*(cos(q2).*sin(q1).*sin(q3) + cos(q3).*sin(q1).*sin(q2)) + 392.0.*sin(q1).*sin(q2).*sin(q3) - 392.0.*cos(q2).*cos(q3).*sin(q1); 425.0.*sin(q2) + 392.0.*cos(q2).*sin(q3) + 392.0.*cos(q3).*sin(q2) + 94.75.*cos(q4).*(cos(q2).*cos(q3) - 1.0.*sin(q2).*sin(q3)) - 94.75.*sin(q4).*(cos(q2).*sin(q3) + cos(q3).*sin(q2)) + 89.2];
%         %     jnt6 = [94.75.*sin(q4).*(cos(q1).*cos(q2).*cos(q3) - 1.0.*cos(q1).*sin(q2).*sin(q3)) - 425.0.*cos(q1).*cos(q2) - 82.5.*sin(q5).*(cos(q4).*(cos(q1).*cos(q2).*cos(q3) - 1.0.*cos(q1).*sin(q2).*sin(q3)) - 1.0.*sin(q4).*(cos(q1).*cos(q2).*sin(q3) + cos(q1).*cos(q3).*sin(q2))) - 82.5.*cos(q5).*sin(q1) - 109.3.*sin(q1) + 94.75.*cos(q4).*(cos(q1).*cos(q2).*sin(q3) + cos(q1).*cos(q3).*sin(q2)) - 392.0.*cos(q1).*cos(q2).*cos(q3) + 392.0.*cos(q1).*sin(q2).*sin(q3); 109.3.*cos(q1) + 82.5.*cos(q1).*cos(q5) - 425.0.*cos(q2).*sin(q1) - 94.75.*sin(q4).*(sin(q1).*sin(q2).*sin(q3) - 1.0.*cos(q2).*cos(q3).*sin(q1)) + 82.5.*sin(q5).*(cos(q4).*(sin(q1).*sin(q2).*sin(q3) - 1.0.*cos(q2).*cos(q3).*sin(q1)) + sin(q4).*(cos(q2).*sin(q1).*sin(q3) + cos(q3).*sin(q1).*sin(q2))) + 94.75.*cos(q4).*(cos(q2).*sin(q1).*sin(q3) + cos(q3).*sin(q1).*sin(q2)) + 392.0.*sin(q1).*sin(q2).*sin(q3) - 392.0.*cos(q2).*cos(q3).*sin(q1); 425.0.*sin(q2) + 392.0.*cos(q2).*sin(q3) + 392.0.*cos(q3).*sin(q2) + 94.75.*cos(q4).*(cos(q2).*cos(q3) - 1.0.*sin(q2).*sin(q3)) - 94.75.*sin(q4).*(cos(q2).*sin(q3) + cos(q3).*sin(q2)) + 82.5.*sin(q5).*(sin(q4).*(cos(q2).*cos(q3) - 1.0.*sin(q2).*sin(q3)) + cos(q4).*(cos(q2).*sin(q3) + cos(q3).*sin(q2))) + 89.2];
%         jnt1 = [0; 0; 89.2];
%         jnt2 = [-425.0.*cos(q1(i)).*cos(q2(i)); -425.0.*cos(q2(i)).*sin(q1(i)); 425.0.*sin(q2(i)) + 89.2];
%         jnt3 = [392.0.*cos(q1(i)).*sin(q2(i)).*sin(q3(i)) - 392.0.*cos(q1(i)).*cos(q2(i)).*cos(q3(i)) - 425.0.*cos(q1(i)).*cos(q2(i)); 392.0.*sin(q1(i)).*sin(q2(i)).*sin(q3(i)) - 425.0.*cos(q2(i)).*sin(q1(i)) - 392.0.*cos(q2(i)).*cos(q3(i)).*sin(q1(i)); 425.0.*sin(q2(i)) + 392.0.*cos(q2(i)).*sin(q3(i)) + 392.0.*cos(q3(i)).*sin(q2(i)) + 89.2];
%         jnt4 = [392.0.*cos(q1(i)).*sin(q2(i)).*sin(q3(i)) - 425.0.*cos(q1(i)).*cos(q2(i)) - 392.0.*cos(q1(i)).*cos(q2(i)).*cos(q3(i)) - 109.3.*sin(q1(i)); 109.3.*cos(q1(i)) - 425.0.*cos(q2(i)).*sin(q1(i)) + 392.0.*sin(q1(i)).*sin(q2(i)).*sin(q3(i)) - 392.0.*cos(q2(i)).*cos(q3(i)).*sin(q1(i)); 425.0.*sin(q2(i)) + 392.0.*cos(q2(i)).*sin(q3(i)) + 392.0.*cos(q3(i)).*sin(q2(i)) + 89.2];
%         jnt5 = [94.75.*sin(q4(i)).*(cos(q1(i)).*cos(q2(i)).*cos(q3(i)) - 1.0.*cos(q1(i)).*sin(q2(i)).*sin(q3(i))) - 425.0.*cos(q1(i)).*cos(q2(i)) - 109.3.*sin(q1(i)) + 94.75.*cos(q4(i)).*(cos(q1(i)).*cos(q2(i)).*sin(q3(i)) + cos(q1(i)).*cos(q3(i)).*sin(q2(i))) - 392.0.*cos(q1(i)).*cos(q2(i)).*cos(q3(i)) + 392.0.*cos(q1(i)).*sin(q2(i)).*sin(q3(i)); 109.3.*cos(q1(i)) - 425.0.*cos(q2(i)).*sin(q1(i)) - 94.75.*sin(q4(i)).*(sin(q1(i)).*sin(q2(i)).*sin(q3(i)) - 1.0.*cos(q2(i)).*cos(q3(i)).*sin(q1(i))) + 94.75.*cos(q4(i)).*(cos(q2(i)).*sin(q1(i)).*sin(q3(i)) + cos(q3(i)).*sin(q1(i)).*sin(q2(i))) + 392.0.*sin(q1(i)).*sin(q2(i)).*sin(q3(i)) - 392.0.*cos(q2(i)).*cos(q3(i)).*sin(q1(i)); 425.0.*sin(q2(i)) + 392.0.*cos(q2(i)).*sin(q3(i)) + 392.0.*cos(q3(i)).*sin(q2(i)) + 94.75.*cos(q4(i)).*(cos(q2(i)).*cos(q3(i)) - 1.0.*sin(q2(i)).*sin(q3(i))) - 94.75.*sin(q4(i)).*(cos(q2(i)).*sin(q3(i)) + cos(q3(i)).*sin(q2(i))) + 89.2];
%         jnt6 = [94.75.*sin(q4(i)).*(cos(q1(i)).*cos(q2(i)).*cos(q3(i)) - 1.0.*cos(q1(i)).*sin(q2(i)).*sin(q3(i))) - 425.0.*cos(q1(i)).*cos(q2(i)) - 82.5.*sin(q5(i)).*(cos(q4(i)).*(cos(q1(i)).*cos(q2(i)).*cos(q3(i)) - 1.0.*cos(q1(i)).*sin(q2(i)).*sin(q3(i))) - 1.0.*sin(q4(i)).*(cos(q1(i)).*cos(q2(i)).*sin(q3(i)) + cos(q1(i)).*cos(q3(i)).*sin(q2(i)))) - 82.5.*cos(q5(i)).*sin(q1(i)) - 109.3.*sin(q1(i)) + 94.75.*cos(q4(i)).*(cos(q1(i)).*cos(q2(i)).*sin(q3(i)) + cos(q1(i)).*cos(q3(i)).*sin(q2(i))) - 392.0.*cos(q1(i)).*cos(q2(i)).*cos(q3(i)) + 392.0.*cos(q1(i)).*sin(q2(i)).*sin(q3(i)); 109.3.*cos(q1(i)) + 82.5.*cos(q1(i)).*cos(q5(i)) - 425.0.*cos(q2(i)).*sin(q1(i)) - 94.75.*sin(q4(i)).*(sin(q1(i)).*sin(q2(i)).*sin(q3(i)) - 1.0.*cos(q2(i)).*cos(q3(i)).*sin(q1(i))) + 82.5.*sin(q5(i)).*(cos(q4(i)).*(sin(q1(i)).*sin(q2(i)).*sin(q3(i)) - 1.0.*cos(q2(i)).*cos(q3(i)).*sin(q1(i))) + sin(q4(i)).*(cos(q2(i)).*sin(q1(i)).*sin(q3(i)) + cos(q3(i)).*sin(q1(i)).*sin(q2(i)))) + 94.75.*cos(q4(i)).*(cos(q2(i)).*sin(q1(i)).*sin(q3(i)) + cos(q3(i)).*sin(q1(i)).*sin(q2(i))) + 392.0.*sin(q1(i)).*sin(q2(i)).*sin(q3(i)) - 392.0.*cos(q2(i)).*cos(q3(i)).*sin(q1(i)); 425.0.*sin(q2(i)) + 392.0.*cos(q2(i)).*sin(q3(i)) + 392.0.*cos(q3(i)).*sin(q2(i)) + 94.75.*cos(q4(i)).*(cos(q2(i)).*cos(q3(i)) - 1.0.*sin(q2(i)).*sin(q3(i))) - 94.75.*sin(q4(i)).*(cos(q2(i)).*sin(q3(i)) + cos(q3(i)).*sin(q2(i))) + 82.5.*sin(q5(i)).*(sin(q4(i)).*(cos(q2(i)).*cos(q3(i)) - 1.0.*sin(q2(i)).*sin(q3(i))) + cos(q4(i)).*(cos(q2(i)).*sin(q3(i)) + cos(q3(i)).*sin(q2(i)))) + 89.2];
%         
%         
%         %Set Joint Positions
%         vrep.simxSetJointPosition(clientID,handle(1),jnt1,vrep.simx_opmode_oneshot);
%         vrep.simxSetJointPosition(clientID,handle(2),jnt2,vrep.simx_opmode_oneshot);
%         vrep.simxSetJointPosition(clientID,handle(3),jnt3,vrep.simx_opmode_oneshot);
%         vrep.simxSetJointPosition(clientID,handle(4),jnt4,vrep.simx_opmode_oneshot);
%         vrep.simxSetJointPosition(clientID,handle(5),jnt5,vrep.simx_opmode_oneshot);
%         vrep.simxSetJointPosition(clientID,handle(6),jnt6,vrep.simx_opmode_oneshot);
%         
%         vrep.simxAddStatusbarMessage(clientID,['time : ' num2str(i)],vrep.simx_opmode_oneshot);
%     end
%     disp('Finished!')
% else
%     disp('Connection Failed!')
% end
% 
% vrep.simxFinish(clientID);