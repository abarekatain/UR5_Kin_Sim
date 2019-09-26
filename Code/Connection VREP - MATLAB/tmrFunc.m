function tmrFunc(hObj, eventdata)

Q = hObj.UserData;
q1 = Q{1};
q2 = Q{2};
q3 = Q{3};
q4 = Q{4};
q5 = Q{5};
q6 = Q{6};
clientID = Q{7};
vrep = Q{8};
t = Q{9};

handle = zeros(6,1);
for i = 1:6
    %Get Joint Handles
    [returnCode,handle(i)]=vrep.simxGetObjectHandle(clientID,['UR5_joint',num2str(i)],vrep.simx_opmode_blocking);
end

i = hObj.TasksExecuted*50;
s = size(q1,1);
if i > s
    stop(t)
else
    vrep.simxSetJointPosition(clientID,handle(1),q1(i),vrep.simx_opmode_oneshot);
    vrep.simxSetJointPosition(clientID,handle(2),q2(i) - pi/2,vrep.simx_opmode_oneshot);
    vrep.simxSetJointPosition(clientID,handle(3),q3(i),vrep.simx_opmode_oneshot);
    vrep.simxSetJointPosition(clientID,handle(4),q4(i) + pi/2,vrep.simx_opmode_oneshot);
    vrep.simxSetJointPosition(clientID,handle(5),q5(i),vrep.simx_opmode_oneshot);
    vrep.simxSetJointPosition(clientID,handle(6),q6(i),vrep.simx_opmode_oneshot);
end

%Set Joint Positions


% vrep.simxAddStatusbarMessage(clientID,['time : ' num2str(i)],vrep.simx_opmode_oneshot);


end