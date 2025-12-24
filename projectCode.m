clc; clear;

load("real_data.mat");

%converting into double arrays
latitudeValues = Position.latitude;
longitudeValues = Position.longitude;
altitudeValues = Position.altitude;

%getting accelerations
accX = Acceleration.X;
accY = Acceleration.Y;
accZ = Acceleration.Z;

%setting origin point as the location of first recorded data
origin = [latitudeValues(1), longitudeValues(1), altitudeValues(1)];

%changing into the coordinates
[xEast, yNorth, zUp] = latlon2local(latitudeValues, longitudeValues, altitudeValues, origin);

%finding position times
DTPos = datetime(Position.Timestamp);
% hours can be ignored cause everything happened in the same hour
minPos = minute(DTPos);
minInSecPos = minPos*60;
secPos = second(DTPos);
timeInSecPos = round(minInSecPos + secPos)- 19; %for some reason it starts from 19

%finding acceleration times
DTAcc = datetime(Acceleration.Timestamp);
% hours can be ignored cause everything happened in the same hour
minAcc = minute(DTAcc);
minInSecAcc = minAcc*60;
secAcc = second(DTAcc);
timeInSecAcc = round(minInSecAcc + secAcc);

numOfDataPos = length(xEast);
numOfDataAcc = length(accX);


%% task 1 - plotting X, Y, and Z coordiantes as functions of time
figure(1)
subplot(3,1,1)
plot(timeInSecPos,xEast, '-r')
title('X trajectory')
xlabel('time (s)')
ylabel('X coordinates (m)')
subplot(3,1,2)
plot(timeInSecPos,yNorth, '-')
title('Y trajectory')
xlabel('time (s)')
ylabel('Y coordinates (m)')
subplot(3,1,3)
plot(timeInSecPos,zUp, '-g')
title('Z trajectory')
xlabel('time (s)')
ylabel('Z coordinates (m)')

%% task 2 - plotting x-y trajectory

figure(2)
plot(xEast, yNorth, '-', 'LineWidth', 1);
title('X-Y trajectory');
xlabel('x');
ylabel('y');

%% task 2 - overlaying trajectory on map
map2 = imread("mapForDSA3.png");

%scaled x and y to overlay on map 
    %make better
xScaled = 500 + 1.68*xEast;
yScaled = 900 - 1.5*yNorth;
 
figure(3)
imshow(map2);
hold on
plot(xScaled, yScaled, '-r', 'LineWidth', 2);
hold off
axis square


%% task 3 - distance traveled

total_distance = 0;
total_distanceZ = 0;

for i = 1: (numOfDataPos - 1) 
    distance = sqrt( (xEast(i)- xEast(i+1)).^2 + (yNorth(i) - yNorth(i+1)).^2 ); %include z
    distanceZ = sqrt( (xEast(i)- xEast(i+1)).^2 + (yNorth(i) - yNorth(i+1)).^2 + (zUp(i) - zUp(i+1)).^2 ); %include z

    total_distance = total_distance + distance;
    total_distanceZ = total_distanceZ + distanceZ;

end 

error = (total_distanceZ - total_distance)/(total_distance)*100;

disp(total_distance);
disp(total_distanceZ);
output = ['The error is: ', num2str(error), '%'];
disp(output)


%% task 4 - finding velocities from position

%unfiltered velocities
initialVelX = (gradient(xEast));
initialVelY = (gradient(yNorth));
initialVelZ = (gradient(zUp));

%plotting unfiltered velocities
figure(4)
subplot(3,1,1)
plot(timeInSecPos,initialVelX, '-r')
title('Non-filtered X Velocity')
xlabel('time')
ylabel('velocity X')
subplot(3,1,2)
plot(timeInSecPos,initialVelY, '-')
title('Non-filtered Y Velocity')
xlabel('time')
ylabel('velocity Y')
subplot(3,1,3)
plot(timeInSecPos,initialVelZ, '-g')
title('Non-filtered Z Velocity')
xlabel('time')
ylabel('Velocity Z')

%% 
%matlab built in filtering

%filtering using matlab filter
windowSize = 5; %becuase it makes sense 
b = (1/windowSize)*ones(1,windowSize);

filtX = filter(b, 1, initialVelX);
filtY = filter(b, 1, initialVelY);
filtZ = filter(b, 1, initialVelZ);

figure(5)
subplot(3,1,1)
plot(timeInSecPos,filtX, '-r')
title('Filtered once X Velocity')
xlabel('time')
ylabel('velocity X')
subplot(3,1,2)
plot(timeInSecPos,filtY, '-')
title('Filtered once Y Velocity')
xlabel('time')
ylabel('velocity Y')
subplot(3,1,3)
plot(timeInSecPos,filtZ, '-g')
title('Filtered once Z Velocity')
xlabel('time')
ylabel('Velocity Z')


%%
%plotting xy at big offset in velocity

%first big offset (for y) - and do for all - comment on the points in the
%map

%625-747,  1239-1334, and 1574-1646
startTime1 = 625;
endTime1 = 747;

startTime2 = 1239;
endTime2 = 1334;

startTime3 = 1574;
endTime3 = 1646;

start1 = 0; start2 =0; start3 = 0;
end1 = 0; end2 = 0; end3 =0;

%getting i equivelant
for i = 1:length(timeInSecPos)
    if (timeInSecPos(i) == startTime1)
        start1 = i;
    elseif (timeInSecPos(i) == endTime1)
        end1 = i;
    elseif (timeInSecPos(i) == startTime2)
        start2 = i;
    elseif (timeInSecPos(i) == endTime2)
        end2 = i;
    elseif (timeInSecPos(i) == startTime3)
        start3 = i;
    elseif (timeInSecPos(i) == endTime3)
        end3 = i;
    end 
end 

xRange1 = xEast(start1:end1, 1);
yRange1 = yNorth(start1:end1, 1);

xRange2 = xEast(start2:end2, 1);
yRange2 = yNorth(start2:end2, 1);

xRange3 = xEast(start3:end3, 1);
yRange3 = yNorth(start3:end3, 1);

figure(6)
plot(xEast, yNorth, '-k', 'LineWidth', 2)
hold on
plot(xRange1, yRange1, '-m', 'LineWidth', 2)
hold on
plot(xRange2, yRange2, '-g', 'LineWidth', 2)
hold on 
plot(xRange3, yRange3, '-r', 'LineWidth', 2)
axis equal
hold off
legend('x-y trajectory', 'Error 1', 'Error 2', 'Error 3')

%%
%manual filtering for values that are still off
maxSpeedXYZ = 9; %considering maximum speed I could have reached is 9ms-1 (still a stretch)

finalFiltX = filtX;
finalFiltY = filtY;
finalFiltZ = filtZ;

for i = 1:numOfDataPos
    %replacing outliers with decided max value
    if (finalFiltX(i) > maxSpeedXYZ)
        finalFiltX(i) = maxSpeedXYZ;
    elseif(finalFiltX(i) < -maxSpeedXYZ)
        finalFiltX(i) = -maxSpeedXYZ;
    end 
    if (finalFiltY(i) > maxSpeedXYZ)
        finalFiltY(i) = maxSpeedXYZ;
    elseif(finalFiltY(i) < -maxSpeedXYZ)
        finalFiltY(i) = -maxSpeedXYZ;
    end 
    if (finalFiltZ(i) > maxSpeedXYZ)
        finalFiltZ(i) = maxSpeedXYZ;
    elseif(finalFiltZ(i) < -maxSpeedXYZ)
        finalFiltZ(i) = -maxSpeedXYZ;
    end 
end

%%
%plotting final filtered values
figure(7)
subplot(3,1,1)
plot(timeInSecPos,finalFiltX, '-r')
title('Final X Velocity')
xlabel('time')
ylabel('velocity X')
subplot(3,1,2)
plot(timeInSecPos,finalFiltY, '-')
title('Final Y Velocity')
xlabel('time')
ylabel('velocity Y')
subplot(3,1,3)
plot(timeInSecPos,finalFiltZ, '-g')
title('Final Z Velocity')
xlabel('time')
ylabel('Velocity Z')


%plotting together with initial values
figure(8)
subplot(3,1,1)
plot(timeInSecPos,initialVelX, '-r')
hold on 
plot(timeInSecPos,finalFiltX, '-b')
title('X Velocity')
xlabel('time')
ylabel('velocity X')
legend('Not filtered', 'Filtered')
hold off
subplot(3,1,2)
plot(timeInSecPos,initialVelY, '-r')
hold on 
plot(timeInSecPos,finalFiltY, 'b-')
title('Y Velocity')
xlabel('time')
ylabel('velocity Y')
legend('Not filtered', 'Filtered')
hold off
subplot(3,1,3)
plot(timeInSecPos,initialVelZ, '-r')
hold on 
plot(timeInSecPos,finalFiltZ, '-b')
title('Z Velocity')
xlabel('time')
ylabel('Velocity Z')
legend('Not filtered', 'Filtered')
hold off

%%
%plotting all three filtering stages
figure(9)
subplot(3,1,1)
plot(timeInSecPos,initialVelX, '-g')
hold on 
plot(timeInSecPos,filtX, '-k')
hold on
plot(timeInSecPos,finalFiltX, '-r')
title('X Velocity')
xlabel('time')
ylabel('velocity X')
legend('Not filtered', 'Filtered once', 'Final values')
hold off
subplot(3,1,2)
plot(timeInSecPos,initialVelY, '-g')
hold on
plot(timeInSecPos,filtY, '-k')
hold on
plot(timeInSecPos,finalFiltY, 'r-')
title('Y Velocity')
xlabel('time')
ylabel('velocity Y')
legend('Not filtered', 'Filtered', 'Final values') 
hold off
subplot(3,1,3)
plot(timeInSecPos,initialVelZ, '-g')
hold on 
plot(timeInSecPos,filtZ, '-k')
hold on
plot(timeInSecPos,finalFiltZ, '-r')
title('Z Velocity')
xlabel('time')
ylabel('Velocity Z')
legend('Not filtered', 'Filtered', 'Final values') 
hold off

%%
%assigning good name
velX = finalFiltX;
velY = finalFiltY;
velZ = finalFiltZ;

%finding total velocity

unfilteredSpeed = sqrt(initialVelX.^2 + initialVelY.^2 + initialVelZ.^2);

speedInitial = sqrt(velX.^2 + velY.^2 + velZ.^2);

%filtering total speed
speed = filter(b, 1, speedInitial);

figure(10)
plot(timeInSecPos, speed, '-r', 'LineWidth',2);
hold on
plot(timeInSecPos, speedInitial, '-b')
hold off
xlabel('Time (s)')
ylabel('Speed (ms-1)')
title('Total speed')
legend('Filtered', 'Not filtered')


%% task 5 - finding velocities from acceleration
clc;

speedFromAccX = zeros(length(accX), 1);
speedFromAccY = zeros(length(accX), 1);
speedFromAccZ = zeros(length(accX), 1);

timeAccFromZero = timeInSecAcc - 19;

for i = 1:numOfDataAcc
    if (i == 1)
        speedFromAccX(i) = 0;
        speedFromAccY(i) = 0;
        speedFromAccZ(i) = 0;
    else 
        dt = (timeAccFromZero(i) - timeAccFromZero(i - 1));
        speedFromAccX(i) = speedFromAccX(i - 1) + (accX(i))*dt;
        speedFromAccY(i) = speedFromAccY(i - 1) + (accY(i))*dt;
        speedFromAccZ(i) = speedFromAccZ(i - 1) + (accZ(i))*dt;
    end 
end

%finding total acceleration 
totalSpeed = sqrt(speedFromAccX.^2 + speedFromAccY.^2 + speedFromAccZ.^2 );

%plotting
figure(11)
subplot(3,1,1)
plot(timeAccFromZero, speedFromAccX, '-')
title('x')
subplot(3,1,2)
plot(timeAccFromZero, speedFromAccY, '-')
title('y')
subplot(3,1,3)
plot(timeAccFromZero, speedFromAccZ, '-')
title('z')

figure(12)
plot(timeAccFromZero, totalSpeed);
title('Total velocity')
xlabel('Time')
ylabel('Velocity')


%big report analysis
%predicate 1: i am not a rocket 
%since i am not a rocket, it is unrealistc
%use other tests as well


%% task 6.1 - automatically detect movement modes 
clc;

%numOfDataPos = 1391
n = 5;
numOfAverageData = floor(length(speed)/n);
averageSpeed = zeros(numOfAverageData, 1);
averageTime = zeros(numOfAverageData, 1); 

%making average every n 
for i = 1:(numOfAverageData) %because 1391 cant be divided by 5
    speedSum = 0;
    timeSum = 0;
    for j = 1:n
        speedSum = speedSum + speed(n*(i - 1) + j);
        timeSum = timeSum + timeInSecPos(n*(i -1) + j);
    end
    averageTime(i) = round(timeSum/n);
    averageSpeed(i) = speedSum/n; %because we know it is every n seconds
end

%manually adding the two last in each case 
averageTime(end) = timeInSecPos(end);
averageSpeed(end) = sum(speed((end-(n-1)):end ))/n;


%making dt array
for i = 1:numOfAverageData - 1
    dt(i) = averageTime(i+1) - averageTime(i); 
end

maxSpeed = max(averageSpeed); 

%making upper values for different speeds
strollThroughTheParkVal = 1.1; 
powerWalkingVal = 2.0;
lateForClassVal = maxSpeed;

%making new arrays 
strollThroughThePark = zeros(numOfAverageData, 2); %speed, time 
powerWalking = zeros(numOfAverageData, 2); 
lateForClass = zeros(numOfAverageData, 2);

% variables to check how many are in each speed mode
s = 0;
p = 0;
l = 0;
        
for i = 1:numOfAverageData
    if ((averageSpeed(i)) <= strollThroughTheParkVal)
        strollThroughThePark(i,1) = averageSpeed(i);
        strollThroughThePark(i,2) = averageTime(i); 
        s = s+1;
    elseif (averageSpeed(i) <= powerWalkingVal)
        powerWalking(i,1) = averageSpeed(i);
        powerWalking(i,2) = averageTime(i);
        p = p+1;
    elseif ((averageSpeed(i) <= maxSpeed))
        lateForClass(i, 1) = averageSpeed(i);
        lateForClass(i, 2) = averageTime(i);
        l = l+1;
    end 
end 

%displaying how many are in each speed mode
disp(s)
disp(p)
disp(l)

% seperating where is what 
% checking times bc speed could potentially be 0
speedModes = zeros(numOfAverageData, 2);

for i = 1: numOfAverageData
    if (averageSpeed(i) <= strollThroughTheParkVal)
        speedModes(i, 1) = 1;
        speedModes(i, 2) = averageTime(i);
    elseif (averageSpeed(i) <= powerWalkingVal)
        speedModes(i) = 2;
        speedModes(i, 2) = averageTime(i);
    else
        speedModes(i) = 3;
        speedModes(i, 2) = averageTime(i);
    end
end

% if there are single measurements in larger groups of different, add them to the different so its not too crowded
for i = 2:(numOfAverageData - 1) % not the first and last one
    if (speedModes(i - 1) == speedModes(i + 1)) %checking if the previous and next one are the same
        if (speedModes(i) ~= speedModes(i -1)) %checking if the current one is not equal to them
            speedModes(i) = speedModes(i - 1); %changing them
        end 
    end 
end 

%% task 6.1 - plotting

%making new speeds without the times = 0 to scatter plot points
j1 = 1;
j2 = 1;
j3 = 1;
strollThroughTheParkPlot = zeros(s, 2); %s is num of points
powerWalkingPlot = zeros(p, 2); %p is num of points
lateForClassPlot = zeros(l, 2); %l is num of points

for i = 1:numOfAverageData
    if(strollThroughThePark(i, 2) ~= 0)
       strollThroughTheParkPlot(j1, 1) = strollThroughThePark(i, 1); 
       strollThroughTheParkPlot(j1, 2) = strollThroughThePark(i, 2); 
       j1 = j1+1;
    end 
    if(powerWalking(i, 2) ~= 0)
       powerWalkingPlot(j2, 1) = powerWalking(i, 1); 
       powerWalkingPlot(j2, 2) = powerWalking(i, 2); 
       j2 = j2+1;
    end 
    if(lateForClass(i, 2) ~= 0)
       lateForClassPlot(j3, 1) = lateForClass(i, 1); 
       lateForClassPlot(j3, 2) = lateForClass(i, 2); 
       j3 = j3+1;
    end
end 

figure(12)
plot(averageTime, averageSpeed, '-k')
hold on
plot(strollThroughTheParkPlot(:, 2), strollThroughTheParkPlot(:, 1), 'ro')
hold on 
plot(powerWalkingPlot(:, 2), powerWalkingPlot(:, 1), 'bo')
hold on 
plot(lateForClassPlot(:, 2), lateForClassPlot(:, 1), 'go')
hold off
title('Different speed modes')
legend('Total Speed', 'Stroll Through The Park', 'Power Walk', 'Late For Class')
xlabel('Time (s)')
ylabel('Speed (ms-1)')


%% task 6 - converting back to 1391 points

%converting back to all data and getting time stamps to plot all data
%we want an array with: time it starts, mode

speedModesBigArray = zeros(numOfDataPos, 3);
speedModesBigArray(:, 2) = timeInSecPos(:);
speedModesBigArray(:, 3) = speed(:);


for i = 1:length(speedModes)
    for j = 1:n
        speedModesBigArray((i - 1)*n + j, 1) = speedModes(i, 1);
    end 
end 

speedModesBigArray (end, 1) = speedModes(end, 1);

s = 0;
p = 0;
l = 0;
        
for i = 1:numOfDataPos
    if (speedModesBigArray(i, 1) == 1)
        s = s+1;
    elseif (speedModesBigArray(i, 1) == 2)
        p = p+1;
    elseif (speedModesBigArray(i, 1) == 3)
        l = l+1;
    end 
end 

%making scatter plot arrays
j1 = 1;
j2 = 1;
j3 = 1;
mode1 = zeros(s, 2); %s is num of points
mode2 = zeros(p, 2); %p is num of points
mode3 = zeros(l, 2); %l is num of points

for i = 1:numOfDataPos
    if(speedModesBigArray(i, 1) == 1)
       mode1(j1, 2) = speedModesBigArray(i, 2); 
       mode1(j1, 1) = speedModesBigArray(i, 3); 
       j1 = j1+1;
    end 
    if(speedModesBigArray(i, 1) == 2)
       mode2(j2, 2) = speedModesBigArray(i, 2); 
       mode2(j2, 1) = speedModesBigArray(i, 3); 
       j2 = j2+1;
    end 
    if(speedModesBigArray(i, 1) == 3)
       mode3(j3, 2) = speedModesBigArray(i, 2); 
       mode3(j3, 1) = speedModesBigArray(i,3); 
       j3 = j3+1;
    end
end 

%plotting big nums
figure(13)
plot(timeInSecPos, speed, '-k')
hold on
plot(mode1(:, 2), mode1(:, 1), 'ro')
hold on 
plot(mode2(:, 2), mode2(:, 1), 'bo')
hold on 
plot(mode3(:, 2), mode3(:, 1), 'go')
hold off
title('Different speed modes')
legend('Total Speed', 'Stroll Through The Park', 'Power Walk', 'Late For Class')
xlabel('Time (s)')
ylabel('Speed (ms-1)')

%comment on not all values existing so avergae isnt fully accurate

% making array with begginings of new modes
%to find size of array
s = 0;
p = 0;
l = 0;

changes = 0;
for i = 2:numOfDataPos
    if (speedModesBigArray(i,1) ~= speedModesBigArray(i-1, 1))
        changes = changes + 1;
        if (speedModesBigArray(i,1) == 1)
            s = s+1;
        elseif (speedModesBigArray(i,1) == 2)
            p = p+1;
        elseif (speedModesBigArray(i,1) == 3)
            l = l+1;
        end 
    end 
end 

speedChangesAndTimes = zeros(changes, 3); %mode, time, i
j = 1;
for i = 2:numOfDataPos 
    if (speedModesBigArray(i,1) ~= speedModesBigArray(i-1, 1))
        speedChangesAndTimes(j, 1) = speedModesBigArray(i, 1);
        speedChangesAndTimes(j, 2) = timeInSecPos(i);
        speedChangesAndTimes(j, 3) = i;
        j = j+1;
    end 
end 


%manually adding first line
speedChangesAndTimes(1, 1) = speedModesBigArray(1, 1);
speedChangesAndTimes(1, 2) = timeInSecPos(1);
speedChangesAndTimes(1, 3) = 1;

%making scatter plot arrays
j1 = 1;
j2 = 1;
j3 = 1;

mode1changes = zeros(s, 3); %s is num of points
mode2changes = zeros(p, 3); %p is num of points
mode3changes = zeros(l, 3); %l is num of points %speed, time, i

for i = 1:changes
    if(speedChangesAndTimes(i, 1) == 1)
       mode1changes(j1, 2) = speedChangesAndTimes(i, 2); 
       mode1changes(j1, 1) = speed(speedChangesAndTimes(i, 3)); 
       mode1changes(j1, 3) = speedChangesAndTimes(i, 3); 
       j1 = j1+1;
    end 
    if(speedChangesAndTimes(i, 1) == 2)
       mode2changes(j2, 2) = speedChangesAndTimes(i, 2); 
       mode2changes(j2, 1) = speed(speedChangesAndTimes(i, 3)); 
       mode2changes(j2, 3) = speedChangesAndTimes(i, 3); 
       j2 = j2+1;
    end 
    if(speedChangesAndTimes(i, 1) == 3)
       mode3changes(j3, 2) = speedChangesAndTimes(i, 2); 
       mode3changes(j3, 1) = speed(speedChangesAndTimes(i,3)); 
       mode3changes(j3, 3) = speedChangesAndTimes(i, 3); 
       j3 = j3+1;
    end
end 


%manually adding last measurement
if (speedModesBigArray(numOfDataPos, 1) == 1)
    mode1changes(end, 1) = speedModesBigArray(end, 3);
    mode1changes(end, 2) = speedModesBigArray(end, 2);
    mode1changes(end, 3) = 1391;
elseif (speedModesBigArray(numOfDataPos, 1) == 2)
    mode2changes(end, 1) = speedModesBigArray(end, 3);
    mode2changes(end, 2) = speedModesBigArray(end, 2);
    mode2changes(end, 3) = 1391;
elseif (speedModesBigArray(numOfDataPos, 1) == 3)
    mode3changes(end, 1) = speedModesBigArray(end, 3);
    mode3changes(end, 2) = speedModesBigArray(end, 2);
    mode3changes(end, 3) = 1391;
end


figure(14)
plot(timeInSecPos, speed, '-k')
hold on
plot(mode1changes(:, 2), mode1changes(:,1), 'ro')
hold on
plot(mode2changes(:, 2), mode2changes(:,1), 'bo')
hold on 
plot(mode3changes(:, 2), mode3changes(:,1), 'go')
hold off
title('Speed with time-stamps for mode changes')
xlabel('Time (s)')
ylabel('Speed (ms-1)')
legend('Speed', 'Speed mode changes to mode 1', 'Speed mode changes to mode 2', 'Speed mode changes to mode 3')
%% task 6.3 - making a linked list of speed modes
clc;

% Declare a variable (initially empty) to store the head node
%head = LinkedListNode.empty;

% Declare a variable (initially empty) to store the previous node in the list,
% to be used to attach the next pointers
nodePrevious = LinkedListNode.empty;

%head = LinkedListNode.empty;

% Traverse the array where we have stored modes and timestamps
for i = 1: length(speedChangesAndTimes)
    i
    % Create a new SpeedMode for this movement mode
    mode = SpeedMode(speedChangesAndTimes(i, 1), speedChangesAndTimes(i, 2));
    % Create a new linked list node, inserting the speed mode
    linkedListNode = LinkedListNode(mode);
    if (i == 1)
	 % If it’s the first node in the list, assign it to the head variable
        head = linkedListNode;
    else
        % Otherwise make the previous list node point to this node 
        nodePrevious.Next = linkedListNode;
    % Update the previous node variable by setting it to the value of 
    % the current node 
    end
    nodePrevious = linkedListNode;
end 

currentNode = head;
disp(head.Data)
disp(head.Next)
while (~isempty(currentNode.Next))
    disp('yes')
    disp(currentNode.Data);
    currentNode = currentNode.Next;
end 


%% task 6.4 - plotting timestamps on x-y coordinates

%converting time stamps into x-y corrdinates
mode1xy = zeros(length(mode1changes), 2);
mode2xy = zeros(length(mode2changes), 2);
mode3xy = zeros(length(mode3changes), 2);

for i = 1:length(mode1changes)
    mode1xy(i, 1) = xEast(mode1changes(i, 3));
    mode1xy(i, 2) = yNorth(mode1changes(i, 3));
end 

for i = 1:length(mode2changes)
    mode2xy(i, 1) = xEast(mode2changes(i, 3));
    mode2xy(i, 2) = yNorth(mode2changes(i, 3));
end 

for i = 1:length(mode3changes) 
    mode3xy(i, 1) = xEast(mode3changes(i, 3));
    mode3xy(i, 2) = yNorth(mode3changes(i, 3));
end 


figure(14)
plot(xEast, yNorth, '-k')
hold on 
plot(mode1xy(:, 1), mode1xy(:, 2), 'ro', 'Linewidth', 3)
hold on 
plot(mode2xy(:, 1), mode2xy(:, 2), 'bo', 'Linewidth', 3)
hold on 
plot(mode3xy(:, 1), mode3xy(:, 2), 'go', 'Linewidth', 3)
title('Trajectory with time-stamps of different speed modes')
xlabel('x')
ylabel('y')
legend('x-y trajectory', 'Speed mode changes to mode 1', 'Speed mode changes to mode 2', 'Speed mode changes to mode 3')



%% task 7 - trajectory prediction
clc;
%the straight line section is from i = 1068 to i = 1163
%finding times for these measuremets
start = 1068;
finish = 1163;

timeStart = timeInSecPos(start);
timeFinish = timeInSecPos(finish);
timeMidPoint = round((timeStart + timeFinish)/2)

%checking if that specific time exists in timeInSecPos

for i = 1:numOfDataPos
    if(timeMidPoint == timeInSecPos(i))
        disp('yes')
        midpoint = i;
    end 
end 

fullLineX = xEast(start:finish, 1);
fullLlineY = yNorth(start:finish, 1);

figure(15)
plot(xEast, yNorth, '-r', 'Linewidth', 2)
hold on
plot(fullLineX, fullLlineY, 'c-', 'Linewidth', 2)
hold on
plot (xEast(midpoint), yNorth(midpoint), 'bo', 'Linewidth', 2)
hold on
plot (xEast(start), yNorth(start), 'mo', 'Linewidth', 2)
hold on
plot (xEast(finish), yNorth(finish), 'ko', 'Linewidth', 2)
hold off
legend('Full trajectory', 'Straight line segment', 'Midpoint', 'Start', 'End')
xlabel('x')
ylabel('y')
title('x-y trajectory with straight line segment')

%using regression in the first half to find besti fit line 
x = xEast(start:midpoint, 1);
y = yNorth(start:midpoint, 1);

%line in the form: y  = a + b*x
a0=[0,0];

fun = @(a) sum((y -a(1) -a(2)*x).^2)
a_optimal = fmincon(fun,a0);

a = a_optimal(1);
b = a_optimal(2);
f = @(x) a + b*x;

figure(16)
plot(xEast, yNorth, '-r', 'Linewidth', 2)
hold on
plot(x, f(x), '-k', 'Linewidth', 2)
hold off
legend('x-y trajectory', 'predicted first half of line')
xlabel('x')
ylabel('y')
title('x-y trajectory and predicted first half of line')

% actual position
xEndTrue = xEast(finish);
yEndTrue = yNorth(finish);

yEndPred = f(xEndTrue);

figure(17)
plot(xEast, yNorth, '-r', 'Linewidth', 2)
hold on 
plot(fullLineX, f(fullLineX), '-k', 'Linewidth', 2)
hold on
plot(xEndTrue, yEndTrue, 'bo', 'Linewidth', 3)
hold on
plot(xEndTrue, yEndPred, 'go', 'Linewidth', 3)
hold off
title('x-y trajectory and predicted end point')
legend('x-y trajectory', 'predicted line', 'true end point', 'predicted end point')
xlabel('x')
ylabel('y')

%finding error between values
error = abs(yEndPred - yEndTrue);
percentageError = error/yEndTrue*100


%% task 8 - data metrics

clc;

%finding max speed and time at which it occurs and average speed
sumSpeed = 0;
maxSpeed = speed(1);
maxTimeSum = 0; %measuring how many times max speed occurs

for i = 1:numOfDataPos
    sumSpeed = sumSpeed + speed(i);
    if (speed(i) > maxSpeed)
        maxTimeSum = 1;
        maxSpeed = speed(i);
    elseif (speed(i) == maxSpeed)
        maxTimeSum = maxTimeSum + 1;
        maxSpeed = speed(i);
    end 
end 

maxTimeSum; %occurs only once
averageSpeed = sumSpeed/numOfDataPos;

%keeping times when max speed occurs
maxSpeedTimes = zeros(maxTimeSum, 1);

%for the general case that it occurs more than once
j = 1;
for i = 1:numOfDataPos
    if (speed(i) == maxSpeed)
        maxSpeedTimes(j) = timeInSecPos(i)
        j = j+1;
    end
end 

%to calculate burned calories: Calories Burned = MET x Body Weight (kg) x Duration of Running (hours) 
%where MET depends on the speed in mph
%my speed in mph is 4.0mph
% Brisk walking (3.5–4 mph): MET = 5
MET = 5;
bodyWeight = 50; %in kg
duration = (timeInSecPos(end) - timeInSecPos(1))/60; %in minutes
%converting into hours
hours = duration/60; %in percentage


caloriesBurned = MET*bodyWeight*hours

%calories of a bigMac
BigMacCal = 257; 

%how many bigMacs burned:
burnedBigMacs = caloriesBurned/BigMacCal %not great...


%% task 3 - getting data from marshagte 
clc; clear;

load("marshgate.mat");
origin = [Position.latitude(1), Position.longitude(1), Position.altitude(1)];
[xEastM, yNorthM, zUpM] = latlon2local(Position.latitude, Position.longitude, Position.altitude, origin);
numOfDataPos = length(xEastM);

DT = datetime(Position.Timestamp);
% hours can be ignored cause everything happened in the same hour
min = minute(DT);
minInSec = min*60;
sec = second(DT);
timeInSec = minInSec + sec;

%calculating building's perimeter

total_distance = 0;

for i = 1: (numOfDataPos - 1) 
    distance = sqrt( (xEastM(i)- xEastM(i+1)).^2 + (yNorthM(i) - yNorthM(i+1)).^2 ); %include z
    total_distance = total_distance + distance;
end 

disp(total_distance)

%plotting data
figure(2)
plot(xEastM, yNorthM, '-', 'LineWidth', 2);
title('X-Y trajectory');
xlabel('x');
ylabel('y');


%%


%% notes 

%(connection losy etc)

