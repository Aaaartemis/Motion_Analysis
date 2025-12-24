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
figure(1)
plot(xEastM, yNorthM, '-', 'LineWidth', 2);
title('X-Y trajectory');
xlabel('x');
ylabel('y');
