%% task 3 - getting data from bridge
clc; clear;

load("bridge.mat");
origin = [Position.latitude(1), Position.longitude(1), Position.altitude(1)];
[xEastB, yNorthB, zUpB] = latlon2local(Position.latitude, Position.longitude, Position.altitude, origin);
numOfDataPos = length(xEastB);

DT = datetime(Position.Timestamp);
% hours can be ignored cause everything happened in the same hour
min = minute(DT);
minInSec = min*60;
sec = second(DT);
timeInSec = minInSec + sec;

figure(1)
plot(xEastB, yNorthB, '-', 'LineWidth', 2);
title('X-Y trajectory');
xlabel('x');
ylabel('y');

%%

%removign first 4 points bc they are abv wrong
x = xEastB(4:end);
y = yNorthB(4:end);
plot(x, y, '-')


total_distance = 0;

for i = 1: (numOfDataPos - 4) 
    distance = sqrt( (x(i)- x(i+1)).^2 + (y(i) - y(i+1)).^2 ); %include z
    total_distance = total_distance + distance;
end 

disp(total_distance);

%plotting data
figure(2)
plot(x, y, '-', 'LineWidth', 2);
title('X-Y trajectory');
xlabel('x');
ylabel('y');

% difference is big and that is excluding the first 4 points which were
% completely off and considering its supposed to be a fully straight line 
