function out = hubenySub(x1_deg, y1_deg, x2_deg, y2_deg)
% hubenySub : Lat,Lon -> [m]
% subroutine hubeny_sub
% out = hubenySub(x1_deg, y1_deg, x2_deg, y2_deg)
%
%
% [ref]


x1 = double(x1_deg) * pi / 180;
y1 = double(y1_deg) * pi / 180;
x2 = double(x2_deg) * pi / 180;
y2 = double(y2_deg) * pi / 180;

dy = y1 - y2;
dx = x1 - x2;
mu = (y1 + y2) / 2;

a = 6378137.000;  % Semi-Major Axis
b = 6356752.314;  % Semi-Minor Axis

e = sqrt((a^2 - b^2) / (a^2));

W = sqrt(1 - e^2 * (sin(mu))^2);

N = a / W;

M = a * (1 - e^2) / W^3;

out = sqrt( double((dy * M)^2 + (dx * N * cos(mu))^2) );

end