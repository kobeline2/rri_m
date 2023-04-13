function [dx, dy] = metricCellsize(XLLCORNER, YLLCORNER, NX, NY, CELLSIZE, utm)
% metricCellsize

% d1:south side length (x1,y1):右上, (x2,y2):右下
x1 = XLLCORNER;
y1 = YLLCORNER;
x2 = XLLCORNER + NX * CELLSIZE;
y2 = YLLCORNER;
if( utm == 0 ); d1 = hubenySub( x1, y1, x2, y2 ); end % 緯度経度→メートル

% d2:north side length (x1,y1):左上, (x2,y2):左下
x1 = XLLCORNER;
y1 = YLLCORNER + NY * CELLSIZE;
x2 = XLLCORNER + NX * CELLSIZE;
y2 = YLLCORNER + NY * CELLSIZE;
if( utm == 0 ); d2 = hubenySub( x1, y1, x2, y2 ); end % 緯度経度→メートル 

% d3:west side length (x1,y1):右上, (x2,y2):左上
x1 = XLLCORNER;
y1 = YLLCORNER;
x2 = XLLCORNER;
y2 = YLLCORNER + NY * CELLSIZE;
if( utm == 0 ); d3 = hubenySub( x1, y1, x2, y2 ); end % 緯度経度→メートル 

% d1:east side length (x1,y1):右下, (x2,y2):左下
x1 = XLLCORNER + NX * CELLSIZE;
y1 = YLLCORNER;
x2 = XLLCORNER + NX * CELLSIZE;
y2 = YLLCORNER + NY * CELLSIZE;
if( utm == 0 ); d4 = hubenySub( x1, y1, x2, y2 ); end % 緯度経度→メートル 

if utm == 1
    dx = CELLSIZE;
    dy = CELLSIZE;
else
    dx = double( (d1 + d2) / 2 / NX);
    dy = double( (d3 + d4) / 2 / NY);
end
disp(['dx [m] : ', num2str(dx),'    dy [m] : ', num2str(dy)] )
end