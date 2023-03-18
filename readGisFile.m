function gisData = readGisFile(fn,NX,NY,XLLCORNER,YLLCORNER,CELLSIZE)
% readGisFile : ASCIIデータ読み取り
% subroutine read_gis_real
% gisData = readGisFile(fn,NX,NY,XLLCORNER,YLLCORNER,CELLSIZE)
% outputの行列はまだ転置されていない（サブルーチン外で転置される）
% NX：列数，NY：行数
%
% [ref]

gisData = readmatrix(fn, 'NumHeaderLines', 6, 'Delimiter', " ");

% data validation
header = readvars(fn, 'Range', 'A1:C5');
header = cellfun(@split, header, 'UniformOutput', false);
xllCornerInFile = str2double(cell2mat(header{3}(2)));
yllCornerInFile = str2double(cell2mat(header{4}(2)));
cellSizeInFile  = str2double(cell2mat(header{5}(2)));
if size(gisData, 2) ~= NX, error("error in gis input data"), end
if size(gisData, 1) ~= NY, error("error in gis input data"), end
if abs(xllCornerInFile - XLLCORNER) > 0.01
    error("error in gis input data");
end
if abs(yllCornerInFile - YLLCORNER) > 0.01
    error("error in gis input data");
end
if abs(cellSizeInFile - CELLSIZE) > 0.01
    error("error in gis input data (real)");
end

end