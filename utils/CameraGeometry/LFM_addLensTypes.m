% oLaF - a flexible 3D reconstruction framework for light field microscopy
% Copyright (c)2017-2020 Anca Stefanoiu and Josue Page

function centersWithTypes = LFM_addLensTypes(lensCentersPx, matrixCenter)

%% Add lenslet types (different focal lengths) for multi-focus arrays 1-2-3
no_rows = size(lensCentersPx,1);
no_cols = size(lensCentersPx,2);

centersWithTypes = zeros(no_rows, no_cols, 3);
centersWithTypes(:,:,1:2) = lensCentersPx;

centersWithTypes(matrixCenter(1), matrixCenter(2), 3) = 1;

center_line_first = 0;
if centersWithTypes(matrixCenter(1),1,1) < centersWithTypes(matrixCenter(1) + 1,1,1)
    center_line_first = 1;
end

% even/odd rows
rows_11 = [fliplr(matrixCenter(1)-2 :-2: 1), matrixCenter(1):  2 : no_rows];
rows_12 = [fliplr(matrixCenter(1)-1 :-2: 1), matrixCenter(1) + 1:  2 : no_rows];

%% populating 1s
cols_11 = [fliplr(matrixCenter(2)-3 :-3: 1), matrixCenter(2):  3 : no_cols];
centersWithTypes(rows_11,cols_11,3) = 1;

if (center_line_first)
    cols_12 = [fliplr(matrixCenter(2)-2 :-3: 1), matrixCenter(2) + 1:  3 : no_cols];
    centersWithTypes(rows_12,cols_12,3) = 1;
else
    cols_12 = [fliplr(matrixCenter(2)-1 :-3: 1), matrixCenter(2) + 2:  3 : no_cols];
    centersWithTypes(rows_12,cols_12,3) = 1;
end

%% populating 2s
cols_21 = [fliplr(matrixCenter(2)-2 :-3: 1), matrixCenter(2) + 1:  3 : no_cols];
centersWithTypes(rows_11,cols_21,3) = 2;

if (center_line_first)
    cols_22 = [fliplr(matrixCenter(2)-1 :-3: 1), matrixCenter(2) + 2:  3 : no_cols];
    centersWithTypes(rows_12,cols_22,3) = 2;
else
    cols_22 = [fliplr(matrixCenter(2):-3: 1), matrixCenter(2):  3 : no_cols];
    centersWithTypes(rows_12,cols_22,3) = 2;
end
%% populating 3s
centersWithTypes(centersWithTypes == 0) = 3;

