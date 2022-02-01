
% simplex algorithm
clc
matrix=[-1,1,-1;0,-1,-2;-1,1,1;1,0,4;-1,-2,0]; %matrix for phase 1

matrix2 = [0 -1 -2;-1 1 -1;-1 1 1;1 0 4;-3 2 2];

matrix3 = [6 5 60;1 2 15;1 0 8;-10 -9 0]; % phase2

matrix4 = [-2 -3 -1;-1 1 0;-3 -6 -2;2 2 1;2 1 2]; % anopther
matrix5 = [0 1 4;-1 1 -1;2 -1 8;-3 -1 0];



updateTableau(matrix5)

function updateTableau(matrix)
%     check the n+1 column for negative values. IF negative, use phase 1
    b = matrix(1:end-1, end);
    if any(b<0)
        [i0,j0,pivotVal] = phaseOne(matrix);
    else
        [i0,j0,pivotVal] = phaseTwo(matrix);
    end
numOfRecursions = 0;
if numOfRecursions ==0
    disp(matrix)
end
fprintf('Pivot value is %f at (%d,%d)\n', (pivotVal), i0,j0)
% start updating the table
 updatedTableu = zeros(size(matrix));
 [i_,j_] = size(matrix);

 for i =1:i_
    for j =1:j_
        if (i == i0) || (j == j0)
                
        else
            updatedTableu(i,j)=matrix(i,j) - ((matrix(i0, j)*matrix(i,j0))/matrix(i0,j0));
        end
    end
 end
updatedTableu(i0,j0) = 1/matrix(i0,j0);
updatedTableu(i0,:) = matrix(i0,:)/matrix(i0,j0);
updatedTableu(:,j0) = matrix(:,j0)/(-1*matrix(i0,j0));
updatedTableu(i0,j0) = 1/matrix(i0,j0); 

format rational
disp(updatedTableu)
% fprintf('Pivot value is %f at (%d,%d)\n', (pivotVal), i0,j0)

disp('**********************')
b_new = updatedTableu(1:end-1, end);
c_new = updatedTableu(end,:);
% perform recursion

if any(b_new<0) || any(c_new<0)
    updateTableau(updatedTableu)
%     numOfRecursions = numOfRecursions + 1;
else
    disp('Optimal value obtained')
end
end


function [i0,j0,pivotVal]=phaseTwo(matrix)
    disp('Perfoming Phase Two')
    %extract the last row, c
    c = matrix(end,:);
    if any(c<0)
       [~,j0] = min(c);
    end
    ratio = matrix(1:end-1,end)./matrix(1:end-1, j0);
    [~, i0] = minpositive(ratio);
    pivotVal = matrix(i0,j0);
end

function [i0,j0, pivotVal] = phaseOne(matrix)
    % extract b
    disp('Perfoming Phase One')
    b = matrix(1:end-1, end);
    i_caret = find(b <0, 1, 'last');
    d = i_caret-1;

    fprintf('Location of i_hat is %d \n', i_caret)

    tmp = matrix(i_caret,1:end-1);
    j0 = min(find(tmp<0,1));
    
    div = matrix(i_caret:end-1,end)./matrix(i_caret:end-1,j0);
    [~, i0 ] = minpositive(div);

    i0 = d + i0;
    
    pivotVal = matrix(i0,j0);  
end


function [mins, idxes] = minpositive(Array, varargin)
  Array(Array<=0) = nan;
  [mins, idxes] = min(Array, varargin{:});
end