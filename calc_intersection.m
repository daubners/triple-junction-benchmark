function [x0,y0] = calc_intersection(x1,y1,x2,y2)
    if length(x1)~=length(y1) || length(x2)~=length(y2)
        error('Input vectors do not have correct size!');
    end

    % Force all inputs to be column vectors.
    x1 = x1(:);
    y1 = y1(:);
    x2 = x2(:);
    y2 = y2(:);

    n1 = length(x1) - 1;
    n2 = length(x2) - 1;
    xy1 = [x1 y1];
    xy2 = [x2 y2];
    dxy1 = diff(xy1);
    dxy2 = diff(xy2);

    [i,j] = find( ...
		mvmin(x1) <= mvmax(x2).' & mvmax(x1) >= mvmin(x2).' & ...
		mvmin(y1) <= mvmax(y2).' & mvmax(y1) >= mvmin(y2).');

    % Initialize matrices.  We'll put the T's and B's in matrices and use them
    % one column at a time.  AA is a 3-D extension of A where we'll use one
    % plane at a time.
    n = length(i);
    T = zeros(4,n);
    AA = zeros(4,4,n);
    AA([1 2],3,:) = -1;
    AA([3 4],4,:) = -1;
    AA([1 3],1,:) = dxy1(i,:).';
    AA([2 4],2,:) = dxy2(j,:).';
    B = -[x1(i) x2(j) y1(i) y2(j)].';

    for k = 1:n
        [L,U] = lu(AA(:,:,k));
        T(:,k) = U\(L\B(:,k));
    end

    % Find where t1 and t2 are between 0 and 1 and return the corresponding
    % x0 and y0 values.
    in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) < 1 & T(2,:) < 1).';
    x0 = T(3,in_range).';
    y0 = T(4,in_range).';
end

function y = mvmin(x)
    % Faster implementation of movmin(x,k) when k = 1.
    y = min(x(1:end-1),x(2:end));
end

function y = mvmax(x)
    % Faster implementation of movmax(x,k) when k = 1.
    y = max(x(1:end-1),x(2:end));
end
