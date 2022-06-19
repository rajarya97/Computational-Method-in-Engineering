method = input('Please select the method to run:\n(1) Linear Spline\n(2) Quadratic Spline\n(3) Natural Cubic Spline\n(4) Not-a-Knot Cubic Spline\n(5) Periodic Cubic Spline\n(6) Clamped Cubic Spline\n');

x = [];
y = [];
X = [];
fileid=fopen('input.txt','r');
data = textscan(fileid, '%*s %*s %*s %*s\n');
data = textscan(fileid, '%f %f');
x(1:numel(data{1}), 1) = data{1};
y(1:numel(data{2}), 1) = data{2};
data = textscan(fileid, '%*s %*s %*s %*s %*s %*s %*s %*s\n');
data1 = textscan(fileid, '%f ');
X(1:numel(data1{1}), 1) = data1{1};
if(method == 6)
    s = zeros(2,1);
    data = textscan(fileid, '%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s\n');
    s(1,1) = fscanf(fileid,'%f', 1);
    s(2,1) = fscanf(fileid,'%f', 1);
end
fclose(fileid);

%points = table2array(x(1));
%eval = table2array(x(2));
points = [x,y];
eval = X;

if (method == 1)
    linearspline(points, eval);
elseif (method == 2)
    quadraticspline(points, eval);
else
    file = fopen('output.txt', 'w');
    
    n = size(points(:, 1));
    n = n(1) - 1;
    mat = cubicspline(points);
    title = '';
    switch (method)
        case 3
            title = 'Natural';
            mat(4*n-1, 3) = 2;
            mat(4*n-1, 4) = 6*points(1, 1);
            mat(4*n, 4*n-1) = 2;
            mat(4*n, 4*n) = 6*points(n+1, 1);
            fprintf(file, 'Natural Spline\n\n');
        case 4
            title = 'Not-a-Knot';
            mat(4*n-1, 4) = 1;
            mat(4*n-1, 8) = -1;
            mat(4*n, 4*n-4) = 1;
            mat(4*n, 4*n) = -1;
            fprintf(file, 'Not-a-Knot Spline\n\n');
        case 5
            title = 'Periodic';
            mat(4*n-1, 2) = 1;
            mat(4*n-1, 3) = 2*points(1, 1);
            mat(4*n-1, 4) = 3*points(1, 1)^2;
            mat(4*n-1, 4*n-2) = -1;
            mat(4*n-1, 4*n-1) = -2*points(n+1, 1);
            mat(4*n-1, 4*n) = -3*points(n+1, 1)^2;
            mat(4*n, 3) = 2;
            mat(4*n, 4) = 6*points(1, 1);
            mat(4*n, 4*n-1) = -2;
            mat(4*n, 4*n) = -6*points(n+1, 1);
            fprintf(file, 'Periodic Spline\n\n');
        case 6
            slopes = s;
            mat(4*n-1, 2) = 1;
            mat(4*n-1, 3) = 2*points(1, 1);
            mat(4*n-1, 4) = 3*points(1, 1)^2;
            mat(4*n-1, 4*n+1) = slopes(1);
            mat(4*n, 4*n-2) = 1;
            mat(4*n, 4*n-1) = 2*points(n+1, 1);
            mat(4*n, 4*n) = 3*points(n+1, 1)^2;
            mat(4*n, 4*n+1) = slopes(2);
            fprintf(file, 'Clamped Spline\n\n');
            
    end
    coeff = inv(mat(:, 1:4*n))*mat(:, 4*n+1);
    coeff = flipud(reshape(coeff, [4, n]));
    for i=1:size(eval)
        for j=1:n
            if (eval(i) >= points(j, 1) && eval(i) <= points(j+1, 1))
                fprintf(file, '%.4f\t%.4f\n', eval(i), polyval(coeff(:, j), eval(i)));
            end
        end
    end
    
    y = [];
    p1 = [];
    for i=1:n
        p = [points(i, 1):0.01:points(i+1,1)];
        p1 = [p1 p];
        y = [y polyval(coeff(:, i), p)];
    end
    plot(p1, y);
end

function linearspline(points, x)

n = size(points(:, 1));
n = n(1) - 1;

mat = zeros(2*n, 2*n);
fun = zeros(2*n, 1);

cnt = 1;
cnt2 = 0;
for i=1:n
    mat(cnt, cnt2*2+1) = 1;
    mat(cnt, cnt2*2+2) = points(i, 1);
    fun(cnt) = points(i, 2);
    cnt = cnt + 1;
    mat(cnt, cnt2*2+1) = 1;
    mat(cnt, cnt2*2+2) = points(i+1, 1);
    fun(cnt) = points(i+1, 2);
    cnt = cnt + 1;
    cnt2 = cnt2 + 1;
end

coeff = inv(mat)*fun;
coeff = flipud(reshape(coeff, [2, n]));

file = fopen('output.txt', 'w');
fprintf(file, 'Linear Spline\n');
for i=1:size(x)
    for j=1:n
        if (x(i) >= points(j, 1) && x(i) <= points(j+1, 1))
            fprintf(file, '%.4f\t%.4f\n', x(i), polyval(coeff(:, j), x(i)));
        end
    end
end
    
y = [];
p1 = [];
for i=1:n
    p = [points(i, 1):0.01:points(i+1,1)];
    p1 = [p1 p];
    y = [y polyval(coeff(:, i), p)];
end
plot(p1, y);

end

function quadraticspline(points, x)

n = size(points(:, 1));
n = n(1) - 1;

mat = zeros(3*n, 3*n);
fun = zeros(3*n, 1);

cnt = 1;
cnt2 = 0;
for i=1:n
    mat(cnt, cnt2*3+1) = 1;
    mat(cnt, cnt2*3+2) = points(i, 1);
    mat(cnt, cnt2*3+3) = points(i, 1)^2;
    fun(cnt) = points(i, 2);
    cnt = cnt + 1;
    mat(cnt, cnt2*3+1) = 1;
    mat(cnt, cnt2*3+2) = points(i+1, 1);
    mat(cnt, cnt2*3+3) = points(i+1, 1)^2;
    fun(cnt) = points(i+1, 2);
    cnt = cnt + 1;
    cnt2 = cnt2 + 1;
end

for i=1:n-1
    mat(cnt, (i-1)*3 + 2) = 1;
    mat(cnt, i*3 + 2) = -1;    
    mat(cnt, (i-1)*3 + 3) = 2*points(i+1);
    mat(cnt, i*3 + 3) = -2*points(i+1);
    cnt = cnt + 1;
end
mat(cnt, 3) = 1;

coeff = inv(mat)*fun;
coeff = flipud(reshape(coeff, [3, n]));

file = fopen('output.txt', 'w');
fprintf(file, 'Quadratic Spline\n');
for i=1:size(x)
    for j=1:n
        if (x(i) >= points(j, 1) && x(i) <= points(j+1, 1))
            fprintf(file, '%.4f\t%.4f\n', x(i), polyval(coeff(:, j), x(i)));
        end
    end
end
    
y = [];
p1 = [];
for i=1:n
    p = [points(i, 1):0.01:points(i+1,1)];
    p1 = [p1 p];
    y = [y polyval(coeff(:, i), p)];
end
plot(p1, y);

end

function mat = cubicspline(points)

n = size(points(:, 1));
n = n(1) - 1;

mat = zeros(4*n, 4*n+1);
fun = zeros(4*n, 1);

cnt = 1;
cnt2 = 0;
for i=1:n
    mat(cnt, cnt2*4+1) = 1;
    mat(cnt, cnt2*4+2) = points(i, 1);
    mat(cnt, cnt2*4+3) = points(i, 1)^2;
    mat(cnt, cnt2*4+4) = points(i, 1)^3;
    mat(cnt, 4*n+1) = points(i, 2);
    cnt = cnt + 1;
    mat(cnt, cnt2*4+1) = 1;
    mat(cnt, cnt2*4+2) = points(i+1, 1);
    mat(cnt, cnt2*4+3) = points(i+1, 1)^2;
    mat(cnt, cnt2*4+4) = points(i+1, 1)^3;
    mat(cnt, 4*n+1) = points(i+1, 2);
    cnt = cnt + 1;
    cnt2 = cnt2 + 1;
end

for i=1:n-1
    mat(cnt, (i-1)*4 + 2) = 1;
    mat(cnt, i*4 + 2) = -1;
    mat(cnt, (i-1)*4 + 3) = 2*points(i+1);
    mat(cnt, i*4 + 3) = -2*points(i+1);
    mat(cnt, (i-1)*4 + 4) = 3*points(i+1)^2;
    mat(cnt, i*4 + 4) = -3*points(i+1)^2;
    cnt = cnt + 1;
end
for i=1:n-1
    mat(cnt, (i-1)*4 + 3) = 2;
    mat(cnt, i*4 + 3) = -2;
    mat(cnt, (i-1)*4 + 4) = 6*points(i+1);
    mat(cnt, i*4 + 4) = -6*points(i+1);
    cnt = cnt + 1;
end
end