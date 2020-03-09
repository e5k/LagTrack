function [ y2 ] = spline( x, y, n, yp1, ypn )



if yp1 > 9.9e29
        y2(1) = 0;
        u(1) = 0;
else
        y2(1) = -0.5;
        u(1) = (3 / (x(2) - x(1))) * ((y(2) - y(1)) / (x(2) - x(1)) - yp1);
end

for i = 2:(n - 1)
        sig = (x(i) - x(i - 1)) / (x(i + 1) - x(i - 1));
        p = sig * y2(i - 1) + 2;
        y2(i) = (sig - 1) / p;
        u(i) = (6 * ((y(i + 1) - y(i)) / (x(i + 1) - x(i)) - (y(i) - y(i - 1)) / (x(i) - x(i - 1))) / (x(i + 1) - x(i - 1)) - sig * u(i - 1)) / p;
end

if ypn > 9.9e29
        qn = 0;
        un = 0;
else
        qn = 0.5;
        un = (3 / (x(n) - x(n - 1))) * (ypn - (y(n) - y(n - 1)) / (x(n) - x(n - 1)));
end

      y2(n) = (un - qn * u(n - 1)) / (qn * y2(n - 1) + 1);
      
for k = (n - 1) : -1 : 1 
        y2(k) = y2(k) * y2(k + 1) + u(k);
end


end

