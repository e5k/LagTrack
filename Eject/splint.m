function [ y ] = splint( xa, ya, y2a, n, x )

klo = 1;
khi = n;

while (khi - klo) > 1
        k = round((khi + klo) / 2); 
        if xa(k) > x 
          khi = k;
        else
          klo = k;
        end
end

h = xa(khi) - xa(klo);

if h == 0
        disp('bad xa input in splint')
        return
end

      a = (xa(khi) - x) / h;
      b = (x - xa(klo)) / h;
      y = a * ya(klo) + b * ya(khi) + ((a ^ 3 - a) * y2a(klo) + (b ^ 3 - b) * y2a(khi)) * (h ^ 2) / 6;

end

