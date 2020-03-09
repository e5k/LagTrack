function [ rk ] = rk( x,z,vx,vz,w,dt,rhoa,cdnow,area,mass,g,rhor )
 
    vwind = sqrt((vx - w) ^ 2 + vz * vz);
    vvxwd = vx - w;
    vxwind = vx - w;
    vvz = vz;
    
    derivs1=derivs(vwind, vvxwd, vvz,w,dt,rhoa,cdnow,area,mass,g,rhor) ;       
    
    vwind = sqrt((vx + 0.5 * derivs1(3)) ^ 2 + (vz + 0.5 * derivs1(4)) ^ 2);
    vvxwd = vx + 0.5 * derivs1(3);
    vvz = vz + 0.5 * derivs1(4);
    
    derivs2=derivs(vwind, vvxwd, vvz,w,dt,rhoa,cdnow,area,mass,g,rhor);
    
    vwind = sqrt((vx + 0.5 * derivs2(3)) ^ 2 + (vz + 0.5 * derivs2(4)) ^ 2);
    vvxwd = vx + 0.5 * derivs2(3);
    vvz = vz + 0.5 * derivs2(4);
    
    derivs3=derivs(vwind, vvxwd, vvz,w,dt,rhoa,cdnow,area,mass,g,rhor);
    
    vwind = sqrt((vx + 0.5 * derivs3(3)) ^ 2 + (vz + 0.5 * derivs3(4)) ^ 2);
    vvxwd = vx + 0.5 * derivs3(3);
    vvz = vz + 0.5 * derivs3(4);
    
    derivs4=derivs(vwind, vvxwd, vvz,w,dt,rhoa,cdnow,area,mass,g,rhor);
    
    r=[derivs1' derivs2' derivs3' derivs4'];
    
    d=zeros(4);
    
    for i=1:4
      d(i) = r(i, 1) / 6 + r(i, 2) / 3 + r(i, 3) / 3 + r(i, 4) / 6;
    end

    rk(1) = x + d(1); %new x
    rk(2) = z + d(2); %new z
    rk(3) = vxwind + d(3) + w; %new vx
    rk(4) = vz + d(4); %new vz

end

