function [ derivs ] = derivs( vwind,vvxwd,vvz,w,dt,rhoa,cdnow,area,mass,g,rhor )
 
    derivs(1) = (vvxwd + w) * dt;
    derivs(2) = vvz * dt;
    derivs(3) = -rhoa * cdnow * area * vwind ^ 2 * vvxwd / (2 * mass * vwind) * dt;
    derivs(4) = (-g * (rhor - rhoa) / rhor - rhoa * cdnow * area * vwind ^ 2 * vvz /(2 * mass * vwind)) * dt;

end

