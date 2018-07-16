function f = funca(v,u,b,I,dT)
%f=0.04*v.^2-v.*3+I+40+b.*100-u;
%f=v+dT*(0.04*v.^2-v.*3+40+100*b-u+I);
f=(0.04*v.^2-v.*3+40+100*b-u+I);
end
