function x0 = newton_raphson_method_v2(init, V, Isc, I0, Rs, Ns, Vt)

clc

x0 = init;
x01 = init;
% outfile = fopen('outfile.txt','w');

while abs(f(x01, V, Isc, I0, Rs, Ns, Vt)) > 1e-10

    x1 = x0 - f(x0, V, Isc, I0, Rs, Ns, Vt)/fprime(x0, V, I0, Rs, Ns, Vt);
%     fprintf(outfile,'%13.11f  %13.11f %13.11f  %13.11f\n',x0,V,abs(f(x0, V, Isc, I0, Rs, Ns, Vt)),abs(fprime(x0, V, I0, Rs, Ns, Vt)));
    x01 = x0;
    x0 = x1;
    
end

% fclose(outfile);

function out = f(in, V, Isc, I0, Rs, Ns, Vt)
out = Isc-(I0*(exp((V+in*Rs)/(Ns*Vt))+exp((V+in*Rs)/(2*Ns*Vt))-2))-in;

function out = fprime(in, V, I0, Rs, Ns, Vt)
out = -I0*((exp((V+in*Rs)/(Ns*Vt))*(Rs/(Ns*Vt)))+(exp((V+in*Rs)/(2*Ns*Vt))*(Rs/(2*Ns*Vt))))-1;