function [Tc,Vt,Isc,Voc,Im,Vm,I0,Rs] = calc(Ns,Np,Iscr,Vocr,Imr,Vmr,coef_Iscr,coef_Vocr,NOCT,Tr,Gr,Ta,G)

clc

%% Constants %%
q = 1.602176565e-19; %electron charge (C)
kb = 1.3806488e-23; %boltzman constant (J/K)

%% Calculations %%
Tc = Ta+((NOCT-20)*(G/800));
Vt = kb/q*(Tc + 273);
Isc = (G/Gr)*(Iscr+(coef_Iscr*(Tc-Tr)));
Voc = Vocr+(coef_Vocr*(Tc-Tr))+(Ns*Vt*log(Isc/Iscr));
Im = (G/Gr)*(Imr+(coef_Iscr*(Tc-Tr)));
I0 = Isc/(exp(Voc/(Ns*Vt))+exp(Voc/(2*Ns*Vt)));
FF0 = ((Voc/(Ns*Vt))-log((Voc/(Ns*Vt))+0.72))/(1+(Voc/(Ns*Vt)));
Vm_approx = Vmr+(coef_Vocr*(Tc-Tr))+(Ns*Vt*log(Isc/Iscr));
Rs_root = (1.21*(FF0^2))+((20/27)*(((Im*Vm_approx)/(Isc*Voc))-FF0));
Rs = (2.7*Voc*((1.1*FF0)-sqrt(Rs_root)))/Isc;
Vm = (Ns*Vt*log(((Isc-Im)/Isc)*(exp(Voc/(Ns*Vt))+exp(Voc/(2*Ns*Vt)))))-(Im*Rs);