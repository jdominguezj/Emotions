%prueba realizada con 2000 us de periodo de muestreo
%115200 puerto
clear,clc
datalb= dlmread('08_M_espontaneo_L_lb.log');
data=dlmread('08_M_espontaneo_L.log');

slb=datalb(:,1).*5/1024;%Linea base conductividad
plb=datalb(:,2).*5/1024;%Linea base pulso
s=data(:,1).*5/1024; %Conductividad
p=data(:,2).*5/1024; %Pulso
% s=data(25001:102500,1).*5/1024; %Conductividad
% p=data(25001:102500,2).*5/1024; %Pulso
Fs=500; %F. sampling
T=1/Fs; %T. sampling
%%Prueba eliminacion 500 muestras cada 6 segundos

%%
Tlb=length(slb);    %Numero de muestras de linea base
Tamp=length(p);     %Numero de muestras de tratamiento

%Vector de tiempo  tratamiento
t=(0:Tamp-1)*T; 
t1=t(1001:end); %Vector de tiempo para conductividad
t2=(0:Tamp-1)*T;%Vector de tiempo para la señal de Temperatura
t=t(51:end);    %Vector de tiempo para ritmo

%Vector de tiempo  linea base
Tbase=(0:Tlb-1)*T;
tglb=Tbase(1001:end);%Vector de tiempo para conductividad
tplb=Tbase(51:end); %Vector de tiempo para ritmo

%Definicion ventana rectangulares 
 
win=hamming(101);   %Tamaño de ventana para Ritmo: 101
win2=hamming(1001); %Tamaño de ventana para Conductividad: 1001

%Filtro pasabanda tipo FIR - Ritmo cardiaco
f1=0.5;
f2=10;

w1=2*f1/Fs;
w2=2*f2/Fs;
B=fir1(100,[w1 w2],'bandpass',win); %Filtro pasabanda tipo FIR
y=filter(B,1,p);
y=y(51:end);        %Eliminacion de primeras 50 muestras
pulsoaf=y;   

%Filtro pasabanda tipo FIR - Ritmo cardiaco linea base
ylb=filter(B,1,plb);
ylb=ylb(51:end);        %Eliminacion de primeras 50 muestras
pulsolb=ylb;   


%Filtrado pasabajos GSR
f3=1.5;
w3=2*f3/Fs;
B2=fir1(1000,w3,'low',win2);
y1=filter(B2,1,s);
y1=y1(1001:end);      %Eliminacion de primeras 1000 muestras
gsr = y1;

%Filtrado linea base GSR
y2=filter(B2,1,slb);
y2=y2(1001:end);      %Eliminacion de primeras 1000 muestras
gsrlb = y2;
gsrlb=gsrlb-mean(gsrlb);

%Segmentacion en ventanas de 5 segundos

 bx = 2500; %Division del total de los clips en ventanas de 5 segundos  
 na = numel(pulsoaf); %Numero de elementos de la señal de Ritmo cardiaco
 c = mat2cell(pulsoaf,diff([0:bx:na-1,na]));%Division de señal en ventanas de 5000 muestras
 c=c(1:end-1);
 n_iter = length(c);
 
 
%Segmentacion en ventanas de 5 segundos (Linea base)

 blb= 2500; %Division del total de los clips en ventanas de 5 segundos  
 nlb = numel(pulsolb); %Numero de elementos de la señal de Ritmo cardiaco
 clb = mat2cell(pulsolb,diff([0:blb:nlb-1,nlb]));%Division de señal en ventanas de 5000 muestras
 clb=clb(1:end-1);
 n_iterlb = length(clb);
 
%% 
%Definicion ventana hamming para FFT
winh=hamming(2500);

%Calculo de la FFT para la linea base del ritmo cardiaco
 for i=1:n_iterlb
 Xlb=fft(clb{i}.*winh,15000);
 bpmlb=abs(Xlb);
 bpmlb(1:30)=0;                 
 bpmlb(54:end)=0;
 [ylb,in]=max(bpmlb);
 hrlb(i)=(in-1)*60*Fs/(length(Xlb)); % Vector de ritmo cardiaco de linea base en BPM
 end
 
plm=mean(hrlb); %Media de linea base de ritmo cardiaco
pls=std(hrlb);  %Desv standar de linea base de ritmo cardiaco

%Calculo de la FFT para la señal de Ritmo
for i=1:n_iter
 X=fft(c{i}.*winh,15000);
 bpm=abs(X);
 bpm(1:30)=0;            
 bpm(54:end)=0;
[yvalue,k]=max(bpm);
 hr(i)=(k-1)*60*Fs/(length(X)); % Vector de ritmo cardiaco en BPM
 hrv(i)=(hr(i)-plm)/pls;         % Calculo de SMD(Standarized mean difference)
end 
%%
%Estimacion de beats/min en el tiempo
[amplitud,muestra,ancho,prominence]=findpeaks(pulsoaf,'MinPeakDistance',250,'MinPeakProminence',0.4,'MaxPeakWidth',400);
interbeat=60./(diff(muestra)/500);
%%
%Estimacion de conductividad promedio cada 5 segundos
nk=250;
bk=arrayfun(@(i) mean(gsr(i:i+nk-1)),1:nk:length(gsr)-nk+1);
% siem=1024+2*s+10000/512-s; %Valor en Komhs
%%
%Prueba 2 
scrb=downsample(gsr,25);
scr=scrb-mean(scrb);
comd=diff(scr);
co=conv(bartlett(20),comd);


scrlb=downsample(gsrlb,25);
scrlb=scrlb-mean(scrlb);
comlb=diff(scrlb);
colb=conv(bartlett(20),comlb);
%%    
figure(1)
subplot(2,2,1)
plot(t1,gsr,'r')

subplot(2,2,2)
plot(tglb,gsrlb,'b')
 
subplot(2,2,3);
plot(t,pulsoaf,'r')

subplot(2,2,4); 
plot(tplb,pulsolb,'b')

figure(2)

findpeaks(co,'MinPeakProminence',0.2*max(co))


figure(3)
plot(co)
%%x
disp('Extraccion de caractersticas')
disp('Conductividad')
fprintf('SCR B.L  Mean %8f .\n',mean(gsrlb));
fprintf('SCR B.L  Standard deviation %8f .\n',std(gsrlb));
fprintf('SCR Mean %8f .\n',mean(gsr));
fprintf('SCR Standard deviation %8f .\n',std(gsr));

fprintf('Heart Rate B.L Mean %8f .\n',mean(hrlb));
fprintf('Heart Rate B.L Standard deviation %8f .\n',std(hrlb));
fprintf('Heart Rate Mean %8f .\n',mean(hr));
fprintf('Heart Rate Standard deviation %8f .\n',std(hr));
fprintf('Heart Rate SMD %8f .\n',mean(hrv));