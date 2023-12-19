function [Vanout,Vanout_new,tout_new,Iout_new,tout,Iout,f,P1,P1_V, THD_percentage, simOut]=partfourierexperiment_pwm_1 (stepsize,Vdc,V,m,frequency,R,L,fs,stoptime,samplefrequency, u_sampletime)

% define a function the left side is the output variables the rigt side is 
% the input variables. Vanout is the output phase volatge, Vanout_new is 
% the output voltage signal in a specific time interval namely 1/frequency, frequency
% is the frequency of the three pahse sinusodial signal, tout is simulation run 
% time,tout_new is time interval of the output phase volatge and pahse current signal from 
% 5*L*R to (5*L*R+(1/frequency)),it is saved as a matrix. f is a matrix which contanis all ticks of frequency achse of
% FFT spectrum analysis. Iout is the output phase current signal, 
% only phase a is analyzied since the only difference between phase a, b, c  
% is a phase delay of -2*pi/3. In other words a b and c is symmertric. P1 is a matrix which contains spectrum of the FFT analysis. 
% THD_percentage represents the total harmonic distortaion of the 
% output current signal. simOut is the Simulink file output, 
% in Simulink file SVPWM_slef.slx the variables are saved in model workspace for 
% further processing. partfourierexperiment_pwm_1 is the name of
% this MATLAB file, Furthermore is is also the name of the involked
% function. stepsize is the 1/samplefrequency which indicates how many
% distance between two neighbouring sampling point. V is the amplitude of
% the three phase sinusoidal signal. frequency represents the frequency of
% the three phase sinusoidal signal.R is the resistance of the electrcial
% machine, L is the inductance of the electrical machine. fs is the
% frequency of the carrier signal it can be sawtooth wave (Bipolar SPWM) or triangular
% wave (unipolar SPWM). Vdc is the DC-link volatge. stoptime is the simulation stoptime of
% the Simulink file. samplefrequency indicates how many sampling pouints are
% in one second time interval. u_sampletime is the sample time of the
% unitdelay block



simOut = sim("PWM_self_modified_PWM.slx");
% simOut is the Simulation output of the Simulink file  PWM_self_modified_PWM.slx
Iout=simOut.Iout(:,1);
Vanout=simOut.Vanout(:,1);
tout=simOut.tout;
% Iout is a matrix which contains amplitude of the three phase output current
% the sequence of the matrix elements is a structure with time
% Vanout is a matrix smimilat to Iout, which contains phase a output current
% tout is the simulation time
t_lr=L/R;

% t_lr is the transient circuit's time constant, since the machine load is 
% inductive, after 5*t_Lr the inductor charging stage will end and the 
% the FFT analysis after transient circuit’s time constant is more accurate

k3=int64((5*t_lr)/stepsize);
k4=int64((1/frequency)/stepsize);
tout_new=tout((k3):(k3+k4-1));
Iout_new=Iout((k3):(k3+k4-1));
Vanout_new=Vanout((k3):(k3+k4-1));
samplesize_new=(1/frequency)/stepsize;
samplesize_1=int64(samplesize_new);

% 5*t_lr represents after how lang time the circuit is in transient
% state,5*t_lr/stepsize is the number of sample pouint before transient state
% uint64 () function makes the data type an integer, since matrix index can
% only be integer. samplesize_new indicates how many sampling points in a
% period of modulating sinusiodal signal are

Y=fft(Iout_new);
Y_V=fft(Vanout_new);

P2 = abs(Y/samplesize_new);
P2_V=abs(Y_V/samplesize_new);
P1 = P2(1:samplesize_new/2+1);
P1_V = P2_V(1:samplesize_new/2+1);
P1(2:end-1) = 2*P1(2:end-1);
P1_V(2:end-1) = 2*P1_V(2:end-1);
Fs=samplefrequency;
f = Fs*(0:(samplesize_new/2))/samplesize_new;
k5=Fs*(samplesize_new/2)/samplesize_new;
k6=(numel(f)-1);
k7=k5/k6;

% FFT analysis orginal code from matlab official website is as follwing

% Fs = 1000;            % Sampling frequency                    
% T = 1/Fs;             % Sampling period       
% L = 1500;             % Length of signal
% t = (0:L-1)*T;        % Time vector
% S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);
% X = S + 2*randn(size(t));
% plot(1000*t(1:50),X(1:50))
% title("Signal Corrupted with Zero-Mean Random Noise")
% xlabel("t (milliseconds)")
% ylabel("X(t)")
% Y = fft(X);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = Fs*(0:(L/2))/L;
% plot(f,P1) 
% title("Single-Sided Amplitude Spectrum of X(t)")
% xlabel("f (Hz)")
% ylabel("|P1(f)|")

% line 59 to 67 use the same methods as line 86 to 90, since the FFT analysis only
% take a period of the output phase current signal or output phase voltage signal, Legth of the signal
% should be euqal to samplesize_new
% numel(f) indicates how many elements are in the matrix f, thus k6=(numel(f)-1)
% indicates how many intervals between those elements are
% k7=k5/k6 indicates the distance between two neighboring ticks in frequency achse 
% the mthod to calculate volatge spectrum P1_V is same as current spectrum
% P1
k8=(frequency/k7);
% k8 indicates the number of  intervals betwenen zero to frequency of 
% of first harmonic, this number of intervals is also the frequency of
% first harmonic to frequency of second harmonic, and so on. The aim to
% calcualte the number of intervals is, P1 is a mtrix whcich save the value
% of FFT spectrum as structure of time, when how many intervals betwenn nth
% harmonic to (n+1) th harmonic is calculated, the spectrum of the harmonic
% can be find with this index.
    k9=k6/k8;

P3=[];
P3(1,1)=P1(1,1);
for n2=1:k9
m1=uint64(1+n2*k8);
P3=[P3,P1(m1,1)];
end
% k9 indicates the number of higher harmonics, k6=(numel(f)-1) indicates the
% number of interval between n th harmonic and (n+1) yh harmonics in frequency achse, k6 divide k8 indicates
% how many k8 should it leap to go to the FFT analysis maximal frequency in
% x achse

% since FFT analysis only take the integer oder of harmonics like 1st 2nd 3rd
% the first element in matrix f is 0, which indicates the zero Hz in
% frequency achse.then after k8 interval is the frequency of 
% the first harmonic, after 2*k8 is the frequency of the second harmonic, and so on...
% after k9*k8 is the frquency of the k9 th harmonic
% P3 is an empty matrix, P3(1,1) is the first elment of P1 which is the 
% spectrum of 0 Hz,using a for loop with P3=[P3,P1(m1,1)]; to get the spectrum form 0Hz , first harmonic, second
% harmonic and so on... to k9 th harmonic and save all the spectrum data in matrix P3

P3_square=P3.^2;
k10=sum(P3_square(:));
k11=k10-P3_square(1,2);
k12=k10-k11;
THD_percentage=sqrt(k11/k12)*100;
end
% make the square of all the data in P3, k10 is the  sum of them. 
% k10 minus the square of spectrum of the first harmonic is k11, k11
% indicates the sum of the square of the 2 nd 3 rd to k9 th harmonic, k12
% is the spectrum of the first harmonic
% according to defination of THD_percentage is THD in percentage is equal to (sqrt(k11/k12))*100