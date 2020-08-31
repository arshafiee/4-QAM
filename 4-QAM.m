%Authors : Hanie Azad & Armin Shafiee
%Student#:  95101033  &   95101785
%emails  : hanieazad1377@yahoo.com & arminshafiee7712@yahoo.com

%-------------------------------Project---------------------------
%-----------------------------------------------------------------

%---------------------------Transmitter---------------------------
%-----------------------------------------------------------------

clear;clc;
%Student Number = 95101785 -------> 4-QAM
num_of_symbols = 1e3;
bits = rand(1,2*num_of_symbols)>=0.5;%generating random bits
bits_reshape=reshape(bits,2,num_of_symbols)';

%--------------------------Generating symbols---------------------

for(j=1:1:num_of_symbols)
   for(i=1:1:2)
       symbols(j,i)=num2str(bits_reshape(j,i));
   end
end
symbols=bin2dec(symbols);
symbols=symbols';

%-------------------------QAM Modulation---------------------------

p = qammod(symbols,4,'gray'); %constalation diagram for 4-QAM
scatterplot(p);
title('constellation of the message signal')
symbol_rate = 4e7;
symbol_period = 1/symbol_rate;
bit_rate = 2*symbol_rate;
bit_period = 1/(2*symbol_rate);
Real = real(p);
Imag = imag(p);
f = bit_rate;
t = symbol_period/100:symbol_period/100:symbol_period;
trans_signal = [];
for(k = 1:1:length(Real))
    yr = Real(k)*cos(2*pi*f*t);% inphase or real component
    yim = Imag(k)*sin(2*pi*f*t);% Quadrature or imagenary component 
    trans_signal(k,:) = yr+yim;
end

%------------------------------Plots-------------------------------

bit = [];
for n=1:1:length(bits)
    if bits(n)==1;
       se=ones(1,100);
    else bits(n)==0;
       se=zeros(1,100);
    end
     bit=[bit se];
end
t1 = 0:bit_period/100:100*10*(bit_period/100)-bit_period/100;
t2 = 0:symbol_period/100:100*5*(symbol_period/100)-symbol_period/100;
figure();
subplot(3,1,1);
plot(t1,bit(1:1000),'lineWidth',2.5);grid on;
axis([ 0 bit_period*10 -.5 1.5]);
ylabel('amplitude(volt)');
xlabel(' time(sec)');
title('transmitting information as digital signal(first 10 bits)');

subplot(3,1,2)
stem(symbols(1:5),'Linewidth',2.0);
title('serial symbol for 4-QAM modulation at transmitter(first 5 symbols)');
xlabel('n(discrete time)');
ylabel(' magnitude');

subplot(3,1,3)
plot(t2,trans_signal(1:500));
axis([ 0 bit_period*10 -1.5 1.5]);
title('waveform for 4-QAM modulation according to symbolic information');
xlabel('time(sec)');
ylabel('amplitude(volt)');
ylim([-2 2]);

trans_signal = trans_signal';
trans_signal = trans_signal(:);
trans_signal = trans_signal';
trans_signal_fft = fftshift(fft(trans_signal));

Fs = symbol_rate*100;
df=Fs/length(trans_signal_fft);
freqvec=-Fs/2+df:df:Fs/2;
figure();
plot(freqvec,abs(trans_signal_fft)/length(trans_signal_fft));
title('Amplitude Spectrum of the Modulated Signal');
xlabel('frequency');
ylabel('Amplitude');

%%
%--------------------------Up-sample & Tx-filter------------------------

Up_sample_factor = 4;
alpha = 0.4; %rolloff factor
Nd = 10; %duration
h_Tx = rcosdesign(alpha,Nd,Up_sample_factor);
h_Tx_fft = fftshift(fft(h_Tx));
Fs = symbol_rate*100;
df = Fs/length(h_Tx);
freqvec = -Fs/2+df:df:Fs/2;
figure();
plot(freqvec,abs(h_Tx_fft/length(h_Tx)));
title('Square root raised-cosine filter(frequency response)');
xlabel('Frequency');
ylabel('Amplitude');
figure();
plot(h_Tx);
title('Square root raised-cosine filter(time domain)');

trans_signal_final = upfirdn(trans_signal,h_Tx,Up_sample_factor);

%-------------------------Truncating------------------------------------

trans_signal_final = trans_signal_final(19:400018);

trans_signal_final_fft = fftshift(fft(trans_signal_final));

eyediagram(p,8,64);%length of blocks = 8

%------------------------Ploting----------------------------------------

Fs = symbol_rate*100;
df = Fs/length(trans_signal_final_fft);
freqvec = -Fs/2+df:df:Fs/2;
figure();
plot(freqvec,abs(trans_signal_final_fft)/length(trans_signal_fft));
title('Amplitude Spectrum of the Transmitted Signal');
xlabel('frequency');
ylabel('Amplitude');

fvtool(h_Tx, 'Analysis', 'impulse')


%%
%--------------------------Transmission over channel---------------------
%------------------------------------------------------------------------

coeff = [-0.4945 - 0.4611i, -0.2813 + 0.8457i, -0.0720 - 0.5260i,...
    0.1201 - 0.4302i, 0.2451 + 0.1999i];
recieved_signal = filter(coeff,1,trans_signal_final);
recieved_signal_fft = fftshift(fft(recieved_signal));
scatterplot(recieved_signal);
title('constellation of Recieved signal(before adding noise)')

%--------------------------Adding noise----------------------------------

SNR = 15;
Noisy_signal = awgn(recieved_signal,SNR,'measured','linear');
Noisy_signal_fft = fftshift(fft(Noisy_signal));
scatterplot(Noisy_signal);
title('constellation of Noisy signal')

%--------------------------Frequency response of the channel-------------

h_channel = impz(coeff , 1 , length(Noisy_signal));
channel_fft = fftshift(fft(h_channel));

%--------------------------Ploting Amplitude of Freq res-----------------

Fs = symbol_rate*100;
df=Fs/length(channel_fft);
freqvec=-Fs/2+df:df:Fs/2;
figure();
plot(freqvec,abs(channel_fft)/length(channel_fft));
title('Amplitude Response of the Channel');
xlabel('frequency');
ylabel('Amplitude');

%--------------------------Ploting Phase of Freq res---------------------

figure();
plot(freqvec,phase(channel_fft));
title('Phase Response of the Channel');
xlabel('frequency');
ylabel('Phase');

%--------------------------Ploting Amp Spec of recieved signal-----------

figure();
plot(freqvec,abs(recieved_signal_fft)/length(recieved_signal_fft));
title('Amplitude Spectrum of the channel output signal(before adding noise)');
xlabel('frequency');
ylabel('Amplitude');

%--------------------------Ploting Amp Spec of Noisy signal-----------

figure();
plot(freqvec,abs(Noisy_signal_fft)/length(Noisy_signal_fft));
title('Amplitude Spectrum of the Noisy signal');
xlabel('frequency');
ylabel('Amplitude');

%%
%---------------------------Reciever-------------------------------------
%------------------------------------------------------------------------

%---------------------------Rx-filter & down-sample----------------------

Y = upfirdn(Noisy_signal,h_Tx,1,Up_sample_factor);
scatterplot(Y);
title('constellation of the down sampled signal(Unequaled Signal)')
eyediagram(Noisy_signal(1:1000),8,64);

%%
%---------------------------LMS-algorithm--------------------------------

N = 1000; % Number of samples              
data = randi([0 1],1,N);% Random signal
d = qammod(data,4);    % 4_QAM modulation
r = filter(coeff,1,data);       % Signal after passing through channel
x = awgn(r, SNR);           % Noisy Signal after channel

%---------------------------LMS-parameters-------------------------------

iterations = 10000; % Number of iterations
eta = 0.0001;  % Step size
order = 30;    % Order of the equalizer
U = zeros(1,order); % Input frame
W = zeros(1,order); % Initial Weigths

%---------------------------Algorithm------------------------------------

for k = 1 : iterations
    for n = 1 : N
        U(1,2:end) = U(1,1:end-1);  % Sliding window
        U(1,1) = x(n);   % Present Input
     
        y = (W)*U'; % Calculating output of LMS
        e = d(n) - y;  % Instantaneous error 
        W = W +  eta * e * U ; % Weight update rule of LMS
        J(k,n) = e * e'; % Instantaneous square error
    end
end

%------------------Calculation of performance parameters-----------------

MJ = mean(J,2); % Mean square error
CS = freqz(coeff); % Channel Spectrum
NF = (0:length(CS)-1)./(length(CS)); % Normalized Frequencies
IMR = -10*log10(real(CS).^2 + imag(CS).^2); % Inverse channel magnitude response 
IPR = -imag(CS)./real(CS); % Inverse channel phase response


ES = freqz(W); % Equalizer Spectrum
EMR = 10*log10(real(ES).^2 + imag(ES).^2); % Equalizer magnitude response
EPR = imag(ES)./real(ES); % Equalizer phase response

TS = ES.*CS; %Total response of Channel and Equalizer
TMR = 10*log10(real(TS).^2 + imag(TS).^2); % Total magnitude response
TPR = imag(TS)./real(TS); % Total phase response

%----------------------Plots---------------------------------------------

figure(); % MSE
plot(10*log10(MJ),'->k','linewidth',2)
hg=legend('MSE');
grid on;
xlabel('iterations');
ylabel('Mean squared error (dB)');
title('Cost function');

figure();
subplot(2,1,1)
plot(NF,IMR,'b','linewidth',2)
hold on
plot(NF,EMR,'--r','linewidth',2)
hg=legend('Inverse Channel','Equalizer');
grid on
xlabel('Normalized Frequency');
ylabel('Magnitude (dB)');
title('Magnitude response');

subplot(2,1,2)
plot(NF, IPR,'g','linewidth',2)
hold on
plot(NF, EPR,'--b','linewidth',2)
hg=legend('Inverse Channel','Equalizer');
grid on
xlabel('Normalized Frequency');
ylabel('Phase shift (rad)');
title('Phase response');

figure();
subplot(2,1,1)
plot(NF,-1*IMR,'b','linewidth',2)
hold on
plot(NF,EMR,'--r','linewidth',2)
hold on
plot(NF,TMR,'--g','linewidth',2)
hg=legend('Channel','Equalizer','Total Res of Equ and Channel');
grid on
xlabel('Normalized Frequency');
ylabel('Magnitude (dB)');
title('Magnitude response');

subplot(2,1,2)
plot(NF,-1*IPR,'b','linewidth',2)
hold on
plot(NF, EPR,'--r','linewidth',2)
hold on
plot(NF, TPR,'--g','linewidth',2)
hg=legend('Channel','Equalizer','Total Res of Equ and Channel');
grid on
xlabel('Normalized Frequency');
ylabel('Phase shift (rad)');
title('Phase response');

%%
%----------------------------Equalizer---------------------------------

Equalized_signal = filter(W,1,Y);
scatterplot(p);
title('constellation of the originally generated symbols')
scatterplot(Equalized_signal);
title('constellation of the Equaled Signal')
scatterplot(Y);
title('constellation of the Unequaled Signal')

%%
%----------------------------Detector----------------------------------

Y = Equalized_signal(6:100005); %correcting delay
m1 = [];
m2 = [];
for n = length(t):length(t):length(Y)
  y1 = cos(2*pi*f*t);% inphase component
  y2 = sin(2*pi*f*t);% quadrature component
  mm1 = y1.*Y((n-(length(t)-1)):n);                                    
  mm2 = y2.*Y((n-(length(t)-1)):n);                                    
  z1 = trapz(t,mm1);% integration
  z2 = trapz(t,mm2) ;% integration
  zz1 = round(2*z1/symbol_period);
  zz2 = round(2*z2/symbol_period);
  m1 = [m1 zz1];
  m2 = [m2 zz2];
end

%---------------------------De-mapping---------------------------------

clear i;
clear j;
for (k = 1:1:length(m1))  
gt(k) = m1(k)+j*m2(k);
end
scatterplot(gt);
title('constellation of the Output Signal')
ax = qamdemod(gt,4,'gray');
bi_in = dec2bin(ax);
[row col] = size(bi_in);
p = 1;
 for(i =1:1:row)
     for(j =1:1:col)
         bits_output(p) = str2num(bi_in(i,j));
         p = p+1;
     end
 end 

disp('First 10 Generated bits :');
disp(bits(1:10));
disp('First 10 Output bits :');
disp(bits_output(1:10));

bit_output = []; 
for n = 1:1:length(x)
    if bits_output(n) == 1;
       se = ones(1,100);
    else bits_output(n) == 0;
        se = zeros(1,100);
    end
     bit_output = [bit_output se];

end

figure();
subplot(2,1,1);
plot(t1,bit(1:1000),'lineWidth',2.5);grid on;
axis([ 0 bit_period*10 -.5 1.5]);
ylabel('amplitude(volt)');
xlabel(' time(sec)');
title('transmitting information as digital signal(first 10 bits)');

subplot(2,1,2);
plot(t1,bit_output(1:1000),'lineWidth',2.5);grid on;
axis([ 0 bit_period*10 -.5 1.5]);
ylabel('amplitude(volt)');
xlabel(' time(sec)');
title('Recieved information as digital signal(first 10 bits)');

Num_of_correct_predictions = sum(bits == bits_output);
Num_of_incorrect_predictions = length(bits) - Num_of_correct_predictions;
SER = (Num_of_incorrect_predictions/length(bits))*100;
accuracy = (Num_of_correct_predictions/length(bits))*100;

fprintf(1,'\nNumber of correct predictions out of %d bits : %d'...
    ,length(bits),Num_of_correct_predictions);
fprintf(1,'\nNumber of incorrect predictions out of %d bits : %d'...
    ,length(bits),Num_of_incorrect_predictions);
fprintf(1,'\nSER(Symbol Error Rate) : %d percent',SER);
fprintf(1,'\naccuracy : %d percent\n',accuracy);

%%
%----------------------Simulation for different SNRs---------------------

clear;clc;
SNR = 1:20;
for counter = 1:length(SNR)
    num_of_symbols = 1e3;
    bits = rand(1,2*num_of_symbols)>=0.5;%generating random bits
    bits_reshape=reshape(bits,2,num_of_symbols)';
    for(j=1:1:num_of_symbols)
        for(i=1:1:2)
            symbols(j,i)=num2str(bits_reshape(j,i));
        end
    end
    symbols=bin2dec(symbols);
    symbols=symbols';

    p = qammod(symbols,4,'gray'); %constalation diagram for 4-QAM
    clear symbols
    symbol_rate = 4e7;
    symbol_period = 1/symbol_rate;
    bit_rate = 2*symbol_rate;
    bit_period = 1/(2*symbol_rate);
    Real = real(p);
    Imag = imag(p);
    f = bit_rate;
    t = symbol_period/100:symbol_period/100:symbol_period;
    trans_signal = [];
    for(k = 1:1:length(Real))
        yr = Real(k)*cos(2*pi*f*t);% inphase or real component
        yim = Imag(k)*sin(2*pi*f*t);% Quadrature or imagenary component 
        trans_signal(k,:) = yr+yim;
    end
    trans_signal = trans_signal';
    trans_signal = trans_signal(:);
    trans_signal = trans_signal';
    trans_signal_fft = fftshift(fft(trans_signal));

    Up_sample_factor = 4;
    alpha = 0.4; %rolloff factor
    Nd = 10; %duration
    h_Tx = rcosdesign(alpha,Nd,Up_sample_factor);
    h_Tx_fft = fftshift(fft(h_Tx));
    trans_signal_final = upfirdn(trans_signal,h_Tx,Up_sample_factor);
    trans_signal_final = trans_signal_final(19:400018);

    coeff = [-0.4945 - 0.4611i, -0.2813 + 0.8457i, -0.0720 - 0.5260i,...
     0.1201 - 0.4302i, 0.2451 + 0.1999i];
    recieved_signal = filter(coeff,1,trans_signal_final);
    Noisy_signal = awgn(recieved_signal,SNR(counter),'measured','linear');
    Y = upfirdn(Noisy_signal,h_Tx,1,Up_sample_factor);
    N = 1000; % Number of samples              
    data = randi([0 1],1,N);% Random signal
    d = qammod(data,4);    % 4_QAM modulation
    r = filter(coeff,1,data);       % Signal after passing through channel
    x = awgn(r, SNR(counter));           % Noisy Signal after channel
    iterations = 10000; % Number of iterations
    eta = 0.0001;  % Step size
    order = 30;    % Order of the equalizer
    U = zeros(1,order); % Input frame
    W = zeros(1,order); % Initial Weigths
    for k = 1 : iterations
        for n = 1 : N
            U(1,2:end) = U(1,1:end-1);  % Sliding window
            U(1,1) = x(n);   % Present Input
     
            y = (W)*U'; % Calculating output of LMS
            e = d(n) - y;  % Instantaneous error 
            W = W +  eta * e * U ; % Weight update rule of LMS
            J(k,n) = e * e'; % Instantaneous square error
        end
    end
    Equalized_signal = filter(W,1,Y);
    Y = Equalized_signal(6:100005); %correcting delay
    m1 = [];
    m2 = [];
    for n = length(t):length(t):length(Y)
        y1 = cos(2*pi*f*t);% inphase component
        y2 = sin(2*pi*f*t);% quadrature component
        mm1 = y1.*Y((n-(length(t)-1)):n);                                    
        mm2 = y2.*Y((n-(length(t)-1)):n);                                    
        z1 = trapz(t,mm1);% integration
        z2 = trapz(t,mm2) ;% integration
        zz1 = round(2*z1/symbol_period);
        zz2 = round(2*z2/symbol_period);
        m1 = [m1 zz1];
        m2 = [m2 zz2];
    end
    clear i;
    clear j;
    for (k = 1:1:length(m1))  
        gt(k) = m1(k)+j*m2(k);
    end

    ax = qamdemod(gt,4,'gray');
    bi_in = dec2bin(ax);
    [row col] = size(bi_in);
    p = 1;
    for(i =1:1:row)
        for(j =1:1:col)
            bits_output(p) = str2num(bi_in(i,j));
            p = p+1;
        end
    end 
    Num_of_correct_predictions = sum(bits == bits_output);
    Num_of_incorrect_predictions = length(bits) - Num_of_correct_predictions;
    SER(counter) = (Num_of_incorrect_predictions/length(bits))*100;
end

figure();
plot(SNR,SER);
title('SER(Symbol Rate Error percent) vs. SNR');
xlabel('SNR');
ylabel('SER');


    















    













