%   _  _   _   _   _   _     _   _ 
%  /  |_   _) / \ / \ |_    | \ /  
%  \_ |_   _) \_/ \_/ |_)   |_/ \_ 
%                                  
%   Phase 3 Modulation with Error correction codes

close all;
clear;


% Check that user has the Communication Toolbox installed.
hasIPT = license('test', 'communication_toolbox');
if ~hasIPT
  % User does not have the toolbox installed.
  message = sprintf('Sorry, but you do not seem to have the Communications Toolbox.\nDo you want to try to continue anyway?');
  reply = questdlg(message, 'Toolbox missing', 'Yes', 'No', 'Yes');
  if strcmpi(reply, 'No')
    return;
  end
end

fc = 10000; % Carrier Frequency
fs = 16*fc; % Sampling Frequency
N = 1024; % Number of bits to transmit
M = 16; % Modulation order (alphabet size or number of points in signal constellation)
k = log2(M); % Number of bits per symbol
data_rate = 1000; % Message data rate = 1000bps
Tb = 1/data_rate; % bit period
Tsym = k*Tb; % symbol period
Ts = 1/fs; % Sampling period
signal_p = 1; % Signal power

%Reed Solomon
B = 8; % Block length
K = 4; % Message length
rsEnc = comm.RSEncoder('BitInput',true,'CodewordLength',B,'MessageLength',K);
rsDec = comm.RSDecoder('BitInput',true,'CodewordLength',B,'MessageLength',K);

%Turbo
intrlvrInd = randperm(N); % Interleaver indices
turboEnc = comm.TurboEncoder('InterleaverIndicesSource','Input port');
turboDec = comm.TurboDecoder('InterleaverIndicesSource','Input port','NumIterations',4);


% % % % % % % % % Generate Message Data % % % % % % % % %

sig_data_raw = randi([0, 1], 1, N);

% % % % % % % % % Encode Message Data % % % % % % % % %

% Cyclic data
sig_data_cyclic = Encode_Cyclic(sig_data_raw);

% Hamming data
sig_data_hamming  = Encode_Hamming(sig_data_raw);

% Reed solomon data
sig_data_rs  = rsEnc(sig_data_raw(:));
sig_data_rs = sig_data_rs';

% Turbo data
sig_data_turbo = turboEnc(sig_data_raw(:),intrlvrInd);
sig_data_turbo = sig_data_turbo';

% % % % % % % % % Modulate Message Data % % % % % % % % %

ook_sig = OOK_Modulate(sig_data_raw,fc,fs,Tb);
ook_sig_cyclic = OOK_Modulate(sig_data_cyclic,fc,fs,Tb);
ook_sig_hamming = OOK_Modulate(sig_data_hamming,fc,fs,Tb);
ook_sig_turbo = OOK_Modulate(sig_data_turbo,fc,fs,Tb);
ook_sig_rs = OOK_Modulate(sig_data_rs,fc,fs,Tb);

bpsk_sig = BPSK_Modulate(sig_data_raw,fc,fs,Tb);
bpsk_sig_cyclic = BPSK_Modulate(sig_data_cyclic,fc,fs,Tb);
bpsk_sig_hamming = BPSK_Modulate(sig_data_hamming,fc,fs,Tb);
bpsk_sig_turbo = BPSK_Modulate(sig_data_turbo,fc,fs,Tb);
bpsk_sig_rs = BPSK_Modulate(sig_data_rs,fc,fs,Tb);

bfsk_sig = BFSK_Modulate(sig_data_raw,fs,2*fc,fc,Tb); % 1 = 20Khz carrier , 0 = 10Khz carrier
bfsk_sig_cyclic = BFSK_Modulate(sig_data_cyclic,fs,2*fc,fc,Tb); % 1 = 20Khz carrier , 0 = 10Khz carrier
bfsk_sig_hamming = BFSK_Modulate(sig_data_hamming,fs,2*fc,fc,Tb); % 1 = 20Khz carrier , 0 = 10Khz carrier
bfsk_sig_turbo = BFSK_Modulate(sig_data_turbo,fs,2*fc,fc,Tb); % 1 = 20Khz carrier , 0 = 10Khz carrier
bfsk_sig_rs = BFSK_Modulate(sig_data_rs,fs,2*fc,fc,Tb); % 1 = 20Khz carrier , 0 = 10Khz carrier

[qam16_sig, ~] = QAM16_Modulate(sig_data_raw,'gray',fc,fs,Tsym);
[qam16_sig_cyclic, ~] = QAM16_Modulate(sig_data_cyclic,'gray',fc,fs,Tsym);
[qam16_sig_hamming, ~] = QAM16_Modulate(sig_data_hamming,'gray',fc,fs,Tsym);
[qam16_sig_turbo, ~] = QAM16_Modulate(sig_data_turbo,'gray',fc,fs,Tsym);
[qam16_sig_rs, ~] = QAM16_Modulate(sig_data_rs,'gray',fc,fs,Tsym);

% % % % % % % % % SNR Loop % % % % % % % % %

max_snr = 0;
snr_inc = 2;
start = -100;

% Use to index error array
index = 1; 

% Store bit error rate for each iteration of the loop
qam16_err_arr = zeros((max_snr*2)/snr_inc,0);
qam16_err_arr_cyclic = zeros((max_snr*2)/snr_inc,0);
qam16_err_arr_hamming = zeros((max_snr*2)/snr_inc,0);
qam16_err_arr_turbo = zeros((max_snr*2)/snr_inc,0);
qam16_err_arr_rs = zeros((max_snr*2)/snr_inc,0);

ook_err_arr = zeros((max_snr*2)/snr_inc,0);
ook_err_arr_cyclic = zeros((max_snr*2)/snr_inc,0);
ook_err_arr_hamming = zeros((max_snr*2)/snr_inc,0);
ook_err_arr_turbo = zeros((max_snr*2)/snr_inc,0);
ook_err_arr_rs = zeros((max_snr*2)/snr_inc,0);

bpsk_err_arr = zeros((max_snr*2)/snr_inc,0);
bpsk_err_arr_cyclic = zeros((max_snr*2)/snr_inc,0);
bpsk_err_arr_hamming = zeros((max_snr*2)/snr_inc,0);
bpsk_err_arr_turbo = zeros((max_snr*2)/snr_inc,0);
bpsk_err_arr_rs = zeros((max_snr*2)/snr_inc,0);

bfsk_err_arr = zeros((max_snr*2)/snr_inc,0);
bfsk_err_arr_cyclic = zeros((max_snr*2)/snr_inc,0);
bfsk_err_arr_hamming = zeros((max_snr*2)/snr_inc,0);
bfsk_err_arr_turbo = zeros((max_snr*2)/snr_inc,0);
bfsk_err_arr_rs = zeros((max_snr*2)/snr_inc,0);

% Main Process Loop
for snr = start:snr_inc:max_snr
    
    % Transmitted Signal passes through AWGN Channel
    qam16_rcv = awgn(qam16_sig,snr,'measured');
    qam16_rcv_cyclic = awgn(qam16_sig_cyclic,snr,'measured');
    qam16_rcv_hamming = awgn(qam16_sig_hamming,snr,'measured');
    qam16_rcv_turbo = awgn(qam16_sig_turbo,snr,'measured');
    qam16_rcv_rs = awgn(qam16_sig_rs,snr,'measured');
    
    ook_rcv = awgn(ook_sig,snr,'measured');
    ook_rcv_cyclic = awgn(ook_sig_cyclic,snr,'measured');
    ook_rcv_hamming = awgn(ook_sig_hamming,snr,'measured');
    ook_rcv_turbo = awgn(ook_sig_turbo,snr,'measured');
    ook_rcv_rs = awgn(ook_sig_rs,snr,'measured');
    
    bpsk_rcv = awgn(bpsk_sig,snr,'measured');
    bpsk_rcv_cyclic = awgn(bpsk_sig_cyclic,snr,'measured');
    bpsk_rcv_hamming = awgn(bpsk_sig_hamming,snr,'measured');
    bpsk_rcv_turbo = awgn(bpsk_sig_turbo,snr,'measured');
    bpsk_rcv_rs = awgn(bpsk_sig_rs,snr,'measured');
    
    bfsk_rcv = awgn(bfsk_sig,snr,'measured');
    bfsk_rcv_cyclic = awgn(bfsk_sig_cyclic,snr,'measured');
    bfsk_rcv_hamming = awgn(bfsk_sig_hamming,snr,'measured');    
    bfsk_rcv_turbo = awgn(bfsk_sig_turbo,snr,'measured');  
    bfsk_rcv_rs = awgn(bfsk_sig_rs,snr,'measured');  

    % Demodulate Signal
    [qam16_demod_sig, ~] = QAM16_Demodulate(qam16_rcv,'gray',fc,fs,Tsym,'bit');
    [qam16_demod_sig_cyclic, ~] = QAM16_Demodulate(qam16_rcv_cyclic,'gray',fc,fs,Tsym,'bit');
    [qam16_demod_sig_hamming, ~] = QAM16_Demodulate(qam16_rcv_hamming,'gray',fc,fs,Tsym,'bit');
    [qam16_demod_sig_turbo, ~] = QAM16_Demodulate(qam16_rcv_turbo,'gray',fc,fs,Tsym,'approxllr');
    [qam16_demod_sig_rs, ~] = QAM16_Demodulate(qam16_rcv_rs,'gray',fc,fs,Tsym,'bit');
    
    ook_demod_sig = OOK_Demodulate(ook_rcv,fc,fs,Tb,'bit');
    ook_demod_sig_cyclic = OOK_Demodulate(ook_rcv_cyclic,fc,fs,Tb,'bit');
    ook_demod_sig_hamming = OOK_Demodulate(ook_rcv_hamming,fc,fs,Tb,'bit');    
    ook_demod_sig_turbo = OOK_Demodulate(ook_rcv_turbo,fc,fs,Tb,'soft');
    ook_demod_sig_rs = OOK_Demodulate(ook_rcv_rs,fc,fs,Tb,'bit');
     
    bpsk_demod_sig = BPSK_Demodulate(bpsk_rcv,fc,fs,Tb,'bit');
    bpsk_demod_sig_cyclic = BPSK_Demodulate(bpsk_rcv_cyclic,fc,fs,Tb,'bit');
    bpsk_demod_sig_hamming = BPSK_Demodulate(bpsk_rcv_hamming,fc,fs,Tb,'bit');
    bpsk_demod_sig_turbo = BPSK_Demodulate(bpsk_rcv_turbo,fc,fs,Tb,'soft');
    bpsk_demod_sig_rs = BPSK_Demodulate(bpsk_rcv_rs,fc,fs,Tb,'bit');
    
    bfsk_demod_sig = BFSK_Demodulate(bfsk_rcv,fs,2*fc,fc,Tb,'bit');  
    bfsk_demod_sig_cyclic = BFSK_Demodulate(bfsk_rcv_cyclic,fs,2*fc,fc,Tb,'bit'); 
    bfsk_demod_sig_hamming = BFSK_Demodulate(bfsk_rcv_hamming,fs,2*fc,fc,Tb,'bit');
    bfsk_demod_sig_turbo = BFSK_Demodulate(bfsk_rcv_turbo,fs,2*fc,fc,Tb,'soft');
    bfsk_demod_sig_rs = BFSK_Demodulate(bfsk_rcv_rs,fs,2*fc,fc,Tb,'bit');
     
%     %Purposely add error
%     qam16_demod_sig=Add_Error1(qam16_demod_sig);
%     qam16_demod_sig_cyclic=Add_Error1(qam16_demod_sig_cyclic);
%     qam16_demod_sig_hamming=Add_Error1(qam16_demod_sig_hamming);
%     qam16_demod_sig_turbo=Add_Error2(qam16_demod_sig_turbo);
%     qam16_demod_sig_rs=Add_Error1(qam16_demod_sig_rs);
%     
%     ook_demod_sig=Add_Error1(ook_demod_sig);
%     ook_demod_sig_cyclic=Add_Error1(ook_demod_sig_cyclic);
%     ook_demod_sig_hamming=Add_Error1(ook_demod_sig_hamming);
%     ook_demod_sig_turbo=Add_Error2(ook_demod_sig_turbo);
%     ook_demod_sig_rs=Add_Error1(ook_demod_sig_rs);
%     
%     bpsk_demod_sig=Add_Error1(bpsk_demod_sig);
%     bpsk_demod_sig_cyclic=Add_Error1(bpsk_demod_sig_cyclic);
%     bpsk_demod_sig_hamming=Add_Error1(bpsk_demod_sig_hamming);
%     bpsk_demod_sig_turbo=Add_Error2(bpsk_demod_sig_turbo);
%     bpsk_demod_sig_rs=Add_Error1(bpsk_demod_sig_rs);
%     
%     bfsk_demod_sig=Add_Error1(bfsk_demod_sig);
%     bfsk_demod_sig_cyclic=Add_Error1(bfsk_demod_sig_cyclic);
%     bfsk_demod_sig_hamming=Add_Error1(bfsk_demod_sig_hamming);
%     bfsk_demod_sig_turbo=Add_Error2(bfsk_demod_sig_turbo);
%     bfsk_demod_sig_rs=Add_Error1(bfsk_demod_sig_rs);
    
    % Decode signal
    qam16_demod_sig_cyclic_dec = Decode_Cyclic(qam16_demod_sig_cyclic);
    qam16_demod_sig_hamming_dec = Decode_Hamming(qam16_demod_sig_hamming);
    qam16_demod_sig_turbo_dec = turboDec(-(qam16_demod_sig_turbo(:)),intrlvrInd);
    qam16_demod_sig_rs_dec = rsDec(qam16_demod_sig_rs(:));
    
    ook_demod_sig_cyclic_dec = Decode_Cyclic(ook_demod_sig_cyclic);
    ook_demod_sig_hamming_dec = Decode_Hamming(ook_demod_sig_hamming);
    ook_demod_sig_turbo_dec = turboDec((ook_demod_sig_turbo(:)),intrlvrInd);
    ook_demod_sig_rs_dec = rsDec(ook_demod_sig_rs(:));
    
    bpsk_demod_sig_cyclic_dec = Decode_Cyclic(bpsk_demod_sig_cyclic);
    bpsk_demod_sig_hamming_dec = Decode_Hamming(bpsk_demod_sig_hamming);
    bpsk_demod_sig_turbo_dec = turboDec((bpsk_demod_sig_turbo(:)),intrlvrInd);
    bpsk_demod_sig_rs_dec = rsDec(bpsk_demod_sig_rs(:));
    
    bfsk_demod_sig_cyclic_dec = Decode_Cyclic(bfsk_demod_sig_cyclic);
    bfsk_demod_sig_hamming_dec = Decode_Hamming(bfsk_demod_sig_hamming);
    bfsk_demod_sig_turbo_dec = turboDec((bfsk_demod_sig_turbo(:)),intrlvrInd);
    bfsk_demod_sig_rs_dec = rsDec(bfsk_demod_sig_rs(:));   
    
    % Check bit error rate and save the values for plotting
    [~, ber] = biterr(sig_data_raw,ook_demod_sig);
    ook_err_arr(index) = ber*100; 
    [~, ber] = biterr(sig_data_raw,ook_demod_sig_cyclic_dec);
    ook_err_arr_cyclic(index) = ber*100;
    [~, ber] = biterr(sig_data_raw,ook_demod_sig_hamming_dec);
    ook_err_arr_hamming(index) = ber*100;
    [~, ber] = biterr(sig_data_raw,ook_demod_sig_turbo_dec');
    ook_err_arr_turbo(index) = ber*100;
    [~, ber] = biterr(sig_data_raw,ook_demod_sig_rs_dec');
    ook_err_arr_rs(index) = ber*100;   

    [~, ber] = biterr(sig_data_raw,bpsk_demod_sig);
    bpsk_err_arr(index) = ber*100;
    [~, ber] = biterr(sig_data_raw,bpsk_demod_sig_cyclic_dec);
    bpsk_err_arr_cyclic(index) = ber*100;
    [~, ber] = biterr(sig_data_raw,bpsk_demod_sig_hamming_dec);
    bpsk_err_arr_hamming(index) = ber*100;
    [~, ber] = biterr(sig_data_raw,bpsk_demod_sig_turbo_dec');
    bpsk_err_arr_turbo(index) = ber*100;
    [~, ber] = biterr(sig_data_raw,bpsk_demod_sig_rs_dec');
    bpsk_err_arr_rs(index) = ber*100;   

    [~, ber] = biterr(sig_data_raw,bfsk_demod_sig);
    bfsk_err_arr(index) = ber*100; 
    [~, ber] = biterr(sig_data_raw,bfsk_demod_sig_cyclic_dec);
    bfsk_err_arr_cyclic(index) = ber*100;
    [~, ber] = biterr(sig_data_raw,bfsk_demod_sig_hamming_dec);
    bfsk_err_arr_hamming(index) = ber*100;
    [~, ber] = biterr(sig_data_raw,bfsk_demod_sig_turbo_dec');
    bfsk_err_arr_turbo(index) = ber*100;
    [~, ber] = biterr(sig_data_raw,bfsk_demod_sig_rs_dec');
    bfsk_err_arr_rs(index) = ber*100; 
    
    
    [~, ber] = biterr(sig_data_raw,qam16_demod_sig);
    qam16_err_arr(index) = ber*100;     
    [~, ber] = biterr(sig_data_raw,qam16_demod_sig_cyclic_dec);
    qam16_err_arr_cyclic(index) = ber*100;
    [~, ber] = biterr(sig_data_raw,qam16_demod_sig_hamming_dec);
    qam16_err_arr_hamming(index) = ber*100;
    [~, ber] = biterr(sig_data_raw,qam16_demod_sig_turbo_dec');
    qam16_err_arr_turbo(index) = ber*100;
    [~, ber] = biterr(sig_data_raw,qam16_demod_sig_rs_dec');
    qam16_err_arr_rs(index) = ber*100; 
    
    
    index = index + 1;
    
    
    
end

% Graph and Plot the result           
figure(1)
snr = start:snr_inc:max_snr;
semilogx(snr,qam16_err_arr_cyclic,'b-');
hold on;
semilogx(snr,qam16_err_arr_hamming,'r-');
hold on;
semilogx(snr,qam16_err_arr_turbo,'g-');
hold on;
semilogx(snr,qam16_err_arr,'m-');
hold on;
semilogx(snr,qam16_err_arr_rs,'k-');
title({'16-Quadrature Amplitude Modulation','SNR to Bit Error Rate'});
ylabel('Bit Error Rate (%)');
xlabel('SNR[dB]');
grid on;
legend('QAM (Cyclic n=12 k=4)','QAM (Hamming n=7 k=4)','QAM (Turbo interations = 4)','QAM (No Encoding)','QAM (Reed Solomon n=8 k=4)');

% Graph and Plot the result           
figure(2)
semilogx(snr,ook_err_arr_cyclic,'b-');
hold on;
semilogx(snr,ook_err_arr_hamming,'r-');
hold on;
semilogx(snr,ook_err_arr_turbo,'g-');
hold on;
semilogx(snr,ook_err_arr,'m-');
hold on;
semilogx(snr,ook_err_arr_rs,'k-');
title({'On-Off Keying','SNR to Bit Error Rate'});
ylabel('Bit Error Rate (%)');
xlabel('SNR[dB]');
grid on;
legend('OOK (Cyclic n=12 k=4)','OOK (Hamming n=7 k=4)','OOK (Turbo interations = 4)','OOK (No Encoding)','OOK (Reed Solomon n=8 k=4)');

% Graph and Plot the result           
figure(3)
semilogx(snr,bpsk_err_arr_cyclic,'b-');
hold on;
semilogx(snr,bpsk_err_arr_hamming,'r-');
hold on;
semilogx(snr,bpsk_err_arr_turbo,'g-');
hold on;
semilogx(snr,bpsk_err_arr,'m-');
hold on;
semilogx(snr,bpsk_err_arr_rs,'k-');
title({'Binary Phase Shift Keying','SNR to Bit Error Rate'});
ylabel('Bit Error Rate (%)');
xlabel('SNR[dB]');
grid on;
legend('BPSK (Cyclic n=12 k=4)','BPSK (Hamming n=7 k=4)','BPSK (Turbo interations = 4)','BPSK (No Encoding)','BPSK (Reed Solomon n=8 k=4)');

% Graph and Plot the result           
figure(4)
semilogx(snr,bfsk_err_arr_cyclic,'b-');
hold on;
semilogx(snr,bfsk_err_arr_hamming,'r-');
hold on;
semilogx(snr,bfsk_err_arr_turbo,'g-');
hold on;
semilogx(snr,bfsk_err_arr,'m-');
hold on;
semilogx(snr,bfsk_err_arr_rs,'k-');
title({'Binary Frequency Shift Keying','SNR to Bit Error Rate'});
ylabel('Bit Error Rate (%)');
xlabel('SNR[dB]');
grid on;
legend('BFSK (Cyclic n=12 k=4)','BFSK (Hamming n=7 k=4)','BFSK (Turbo interations = 4)','BFSK (No Encoding)','BFSK (Reed Solomon n=8 k=4)');


% Phase 3 plots


%%% OOK

% Set correct data before plotting

modDataC = ook_sig_cyclic;
modDataH = ook_sig_hamming;
modDataR = ook_sig_rs;
modDataT = ook_sig_turbo;

modDataNoiseC = ook_rcv_cyclic;
modDataNoiseH = ook_rcv_hamming;
modDataNoiseR = ook_rcv_rs;
modDataNoiseT = ook_rcv_turbo;

decDataC = ook_demod_sig_cyclic_dec;
decDataH = ook_demod_sig_hamming_dec;
decDataR = ook_demod_sig_rs_dec;
decDataT = ook_demod_sig_turbo_dec;

% Plot message signal
figure(5);
t_in = 0:1:N;
subplot(5,4,1)
stairs(t_in(1:9),[sig_data_raw(1:8),sig_data_raw(8)]);
axis([0 8 -0.5 1.5])
title("Message Signal");

subplot(5,4,2)
stairs(t_in(1:9),[sig_data_raw(1:8),sig_data_raw(8)]);
axis([0 8 -0.5 1.5])
title("Message Signal");

subplot(5,4,3)
stairs(t_in(1:9),[sig_data_raw(1:8),sig_data_raw(8)]);
axis([0 8 -0.5 1.5])
title("Message Signal");

subplot(5,4,4)
stairs(t_in(1:9),[sig_data_raw(1:8),sig_data_raw(8)]);
axis([0 8 -0.5 1.5])
title("Message Signal");

%----------------------------
% Plot all the encoded signals

subplot(5,4,5)
stairs(t_in(1:9),[sig_data_cyclic(1:8),sig_data_cyclic(8)]);
axis([0 8 -0.5 1.5])
title("Cyclic Encoded Signal");

subplot(5,4,6)
stairs(t_in(1:9),[sig_data_hamming(1:8),sig_data_hamming(8)]);
axis([0 8 -0.5 1.5])
title("Hamming Encoded Signal");

subplot(5,4,7)
stairs(t_in(1:9),[sig_data_rs(1:8),sig_data_rs(8)]);
axis([0 8 -0.5 1.5])
title("Reed Solomon Encoded Signal");

subplot(5,4,8)
stairs(t_in(1:9),[sig_data_turbo(1:8),sig_data_turbo(8)]);
axis([0 8 -0.5 1.5])
title("Turbo Encoded Signal");

%----------------------------
% Plot all the modulated signals

t=0:Ts:(2*Tsym)-Ts; % Duration of Message Signal

subplot(5,4,9)
plot(modDataC(1:length(t)),'r-');
title("Modulated Signal");
xlim([0,length(t)]);

subplot(5,4,10)
plot(modDataH(1:length(t)),'r-');
title("Modulated Signal");
xlim([0,length(t)]);

subplot(5,4,11)
plot(modDataR(1:length(t)),'r-');
title("Modulated Signal");
xlim([0,length(t)]);

subplot(5,4,12)
plot(modDataT(1:length(t)),'r-');
title("Modulated Signal");
xlim([0,length(t)]);

%----------------------------
% Plot all the modulated signals + noise

subplot(5,4,13)
plot(modDataNoiseC(1:length(t)),'m-');
title("Modulated Signal + Noise");
xlim([0,length(t)]);

subplot(5,4,14)
plot(modDataNoiseH(1:length(t)),'m-');
title(" Modulated Signal + Noise");
xlim([0,length(t)]);

subplot(5,4,15)
plot(modDataNoiseR(1:length(t)),'m-');
title("Modulated Signal + Noise");
xlim([0,length(t)]);

subplot(5,4,16)
plot(modDataNoiseT(1:length(t)),'m-');
title("Modulated Signal + Noise");
xlim([0,length(t)]);

%----------------------------
% Plot decoded signal

subplot(5,4,17)
stairs(t_in(1:9),[decDataC(1:8),decDataC(8)]);
axis([0 8 -0.5 1.5])
title("Demodulated and Decoded Signal");

subplot(5,4,18)
stairs(t_in(1:9),[decDataH(1:8),decDataH(8)]);
axis([0 8 -0.5 1.5])
title("Demodulated and Decoded Signal");

subplot(5,4,19)
stairs(t_in(1:9),[decDataR(1:8)',decDataR(8)']);
axis([0 8 -0.5 1.5])
title("Demodulated and Decoded Signal");

subplot(5,4,20)
stairs(t_in(1:9),[decDataT(1:8)',decDataT(8)']);
axis([0 8 -0.5 1.5])
title("Demodulated and Decoded Signal");

sgtitle('On-Off Keying');

%%% End of OOK Plots

%%% BPSK

% Set correct data before plotting

modDataC = bpsk_sig_cyclic;
modDataH = bpsk_sig_hamming;
modDataR = bpsk_sig_rs;
modDataT = bpsk_sig_turbo;

modDataNoiseC = bpsk_rcv_cyclic;
modDataNoiseH = bpsk_rcv_hamming;
modDataNoiseR = bpsk_rcv_rs;
modDataNoiseT = bpsk_rcv_turbo;

decDataC = bpsk_demod_sig_cyclic_dec;
decDataH = bpsk_demod_sig_hamming_dec;
decDataR = bpsk_demod_sig_rs_dec;
decDataT = bpsk_demod_sig_turbo_dec;

% Plot message signal
figure(6);
t_in = 0:1:N;
subplot(5,4,1)
stairs(t_in(1:9),[sig_data_raw(1:8),sig_data_raw(8)]);
axis([0 8 -0.5 1.5])
title("Message Signal");

subplot(5,4,2)
stairs(t_in(1:9),[sig_data_raw(1:8),sig_data_raw(8)]);
axis([0 8 -0.5 1.5])
title("Message Signal");

subplot(5,4,3)
stairs(t_in(1:9),[sig_data_raw(1:8),sig_data_raw(8)]);
axis([0 8 -0.5 1.5])
title("Message Signal");

subplot(5,4,4)
stairs(t_in(1:9),[sig_data_raw(1:8),sig_data_raw(8)]);
axis([0 8 -0.5 1.5])
title("Message Signal");

%----------------------------
% Plot all the encoded signals

subplot(5,4,5)
stairs(t_in(1:9),[sig_data_cyclic(1:8),sig_data_cyclic(8)]);
axis([0 8 -0.5 1.5])
title("Cyclic Encoded Signal");

subplot(5,4,6)
stairs(t_in(1:9),[sig_data_hamming(1:8),sig_data_hamming(8)]);
axis([0 8 -0.5 1.5])
title("Hamming Encoded Signal");

subplot(5,4,7)
stairs(t_in(1:9),[sig_data_rs(1:8),sig_data_rs(8)]);
axis([0 8 -0.5 1.5])
title("Reed Solomon Encoded Signal");

subplot(5,4,8)
stairs(t_in(1:9),[sig_data_turbo(1:8),sig_data_turbo(8)]);
axis([0 8 -0.5 1.5])
title("Turbo Encoded Signal");

%----------------------------
% Plot all the modulated signals

t=0:Ts:(2*Tsym)-Ts; % Duration of Message Signal

subplot(5,4,9)
plot(modDataC(1:length(t)),'r-');
title("Modulated Signal");
xlim([0,length(t)]);

subplot(5,4,10)
plot(modDataH(1:length(t)),'r-');
title("Modulated Signal");
xlim([0,length(t)]);

subplot(5,4,11)
plot(modDataR(1:length(t)),'r-');
title("Modulated Signal");
xlim([0,length(t)]);

subplot(5,4,12)
plot(modDataT(1:length(t)),'r-');
title("Modulated Signal");
xlim([0,length(t)]);

%----------------------------
% Plot all the modulated signals + noise

subplot(5,4,13)
plot(modDataNoiseC(1:length(t)),'m-');
title("Modulated Signal + Noise");
xlim([0,length(t)]);

subplot(5,4,14)
plot(modDataNoiseH(1:length(t)),'m-');
title(" Modulated Signal + Noise");
xlim([0,length(t)]);

subplot(5,4,15)
plot(modDataNoiseR(1:length(t)),'m-');
title("Modulated Signal + Noise");
xlim([0,length(t)]);

subplot(5,4,16)
plot(modDataNoiseT(1:length(t)),'m-');
title("Modulated Signal + Noise");
xlim([0,length(t)]);

%----------------------------
% Plot decoded signal

subplot(5,4,17)
stairs(t_in(1:9),[decDataC(1:8),decDataC(8)]);
axis([0 8 -0.5 1.5])
title("Demodulated and Decoded Signal");

subplot(5,4,18)
stairs(t_in(1:9),[decDataH(1:8),decDataH(8)]);
axis([0 8 -0.5 1.5])
title("Demodulated and Decoded Signal");

subplot(5,4,19)
stairs(t_in(1:9),[decDataR(1:8)',decDataR(8)']);
axis([0 8 -0.5 1.5])
title("Demodulated and Decoded Signal");

subplot(5,4,20)
stairs(t_in(1:9),[decDataT(1:8)',decDataT(8)']);
axis([0 8 -0.5 1.5])
title("Demodulated and Decoded Signal");

sgtitle('Binary Phase Shift Keying');

%%% End of BPSK Plots

%%% BFSK

% Set correct data before plotting

modDataC = bfsk_sig_cyclic;
modDataH = bfsk_sig_hamming;
modDataR = bfsk_sig_rs;
modDataT = bfsk_sig_turbo;

modDataNoiseC = bfsk_rcv_cyclic;
modDataNoiseH = bfsk_rcv_hamming;
modDataNoiseR = bfsk_rcv_rs;
modDataNoiseT = bfsk_rcv_turbo;

decDataC = bfsk_demod_sig_cyclic_dec;
decDataH = bfsk_demod_sig_hamming_dec;
decDataR = bfsk_demod_sig_rs_dec;
decDataT = bfsk_demod_sig_turbo_dec;

% Plot message signal
figure(7);
t_in = 0:1:N;
subplot(5,4,1)
stairs(t_in(1:9),[sig_data_raw(1:8),sig_data_raw(8)]);
axis([0 8 -0.5 1.5])
title("Message Signal");

subplot(5,4,2)
stairs(t_in(1:9),[sig_data_raw(1:8),sig_data_raw(8)]);
axis([0 8 -0.5 1.5])
title("Message Signal");

subplot(5,4,3)
stairs(t_in(1:9),[sig_data_raw(1:8),sig_data_raw(8)]);
axis([0 8 -0.5 1.5])
title("Message Signal");

subplot(5,4,4)
stairs(t_in(1:9),[sig_data_raw(1:8),sig_data_raw(8)]);
axis([0 8 -0.5 1.5])
title("Message Signal");

%----------------------------
% Plot all the encoded signals

subplot(5,4,5)
stairs(t_in(1:9),[sig_data_cyclic(1:8),sig_data_cyclic(8)]);
axis([0 8 -0.5 1.5])
title("Cyclic Encoded Signal");

subplot(5,4,6)
stairs(t_in(1:9),[sig_data_hamming(1:8),sig_data_hamming(8)]);
axis([0 8 -0.5 1.5])
title("Hamming Encoded Signal");

subplot(5,4,7)
stairs(t_in(1:9),[sig_data_rs(1:8),sig_data_rs(8)]);
axis([0 8 -0.5 1.5])
title("Reed Solomon Encoded Signal");

subplot(5,4,8)
stairs(t_in(1:9),[sig_data_turbo(1:8),sig_data_turbo(8)]);
axis([0 8 -0.5 1.5])
title("Turbo Encoded Signal");

%----------------------------
% Plot all the modulated signals

t=0:Ts:(2*Tsym)-Ts; % Duration of Message Signal

subplot(5,4,9)
plot(modDataC(1:length(t)),'r-');
title("Modulated Signal");
xlim([0,length(t)]);

subplot(5,4,10)
plot(modDataH(1:length(t)),'r-');
title("Modulated Signal");
xlim([0,length(t)]);

subplot(5,4,11)
plot(modDataR(1:length(t)),'r-');
title("Modulated Signal");
xlim([0,length(t)]);

subplot(5,4,12)
plot(modDataT(1:length(t)),'r-');
title("Modulated Signal");
xlim([0,length(t)]);

%----------------------------
% Plot all the modulated signals + noise

subplot(5,4,13)
plot(modDataNoiseC(1:length(t)),'m-');
title("Modulated Signal + Noise");
xlim([0,length(t)]);

subplot(5,4,14)
plot(modDataNoiseH(1:length(t)),'m-');
title(" Modulated Signal + Noise");
xlim([0,length(t)]);

subplot(5,4,15)
plot(modDataNoiseR(1:length(t)),'m-');
title("Modulated Signal + Noise");
xlim([0,length(t)]);

subplot(5,4,16)
plot(modDataNoiseT(1:length(t)),'m-');
title("Modulated Signal + Noise");
xlim([0,length(t)]);

%----------------------------
% Plot decoded signal

subplot(5,4,17)
stairs(t_in(1:9),[decDataC(1:8),decDataC(8)]);
axis([0 8 -0.5 1.5])
title("Demodulated and Decoded Signal");

subplot(5,4,18)
stairs(t_in(1:9),[decDataH(1:8),decDataH(8)]);
axis([0 8 -0.5 1.5])
title("Demodulated and Decoded Signal");

subplot(5,4,19)
stairs(t_in(1:9),[decDataR(1:8)',decDataR(8)']);
axis([0 8 -0.5 1.5])
title("Demodulated and Decoded Signal");

subplot(5,4,20)
stairs(t_in(1:9),[decDataT(1:8)',decDataT(8)']);
axis([0 8 -0.5 1.5])
title("Demodulated and Decoded Signal");

sgtitle('Binary Frequency Shift Keying');

%%% End of BFSK Plots

%%% 16-QAM

% Set correct data before plotting

modDataC = qam16_sig_cyclic;
modDataH = qam16_sig_hamming;
modDataR = qam16_sig_rs;
modDataT = qam16_sig_turbo;

modDataNoiseC = qam16_rcv_cyclic;
modDataNoiseH = qam16_rcv_hamming;
modDataNoiseR = qam16_rcv_rs;
modDataNoiseT = qam16_rcv_turbo;

decDataC = qam16_demod_sig_cyclic_dec;
decDataH = qam16_demod_sig_hamming_dec;
decDataR = qam16_demod_sig_rs_dec;
decDataT = qam16_demod_sig_turbo_dec;

% Plot message signal
figure(8);
t_in = 0:1:N;
subplot(5,4,1)
stairs(t_in(1:9),[sig_data_raw(1:8),sig_data_raw(8)]);
axis([0 8 -0.5 1.5])
title("Message Signal");

subplot(5,4,2)
stairs(t_in(1:9),[sig_data_raw(1:8),sig_data_raw(8)]);
axis([0 8 -0.5 1.5])
title("Message Signal");

subplot(5,4,3)
stairs(t_in(1:9),[sig_data_raw(1:8),sig_data_raw(8)]);
axis([0 8 -0.5 1.5])
title("Message Signal");

subplot(5,4,4)
stairs(t_in(1:9),[sig_data_raw(1:8),sig_data_raw(8)]);
axis([0 8 -0.5 1.5])
title("Message Signal");

%----------------------------
% Plot all the encoded signals

subplot(5,4,5)
stairs(t_in(1:9),[sig_data_cyclic(1:8),sig_data_cyclic(8)]);
axis([0 8 -0.5 1.5])
title("Cyclic Encoded Signal");

subplot(5,4,6)
stairs(t_in(1:9),[sig_data_hamming(1:8),sig_data_hamming(8)]);
axis([0 8 -0.5 1.5])
title("Hamming Encoded Signal");

subplot(5,4,7)
stairs(t_in(1:9),[sig_data_rs(1:8),sig_data_rs(8)]);
axis([0 8 -0.5 1.5])
title("Reed Solomon Encoded Signal");

subplot(5,4,8)
stairs(t_in(1:9),[sig_data_turbo(1:8),sig_data_turbo(8)]);
axis([0 8 -0.5 1.5])
title("Turbo Encoded Signal");

%----------------------------
% Plot all the modulated signals

t=0:Ts:(2*Tsym)-Ts; % Duration of Message Signal

subplot(5,4,9)
plot(modDataC(1:length(t)),'r-');
title("Modulated Signal");
xlim([0,length(t)]);

subplot(5,4,10)
plot(modDataH(1:length(t)),'r-');
title("Modulated Signal");
xlim([0,length(t)]);

subplot(5,4,11)
plot(modDataR(1:length(t)),'r-');
title("Modulated Signal");
xlim([0,length(t)]);

subplot(5,4,12)
plot(modDataT(1:length(t)),'r-');
title("Modulated Signal");
xlim([0,length(t)]);

%----------------------------
% Plot all the modulated signals + noise

subplot(5,4,13)
plot(modDataNoiseC(1:length(t)),'m-');
title("Modulated Signal + Noise");
xlim([0,length(t)]);

subplot(5,4,14)
plot(modDataNoiseH(1:length(t)),'m-');
title(" Modulated Signal + Noise");
xlim([0,length(t)]);

subplot(5,4,15)
plot(modDataNoiseR(1:length(t)),'m-');
title("Modulated Signal + Noise");
xlim([0,length(t)]);

subplot(5,4,16)
plot(modDataNoiseT(1:length(t)),'m-');
title("Modulated Signal + Noise");
xlim([0,length(t)]);

%----------------------------
% Plot decoded signal

subplot(5,4,17)
stairs(t_in(1:9),[decDataC(1:8),decDataC(8)]);
axis([0 8 -0.5 1.5])
title("Demodulated and Decoded Signal");

subplot(5,4,18)
stairs(t_in(1:9),[decDataH(1:8),decDataH(8)]);
axis([0 8 -0.5 1.5])
title("Demodulated and Decoded Signal");

subplot(5,4,19)
stairs(t_in(1:9),[decDataR(1:8)',decDataR(8)']);
axis([0 8 -0.5 1.5])
title("Demodulated and Decoded Signal");

subplot(5,4,20)
stairs(t_in(1:9),[decDataT(1:8)',decDataT(8)']);
axis([0 8 -0.5 1.5])
title("Demodulated and Decoded Signal");

sgtitle('16-Quadrature Amplitude Modulation');

%%% End of 16-QAM Plots
