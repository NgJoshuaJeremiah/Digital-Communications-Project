%   _  _   _   _   _   _     _   _ 
%  /  |_   _) / \ / \ |_    | \ /  
%  \_ |_   _) \_/ \_/ |_)   |_/ \_ 
%                                  
%   Phase 2 Different Modulations

clear all;


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
N = 3096; % Number of bits to transmit
M = 16; % Modulation order (alphabet size or number of points in signal constellation)
k = log2(M); % Number of bits per symbol
data_rate = 1000; % Message data rate = 1000bps
Tb = 1/data_rate; % bit period
Tsym = k*Tb; % symbol period
Ts = 1/fs; % Sampling period
signal_p = 1; % Signal power

% % % % % % % % % Generate Message Data % % % % % % % % %

sig_data_raw = randi([0, 1], 1, N);

% % % % % % % % % Modulate Message Data % % % % % % % % %

ook_sig = OOK_Modulate(sig_data_raw,fc,fs,Tb);
bpsk_sig = BPSK_Modulate(sig_data_raw,fc,fs,Tb);
bfsk_sig = BFSK_Modulate(sig_data_raw,fs,2*fc,fc,Tb); % 1 = 20Khz carrier , 0 = 10Khz carrier
[qam16_sig, complex_env] = QAM16_Modulate(sig_data_raw,'gray',fc,fs,Tsym);

% % % % % % % % % SNR Loop % % % % % % % % %

max_snr = 0;
snr_inc = 1;
start = -100;

% Use to index error array
index = 1; 

% Store bit error rate for each iteration of the loop
qam16_err_arr = zeros((max_snr*2)/snr_inc,0);
ook_err_arr = zeros((max_snr*2)/snr_inc,0);
bpsk_err_arr = zeros((max_snr*2)/snr_inc,0);
bfsk_err_arr = zeros((max_snr*2)/snr_inc,0);

% Main Process Loop
for snr = start:snr_inc:max_snr
    
    % Transmitted Signal passes through AWGN Channel
    qam16_rcv = awgn(qam16_sig,snr,'measured');
    ook_rcv = awgn(ook_sig,snr,'measured');
    bpsk_rcv = awgn(bpsk_sig,snr,'measured');
    bfsk_rcv = awgn(bfsk_sig,snr,'measured');
    
    % Demodulate Signal
    [qam16_demod_sig, complex_demod_sig] = QAM16_Demodulate(qam16_rcv,'gray',fc,fs,Tsym,'bit');
    ook_demod_sig = OOK_Demodulate(ook_rcv,fc,fs,Tb);
    bpsk_demod_sig = BPSK_Demodulate(bpsk_rcv,fc,fs,Tb,'bit');
    bfsk_demod_sig = BFSK_Demodulate(bfsk_rcv,fs,2*fc,fc,Tb);   
    
    % Check bit error rate and save the values for plotting
    [~, ber] = biterr(sig_data_raw,ook_demod_sig);
    ook_err_arr(index) = ber*100; 
    [~, ber] = biterr(sig_data_raw,bpsk_demod_sig);
    bpsk_err_arr(index) = ber*100;
    [~, ber] = biterr(sig_data_raw,qam16_demod_sig);
    qam16_err_arr(index) = ber*100;  
    [~, ber] = biterr(sig_data_raw,bfsk_demod_sig);
    bfsk_err_arr(index) = ber*100; 
    
    index = index + 1;
    
    
    
end

% For plotting generate the SNR values
snr = start:snr_inc:max_snr;

%Graph and Plot the result           
figure(1)
semilogx(snr,qam16_err_arr,'b-');
hold on;
semilogx(snr,ook_err_arr,'r-');
hold on;
semilogx(snr,bpsk_err_arr,'g-');
hold on;
semilogx(snr,bfsk_err_arr,'m-');
title('SNR to Bit Error Rate');
ylabel('Bit Error Rate (%)');
xlabel('SNR[dB]');
grid on;
legend('16-QAM','OOK','BPSK','BFSK');



%%% Phase 2 plots 

%%% OOK


% Plot message signal
figure(2);
t_in = 0:1:N;
subplot(5,1,1)
stairs(t_in(1:17),[sig_data_raw(1:16),sig_data_raw(16)]);
title("Message Signal (Showing 16 bits or 4 symbols)");
axis([0 16 -0.5 1.5])

% Plot modulated signal

t=0:Ts:(4*Tsym)-Ts; % Duration of Message Signal

subplot(4,1,2)
plot(ook_sig(1:length(t)),'r-');
title("OOK Modulated Signal");
xlim([0,length(t)]);

subplot(4,1,3)
plot(ook_rcv(1:length(t)),'b-');
title("OOK Modulated Signal with Noise");
xlim([0,length(t)]);


% Plot demodulated signal

subplot(4,1,4)
stairs(t_in,[ook_demod_sig,ook_demod_sig(end)]);
title("OOK Demodulated signal");
axis([0 16 -0.5 1.5])


%%% BPSK


% Plot message signal
figure(3);
t_in = 0:1:N;
subplot(5,1,1)
stairs(t_in(1:17),[sig_data_raw(1:16),sig_data_raw(16)]);
title("Message Signal (Showing 16 bits or 4 symbols)");
axis([0 16 -0.5 1.5])

% Plot modulated signal

t=0:Ts:(4*Tsym)-Ts; % Duration of Message Signal

subplot(4,1,2)
plot(bpsk_sig(1:length(t)),'r-');
title("BPSK Modulated Signal");
xlim([0,length(t)]);

subplot(4,1,3)
plot(bpsk_rcv(1:length(t)),'b-');
title("BPSK Modulated Signal with Noise");
xlim([0,length(t)]);


% Plot demodulated signal

subplot(4,1,4)
stairs(t_in,[bpsk_demod_sig,bpsk_demod_sig(end)]);
title("BPSK Demodulated signal");
axis([0 16 -0.5 1.5])





%%% 16-QAM 


% Plot message signal
figure(4);
t_in = 0:1:N;
subplot(4,1,1)
stairs(t_in(1:17),[sig_data_raw(1:16),sig_data_raw(16)]);
title("Message Signal (Showing 16 bits or 4 symbols)");
axis([0 16 -0.5 1.5])

clear i;
for i=4:4:16
  xline(i,'--',{'Symbol','Divider'});
end

% Plot modulated signal

t=0:Ts:(4*Tsym)-Ts; % Duration of Message Signal

subplot(4,1,2)
plot(qam16_sig(1:length(t)),'r-');
title("16-QAM Modulated Signal");
xlim([0,length(t)]);

clear i;
for i=640:640:2560
  xline(i,'--',{'Symbol','Divider'});
end

subplot(4,1,3)
plot(qam16_rcv(1:length(t)),'b-');
title("16-QAM Modulated Signal with Noise");
xlim([0,length(t)]);

clear i;
for i=640:640:2560
  xline(i,'--',{'Symbol','Divider'});
end

% Plot demodulated signal

subplot(4,1,4)
stairs(t_in,[qam16_demod_sig,qam16_demod_sig(end)]);
title("16-QAM Demodulated signal");
axis([0 16 -0.5 1.5])

clear i;
for i=4:4:16
  xline(i,'--',{'Symbol','Divider'});
end

%%% BFSK


% Plot message signal
figure(5);
t_in = 0:1:N;
subplot(4,1,1)
stairs(t_in(1:17),[sig_data_raw(1:16),sig_data_raw(16)]);
title("Message Signal (Showing 16 bits or 4 symbols)");
axis([0 16 -0.5 1.5])

% Plot modulated signal

t=0:Ts:(4*Tsym)-Ts; % Duration of Message Signal

subplot(4,1,2)
plot(bfsk_sig(1:length(t)),'r-');
title("BFSK Modulated Signal");
xlim([0,length(t)]);

subplot(4,1,3)
plot(bfsk_rcv(1:length(t)),'b-');
title("BFSK Modulated Signal with Noise");
xlim([0,length(t)]);

% Plot demodulated signal

subplot(4,1,4)
stairs(t_in,[bfsk_demod_sig,bfsk_demod_sig(end)]);
title("BFSK Demodulated signal");
axis([0 16 -0.5 1.5])



sPlotFig = scatterplot(complex_demod_sig,1,0,'g.');
hold on
scatterplot(complex_env,1,0,'k*',sPlotFig)



