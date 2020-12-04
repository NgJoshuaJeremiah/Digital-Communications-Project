
%   _  _   _   _   _   _     _   _ 
%  /  |_   _) / \ / \ |_    | \ /  
%  \_ |_   _) \_/ \_/ |_)   |_/ \_ 
%                                  
%   Phase 1

close all; 
clc;

N = 1024; % No of bits

max_snr = 50;
snr_inc = 5;
start = 0;

threshold = 0;
signal_p = 1; % Signal power assume 1 unit

sig_data_raw = randi([0, 1], [1, N]); % Generate binary digits and store in array 
sig_data = sig_data_raw.*2-1; % convert to +- 1


i = 1;


err_arr = start:snr_inc:max_snr;

% snr is in dB
for snr = start:snr_inc:max_snr
    
  noise_v = signal_p /(10^(snr/10)); 
  noise = sqrt(noise_v/2).*randn(1,N); %Assume two sided white noise
  rcv_sig = sig_data+noise;
  
  % Convert rcv signal based on threshold
  rcv_sig(rcv_sig>=threshold)=1; 
  rcv_sig(rcv_sig<threshold)=0;

  
  %Calculate number of bit errors
  
  err = 0;
  
  for x = 1:1:N
      if rcv_sig(x) ~= sig_data_raw(x) 
        err=err+1;
      end 
  end
  
  % Save the error values for usage later
  err_arr(i) = (err/N)*100;
  i=i+1;
  
  fprintf("SNR: %f     Bit Error Rate: %f\n",snr,(err/N));
  
end


% For plotting generate the SNR values
snr = start:snr_inc:max_snr;

%Graph and Plot the result           
figure(1)
semilogx(snr,err_arr,'bo-');
title('SNR to Bit Error Rate');
ylabel('Bit Error Rate (%)');
xlabel('SNR[dB]');
grid on;
hold on;
