function out = BPSK_Demodulate(in,fc,fs,symbol_duration,output_type)

    Ts = 1/fs;                          % Sampling period
    t=0:Ts:symbol_duration-Ts;          % No of samples for 1 duration of one symboll
    
    r_cos_carrier = 2*cos(2*pi*fc*t);   % Receiver oscillator (Phase locked to transmitter)
    
    I1 = [];
    
    for i=1:length(t):length(in);
        I1=[I1 r_cos_carrier.*in(1,i:i+(length(t)-1))]; % for I-channel
    end
    
    %LPF
    [b,a] = butter(6,0.2,'low');        % Butterworth filter of 6th order with cut off freq 1.6Khz
    FI1 = filtfilt(b,a,I1);

    %%%% Thresholding and Phase Detection %%%%

    out = [];
    
    clear i;
    for i=1:length(t):length(in);       % Calculating Average values for I-channel 
            FI2=FI1(1,i:i+(length(t)-1));
            out=[out mean(FI2)];        
    end
    
    
    %Decision
    threshold = 0;
    
    if strcmp(output_type,'soft')  % Instead of making a hard decision send out soft bits
    else
        out(out>=threshold)=1; 
        out(out<threshold)=0;        
    end
   
end