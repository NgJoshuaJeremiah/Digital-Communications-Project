function [out, demod_sig] = QAM16_Demodulate(in,type,fc,fs,symbol_duration,output_type)

    M = 16;                         % Modulation order (alphabet size or number of points in signal constellation)
    Ts = 1/fs; 
    t=0:Ts:symbol_duration-Ts;      % Duration of Message Signal
    
    %%%% Coherant Demodulation %%%%
    I1=[];  % Array to store in-phase component
    Q1=[];  % Array to store quadrature component

    r_cos_carrier = 2*cos(2*pi*fc*t);
    r_sin_carrier = 2*sin(2*pi*fc*t);
    
    for i=1:length(t):length(in);
        I1=[I1 r_cos_carrier.*in(1,i:i+(length(t)-1))]; % for I-channel
        Q1=[Q1 r_sin_carrier.*in(1,i:i+(length(t)-1))]; % for Q-channel
    end
    
    %%%% Low-Pass Filtering %%%%

    wn=fc/(fs/2);                 % Cut-off Frequency is 2Khz
    [b1,a1]=butter(6,wn,'low');   % ButterWorth filter of Order 4
    
    FI1=filter(b1,a1,I1);          
    FQ1=filter(b1,a1,Q1);
    
    %%%% Thresholding and Phase Detection %%%%
    Iv = [];
    Qv = [];
    
    
    clear i;
    for i=1:length(t):length(in);           % Calculating Average values for I-channel 
            FI2=FI1(1,i:i+(length(t)-1));
            Iv=[Iv mean(FI2)];        
    end
    
    clear i;
    for i=1:length(t):length(in);           % Calculating Average values for Q-channel
            FQ2=FQ1(1,i:i+(length(t)-1));
            Qv=[Qv mean(FQ2)];       
    end
    
    
    demod_sig = complex(Iv,Qv);
    out = qamdemod(demod_sig(:),M,type,'OutputType',output_type);
    out = out(:)';

end