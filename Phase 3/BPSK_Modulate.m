function out = BPSK_Modulate(in,fc,fs,symbol_duration)

    in = in.*2-1;                     % BPSK is basically ASK with in-phase component from of -1 and 1
    
    Ts = 1/fs;                      % Sampling period
    t=0:Ts:symbol_duration-Ts;      % Duration of Message Signal
    
    out=[];

    for i = 1:length(in) 
            s1=in(i)*cos(2*pi*fc*t);
            out=[out s1];
    end 
    
end



