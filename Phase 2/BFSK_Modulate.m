function out = BFSK_Modulate(in,fs,f1,f2,symbol_duration)

    Ts = 1/fs;                      % Sampling period
    t=0:Ts:symbol_duration-Ts;      % No of samples for 1 duration of one symbol

    c1=cos(2*pi*f1*t);              % Frequency 1 use for bit 1
    c2=cos(2*pi*f2*t);              % Frequency 2 use for bit 0

    out = [];
    
    for i =1:length(in) 
        if in(i)==1
            out = [out c1];
        else
            out = [out c2];
        end
    end

end