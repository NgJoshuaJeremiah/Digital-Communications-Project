function out = OOK_Modulate(in,fc,fs,symbol_duration)

    Ts = 1/fs;                      % Sampling period
    t=0:Ts:symbol_duration-Ts;      % Duration of Message Signal
    
    out=[];

    for i = 1:length(in) 
            s1=in(i)*cos(2*pi*fc*t);
            out=[out s1];
    end 
    
end