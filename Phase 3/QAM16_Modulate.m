function [out, complex_env] = QAM16_Modulate(in,type,fc,fs,symbol_duration)
    

    M = 16;                         % Modulation order (alphabet size or number of points in signal constellation)
    k = log2(M);
    Ts = 1/fs;                      % Sampling period
    t=0:Ts:symbol_duration-Ts;      % Duration of Message Signal
    
    
    sig_data = reshape(in,length(in)/k,k);
    complex_env = qammod(sig_data(:),M,type,'InputType','bit');
    
    
    I = real(complex_env);          % Get in-phase component
    Q = imag(complex_env);          % Get quadrature component
    
    out=[];

    for i = 1:length(complex_env) 
            s1=I(i)*cos(2*pi*fc*t);
            s2=Q(i)*sin(2*pi*fc*t);
            out=[out (s1+s2)];
    end 

end
