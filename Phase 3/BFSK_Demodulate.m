function out = BFSK_Demodulate(in,fs,f1,f2,symbol_duration,output_type) 

    Ts = 1/fs;                      % Sampling period
    t=0:Ts:symbol_duration-Ts;      % Duration of Message Signal

    c1=cos(2*pi*f1*t);              % Frequency 1 use for bit 1 (Phase locked to transmitter)
    c2=cos(2*pi*f2*t);              % Frequency 2 use for bit 0 (Phase locked to transmitter)

    y1=[];
    y2=[];
    
    
    for i=1:length(t):length(in)
        y1=[y1 c1.*in(1,i:i+(length(t)-1))]; 
        y2=[y2 c2.*in(1,i:i+(length(t)-1))];
    end
    
    
    
    
    int_op1 = [];
    int_op2 = [];
    
    
    S = (length(t)-1);
    
    %First correlator
    for i=1:length(t):length(in)
        int_o=(1/S)*trapz(y1(1,i : i+S));   % Integrate over symbol period
        int_op1=[int_op1 int_o];
    end

    %Second correlator
    for i=1:length(t):length(in)
        int_o=(1/S)*trapz(y2(1,i : i+S));   % Integrate over symbol period
        int_op2=[int_op2 int_o];
    end

    out = [];
    out = int_op1-int_op2;

    %Decision
    threshold = 0;
    if strcmp(output_type,'soft')  % Instead of making a hard decision send out soft bits
    else
        out(out>=threshold)=1; 
        out(out<threshold)=0;        
    end

end