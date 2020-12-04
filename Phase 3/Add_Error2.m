function out = Add_Error2(in)

    mask=randsrc(1,length(in),[-10 -20 10 20 0; 0.05 0.05 0.05 0.05 0.80]);
    out=mask+in;
    
end