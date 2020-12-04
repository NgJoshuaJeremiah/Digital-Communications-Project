function out = Add_Error1(in)

    mask=randsrc(1,length(in),[1 0; 0.2 0.8]);
    out=xor(in,mask);
    
end