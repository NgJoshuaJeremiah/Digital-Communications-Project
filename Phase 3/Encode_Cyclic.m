function encData_cyclic = Encode_Cyclic(s)
    n=12; %hamming: 2^k-1
    k=4; %hamming: 2^k-1-k
    
    encData_cyclic = encode(s,n,k,'cyclic/binary'); %Encode signal
    
end