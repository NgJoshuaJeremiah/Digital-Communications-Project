function encData_hamming = Encode_Hamming(s)
    n=7; %hamming: 2^k-1
    k=4; %hamming: 2^k-1-k
    
    encData_hamming = encode(s,n,k,'hamming/binary'); %Encode signal
    
end