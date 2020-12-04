function decData_hamming = Decode_Hamming(s)
    n=7; %hamming: 2^k-1
    k=4; %hamming: 2^k-1-k
    decData_hamming = decode(s,n,k,'hamming/binary');
    
end