function decData_cyclic = Decode_Cyclic(s)
    n=12; %hamming: 2^k-1
    k=4; %hamming: 2^k-1-k
    decData_cyclic = decode(s,n,k,'cyclic/binary');
    
end