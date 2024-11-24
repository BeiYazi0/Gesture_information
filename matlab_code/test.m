Range_Number=3;
Doppler_Number=1;
Tx_Number=2;
Rx_Number=2;
a=[5,5;0,0;0,0];
range_profile(:,:,1)=a;
range_profile(:,:,2)=a;
fft1d_sub =zeros(Range_Number,Doppler_Number*Tx_Number,Rx_Number);
for n=1:Rx_Number
    avg = sum(range_profile(:,:,n),2)/Doppler_Number/Tx_Number;
    for m=1:Doppler_Number*Tx_Number
        fft1d_sub(:,m,n) = range_profile(:,m,n)-avg;
    end
end
range_profile =fft1d_sub;