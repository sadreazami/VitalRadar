%----------------- Function for creating the surrogates -------------------
function surr = createsurr(sig,NS)

L=length(sig); surr=zeros(NS,L); sig=sig(:)';


ll=ceil((L+1)/2)-1; ftsig=fft(sig,L);
for sn=1:NS
    surr(sn,1)=ftsig(1); randph=2*pi*rand(1,ll-1);
    surr(sn,2:ll)=ftsig(2:ll).*exp(1i*randph);
    surr(sn,2+L-ll:L)=conj(fliplr(surr(sn,2:ll)));
    
    surr(sn,:)=real(ifft(surr(sn,:),L));
end

end
