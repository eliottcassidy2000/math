# Longest repeated window and maximum-order complexity (shortest feedback shift register with ANY feedback function)
import sys, numpy as np
raw=np.frombuffer(open(sys.argv[1],'rb').read(),dtype=np.uint8); N=len(raw)
def windows(k):
    w=np.zeros(N-k+1,dtype=np.uint64)
    for i in range(k): w=(w<<np.uint64(1))|raw[i:N-k+1+i].astype(np.uint64)
    return w
def has_repeat(k):
    w=windows(k); u=np.unique(w); return len(u)<len(w)
def has_ambiguous(k):  # some k-window occurs twice with different next bits
    w=windows(k)[:-1]; nxt=raw[k:].astype(np.uint64)
    keys=(w<<np.uint64(1))|nxt
    u_w=np.unique(w); u_k=np.unique(keys)
    return len(u_k)>len(u_w)
lo,hi=1,63
while lo<hi:
    mid=(lo+hi+1)//2
    if has_repeat(mid): lo=mid
    else: hi=mid-1
lam=lo
lo,hi=1,63
while lo<hi:
    mid=(lo+hi+1)//2
    if has_ambiguous(mid): lo=mid
    else: hi=mid-1
print(f"N={N}: longest repeated window = {lam} bits; longest ambiguous window = {lo} bits; maximum-order complexity = {lo+1}; 2*log2(N) = {2*np.log2(N):.1f}")
