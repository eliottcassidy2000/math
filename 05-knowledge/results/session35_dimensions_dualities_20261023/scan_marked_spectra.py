import numpy as np,json,time
from collections import defaultdict

def ham(out,n):
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for v in range(n):
            a=dp[mask][v]
            if not a:continue
            nxt=out[v]&~mask
            while nxt:
                b=nxt&-nxt; nxt-=b; w=b.bit_length()-1
                dp[mask|b][w]+=a
    return sum(dp[-1])

def run(n):
    free=[(i,j) for i in range(n) for j in range(i+2,n)]
    seen={}; groups=defaultdict(set); witness=None; start=time.time()
    for mask in range(1<<len(free)):
        M=np.zeros((n,n),dtype=np.int64)
        for i in range(n-1):M[i,i+1]=1;M[i+1,i]=-1
        for b,(i,j) in enumerate(free):M[i,j]=1 if mask>>b&1 else -1;M[j,i]=-M[i,j]
        power=np.eye(n,dtype=np.int64); sig=[]
        for k in range(1,n+1):
            power=power@M;sig.extend((int(np.trace(power)),int(power.sum())))
        key=tuple(sig)
        out=[sum(1<<j for j in range(n) if M[i,j]>0) for i in range(n)]
        h=ham(out,n);groups[key].add(h)
        if key in seen and seen[key][1]!=h and witness is None:
            witness={'first_mask':seen[key][0],'second_mask':mask,'H':[seen[key][1],h],'signature':sig,'out_second':out}
        else:seen.setdefault(key,(mask,h))
    if witness:
        mask=witness['first_mask']; M=np.zeros((n,n),dtype=np.int64)
        for i in range(n-1):M[i,i+1]=1;M[i+1,i]=-1
        for b,(i,j) in enumerate(free):M[i,j]=1 if mask>>b&1 else -1;M[j,i]=-M[i,j]
        witness['out_first']=[sum(1<<j for j in range(n) if M[i,j]>0) for i in range(n)]
    return {'n':n,'presentations':1<<len(free),'marked_spectral_classes':len(groups),'H_ambiguous_classes':sum(len(v)>1 for v in groups.values()),'witness':witness}
print(json.dumps([run(n) for n in range(3,8)],indent=2))
