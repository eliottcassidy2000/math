"""Independent unordered-set-partition profiles of supplied literal witnesses."""
import json
from pathlib import Path
from functools import lru_cache
from tournament_kernel import *

def profile(a):
    n=len(a);N=1<<n;dp=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
    for i in range(n):dp[1<<i][i]=1
    for s in range(1,N):
        H[s]=sum(dp[s])
        for j in range(n):
            if dp[s][j]:
                for k in range(n):
                    if not(s>>k&1) and a[j][k]:dp[s|1<<k][k]+=dp[s][j]
    @lru_cache(None)
    def covers(s):
        if not s:return (1,)
        least=s&-s;other=s^least;ans=[0]*(s.bit_count()+1);t=other
        while True:
            block=t|least;rest=s^block
            for k,value in enumerate(covers(rest)):ans[k+1]+=H[block]*value
            if t==0:break
            t=(t-1)&other
        return tuple(ans)
    return covers(N-1)[1:]

def signature(p,M):
    precisions={3:[3,3],4:[4,2],8:[8,4,2],12:[12,6]}[M]
    return tuple((p[i] if i<len(p) else 0)%m for i,m in enumerate(precisions))

def run():
    payload=json.loads(Path(__file__).with_name('width_witnesses.json').read_text());witnesses=0
    classes={}
    for key,rows in payload['states'].items():
        M=int(key);states=set()
        for r in rows:
            masks=r['adjacency'];n=len(masks);a=tuple(tuple(masks[i]>>j&1 for j in range(n)) for i in range(n))
            p=profile(a);sig=signature(p,M)
            require(sig==tuple(r['state']),'literal state witness')
            require(p[0]==r['H'],'literal H witness');states.add(sig);witnesses+=1
        classes[key]=sorted(states)
    require({int(k):len(v) for k,v in classes.items()}=={3:5,4:3,8:14,12:11},'state counts')
    checked=0
    for n,bank in class_bank(7):
        if not n:continue
        for code in bank:
            p=profile(matrix(n,code));H=p[0];p2=p[1] if n>=2 else 0;p3=p[2] if n>=3 else 0
            require((p2%3)==(0 if n%2 else H%3),'reversal residue')
            require(p2%2==int(n>=2),'p2 parity')
            require(p3%2==int(n>=3 and n%2==1),'p3 parity')
            if n>=4 and n%2==0:require(p2%4==(H+2)%4,'even OCF residue')
            for M in (3,4,8,12):require(signature(p,M) in classes[str(M)],'upper-state classification')
            checked+=1
    G={(1,0)}
    while True:
        old=set(G)
        for h,p in old:G.add(((3*h+2*p)%12,(3*h+8*p)%6))
        for h,p in old:
            for g,q in old:G.add((h*g%12,(h*g+2*h*q+2*p*g+2*p*q)%6))
        if G==old:break
    require(G=={(1,0),(1,1),(1,3),(3,3),(5,3),(5,5),(9,3)},'restricted grammar closure')
    print(json.dumps({'scope':'FINITE-EXACT independently evaluated literal witnesses and all classes n<=7; universal upper identities reviewed mathematically',
                      'source_witness_sha256':payload['source_sha256'],'literal_witnesses':witnesses,
                      'all_class_profile_checks':checked,'exact_context_classes':classes,
                      'restricted_G_mod12_states':sorted(G)},sort_keys=True,indent=2))
if __name__=='__main__':run()
