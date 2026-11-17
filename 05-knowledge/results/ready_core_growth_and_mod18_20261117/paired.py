from replay import *
from collections import defaultdict

def near(g, Cs, cap):
    best=defaultdict(dict)
    for C in Cs:
        for l in range(cap+1):
            for r in range(cap-l+1):
                D=rows_ord([(0,)]*l+list(C)+[(0,)]*r)
                G,K,V=arrays(D)
                assert G==g
                n=len(D);N=1<<n;rows=[sum(K[i*n:i*n+n]) for i in range(n)]
                for U in range(N):
                    ids=[i for i in range(n) if U>>i&1]
                    k=sum(K[i*n+j] for i in ids for j in ids)
                    if abs(k-g)>1:continue
                    j=sum(rows[i] for i in ids)-k
                    E=sum(V[i] for i in ids)
                    O=g+j+k
                    if O-E!=g or E==0:continue
                    v=2*E;h=g+2*j
                    c=2*(g-k)
                    t=(v&-v).bit_length()-1;b=v>>t
                    z=best[b].get(c)
                    if z is None or t<z[0]:best[b][c]=(t,h,D,U)
    for a,d in sorted(best.items()):
        if 0 in d:
            print(a,{c:(t,h,len(D)+1) for c,(t,h,D,U) in d.items()},flush=True)
            for c,(t,h,D,U) in d.items():
                if c:print('masks',c,D,U)
    return best
if __name__=='__main__':
    k=near(5,[((2,4,9,3),)],6)
    assert {a:sorted(d) for a,d in k.items() if 0 in d}=={3:[0],5:[-2,0],9:[0],15:[0],19:[0,2],21:[0],23:[0],25:[0],71:[0]}
    assert k[19][2][0]==1 and k[19][0][0]==2
    banks=[(62,36,104,114,102,64,3),(242,244,185,243,160,192,20,64)]
    def recode(A,r):
        js=[i for i in range(len(A)) if i!=r]
        D=tuple(sum((A[i]>>j&1)<<p for p,j in enumerate(js)) for i in js)
        U=sum((A[r]>>i&1)<<p for p,i in enumerate(js))
        return literal(D,U)
    def add(A,r,b):
        n=len(A)
        return tuple(x+((i!=r or not b)<<n) for i,x in enumerate(A))+((1<<r) if b else 0,)
    for b,h in zip(banks,(71,73)):
        assert recode(b,6)==(h,5,38,43)
    aa,bb=banks
    for _ in range(4):aa=add(aa,6,1);bb=add(bb,6,1)
    assert recode(aa,6)==(1211,5,608,613)
    assert recode(bb,6)==(1213,5,608,613)
    assert 1211==7*173
    assert all(v%d for v in (173,1213) for d in range(2,int(v**.5)+1))
    print("literal 1211/1213 switch verified with identical continuation")
