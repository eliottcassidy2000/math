from itertools import product

def p(t):
    if not t:
        return 0
    u,b=divmod(t,2)
    return p(u)+(u+(u&1) if b==0 else 0)

def cut(w):
    N,S,d=1,0,0
    for b in map(int,w):
        N,S,d=N+(1-b)*(S+d),2*S+b,b
    return N,S,d

def literal(w):
    m=len(w); n=m+1; edges=[0]*n
    for i in range(m):
        edges[i]=((1<<m)-1)^((1<<(i+1))-1)
        if w[i]=='0': edges[i]|=1<<m
        else: edges[m]|=1<<i
    D=[[0]*n for _ in range(1<<n)]
    for i in range(n): D[1<<i][i]=1
    for s in range(1,1<<n):
        for i,x in enumerate(D[s]):
            if not x: continue
            cand=edges[i]&~s
            while cand:
                b=cand&-cand; cand-=b
                D[s|b][b.bit_length()-1]+=x
    return sum(D[-1])

for e,cutoff,s,A,t0 in [(7,4,'1110',27,8),(13,5,'10110',57,18),(21,6,'11100',75,24)]:
    for m in range(2,cutoff+1):
        for w in product('01',repeat=m-2):
            w='1'+''.join(w)+'0'
            assert literal(w)==cut(w)[0]!=e
    for r in range(1,8):
        w='11'+'0'*r+s
        assert cut(w)[0]==e+A*2**r
        if r<=4: assert literal(w)==e+A*2**r
    assert p(t0)==e+1 and all(p(t)!=e+1 for t in range(0,t0,2))
    assert A==3*(t0+1) and (A-2)//e<5
print('three literal rays, small forbidden spectra, and sharp normalized odd coefficients passed')

for e in range(1,1024,2):
    b=max(1,e.bit_length())
    c=1+2**b
    js=[j for j in range(1,b+1) if ((c-e)>>j)&1]
    W='1'+'0'*b
    assert cut(W)[0]==c
    for r in (1,2,5,13):
        z=''.join('1'+'0'*j+'1'+'0'*r for j in js)+W
        h,s,_=cut(z)
        assert (h-e)%2**r==0 and (s-int(W,2))%2**r==0
        if len(z)<=11: assert h==literal(z)
print('general compiler: count and full-address congruences passed')
