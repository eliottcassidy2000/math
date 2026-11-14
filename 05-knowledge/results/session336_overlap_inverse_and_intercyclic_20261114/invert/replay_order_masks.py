from itertools import combinations
def read(w):
    m=len(w); n=m+2; x,y=m,m+1
    o=[sum(1<<j for j in range(i+1,m)) for i in range(m)]+[1<<y,0]
    for i,d in enumerate(map(int,w)):
        for h in (x,y):
            if d>>(y-h)&1:o[h]|=1<<i
            else:o[i]|=1<<h
    D=[[0]*n for _ in range(1<<n)]
    for v in range(n):D[1<<v][v]=1
    for mask in range(1,1<<n):
       for v,k in enumerate(D[mask]):
        if not k:continue
        r=o[v]&~mask
        while r:
          l=r&-r;r^=l;D[mask|l][l.bit_length()-1]+=k
    return o,sum(D[-1])
for w,N,(i,j,k,l) in [
    ('1030221',385,(0,2,5,6)),('20302212',573,(0,2,6,7))]:
    o,h=read(w);x=len(w);y=x+1
    assert h==N
    def arc(a,b):return bool(o[a]&1<<b)
    assert any(arc(x,b) and arc(a,x) and arc(y,d) and arc(c,y)
        for b,d in [(i,j),(j,i)] for a,c in [(k,l),(l,k)])
    print(N,w,'alpha',j-i-1+l-k-1,'beta',k-j-1)
