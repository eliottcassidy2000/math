"""Independent subset-path oracle for the exact tape readout."""
from fractions import Fraction as F
from itertools import permutations, product
from random import Random
from math import prod

def tournament(n, code):
    a=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate((i,j) for i in range(n) for j in range(i+1,n)):
        a[i][j]=code>>k&1; a[j][i]=1-a[i][j]
    return a

def path_dp(a):
    n=len(a); dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for s in range(1,1<<n):
        for u,w in enumerate(dp[s]):
            if w:
                for v in range(n):
                    if not s>>v&1 and a[u][v]: dp[s|1<<v][v]+=w
    return sum(dp[-1])

def block(a,V):
    n=len(a); f=(1<<n)-1
    starts=[[0]*n for _ in range(1<<n)]; ends=[[0]*n for _ in range(1<<n)]
    for k in range(1,n+1):
        for w in permutations(range(n),k):
            if all(a[u][v] for u,v in zip(w,w[1:])):
                mask=sum(1<<v for v in w)
                starts[mask][w[0]]+=1; ends[mask][w[-1]]+=1
    total=[1]+[sum(row) for row in starts[1:]]
    sv=[sum(row[v] for v in V) for row in starts]
    ev=[sum(row[v] for v in V) for row in ends]
    return total[f], sum(total[s]*total[f^s] for s in range(f+1)),sum(ev[s]*total[f^s] for s in range(f+1)),sum(ev[s]*sv[f^s] for s in range(f+1)),sum(total[s]*sv[f^s] for s in range(f+1)),ev[f],sv[f]

def lift(tape,marks):
    labels=[(i,v) for i,a in enumerate(tape) for v in range(len(a))]
    n=len(labels); a=[[0]*(n+2) for _ in range(n+2)]
    for u,(i,v) in enumerate(labels):
        for w,(j,z) in enumerate(labels):
            a[u][w]=int(i<j) if i!=j else tape[i][v][z]
        a[u][n]=1
        a[u][n+1]=int(v in marks[i]); a[n+1][u]=1-a[u][n+1]
    a[n][n+1]=1
    return a

def predict(tape,marks):
    hs,qs,es,js,ls,ts,ss=zip(*(block(a,V) for a,V in zip(tape,marks)))
    r=[F(q,h) for h,q in zip(hs,qs)]
    B=1; correction=0
    for i,h in enumerate(hs):
        correction += F(es[i]*B-js[i]-ls[i]*prod(r[i+1:]),h)
        if i+1<len(hs): correction-=F(ts[i]*ss[i+1],h*hs[i+1])
        B=r[i]*B-F(ls[i],h)
    return prod(hs)+prod(qs)+prod(hs)*correction

def main():
    rng=Random(1881); g=tournament(1,0)
    cyclic=next(tournament(3,c) for c in range(8) if path_dp(tournament(3,c))==3)
    best=0; trials=0
    for d in range(1,5):
        tape=[g]+[cyclic]*d+[g]
        for word in product(range(4), repeat=d):
            marks=[set()]+[set(range(s)) for s in word]+[set()]
            if d<=3 or rng.random()<.12:
                v=path_dp(lift(tape,marks))
                assert v==predict(tape,marks) and v>=2*12**d
                if d>=2:
                    div=3**(d-2)
                    assert v%div==0 and (v//div+sum(a*b for a,b in zip(word,word[1:])))%3==0
                trials+=1
    for _ in range(60):
        tape=[g]+[tournament(n,rng.randrange(1<<(n*(n-1)//2))) for n in [rng.randrange(1,4) for _ in range(3)]]+[g]
        marks=[set()]+[{v for v in range(len(a)) if rng.randrange(2)} for a in tape[1:-1]]+[set()]
        assert path_dp(lift(tape,marks))==predict(tape,marks); trials+=1
        bands=[]
        for i in range(1,4):
            values=[]
            for s in range(1,1<<len(tape[i])):
                masks=[set() for _ in tape]; masks[i]={v for v in range(len(tape[i])) if s>>v&1}
                values.append(predict(tape,masks))
            bands.append((min(values),max(values)))
        assert all(bands[i][1]<bands[i+1][0] for i in range(2))
    edge=tournament(2,1)
    for l, delta in [(1,6),(2,-4),(3,0)]:
        before=lift([g,cyclic,g,edge,g],[set(),set(range(l)),set(),set(),set()])
        V={1} if l==1 else {0} if l==2 else {0,1}
        after=lift([g,edge,g,cyclic,g],[set(),V,set(),set(),set()])
        assert path_dp(after)-path_dp(before)==delta
    print(f'{trials} independent subset-path comparisons; 120 adjacent-band tests')
if __name__=='__main__': main()
