"""Exact CRT witness for one fixed physical gap-two rail."""
from math import gcd

T=(1<<60)*((1<<61)-1)-1
A=522665
C=-4643

def H(n):
    v=A*4**n+174*n+C
    assert v%9==0
    return v//9

def Hmod(n,m):
    v=(A*pow(4,n,9*m)+174*n+C)%(9*m)
    assert v%9==0
    return v//9

def crt(r,m,s,n):
    g=gcd(m,n)
    assert (s-r)%g==0, (r,m,s,n)
    q=n//g
    t=((s-r)//g*pow(m//g,-1,q))%q if q>1 else 0
    return (r+m*t)%(m*q),m*q

def primes(p):
    d=2
    while d*d<=p:
        if p%d==0:
            yield d
            while p%d==0:p//=d
        d+=1
    if p>1:yield p

def order4(p):
    q=p-1
    for d in primes(q):
        while q%d==0 and pow(4,q//d,p)==1:q//=d
    return q

def lift(n,step,p,mod):
    winners=[n+a*step for a in range(p) if Hmod(n+a*step,mod)==T%mod]
    assert len(winners)==1,(n,step,p,mod,winners)
    return winners[0]

def witness(prec, free=()):
    # Extra base precision for the coherence of the exceptional torsion.
    s2=max(2,prec.get(2,0));s3=max(1,prec.get(3,0))
    s7=max(1,prec.get(7,0));s29=max(1,prec.get(29,0))
    L=2**(s2-1)
    r=((9*T-C)//2*pow(87,-1,L))%L
    n=15
    for j in range(1,s3):n=lift(n,3**j,3,3**(j+1))
    r,L=crt(r,L,n,3**s3)
    n=15
    for j in range(s7):n=lift(n,3*7**j,7,7**(j+1))
    r,L=crt(r,L,n,7**s7)
    assert r%14==1
    n=15
    for j in range(s29-1):n=lift(n,14*29**j,29,29**(j+2))
    r,L=crt(r,L,n,29**(s29-1))
    for free_r,free_mod in free:r,L=crt(r,L,free_r,free_mod)
    for p,s in sorted(prec.items()):
        if p in (2,3,7,29):continue
        q=order4(p);n=r%q
        assert L%p!=0
        for j in range(s):n=lift(n,q*p**j,p,p**(j+1))
        r,L=crt(r,L,n,q*p**s)
    n=r
    if n<max(1,s2):n+=((max(1,s2)-n+L-1)//L)*L
    for p,s in prec.items():assert Hmod(n,p**s)==T%(p**s)
    return n,L

def p(y):
    z=0
    while y:
        y,b=divmod(y,2)
        if b==0:z+=y+(y&1)
    return z

assert p(73)==32
for n in range(1,13):
    y=int(bin(73)[2:]+"10"*n,2)
    assert H(n)==493+769*y+29*p(y)
assert H(52)<T<H(53)
assert [(H(15)%q,T%q) for q in (3,7,29)]==[(0,0),(0,0),(2,2)]
assert ((H(18)-H(15))//3)%3==2 and ((H(29)-H(15))//29)%29==2
Ogg=(2,3,5,7,11,13,17,19,23,29,31,41,47,59,71)
for name,prec in [
    ("Ogg",dict.fromkeys(Ogg,1)),
    ("Ogg squared",dict.fromkeys(Ogg,2)),
    ("mixed high",{2:60,3:20,5:17,7:10,29:15,59:5,211:7,337:3,593:2})
]:
    n,L=witness(prec, ((15,29),) if name=="Ogg" else ())
    print(name,"blocks",n,"arithmetic period",L,"vertices",16+2*n)
print("Exact non-hit bracket:",H(52),T,H(53))

# Optional independent literal subset-DP fixture at the first width.
if __name__ == "__main__":
    import sys
    if "--literal" in sys.argv:
        from functools import lru_cache
        def tournament(x,y,upper):
            rows=[4,17,34,7,13,27] if upper else [4,17,42,3,13,27]
            l=[int(i) for i in bin(x)[2:]]; r=[1-int(i) for i in bin(y)[:1:-1]]
            tags=l+[(22>>i)&1 for i in range(6)]+r
            n=len(tags)+1
            final=[]
            for i in range(n-1):
                row=0
                for j in range(n-1):
                    if len(l)<=i<len(l)+6 and len(l)<=j<len(l)+6:
                        b=(rows[i-len(l)]>>(j-len(l)))&1
                    else:
                        b=i<j
                    row|=b<<j
                final.append(row|((1-tags[i])<<(n-1)))
            return tuple(final+[sum(a<<i for i,a in enumerate(tags))])
        def calc(A):
            @lru_cache(None)
            def F(S,i):
                T=S^(1<<i)
                return 1 if not T else sum(F(T,j) for j in range(len(A)) if T>>j&1 and A[j]>>i&1)
            return sum(F((1<<len(A))-1,i) for i in range(len(A)))
        a=tournament(2,294,False)
        b=tournament(2,294,True)
        assert len(a)==len(b)==18
        assert [x^y for x,y in zip(a,b)]==[0]*4+[1<<5,1<<4]+[0]*12
        assert (calc(a),calc(b))==(H(1),H(1)+2)
        print("literal order-18 marked pair PASS")

from itertools import product
N1=[]
for a,b,c in product((1,3),(0,3,6),range(29)):
    N1.append(witness(dict.fromkeys(Ogg,1),((a,4),(b,9),(c,29))))
N2=[witness(dict.fromkeys(Ogg,2),((a,4),)) for a in (1,3)]
print("Shortest Ogg by full free-digit CRT enumeration:",min(N1),"classes",len(set(N1)))
print("Shortest squared Ogg:",min(N2),"classes",len(set(N2)))


def star_rows(core, tags, x, y):
    """Ordinal left flank, marked core, right flank, and last vertex as apex."""
    left = bin(x)[2:] if x else ''
    right = bin(y)[2:][::-1].translate(str.maketrans('01', '10')) if y else ''
    tape = left + ''.join(str((tags >> i) & 1) for i in range(len(core))) + right
    n = len(tape)
    mask = (1 << n) - 1
    rows = [mask ^ ((1 << (i + 1)) - 1) for i in range(n)]
    width = len(core)
    core_mask = ((1 << width) - 1) << len(left)
    for i, row in enumerate(core):
        index = len(left) + i
        rows[index] = (rows[index] & ~core_mask) | (row << len(left))
    apex = 0
    for i, bit in enumerate(tape):
        if bit == '1':
            apex |= 1 << i
        else:
            rows[i] |= 1 << n
    return tuple(rows + [apex])

def subpaths(rows, starts=None):
    n = len(rows)
    size = 1 << n
    end = [[0] * n for _ in range(size)]
    for i in range(n):
        if starts is None or starts >> i & 1:
            end[1 << i][i] = 1
    for subset in range(1, size):
        for i, count in enumerate(end[subset]):
            if count:
                choices = rows[i] & ~subset
                while choices:
                    bit = choices & -choices
                    end[subset | bit][bit.bit_length() - 1] += count
                    choices ^= bit
    return [sum(row) for row in end]

def profile(core, tags):
    n = len(core)
    full = (1 << n) - 1
    paths = subpaths(core)
    paths[0] = 1
    h = paths[-1]
    q = sum(paths[s] * paths[full ^ s] for s in range(full + 1))

    def eo(graph, marked):
        rooted = subpaths(graph, marked)
        return tuple(sum(rooted[s] * paths[full ^ s] for s in range(1, full + 1)
                         if s.bit_count() % 2 == parity) for parity in (0, 1))

    reverse = tuple(row ^ full ^ (1 << i) for i, row in enumerate(core))
    k = subpaths(star_rows(core, tags, 0, 0))[-1]
    return h, q, eo(core, tags), eo(reverse, full ^ tags), k

def strong(rows):
    target=(1<<len(rows))-1
    for u in range(len(rows)):
        mask=1<<u
        while True:
            new=mask
            for i,r in enumerate(rows):
                if mask>>i&1:new|=r
            if new==mask:break
            mask=new
        if mask!=target:return False
    return True
low=(4,17,42,3,13,27)
high=(4,17,34,7,13,27)
assert profile(low,22)==(29,320,(58,71),(72,90),111)
assert profile(high,22)==(29,320,(59,70),(73,89),113)
assert [n for n in range(3) if Hmod(n,3)==T%3]==[0]
assert [n for n in range(14) if Hmod(n,29)==T%29]==[1]
assert [order4(p) for p in (5,7,11,13,17,19,23,29,31,41,47,59,71)]==[2,3,5,6,4,9,11,14,5,10,23,29,35]
n=1
y=int(bin(73)[2:]+"10"*n,2)
pair=[star_rows(c,22,2,y) for c in (low,high)]
assert all(map(strong,pair))
assert sum((x^y).bit_count() for x,y in zip(*pair))==2
print("Core shape and literal physical arc replayed.")
