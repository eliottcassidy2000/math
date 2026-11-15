import itertools, random
from fractions import Fraction

def H_C(rows):
    n = len(rows)
    H = C = 0
    def visit(v,last,used):
        nonlocal H,C
        if used == (1<<n)-1:
            H += 1
            if v==0 and rows[last]&1:
                C += 1
            return
        rem=rows[last]&~used
        while rem:
            lo=rem&-rem
            rem^=lo
            visit(v,lo.bit_length()-1,used|lo)
    for v in range(n): visit(v,v,1<<v)
    return H,C

def make(a,b,ata,btb,between):
    n=a+b-1
    rows=[0]*n
    B=[0]+list(range(a,n))
    amask={i:k for k,i in enumerate(itertools.combinations(range(a),2))}
    bmask={i:k for k,i in enumerate(itertools.combinations(B,2))}
    cross={i:k for k,i in enumerate(itertools.product(range(1,a),range(a,n)))}
    for i,j in itertools.combinations(range(n),2):
        if (i,j) in amask: d=ata>>amask[i,j]&1
        elif (i,j) in bmask: d=btb>>bmask[i,j]&1
        else: d=between>>cross[i,j]&1
        if d: rows[i]|=1<<j
        else: rows[j]|=1<<i
    return rows

def single(n,code):
    return H_C(make(n,1,code,0,0))
rng=random.Random(20261115)
checks=0
least=None
for a,b in [(3,3),(3,4),(3,5),(4,4),(4,5),(5,5)]:
    As=[(w,single(a,w)) for w in range(2**(a*(a-1)//2)) if single(a,w)[1]]
    Bs=[(w,single(b,w)) for w in range(2**(b*(b-1)//2)) if single(b,w)[1]]
    if a+b<=7:
        cases=itertools.product(As,Bs,range(2**((a-1)*(b-1))))
    else:
        cases=((rng.choice(As),rng.choice(Bs),rng.randrange(2**((a-1)*(b-1)))) for _ in range(1200))
    for (x,(hx,cx)),(y,(hy,cy)),z in cases:
        R=make(a,b,x,y,z)
        h,c=H_C(R)
        if h-3*c<hx*hy-3*cx*cy: raise RuntimeError((a,b,R,h,c,hx,cx,hy,cy))
        if c>cx*cy:
            value=Fraction(h-hx*hy,c-cx*cy)
            item=(value,a,b,R,h,c,hx,cx,hy,cy)
            if least is None or item<least: least=item
        checks+=1
print('Strong split controls:',checks)
print('Smallest observed overlap-surplus ratio with added cycles:',least)
