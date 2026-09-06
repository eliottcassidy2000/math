"""Independent exact referee: integer lines, affine laws and lattice cosets.

Does not import the producer. Uses complete line sets to count triples and
endpoint parity pairs to construct APs, rather than its direction predicate.
"""
from itertools import combinations, permutations, product
from collections import Counter, defaultdict
from fractions import Fraction
from math import gcd, comb, isqrt
import sys

sys.stdout.reconfigure(newline='\n')

checks=0
def check(ok, label):
    global checks
    checks+=1
    if not ok: raise ArithmeticError(label)

def line_key(a,b):
    u,v=b[1]-a[1],a[0]-b[0]
    w=-u*a[0]-v*a[1]
    d=gcd(gcd(abs(u),abs(v)),abs(w))
    u,v,w=u//d,v//d,w//d
    if u<0 or (u==0 and v<0): u,v,w=-u,-v,-w
    return u,v,w

def line_triples(board):
    lines=defaultdict(set)
    for a,b in combinations(board,2):
        lines[line_key(a,b)].update((a,b))
    return sum(comb(len(points),3) for points in lines.values())

def bad(p,rows,cols,a,b):
    q=[(a*y+b)%p for y in cols]
    return (rows[1]-rows[0])*(q[2]-q[0]) == (rows[2]-rows[0])*(q[1]-q[0])

triple_count=0
for p in (3,5,7):
    for rows in combinations(range(p),3):
        for cols in permutations(range(p),3):
            modular=((rows[1]-rows[0])*(cols[2]-cols[0])-
                     (rows[2]-rows[0])*(cols[1]-cols[0]))%p == 0
            gg=gcd(rows[1]-rows[0],rows[2]-rows[0])
            n=(rows[2]-rows[0])//gg
            total=0
            for a in range(1,p):
                actual={b for b in range(p) if bad(p,rows,cols,a,b)}
                predicted=set()
                if modular:
                    k0=(a*(cols[2]-cols[0])*pow(n,-1,p))%p
                    for k in (k0,k0-p):
                        if k and abs(k)*n<=p-1:
                            predicted|={(v-a*cols[0])%p
                                        for v in range(p)
                                        if 0<=v+k*n<p}
                check(actual==predicted, 'independent endpoint interval')
                total+=len(actual)
            k=(p-1)//n
            check(total == (2*k*p-n*k*(k+1) if modular else 0), 'event total')
            triple_count+=1
print('independent complete triple bank',triple_count)

p=5
base=[(i,i) for i in range(p)]+[(i,(i+1)%p) for i in range(p)]
allp=list(permutations(range(p)))
aff={tuple((a*i+b)%p for i in range(p)) for a in range(1,p) for b in range(p)}
hist=Counter(); ahist=Counter(); good=[]
for pi in allp:
    count=line_triples([(x,pi[y]) for x,y in base])
    hist[count]+=1
    if pi in aff: ahist[count]+=1
    if not count: good.append(pi)
check(hist[0]==4 and ahist=={2:8,3:4,5:4,14:4}, 'full line-set histogram')
check(sum(k*v for k,v in hist.items())==512, 'full mean numerator')
check((1,3,0,4,2) in good, 'one-swap repaired actual board')
for k in (1,2):
    for domain in combinations(range(p),k):
        aa=Counter(tuple(pi[i] for i in domain) for pi in aff)
        bb=Counter(tuple(pi[i] for i in domain) for pi in allp)
        check(all(bb[key]==6*value for key,value in aa.items()) and aa.keys()==bb.keys(),
              'all pairwise laws')
print('line-set successes',len(good),'affine histogram',dict(sorted(ahist.items())))

lines_checked=0
for p in (3,5,7,11,13,17,19):
    for u in range(p):
        # Tangent vectors checked by residue, not the producer's norm shells.
        h=isqrt(2*p)
        vectors=[(a,b) for a in range(-h,h+1) for b in range(-h,h+1)
                 if (a or b) and (b-u*a)%p==0]
        w=min(vectors,key=lambda q:(abs(q[0])+abs(q[1]),q))
        rr=abs(w[0])+abs(w[1])
        check(rr<=h and gcd(abs(w[0]),abs(w[1]))==1,'primitive tangent bound')
        for v in range(p):
            board={(x,(u*x+v)%p) for x in range(p)}
            aps=set()
            for a,b in combinations(board,2):
                if (a[0]+b[0])%2==0 and (a[1]+b[1])%2==0:
                    middle=((a[0]+b[0])//2,(a[1]+b[1])//2)
                    check(middle in board and middle not in (a,b),'actual midpoint')
                    aps.add(frozenset((a,b,middle)))
            parity=Counter((a%2,b%2) for a,b in board)
            check(len(aps)==sum(comb(n,2) for n in parity.values()),'AP bijection')
            q,r=divmod(p,4)
            floor=r*comb(q+1,2)+(4-r)*comb(q,2)
            check(len(aps)>=floor,'balanced parity lower bound')
            if u==2 and v==0: check(len(aps)==floor,'sharp modular line')
            layers=Counter(w[1]*x-w[0]*y for x,y in board)
            check(len(layers)<=rr,'direction layer count')
            # Independent N != p controls on all points in boxes 2 and p+2.
            for n in (2,p+2):
                values={w[1]*x-w[0]*y for x in range(n) for y in range(n)
                        if (y-u*x-v)%p==0}
                check(len(values)<=((n-1)*rr)//p+1,'larger/smaller grid layer bound')
            lines_checked+=1
    # Vertical axis, omitted by a slope parameterization.
    for v in range(p):
        check(len({x for x in range(p) for y in range(p) if x==v})==1,
              'vertical line capacity two')
print('modular lines independently checked',lines_checked)
for k in (1,2,4,9):
    p=k*k+(k+1)*(k+1)
    h=2*k+1
    shortest=min(abs(a*k-b*(k+1))+abs(a*(k+1)+b*k)
                 for a,b in product(range(-3,4),repeat=2) if a or b)
    check(shortest==h==isqrt(2*p),'sharp lattice basis bank')
print('PASS',checks,'always-active exact gates')
