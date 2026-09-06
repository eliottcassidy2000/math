"""Exact affine-column orbit decoder and modular-line capacity controls.

Standalone standard-library verifier. All collinearity tests are literal
integer determinants unless explicitly named modular. No producer imported.
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from math import gcd, isqrt, comb
import sys

sys.stdout.reconfigure(newline="\n")
GATES = 0


def need(ok, label):
    global GATES
    if not ok:
        raise ArithmeticError(label)
    GATES += 1


def determinant(triple):
    (x,y),(u,v),(a,b) = triple
    return (u-x)*(b-y)-(v-y)*(a-x)


def count_bad(board):
    return sum(determinant(t)==0 for t in combinations(board,3))


def directed_bad(board):
    """Independent fixed-point primitive-direction repeated-ray test."""
    for x,y in board:
        seen = set()
        for a,b in board:
            dx,dy = a-x,b-y
            if dx==dy==0:
                continue
            d = gcd(abs(dx),abs(dy))
            dx,dy = dx//d,dy//d
            if dx<0 or (dx==0 and dy<0):
                dx,dy = -dx,-dy
            if (dx,dy) in seen:
                return True
            seen.add((dx,dy))
    return False


def transform(board,p,a,b):
    return tuple((x,(a*y+b)%p) for x,y in board)


def triple_events(triple,p):
    """Exact forbidden translation masks indexed by nonzero multiplier.

    Domain: distinct rows and columns. Nonmodular triples contribute none.
    The dictionary is generated only over the possible small signed slopes.
    """
    triple = sorted(triple)
    (x0,y0),(x1,y1),(x2,y2) = triple
    if len({x0,x1,x2})<3 or len({y0,y1,y2})<3:
        raise ValueError("triple requires three distinct rows and columns")
    if determinant(triple)%p:
        return {}
    g = gcd(x1-x0,x2-x0)
    m,n = (x1-x0)//g,(x2-x0)//g
    slope = ((y1-y0)*pow(x1-x0,-1,p))%p
    scale = pow(slope*g,-1,p)
    events = {}
    for k in range(-((p-1)//n),(p-1)//n+1):
        if not k:
            continue
        a = k*scale%p
        lo,hi = max(0,-k*n),min(p-1,p-1-k*n)
        mask = sum(1<<((target-a*y0)%p) for target in range(lo,hi+1))
        need(a not in events, "signed slope representatives cannot collide")
        need(mask.bit_count()==p-abs(k)*n, "exact cyclic interval cardinality")
        events[a] = mask
    return events


def probability(triple,p):
    triple = sorted(triple)
    if determinant(triple)%p:
        return Q(0)
    g = gcd(triple[1][0]-triple[0][0],triple[2][0]-triple[0][0])
    n = (triple[2][0]-triple[0][0])//g
    K = (p-1)//n
    return Q(2*K*p-n*K*(K+1),p*(p-1))


def compile_orbit(board,p):
    """Every good affine-column pair, as an exact finite iff predicate."""
    board = tuple(board)
    need(len(set(board))==len(board), "board points distinct")
    need(all(0<=x<p and 0<=y<p for x,y in board), "board lies in prime square")
    full = (1<<p)-1
    masks = {a:0 for a in range(1,p)}
    if max(Counter(x for x,y in board).values(),default=0)>=3 or max(
            Counter(y for x,y in board).values(),default=0)>=3:
        return {a:full for a in masks}
    for triple in combinations(board,3):
        if len({x for x,y in triple})<3 or len({y for x,y in triple})<3:
            continue
        for a,mask in triple_events(triple,p).items():
            masks[a] |= mask
    return masks


def literal_masks(board,p):
    return {a:sum((1<<b) for b in range(p) if count_bad(transform(board,p,a,b)))
            for a in range(1,p)}


# Exhaust every three-distinct-row/column triple at p=3,5,7. This includes
# every nonmodular rejection case, not only the certificate's positive class.
triple_total = 0
for p in (3,5,7):
    for rows in combinations(range(p),3):
        for cols in permutations(range(p),3):
            triple = tuple(zip(rows,cols))
            got = triple_events(triple,p)
            literal = literal_masks(triple,p)
            need(all(got.get(a,0)==literal[a] for a in literal),
                 "complete literal triple orbit versus interval compiler")
            count = sum(mask.bit_count() for mask in literal.values())
            need(Q(count,p*(p-1))==probability(triple,p),
                 "exact primitive-span probability versus literal orbit")
            triple_total += 1
print("ALL NONAXIS TRIPLES p=3,5,7", triple_total)

# Complete row-triple banks at larger primes, retaining every gcd/span.
row_total = 0
for p in (11,13,17,19):
    for rows in combinations(range(p),3):
        triple = tuple((x,x) for x in rows)
        got = triple_events(triple,p)
        literal = literal_masks(triple,p)
        need(all(got.get(a,0)==literal[a] for a in literal),
             "larger-prime full row-span bank")
        need(Q(sum(mask.bit_count() for mask in literal.values()),p*(p-1))==probability(triple,p),
             "larger-prime triangular probability")
        row_total += 1
print("ADDITIONAL COMPLETE ROW-SPAN BANK", row_total)
need(probability(((0,0),(1,1),(3,3)),5)==Q(1,5), "primitive span three")
need(probability(((0,0),(2,2),(4,4)),5)==Q(2,5), "primitive span two")
need(determinant(((0,0),(1,2),(3,1)))==-5, "modular zero need not be integer zero")

# Complete S5 column orbit of the fixed cycle. Independent direction tests
# accompany literal triple counts on every board.
p = 5
base = tuple((r,c) for r in range(p) for c in (r,(r+1)%p))
all_perms = list(permutations(range(p)))
affine_perms = {tuple((a*c+b)%p for c in range(p)) for a in range(1,p) for b in range(p)}
need(len(affine_perms)==20, "all distinct affine columns")
hist = Counter(); aff_hist = Counter(); successes = []
for pi in all_perms:
    board = tuple((r,pi[c]) for r,c in base)
    bad = count_bad(board)
    need((bad>0)==directed_bad(board), "independent geometry of every S5 board")
    hist[bad] += 1
    if pi in affine_perms:
        aff_hist[bad] += 1
    if bad==0:
        successes.append(pi)
    need(compile_orbit(board,p)==literal_masks(board,p),
         "full affine iff decoder on every fixed-row S5 board")
need(hist[0]==4 and aff_hist[0]==0 and min(aff_hist)==2, "success destroyed by affine restriction")
need(Q(sum(k*v for k,v in hist.items()),120)==Q(64,15), "full-column exact triple mean")
need(Q(sum(k*v for k,v in aff_hist.items()),20)==Q(26,5), "affine-column exact triple mean")
print("S5 success",hist[0],"/120; affine success",aff_hist[0],"/20")
print("S5 triple mean",Q(64,15),"; affine triple mean",Q(26,5))
print("AFFINE INTEGER-TRIPLE HISTOGRAM",dict(sorted(aff_hist.items())))

# Equal one- and two-column cylinder laws, before taking a nonlinear target.
for col in range(p):
    A=Counter(pi[col] for pi in affine_perms); B=Counter(pi[col] for pi in all_perms)
    need(all(Q(A[target],20)==Q(B[target],120)==Q(1,5) for target in range(p)),
         "all one-column cylinders agree")
for col0,col1 in permutations(range(p),2):
    A=Counter((pi[col0],pi[col1]) for pi in affine_perms)
    B=Counter((pi[col0],pi[col1]) for pi in all_perms)
    need(all(Q(A[target],20)==Q(B[target],120)==Q(1,20)
             for target in permutations(range(p),2)), "all ordered two-column cylinders agree")
need({tuple((a*c+b)%3 for c in range(3)) for a in (1,2) for b in range(3)}==set(permutations(range(3))),
     "three is too small for an affine/full-law separation")

# The six sharply-two-transitive orbit representatives keep first two images
# fixed. One nonaffine transposition already unlocks a successful affine orbit.
orbit_rows = []
unused = set(all_perms)
while unused:
    pi = min(unused)
    orbit = {tuple((a*pi[c]+b)%p for c in range(p)) for a in range(1,p) for b in range(p)}
    need(len(orbit)==20 and orbit<=unused, "disjoint affine quotient orbit")
    unused -= orbit
    good = sum(partner in successes for partner in orbit)
    orbit_rows.append((pi,good))
need([good for pi,good in orbit_rows]==[0,2,0,2,0,0], "complete six-orbit success inventory")
tau = (0,1,2,4,3)
repaired = tuple((2*tau[c]+1)%5 for c in range(5))
need(repaired==(1,3,0,4,2) and repaired in successes,
     "one nonaffine transposition plus affine map repairs the board")
print("SIX AFFINE QUOTIENT ORBITS",orbit_rows)
print("ONE-SWAP REPAIR",repaired)


def prime(n):
    return n>=2 and all(n%d for d in range(2,isqrt(n)+1))


def shortest_vector(p,slope):
    for r in range(1,isqrt(2*p)+1):
        for dx in range(-r,r+1):
            dyabs = r-abs(dx)
            for dy in sorted({-dyabs,dyabs}):
                if (dx>0 or (dx==0 and dy>0)) and (dy-slope*dx)%p==0:
                    return dx,dy
    raise ArithmeticError("pigeonhole short-vector bound failed")


# Tangent lattices and actual Euclidean layer counts. These checks do not
# establish a universal geometric theorem by sampling; the proof is separate.
lattice_cases=0; parity_cases=0
for p in range(3,44):
    if not prime(p):
        continue
    ap_min=None
    q4,r4=divmod(p,4)
    ap_floor=r4*comb(q4+1,2)+(4-r4)*comb(q4,2)
    for slope in range(p):
        dx,dy=shortest_vector(p,slope)
        r=abs(dx)+abs(dy)
        need(gcd(abs(dx),abs(dy))==1 and r<=isqrt(2*p)<p,
             "primitive shortest tangent vector")
        for intercept in range(p):
            line=tuple((x,(slope*x+intercept)%p) for x in range(p))
            layers=Counter(dy*x-dx*y for x,y in line)
            need(len(layers)<=r, "actual Euclidean layers bounded by tangent L1")
            need(len({q%p for q in layers})==1, "all layer labels retain one modular residue")
            capacity=sum(min(2,count) for count in layers.values())
            need(capacity<=2*r, "complete modular line directional capacity")
            parity_sizes=Counter((x%2,y%2) for x,y in line)
            ap_count=sum(comb(size,2) for size in parity_sizes.values())
            need(ap_count>=ap_floor,"balanced parity AP lower bound")
            if slope:
                ap_min=ap_count if ap_min is None else min(ap_min,ap_count)
            if slope==2 and intercept==0:
                need(ap_count==ap_floor,"slope-two realizes balanced parity equality")
            if p<=13:
                literal_ap=sum(determinant(triple)==0 and
                    (triple[0][0]+triple[2][0]==2*triple[1][0]) and
                    (triple[0][1]+triple[2][1]==2*triple[1][1])
                    for triple in combinations(line,3))
                need(literal_ap==ap_count,"literal equal-step triples versus parity pairs")
            if p>=5:
                parity={};mid=None
                for a in line[:5]:
                    key=(a[0]%2,a[1]%2)
                    if key in parity:
                        b=parity[key];mid=((a[0]+b[0])//2,(a[1]+b[1])//2)
                        need(mid in line and mid not in (a,b) and determinant((a,b,mid))==0,
                             "full modular line midpoint hostile")
                        break
                    parity[key]=a
                need(mid is not None,"five points force a repeated parity class")
                parity_cases+=1
            # Exhaust subsets only in the two smallest odd primes.
            if p in (3,5):
                for mask in range(1<<p):
                    subset=tuple(line[i] for i in range(p) if mask>>i&1)
                    if count_bad(subset)==0:
                        need(len(subset)<=capacity,"every small safe subset obeys layer capacity")
            lattice_cases+=1
    need(ap_min==ap_floor,"exact affine modular-line AP minimum")
print("LATTICE/INTERCEPT CASES",lattice_cases,"; PARITY HOSTILES",parity_cases)
print("EXACT MODULAR-LINE AP MINIMUM: p=4q+r gives r*C(q+1,2)+(4-r)*C(q,2).")
for k in (1,2,4,9):
    p=k*k+(k+1)*(k+1)
    need(prime(p),"declared sharp-direction prime")
    slope=(k+1)*pow(k,-1,p)%p
    dx,dy=shortest_vector(p,slope)
    need(abs(dx)+abs(dy)==2*k+1==isqrt(2*p),"short-direction bound attained")
    need(k*k+(k+1)*(k+1)==p and 2*p>(2*k+1)**2,
         "orthogonal-basis proof controls")
print("SHARP SHORT-DIRECTION CONTROLS p=5,13,41,181; no prime-infinitude assertion.")
print("PASS",GATES,"always-active exact gates; no independence or extremal-asymptotic inference.")
