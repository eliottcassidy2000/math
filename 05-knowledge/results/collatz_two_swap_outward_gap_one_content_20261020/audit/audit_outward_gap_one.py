#!/usr/bin/env python3
"""Independent exact audit: no imports of producer or repository mathematics.

Physical paths: all primitive 2<=p<=120 and 1<=r<p, built from ceiling
positions and cyclic 0110 patterns. Short-lattice resultants: every legal
parameter through p<=60, using QQ denominator clearing and SymPy's resultant
rather than a hand-built Sylvester matrix. Arbitrary-row coprimality: every
independent pair among sup-height<=2 rows for p<=12. All checks remain active
under python -O. Runtime/dependency metadata are separated from mathematical
scope. The output is finite evidence, not an all-length exclusion.
"""
from fractions import Fraction
from hashlib import sha256
from itertools import product, combinations
from math import gcd
from pathlib import Path
import json
import sympy as S

X = S.Symbol('X')
CHECKS = 0

def require(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ArithmeticError(label)


def encode(obj):
    return json.dumps(obj, sort_keys=True, separators=(',', ':')).encode() + b'\n'


def source_word(p, q):
    # Inverse floor jumps, with labels j=1,...,p rather than j=0,...,p-1.
    positions = [((j*q+p-1)//p)-1 for j in range(1, p+1)]
    word = [int(i in set(positions)) for i in range(q)]
    labels = {pos: (-j*(q-p)) % p for j, pos in enumerate(positions, 1)}
    require(sum(word) == p, 'source cardinality')
    return word, labels


def physical_paths(word):
    q = len(word)
    paths = []
    for i in range(q):
        slots = tuple((i+k) % q for k in range(4))
        if len(set(slots)) != 4:
            continue
        if tuple(word[s] for s in slots) == (0, 1, 1, 0):
            paths.append((slots[0], slots[2]))
    return paths


def apply(word, edges):
    target = word[:]
    for e in edges:
        f = (e+1) % len(word)
        require(target[e] != target[f], 'actual swap')
        target[e], target[f] = target[f], target[e]
    return target


def carry(word):
    # Reverse scan and increment the suffix multiplier only on a 1.
    suffix_power = 1
    ans = 0
    for i in range(len(word)-1, -1, -1):
        if word[i]:
            ans += (1 << i) * suffix_power
            suffix_power *= 3
    return ans


def affine_check(word, B, D):
    # This independent rational affine recurrence checks the actual fixed point.
    initial = Fraction(B, D)
    y = initial
    for bit in word:
        y = (3*y+1)/2 if bit else y/2
    require(y == initial, 'rational branch orbit does not close')


def canonical_root(p, q, modulus):
    v = pow(p, -1, q)
    t = (1-p*v)//q
    z = pow(2, v, modulus) * pow(3, t, modulus) % modulus
    require(pow(z, p, modulus) == 2 % modulus, 'root p clock')
    require(pow(z, q, modulus) == 3 % modulus, 'root q clock')
    require(gcd(z*(z-1), modulus) == 1, 'root and root-minus-one units')
    return z


def independent(v, w):
    return any(v[i]*w[j] != v[j]*w[i] for i,j in combinations(range(3), 2))


def orient(v):
    b,c,m = v
    # Orient by c, then m, then b, unlike producer's b,c,m orientation.
    first = next(t for t in (c,m,b) if t)
    return v if first > 0 else tuple(-x for x in v)


def all_short_rows(p,r,u,H):
    rows = set()
    for b,c in product(range(-H,H+1), repeat=2):
        residue = (b*r+c*u) % p
        for m in range(-H,H+1):
            v = (b,c,m)
            if any(v) and m % p == residue:
                rows.add(orient(v))
    return sorted(rows, key=lambda v:(max(map(abs,v)), abs(v[2])+v[1], v[1], v[2], v[0]))


def short_pair(p,r,u):
    first = None
    for H in range(1,p+1):
        rows = all_short_rows(p,r,u,H)
        if rows and first is None:
            first = rows[0]
        if first is not None:
            for v in rows:
                if independent(first,v):
                    h1 = max(map(abs,first)); h2 = max(map(abs,v))
                    require((h1*h2)**3 <= p*p, 'successive minima product')
                    return first,v
    raise ArithmeticError('short pair absent')


def polynomial(p,r,u,v,legal=True):
    b,c,m = orient(v)
    a,rem = divmod(m-b*r-c*u,p)
    require(rem == 0, 'nonrelation row')
    e = a-b+2*c
    constant = S.Rational(-1)**c*S.Rational(2)**e*S.Rational(3)**b
    # Rational-coefficient creation then exact denominator clearing is
    # independent of producer's coefficient/binomial-array construction.
    rational = S.Poly(X**max(m,0)*(X-1)**c-constant*X**max(-m,0), X, domain=S.QQ)
    denominator, f = rational.clear_denoms(convert=True)
    require(not f.is_zero, 'zero polynomial')
    if legal:
        H = max(map(abs,v))
        require(abs(a)<=H, 'a height bound')
        require(f.eval(0)!=0 and f.eval(1)!=0, 'root zero or one')
        require((f.degree()==0)==(m==0 and c==0), 'constant iff')
        require(f.degree()<=2*H, 'degree bound')
        require(sum(abs(int(t)) for t in f.all_coeffs())<=22**H, 'height 22 bound')
    return f


def pair_check(p,r,u,v,w,z,g):
    require(independent(v,w), 'dependent rows passed to pair test')
    f = polynomial(p,r,u,v); h = polynomial(p,r,u,w)
    require(max(f.degree(),h.degree()) > 0, 'two constants')
    R = int(f.resultant(h))
    require(R!=0, 'vanishing resultant')
    require(R%g==0, 'full prime-power content does not divide resultant')
    require(int(f.eval(z))%g==0 and int(h.eval(z))%g==0, 'modular polynomial gate')
    require(f.gcd(h).degree()==0, 'polynomial common factor')
    h1=max(map(abs,v)); h2=max(map(abs,w))
    l1=sum(abs(int(t)) for t in f.all_coeffs()); l2=sum(abs(int(t)) for t in h.all_coeffs())
    require(abs(R)<=l1**h.degree()*l2**f.degree(), 'Hadamard bound')
    require(abs(R)<=22**(4*h1*h2), 'resultant exponent bound')
    return R, [int(t) for t in f.all_coeffs()], [int(t) for t in h.all_coeffs()]


def main():
    path_hash = sha256(); resultant_hash = sha256(); arbitrary_hash = sha256()
    shapes=pairs=cut=positive=negative=small_pairs=resultants=arbitrary_pairs=0
    max_g = None; controls=[]; hits=[]; surplus=[]
    for p in range(2,121):
        for r in range(1,p):
            if gcd(p,r)!=1: continue
            q=p+r; shapes+=1; source,labels=source_word(p,q)
            actual=physical_paths(source)
            by_u={labels[right]: (left,right) for left,right in actual}
            require(len(by_u)==len(actual), 'duplicate residue label')
            predicted=set(range(max(0,p-2*r),min(r,p-r))) if p>=3 else set()
            require(set(by_u)==predicted, 'residue versus physical universe mismatch')
            D=(1<<q)-3**p
            require(D!=0 and gcd(abs(D),6)==1, 'clock units')
            z=canonical_root(p,q,abs(D)) if abs(D)>1 else 0
            for u,edges in sorted(by_u.items()):
                left,right=edges
                require(labels[(left+1)%q]==u+r, 'pair source residues')
                require(u<p-r and 3*r>p and r+u<=p-1, 'legal inequalities')
                target=apply(source,edges)
                require(target==apply(source,edges[::-1]), 'orders not commuting')
                B=carry(target); g=gcd(abs(D),B)
                gate=4+pow(z,u,abs(D))*(z-1) if abs(D)>1 else 0
                require(gcd(abs(D),gate)==g, 'full content iff gate')
                row=dict(p=p,r=r,u=u,D=D,B=B,g=g,edges=list(edges))
                path_hash.update(encode(row)); pairs+=1; cut+=int(q-1 in edges)
                positive+=int(D>0); negative+=int(D<0)
                if p<=20: affine_check(target,B,D)
                if max_g is None or g>max_g['g']: max_g=row
                if g==abs(D): hits.append(row)
                factors=S.factorint(g) if g>1 else {}
                if any(e>1 for e in factors.values()) and len(controls)<8:
                    controls.append(dict(row,factorization={str(ell):int(e) for ell,e in factors.items()}))
                if p<=60:
                    small_pairs+=1
                    v,w=short_pair(p,r,u)
                    R,f,h=pair_check(p,r,u,v,w,z,g)
                    resultants+=1
                    resultant_hash.update(encode(dict(row,v=v,w=w,f=f,h=h,R=R)))
                    if gcd(abs(D),abs(R))>g and len(surplus)<5:
                        surplus.append(dict(row,R=R,resultant_clock_gcd=gcd(abs(D),abs(R)),v=v,w=w))
                if p<=12:
                    rows=all_short_rows(p,r,u,2)
                    polys={v:polynomial(p,r,u,v) for v in rows}
                    for v,w in combinations(rows,2):
                        if not independent(v,w): continue
                        f,h=polys[v],polys[w]
                        require(f.gcd(h).degree()==0, 'arbitrary row common root')
                        R=int(f.resultant(h))
                        require(R!=0 and R%g==0, 'arbitrary full-content resultant')
                        arbitrary_pairs+=1
                        arbitrary_hash.update(encode([p,r,u,v,w,R]))
    require(not hits, 'integral hit')
    require(small_pairs==7367, 'producer pair count')
    # Dependent projections, nonzero constants, and irrelevant resultant factor.
    degeneracies=[]
    for p,r,u in [(3,2,0),(5,2,1),(5,3,0),(26,17,5)]:
        D=(1<<(p+r))-3**p; source,labels=source_word(p,p+r)
        pairs_u={labels[right]:(left,right) for left,right in physical_paths(source)}
        B=carry(apply(source,pairs_u[u]));g=gcd(abs(D),B);z=canonical_root(p,p+r,abs(D))
        v=(0,1,u); w=(p,1,u)
        require(independent(v,w) and v[1:]==w[1:], 'dependent projection control')
        R,f,h=pair_check(p,r,u,v,w,z,g)
        c=(p,0,0)
        Rc,fc,hc=pair_check(p,r,u,c,v,z,g)
        require(polynomial(p,r,u,c).degree()==0,'nonzero constant control')
        degeneracies.append(dict(p=p,r=r,u=u,g=g,dependent_projection_resultant=R,constant_resultant=Rc))
    # Canonical unrestricted-smooth counterexample, kept outside legal checks.
    illegal_v=(1,0,1);illegal_w=(0,2,0)
    f=polynomial(3,1,9,illegal_v,False);h=polynomial(3,1,9,illegal_w,False)
    require(f.eval(S.Rational(3,2))==0 and h.eval(S.Rational(3,2))==0, 'illegal hostile root')
    require(f.resultant(h)==0, 'illegal hostile resultant')
    # Exact height corners and rational archimedean separators.
    corners=[]
    for a,b,c in product((-1,1),(-1,1),(0,1)):
        e=a-b+2*c
        corners.append(2**(max(-e,0)+c)*3**max(-b,0)+2**max(e,0)*3**max(b,0))
    require(corners==[4,10,7,5,7,22,4,14], 'coefficient corner list')
    require(Fraction(3,2)**3>3 and Fraction(5,2)**3<16, 'norm comparison separators')
    result=dict(status='FINITE-EXACT independent audit; not an all-length exclusion',
        script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),sympy_version=S.__version__,
        physical_universe='all coprime 2<=p<=120, 1<=r<p, cyclic four-distinct-position 0110 patterns',
        shapes=shapes,unordered_paths=pairs,ordered_paths=2*pairs,ordered_cut_crossing_paths=2*cut,
        positive_clock_unordered_paths=positive,negative_clock_unordered_paths=negative,
        short_resultant_universe='every legal outward gap-one parameter with p<=60',
        short_resultants=resultants,
        arbitrary_row_universe='all independent pairs of lattice rows of sup-height<=2 for every legal p<=12 parameter, oriented by c,m,b',
        arbitrary_row_pairs=arbitrary_pairs,integer_hits=hits,maximum_content=max_g,
        prime_power_controls=controls,resultant_surplus_controls=surplus,
        degenerate_controls=degeneracies,
        illegal_smooth_common_root=dict(p=3,r=1,u=9,root='3/2',polynomials=[str(f.as_expr()),str(h.as_expr())],resultant=0),
        coefficient_corners=corners,physical_trace_sha256=path_hash.hexdigest(),
        short_resultant_trace_sha256=resultant_hash.hexdigest(),arbitrary_row_trace_sha256=arbitrary_hash.hexdigest(),
        checks=CHECKS)
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
