"""All 4+6+6 ordinal core replacements and all tags at two fixed E61 flank pairs."""
from collections import Counter
from functools import lru_cache
from itertools import product
from math import factorial
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "chanpak_decimal_kernel_20261122"))
from fermat_receiver import A, D, EDGE, TARGET, flipped, join, literal_receiver, nontags, p, tagged_data
from e15_pair import subpaths
from near_witnesses import WINNING_SLICES


def catalogue(n):
    source = HERE.parents[2] / "session212" / "cycle-levels" / f"classes{n}.rows"
    return [tuple(map(int,line.split())) for line in source.read_text().splitlines()]


@lru_cache(None)
def pair_slots(rows):
    paths = subpaths(rows)
    paths[0] = 1
    return sum(v*w for v,w in zip(paths,reversed(paths)))


@lru_cache(None)
def profiles(rows):
    return [tagged_data(rows,s) for s in range(1 << len(rows))]


def scan(specifications, outer, label):
    offsets = [Counter() for row in outer]
    unit = [0]*len(outer)
    shifted = [set() for row in outer]
    total = 0
    for left, mid, right, edge in specifications:
        upper_mid = flipped(mid,*edge)
        # Middle tag positions are held fixed when the central edge reverses.
        for l in profiles(left):
            for m, mm in zip(profiles(mid),profiles(upper_mid)):
                low_prefix, high_prefix = join(l,m), join(l,mm)
                for v in profiles(right):
                    a,b = join(low_prefix,v),join(high_prefix,v)
                    assert a[1]==b[1]==2293760
                    total+=1
                    for j,(x,y,u,base,original) in enumerate(outer):
                        # Both chosen flank pairs are even. base = q*x*y - TARGET.
                        alpha=base+a[0]*u+a[2]*y+a[3]*x+a[6]
                        beta=base+b[0]*u+b[2]*y+b[3]*x+b[6]
                        if abs(alpha-beta)!=2:
                            continue
                        unit[j]+=1
                        delta=min(alpha,beta)
                        offsets[j][delta]+=1
                        # Constant-only repair after freezing all three growing coefficients.
                        lo = a if alpha<beta else b
                        if (lo[0],lo[2],lo[3]) == (2465,1344270,712944):
                            shifted[j].add(delta-original)
    print(label, "triple/arc choices",len(specifications),"tagged presentations",total)
    for j in range(len(outer)):
        best=sorted(offsets[j],key=abs)[:3]
        print("address",j,"unit presentations",unit[j],"unique midpoints",len(offsets[j]),
              "closest offsets and multiplicities",[(o,offsets[j][o]) for o in best],
              "frozen-coefficient shifts",sorted(shifted[j]))
    return offsets


def unmarked(triples, outer):
    in_band = [set() for _ in outer]
    for left,mid,right in triples:
        for a in profiles(left):
            for d in profiles(mid):
                first=join(a,d)
                for v in profiles(right):
                    h,q,c,t,e,r,k=join(first,v)
                    for j,(x,y,u,base,original) in enumerate(outer):
                        residual=base+h*u+c*y+t*x+k
                        if abs(residual)<=10**10:
                            in_band[j].add(residual)
    print("unmarked core types",len(triples),"tagged presentations",len(triples)*65536)
    for j,band in enumerate(in_band):
        print("address",j,"all unmarked H-(E-1) within 10^10",sorted(band))
        assert not ({0,2} & band)
    return in_band


def main():
    central=(1,4)
    assert central==EDGE
    original_q=2293760
    nr,tags = literal_receiver(nontags(4))
    outer=[]
    for *_,(x,y,order,delta) in WINNING_SLICES:
        assert x%2==y%2==0
        outer.append((x,y,p(x)+p(y),original_q*x*y-TARGET,delta))
    original = scan([(A,D,nr,central)], outer, "original 16 rows")
    fours, sixes = catalogue(4),catalogue(6)
    assert (len(fours),len(sixes))==(4,56)
    maxq = max(map(pair_slots,fours))*max(map(pair_slots,sixes))**2
    # Two endpoint charts whose q differs by at least one differ at height by
    # more than x*y minus twice this universal cap on the remainder.
    for x,y,_,_,delta in outer:
        bound = factorial(17)+2*maxq*(x+y+1)
        assert x*y-2*bound>10**29+abs(delta)+2
    specs=[]
    bank=Counter()
    triples=[]
    for left,mid,right in product(fours,sixes,sixes):
        q0, q1, q2 = map(pair_slots,(left,mid,right))
        if q0*q1*q2!=original_q:
            continue
        triples.append((left,mid,right))
        for i in range(6):
            for j in range(i):
                if pair_slots(flipped(mid,i,j))!=q1:
                    continue
                specs.append((left,mid,right,(i,j)))
                bank[(q0,q1,q2)]+=1
    print("Q cap",maxq,"slot-factor triples",sorted(bank.items()))
    full=scan(specs,outer,"all isomorphism representatives")
    in_band=unmarked(triples,outer)
    assert len(triples)==16 and len(specs)==30 and maxq==5971968
    for j,(_,_,_,_,d) in enumerate(outer):
        assert d in original[j] and d in full[j]
        assert in_band[j]=={d,d+2}
        assert sum(original[j].values())==len(original[j])==64
        assert sum(full[j].values())==320 and len(full[j])==160
        assert sorted(full[j],key=abs)[:2]==[d,(10420831809831200142,6135804521871654884)[j]]


if __name__=="__main__":
    main()
