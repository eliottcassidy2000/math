"""Print or invert the exact Hamilton-path profiles on the fixed Heegner-prime frame."""
from argparse import ArgumentParser
from functools import lru_cache
from itertools import combinations
from math import comb, lcm, prod
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"chanpak_ternary_buckets_20261122"))
from target_flip import height

FRAME=(3,7,11,19,43,67,163)
PENTAD=(11,19,43,67,163)
FERMAT=(3,5,17,257,65537)
PRIMES=PENTAD+(37,)
# Mask of the extra singleton toward the frame; mask of the module toward the
# eight outside vertices; prime index of the arc to reverse at the singleton.
EDGES=(
    (11,1,5,(309,134,0,0,0)), (11,5,5,(411,238,0,0,0)),
    (55,63,6,(497,351,72,0,0)), (69,223,2,(627,297,0,0,0)),
    (73,3,3,(883,732,178,0,0)), (119,254,6,(1175,894,0,0,0)),
    (78,141,2,(1329,1697,856,144,0)),
)


def legendre_mask(q, primes):
    return sum(1<<i for i,p in enumerate(primes) if pow(q,(p-1)//2,p)==1)


def tournament(primes):
    return tuple(legendre_mask(p,primes) for p in primes)


def host(a,b,width):
    rows=list(tournament(FRAME))+[a]
    for i in range(7):
        if not (a>>i)&1:
            rows[i]|=1<<7
    for i in range(8):
        if not (b>>i)&1:
            rows[i]|=((1<<width)-1)<<8
    rows.extend(b|(((1<<(width-i-1))-1)<<(9+i)) for i in range(width))
    return rows


def powers(slots):
    return tuple(sum((-1)**(k-j)*comb(k,j)*slots[k-1] for k in range(j,6))
                 for j in range(1,6))


def value(a,s,mod=None):
    return sum(c*pow(j,s,mod) for j,c in enumerate(a,1)) if mod else sum(
        c*j**s for j,c in enumerate(a,1))


@lru_cache(None)
def derivative_rows(p):
    fermat=[(pow(j,p-1,p*p)-1)//p for j in range(2,6)]
    return tuple(tuple(q*pow(j,r,p)%p for j,q in enumerate(fermat,2)) for r in range(p-1))


def derivative(a,row,p):
    return (a[1]*row[0]+a[2]*row[1]+a[3]*row[2]+a[4]*row[3])%p


def widths(a,target,primes=PRIMES):
    """All positive widths on this growing profile having the given target height."""
    if target<=0:
        return set()
    limit=target.bit_length()+5
    hits=set()
    for p in primes:
        rows=derivative_rows(p)
        for r in range(p-1):
            if value(a,r,p)%p!=target%p:
                continue
            d=derivative(a,rows[r],p)
            if not d:
                continue
            s,period,modulus=r,p-1,p
            inverse=pow(d,-1,p)
            while s<=limit and period<=limit:
                next_mod=modulus*p
                diff=(target-value(a,s,next_mod))%next_mod
                s+=period*((diff//modulus*inverse)%p)
                modulus=next_mod
                period*=p
            if 1<=s<=limit and value(a,s)==target:
                hits.add(s)
    return hits


def singular(a,primes):
    period=1
    phases=[0]
    for p in primes:
        rows=derivative_rows(p)
        bad={r for r,row in enumerate(rows) if derivative(a,row,p)==0}
        if not bad:
            return ()
        new_period=lcm(period,p-1)
        phases=[r for t in phases for r in range(t,new_period,period) if r%(p-1) in bad]
        if not phases:
            return ()
        period=new_period
    return tuple(phases)


def read_bank(path):
    bank={}
    for line in path.read_text().splitlines():
        a,b,*p=map(int,line.split())
        if (a,b) in bank or len(p)!=6:
            raise ValueError("unexpected coefficient record")
        bank[a,b]=tuple(p)
    assert len(bank)==128*256
    return bank


def census(bank):
    states=set(bank.values())
    profiles={p[1:] for p in states}
    heights={p[1] for p in states}
    assert len(profiles)==29006 and len(heights)==1214 and max(heights)==2995
    assert 403 not in heights and 2299 not in heights and all(h in heights for h in range(405,2298,2))
    assert len({h%36 for h in heights})==18
    mids=set()
    persistent=set()
    for (a,b),p in bank.items():
        for i in range(15):
            neighbor=(a^(1<<i),b) if i<7 else (a,b^(1<<(i-7)))
            q=bank[neighbor]
            if q[1]!=p[1]+2:
                continue
            mids.add(p[1]+1)
            if i<7 and q[2:]==p[2:]:
                persistent.add(p[1:])
    assert len(mids)==441 and len(persistent)==11
    assert {p for p in persistent if any(p[1:])}=={p for a,b,i,p in EDGES}
    print("heights, midpoints, profiles, growing",len(heights),len(mids),len(profiles),
          sum(any(p[1:]) for p in profiles))
    exceptions={}
    uniform=0
    for slots in profiles:
        if not any(slots[1:]):
            continue
        a=powers(slots)
        if any(all(derivative(a,row,p) for row in derivative_rows(p)) for p in PENTAD):
            uniform+=1
        phases=singular(a,PENTAD)
        if phases:
            exceptions[a]=phases
            assert not singular(a,PRIMES)
    assert (uniform,len(exceptions),sum(map(len,exceptions.values())))==(26209,9,11)
    fail7=[(a,s) for a,ss in exceptions.items() for s in ss
           if not derivative(a,derivative_rows(7)[s%6],7)]
    assert fail7==[((285,-638,-130,612,48),34064)]
    print("pentad uniform charts",uniform,"joint exceptions",len(exceptions),11,"after 37",0)
    for a,b,i,p in EDGES:
        assert bank[a,b][1:]==p
        assert bank[a^(1<<i),b][1:]==(p[0]+2,)+p[1:]


def anchors():
    words={legendre_mask(prod(word),PENTAD) for r in range(6) for word in combinations(FERMAT,r)}
    assert len(words)==32 and prod(FERMAT)==2**32-1
    assert height(tournament(FRAME+(2683,2311)))==743
    assert height(tournament(FRAME+(2683,751)))==745
    assert legendre_mask(2683,FRAME)==27 and legendre_mask(2311,FRAME+(2683,))==167
    assert legendre_mask(751,FRAME+(2683,))==39
    assert tuple(tournament(FRAME))==(4,9,26,33,75,23,47)
    # The ramified pair of fields with discriminants -4,-8 has two completions.
    assert height(host(127,2,1))==161 and height(host(127,130,1))==983
    ordinary=(1,2,5,8,20)
    platonic=(4,8,6,20,12)
    for inputs,total in ((ordinary,4),(platonic,8)):
        signs={legendre_mask(prod(word),PENTAD) for r in range(6) for word in combinations(inputs,r)}
        assert len(signs)==total
    # A copy of the prime 7 is blown up while preserving its two old outneighbors.
    frame=tuple(p for p in FRAME if p!=7)
    for s in range(4):
        word=frame+(7,)*s
        rows=[sum(1<<j for j,q in enumerate(word)
                  if (p==q and i<j) or (p!=q and pow(p,(q-1)//2,q)==1))
              for i,p in enumerate(word)]
        assert height(rows)==(10*3**s+21*2**s-1 if s else 29)
    survivors=[]
    for a,b,i,p in EDGES:
        coef=powers(p)
        good=19 if all(derivative(coef,row,19) for row in derivative_rows(19)) else 67
        assert all(derivative(coef,row,good) for row in derivative_rows(good))
        for s in (1,2,4,63,670):
            assert s in widths(coef,value(coef,s),(good,))
        for s in (1,2):
            assert height(host(a,b,s))==value(coef,s)
            assert height(host(a^(1<<i),b,s))==value(coef,s)+2
        for s in range(2,8):
            if (value(coef,s)+1)%36==28:
                survivors.append((a,b,s))
    assert survivors==[(11,5,5),(78,141,5)]


def main():
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("--bank",type=Path,help="output of heegner_slots")
    parser.add_argument("--census",action="store_true")
    targets=parser.add_mutually_exclusive_group()
    targets.add_argument("--midpoint",type=int,help="invert all seven growing permanent edges")
    targets.add_argument("--height",type=int,help="invert all growing height profiles in the bank")
    args=parser.parse_args()
    anchors()
    if (args.census or args.height is not None) and args.bank is None:
        parser.error("--census and --height need --bank")
    bank=read_bank(args.bank) if args.bank else {}
    if args.census:
        census(bank)
    if args.height is not None:
        first={p[1:]:(a,b) for (a,b),p in bank.items() if any(p[2:])}
        matched=False
        for slots,(a,b) in first.items():
            for s in sorted(widths(powers(slots),args.height)):
                print("height",args.height,"fixed mask",a,"module mask",b,"width",s)
                matched=True
        if not matched:
            print("No growing Heegner-frame chart has this height.")
    if args.midpoint is not None:
        matched=False
        for a,b,i,p in EDGES:
            for s in sorted(widths(powers(p),args.midpoint-1,(19,67))):
                print("midpoint",args.midpoint,"fixed mask",a,"module mask",b,"width",s,
                      "flip fixed singleton and prime",FRAME[i])
                if s<=8:
                    print("lower rows",*host(a,b,s))
                    print("upper rows",*host(a^(1<<i),b,s))
                matched=True
        if not matched:
            print("No growing permanent Heegner-frame edge from this bank crosses this midpoint.")


if __name__=="__main__":
    main()
