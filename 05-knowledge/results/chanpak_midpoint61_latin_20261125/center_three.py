"""Read the paired cap-two profiles. Solve (2,1,1) speed identities exactly."""
from itertools import permutations
from collections import Counter
from pathlib import Path
import sys

EXPONENT=61
MID=(1<<(EXPONENT-1))*((1<<EXPONENT)-1)
HALF=(EXPONENT-1)//2


def weighted(x, y, a, b, total=0):
    low=min(a,b,0)
    return (x<<(a-low))+(y<<(b-low)) == total<<-low


def height(poly, widths):
    return sum(c << sum(w for k,w in enumerate(widths) if j>>k&1)
               for j,c in enumerate(poly))


def mode211(poly):
    if poly[0] or poly[7]<=0 or poly[7]&(poly[7]-1):
        return
    offsets=2-(poly[7].bit_length()-1)
    for i in range(3):
        j,k=[q for q in range(3) if q!=i]
        x,y,z=1<<i,1<<j,1<<k
        for a in range(1-2*HALF,offsets+2*HALF-1):
            if not weighted(poly[x],poly[y|z],a,offsets-a,-2):
                continue
            for b in range(1-HALF, offsets-a+HALF):
                c=offsets-a-b
                if not weighted(poly[y],poly[z],b,c):
                    continue
                if not weighted(poly[x|y],poly[x|z],b,c):
                    continue
                d=[0]*3
                d[i],d[j],d[k]=2*HALF+a,HALF+b,HALF+c
                yield tuple(d)


def min_order_solutions(poly):
    target=2*MID
    max_w=target.bit_length()-1
    def limit(d):
        return max_w if d==0 else min(max_w,(d&-d).bit_length()-1)
    small_limit=limit(target-poly[0])
    for i,j,k in permutations(range(3)):
        a,b,c=1<<i,1<<j,1<<k
        for x in range(1,small_limit+1):
            constant=poly[0]+(poly[a]<<x)
            slope=poly[b]+(poly[a|b]<<x)
            end=poly[c]+(poly[a|c]<<x)
            mixed=poly[b|c]+(poly[7]<<x)
            high=limit(target-constant)
            for y in range(x+(i>j),high+1):
                resid=target-constant-(slope<<y)
                odd_part=end+(mixed<<y)
                if odd_part<=0 or resid<=0:
                    continue
                # In an honest source the last port participates in a returning
                # word even when the other ports have width one.
                val,rem=divmod(resid,odd_part)
                if rem or val&(val-1):
                    continue
                z=val.bit_length()-1
                if y+(j>k) <= z <=max_w:
                    widths=[0]*3
                    widths[i],widths[j],widths[k]=x,y,z
                    yield tuple(widths)


def scan(paths, unrestricted=False):
    unique=set()
    best=None
    hits=0
    signed_heights = set()
    shapes = Counter()
    profiles = set()
    for path in paths:
        for line in Path(path).read_text().splitlines():
            s,*data=line.split()
            data=tuple(map(int,data))
            f,g=data[5:13],data[13:]
            if (f,g) in unique:
                continue
            unique.add((f,g))
            poly=[a+b for a,b in zip(f,g)]
            candidates = min_order_solutions(poly) if unrestricted else mode211(poly)
            for widths in candidates:
                h1,h2=height(f,widths),height(g,widths)
                assert h1+h2==2*MID
                delta=abs(h1-h2)//2
                hits+=1
                signed_heights.add(tuple(sorted((h1,h2))))
                shapes[tuple(sorted(widths))] += 1
                norm=[]
                for perm in permutations(range(3)):
                    ordering=[sum(1<<perm[k] for k in range(3) if j>>k&1)
                              for j in range(8)]
                    norm.append(tuple(sorted(tuple(h[j] for j in ordering) for h in (f,g))))
                profiles.add(min(norm))
                score=delta,len(s),sum(widths)
                if best is None or score<best[0]:
                    best=(score,s,data[:5],f,g,widths,h1,h2)
    print("distinct input pairs",len(unique),"hits",hits,"best",best)
    print("distinct endpoint pairs",len(signed_heights),"width multisets",len(shapes),
          "profiles up to port permutation and flip",len(profiles))
    if 0<hits<=1000:
        print("profiles",sorted(profiles),"endpoint pairs",sorted(signed_heights))


def direct_witnesses():
    source_dir=Path(__file__).resolve().parent.parent/"chanpak_support_birth_grammar_20261122"
    sys.path.insert(0,str(source_dir))
    from check_perfect_port import H, inflate, power_coeffs, signature, strong
    cases=[
        ([462,508,504,368,481,449,128,264,64],(1,2,6),(5,8),(1,1,2),
         ((-1,0,0,-5,-1,0,0,26),(1,0,0,-1,-1,0,0,38)),(391,601),True),
        ([922,1020,1017,1008,608,961,897,784,16,256],(1,2,3),(8,4),(2,1,1),
         ((39,-19,-19,19,49,0,-49,49),(9,-5,-5,5,15,0,-15,15)),(763,229),False),
    ]
    for rows,ports,(a,b),widths,polys,counts,other_strong in cases:
        assert strong(rows) and rows[a]>>b&1
        for expected,count in zip(polys,counts):
            assert power_coeffs(signature(rows,ports))==expected
            assert H(inflate(rows,ports,widths))==height(expected,widths)==count
            rows[a]^=1<<b
            rows[b]^=1<<a
            if expected==polys[0]:
                assert strong(rows)==other_strong
        assert sum(counts)==2*496
    print("Both one-arc pairs and their independent literal exponent-five counts verified.")


if __name__=="__main__":
    if sys.argv[1:2]==["--witness"]:
        direct_witnesses()
    elif sys.argv[1:2]==["--arbitrary"]:
        scan(sys.argv[2:],True)
    else:
        scan(sys.argv[1:])
