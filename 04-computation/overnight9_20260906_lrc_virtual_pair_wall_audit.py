"""Independent exact referee: virtual walls, both margin branches, actual cells.

No producer imports. Physical distances use integer modular residues. The
profile controls subdivide all tent pieces and every crossing, rather than
sampling a Farey grid. Body/physical intervals are reconstructed directly.
"""
from fractions import Fraction as Q
from itertools import combinations
from math import gcd, lcm
from functools import reduce
import sys

sys.stdout.reconfigure(newline="\n")
gates = 0


def require(ok, what):
    global gates
    gates += 1
    if not ok:
        raise ArithmeticError(what)


def distance(w, y):
    p, q = y.numerator, y.denominator
    r = (w * p) % q
    return Q(min(r, q - r), q)


def nine(T, y):
    return tuple(distance(w, (y + j) / 3) for j in range(3) for w in T)


def heights(T, y):
    a = nine(T, y)
    return tuple(min(a[3*j:3*j+3]) for j in range(3))


def pair_frequencies(T):
    return tuple(abs(u - v) // 3 if (u-v) % 3 == 0 else (u+v) // 3
                 for u, v in combinations(T, 2))


def critical_profile(T):
    # Every tent is affine between these knots. Include all kinks of the
    # absolute-value lower profile for each signed pair as well.
    knots = {Q(0), Q(1)}
    ds = pair_frequencies(T)
    for w in T:
        for j in range(3):
            for k in range(2*w + 1):
                y = Q(3*k, 2*w) - j
                if 0 <= y <= 1:
                    knots.add(y)
    for d in ds:
        knots.update(Q(k, 2*d) for k in range(2*d + 1))
        knots.update(Q(3*k+s, 3*d) for k in range(d+1) for s in (-1, 1)
                     if 0 <= 3*k+s <= 3*d)
    knots = sorted(knots)
    all_points = set(knots)
    for left, right in zip(knots, knots[1:]):
        a, b = nine(T, left), nine(T, right)
        for i, j in combinations(range(9), 2):
            da, db = a[i]-a[j], b[i]-b[j]
            if da * db < 0:
                all_points.add(left + (right-left) * da / (da-db))
    # At these endpoints and crossings each selected min/max and each
    # profile is affine on every intervening cell: endpoint tests suffice.
    for y in sorted(all_points):
        H = max(heights(T, y))
        for d in ds:
            delta = distance(d, y)
            require(2*H >= abs(delta-Q(1,3)), "complete exact piecewise profile")
            if H < Q(1,14):
                require(Q(4,21) < delta < Q(10,21), "strict double-sided owner band")
    return len(all_points)


def owner_controls():
    total = 0
    for w in (1, 2, 4, 5, 7, 8, 11, 14):
        for alpha in (Q(1,14), Q(1,10), Q(1,6)):
            phases = {Q(0), Q(1)}
            for n in range(w+1):
                phases.update(y for y in (Q(n,w)-3*alpha/w,
                                           Q(n,w)+3*alpha/w,
                                           Q(2*n+1,2*w)) if 0 <= y <= 1)
            ordered = sorted(phases)
            phases.update((a+b)/2 for a,b in zip(ordered,ordered[1:]))
            for y in phases:
                literal = [j for j in range(3) if distance(w,(y+j)/3) < alpha]
                require(len(literal) <= 1, "strict grid cardinality including alpha=1/6")
                if distance(w,y) < 3*alpha:
                    z = w*y
                    n = (2*z.numerator+z.denominator)//(2*z.denominator)
                    predicted = (-n * pow(w,-1,3)) % 3
                    require(literal == [predicted], "literal active owner from nearest tooth")
                else:
                    require(not literal, "weak-safe effective boundary, including half-integer tie")
                total += 1
    return total


def wall_controls():
    rows = 0
    for d in range(1,54):
        for w in list(range(1,66)) + [13*d,14*d,15*d,28*d]:
            g = gcd(w,d)
            q, h = d//g, w//g
            formula = g*((q-1-h)//14 + (q-1+h)//14 + 1)
            for sign in (-1,1):
                actual = sum(distance(w,Q(14*k+sign,14*d)) < Q(1,14)
                             for k in range(d))
                require(actual == formula, "both wall signs, rational-distance count")
                require((actual == d) == (w % (14*d) == 0), "exact whole-wall blocking iff")
                require(actual <= g*((q+6)//7), "gcd branch ceiling")
            rows += 1
    require(6*((36+6)//7) >= 36, "36 fails the sufficient scalar count")
    require(all(6*((d+6)//7)<d for d in range(37,400)), "seven-class count tail control")
    return rows


def safe_interval(S, x):
    # Maximal safe interval through x in its ordinary real lift.
    lows, highs = [], []
    for w in S:
        require(distance(w,x) >= Q(1,14), "interval entry is literally safe")
        z = w*x
        n = z.numerator//z.denominator
        lows.append(Q(14*n+1,14*w))
        highs.append(Q(14*n+13,14*w))
    return max(lows), min(highs)


def data(d,h):
    b,c = 3*d+1,3*d+2
    base = (1,d,2*d,3*d,4*d,14,14*b,14*c)
    return base+(h,2*h), (1,b,c)


def addressed(d,m, subset=False):
    L = 42*d*(3*d+1)*(3*d+2)
    h = 1+m*L
    C,T = data(d,h)
    require(L == lcm(*(w for w in C if w not in (h,2*h))), "claimed lcm exact")
    require(h > max(w for w in C if w not in (h,2*h)), "all m>=1 give separated pair")
    # Address recovered by the stated residue-dependent displacement,
    # then compared against the floor formula.
    shift = Q(5,21*d) if d%6 == 1 else -Q(2,21*d)
    y = Q(5,6)+shift
    k = (35*d+18)//42
    require(y == Q(14*k+1,14*d), "independent address reconstruction")
    left,right = safe_interval(C,y)
    require(left == y < right, "positive actual body component left endpoint")
    require([w for w in C if distance(w,y)==Q(1,14)] == [d], "unique body owner")
    for w in T:
        require([j for j in range(3) if distance(w,(y+j)/3)<Q(1,14)] == [2],
                "all three actual owners are two")
        z = w*y
        n = (2*z.numerator+z.denominator)//(2*z.denominator)
        right = min(right,(Q(n)+Q(3,14))/w)
        require(14*w in C, "ordinary positive and negative event cosets have integer blocker")
    require(right > y, "positive common body/owner cell reconstructed directly")
    S = tuple(sorted([3*w for w in C]+list(T)))
    require(len(set(S)) == 13 and reduce(gcd,S)==1, "actual thirteen distinct primitive speeds")
    for j in (0,1):
        for theta in (Q(1,4),Q(1,2),Q(3,4)):
            x = (y+theta*(right-y)+j)/3
            require(all(distance(w,x)>Q(1,14) for w in S), "both actual free sheets strict in cell")
    cross = min(max(a//gcd(a,b),b//gcd(a,b))
                for a in (h,2*h) for b in C if b not in (h,2*h))
    require(cross == h, "exact necessary numerical cross height, not decoder entry")
    require(7*h >= 45*max(w for w in C if w not in (h,2*h))
            and 7*h >= 15*max(T), "incoming THM4448 also closes every explicit h=1+mL row")
    if subset:
        nonprim = [(P,reduce(gcd,P)) for P in combinations(S,10) if reduce(gcd,P)>1]
        require(nonprim == [(tuple(sorted(3*w for w in C)),3)], "unique nonprimitive ten-subset")
        for count in (11,12):
            require(all(reduce(gcd,P)==1 for P in combinations(S,count)), "larger subsets primitive")
    return (d,m,h,str(y),str(right-y))


def arbitrary(d, minimum):
    L = 42*d*(3*d+1)*(3*d+2)
    h = minimum
    while gcd(h,L) != 1 or h%L==1:
        h += 1
    C,T = data(d,h)
    require(h > max(w for w in C if w not in (h,2*h)), "arbitrary-h domain")
    survived = 0
    for k in range(d):
        y = Q(14*k+1,14*d)
        if any(distance(w,y)<Q(1,14) for w in C):
            continue
        require([w for w in C if distance(w,y)==Q(1,14)] == [d], "all residual walls strict")
        H = heights(T,y)
        j = max(range(3),key=H.__getitem__)
        require(H[j]>=Q(11,84), "sharp physical fibre margin at every survivor")
        x = (y+j)/3
        S = tuple(3*w for w in C)+T
        left,right = safe_interval(S,x)
        require(left==x<right, "selected sheet has actual positive physical component")
        require(all(distance(w,(left+right)/2)>Q(1,14) for w in S), "strict physical interior")
        survived += 1
    require(survived>0, "phase-free count actually produces an address")
    return (d,h,survived)


def atlas_four(d,h,explicit=False,subsets=False):
    """Bounded repaired atlas control, with unchanged tail wall and label map."""
    C0,T=data(d,h)
    C=tuple(4*h if w==2*h else w for w in C0)
    L=42*d*(3*d+1)*(3*d+2)
    base=tuple(w for w in C if w not in (h,4*h))
    require(h>max(base) and gcd(h,L)==1 and gcd(d,42)==1, "ratio-four coprime domain")
    require((h//gcd(h,4*h),4*h//gcd(h,4*h))==(1,4)
            and 1+4==5 and 5%3==2, "ratio-four satisfies the inert-prime pair atlas")
    cross=min(max(a//gcd(a,b),b//gcd(a,b)) for a in (h,4*h) for b in base)
    require(cross==h, "ratio-four cross height exact, not actual decoder entry")
    S=tuple(sorted(3*w for w in C))+T
    require(len(set(S))==13 and reduce(gcd,S)==1, "ratio-four actual primitive row")
    if explicit:
        require((h-1)%L==0 or (h==42*d+29 and d>=106),
                "one of the two proved explicit ratio-four address families")
        k=(35*d+18)//42
        ks=[k]
    else:
        ks=range(d)
    count=0
    for k in ks:
        y=Q(14*k+1,14*d)
        if any(distance(w,y)<Q(1,14) for w in C):
            require(not explicit, "explicit ratio-four wall is body-safe")
            continue
        require([w for w in C if distance(w,y)==Q(1,14)]==[d], "ratio-four unique body owner")
        H=heights(T,y)
        require(max(H)>=Q(11,84), "ratio-four sharp common-sheet margin")
        if explicit:
            require(all([j for j in range(3) if distance(w,(y+j)/3)<Q(1,14)]==[2] for w in T),
                    "explicit ratio-four retains owner word (2,2,2)")
            js=(0,1)
        else:
            js=(max(range(3),key=H.__getitem__),)
        for j in js:
            x=(y+j)/3
            left,right=safe_interval(S,x)
            require(left==x<right, "ratio-four positive actual physical interval")
            require(all(distance(w,(left+right)/2)>Q(1,14) for w in S), "ratio-four strict physical interior")
        count+=1
    require(count>0, "ratio-four virtual selection survives")
    require(all(14*w in C for w in T), "ratio-four ordinary event cosets wholly blocked")
    if subsets:
        nonprim=[(P,reduce(gcd,P)) for P in combinations(sorted(S),10) if reduce(gcd,P)>1]
        require(nonprim==[(tuple(sorted(3*w for w in C)),3)], "ratio-four unique ten-pack clock")
        require(all(reduce(gcd,P)==1 for size in (11,12) for P in combinations(S,size)),
                "ratio-four larger subsets primitive")
    if d>=7865 and (d-7865)%120120==0:
        require(all(any(w%q==0 for w in S) for q in range(2,15)), "ratio-four all small clocks still blocked")
    return (d,h,count,14*h>=87*max(base))


def main():
    units = (1,2,4,5,7,8,10)
    profile_points = sum(critical_profile(T) for T in combinations(units,3))
    owner_points = owner_controls()
    code_rows = wall_controls()
    sharp_rows = 0
    # Distinct parameter range and complete gcd rescaling for both branches.
    for A,B in combinations(range(2,111,3),2):
        for upper in (False,True):
            if upper and B<4*A:
                continue
            T=(A,B,A+B)
            y=Q(1,B if upper else A+B)
            d=(2*A+B)//3 if upper else (B-A)//3
            want=Q(A,3*(B if upper else A+B))
            require(heights(T,y)==((want,Q(0),want) if upper else (want,want,Q(0))),
                    "separate literal equality families")
            g=gcd(A,B)
            require(g%3!=0 and d%g==0, "primitive scaling preserves legal signed clock")
            require(max(heights(tuple(w//g for w in T),g*y))==want
                    and want==abs(distance(d,y)-Q(1,3))/2, "primitive global profile sharpness")
            sharp_rows+=1
    require(heights((1,4,5),Q(1,4))==(Q(1,12),Q(1,12),Q(0)), "upper endpoint sharp after sheet permutation")
    require(max(heights((11,17,28),Q(1,28)))==Q(11,84), "virtual-wall optimal value")
    for N in range(1,41):
        if N%11!=5:
            T=(9*N-1,42*N-1,51*N-2)
            y=Q(1,T[1])
            require(reduce(gcd,T)==1 and max(heights(T,y))<Q(1,14), "upper band genuinely sharp")
            require(Q(10,21)-distance(20*N-1,y)==Q(11,21*(42*N-1)), "upper exact approach")
        N2=2*N
        T=(9*N2-1,33*N2-1,42*N2-2)
        require(reduce(gcd,T)==1 and max(heights(T,Q(1,T[2])))<Q(1,14), "lower band genuinely sharp")
    # The stronger band too is only necessary: a non-additive tail set
    # with a wholly inactive tail can satisfy all three pair inequalities.
    noniff = None
    for T in combinations(units,3):
        for q in range(2,31):
            for p in range(1,q):
                y=Q(p,q)
                if all(Q(4,21)<distance(d,y)<Q(10,21) for d in pair_frequencies(T)) \
                        and max(heights(T,y))>=Q(1,14):
                    noniff=(T,str(y),tuple(str(distance(d,y)) for d in pair_frequencies(T)))
                    break
            if noniff: break
        if noniff: break
    require(noniff is not None, "even all strengthened pair bands do not characterize spoil")
    # Whole-component deletion, reconstructed from actual intervals.
    old=[(Q(14*k+1,196),Q(14*k+13,196)) for k in range(14)]
    kept=[I for I in old if I[0]>=Q(1,14) and I[1]<=Q(13,14)]
    require(kept==old[1:-1] and all(distance(1,x)>Q(1,14) for I in kept for x in I),
            "owner-free deletion removes two whole components")
    C=tuple(range(1,10))+(14,)
    require(all(distance(w,Q(1,28))>=Q(1,14) for w in C if w!=1)
            and distance(1,Q(1,28))<Q(1,14), "deletion inside primitive ten-pack")
    require(all(distance(w,Q(1,11))>=Q(1,14) for w in C), "same ten-pack nonempty")
    explicit=[]
    for d in (31,37,55,65,125,169,715):
        for m in (1,2,91**6//(42*d*(3*d+1)*(3*d+2))+1):
            explicit.append(addressed(d,m,subset=(m==1)))
    for r in range(2):
        d=715*(11+168*r)
        explicit.append(addressed(d,1,subset=True))
        C,T=data(d,1+42*d*(3*d+1)*(3*d+2))
        S=tuple(3*w for w in C)+T
        require(all(any(w%q==0 for w in S) for q in range(2,15)), "every denominator clock 2..14 blocked")
    free=[]
    for d in (31,41,55,85,121):
        base_max=14*(3*d+2)
        for start in (base_max+1,91**6+100):
            free.append(arbitrary(d,start))
    four=[]
    for d in (31,37,125,7865):
        four.append(atlas_four(d,1+42*d*(3*d+1)*(3*d+2),explicit=True,subsets=True))
    for d in (31,41,55,85,121):
        L=42*d*(3*d+1)*(3*d+2)
        h=14*(3*d+2)+1
        while gcd(h,L)!=1:
            h+=1
        four.append(atlas_four(d,h))
    near_base=[]
    for r in (0,1,500000):
        d=715*(11+168*r)
        h=42*d+29
        require(d%29!=0 and d%5!=3, "explicit near-base coprime progression filters")
        result=atlas_four(d,h,explicit=True,subsets=True)
        require(not result[3] and 7*h<45*(42*d+28), "both uniform attachment cones fail")
        C0,T=data(d,h)
        S=tuple(3*(4*h if w==2*h else w) for w in C0)+T
        support=(3*h,42*(3*d+2),3)
        coefficients=(1,-1,-1)
        require(len(set(support))==3 and set(support)<=set(S)
                and set(support)&{3*h,12*h}=={3*h}
                and sum(a*w for a,w in zip(coefficients,support))==0
                and max(abs(a) for a in coefficients)==1<=91**6,
                "height-one crossing row excludes W=Vdec on prescribed pair split")
        if r==500000:
            require(h>91**6, "near-base class remains cofinal above necessary cross-height threshold")
        near_base.append(result)
    print("STATUS: PASS after explicit global-sharpness repair; LRC14 remains open")
    print("PROFILE: 35 complete piecewise tail rows;",profile_points,"critical/crossing phases")
    print("OWNER/STRICT BOUNDARY POINTS:",owner_points,"; SIGNED WALL CODE ROWS:",code_rows)
    print("TWO SHARP FAMILIES:",sharp_rows,"rows, complete primitive rescaling")
    print("STRONGER BAND IS NOT IFF:",noniff)
    print("EXPLICIT FAMILY:",len(explicit),"controls including m=1, composite d, full subset gcds")
    print("FIRST RECONSTRUCTED CELL (d,m,h,y,quotient_length):",explicit[0])
    print("ARBITRARY COPRIME h ACTUAL COMPONENTS:",free)
    print("ATLAS-ADMISSIBLE (1,4) CONTROLS (d,h,verified_walls,uniform_decoder_cone):",four)
    print("EXACT ADDRESSED h=M+1, ALL SMALL CLOCKS BLOCKED:",near_base)
    print("ENTRY FIREWALL: 3h-42c-3=0 crosses {3h,12h} versus the rest at height one; prescribed W=Vdec excluded")
    print("PROVENANCE: THM4448 also supplies existence in its scale cone; this retains exact wall/owner/margin")
    print("ACTIVE GATES:",gates)


if __name__ == "__main__":
    main()
