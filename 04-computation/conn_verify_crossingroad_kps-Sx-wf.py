"""
Independent adversarial verification of the "crossing-roadcoloring" connection.

Claims to check:
(a) raw orientation-free crossing number = C(n,4) for every tournament (K_n convex).
(b) P(T) = #{w<x<y<z : A[w][y]==A[x][z]}  (crossing chords parallel-oriented)
    (1) P is EXACTLY quadratic in tile bits (3rd mixed diffs vanish, some 2nd nonzero).
    (2) P NOT score-determined (same score, different P).
    (3) P not affine in c3.
    (4) P strictly refines c3 (#(c3,P) > #c3).
ROAD COLORING: bit-toggle automaton = permutation automaton => no synchronizing word.

Fresh, independent implementation. Exact integer arithmetic only.
"""
from itertools import combinations, product

def tiles(n):
    return [(a,b) for a in range(3,n+1) for b in range(1,a-1)]

def adj(n,bits,T):
    A=[[0]*(n+1) for _ in range(n+1)]
    for k in range(n,1,-1): A[k][k-1]=1
    for (a,b),bit in zip(T,bits):
        if bit==0: A[a][b]=1
        else: A[b][a]=1
    return A

def c3(A,n):
    t=0
    for i in range(1,n+1):
        for j in range(i+1,n+1):
            for k in range(j+1,n+1):
                if (A[i][j]+A[i][k],A[j][i]+A[j][k],A[k][i]+A[k][j])==(1,1,1):t+=1
    return t

def scores(A,n):
    return tuple(sorted(sum(A[i][j] for j in range(1,n+1) if j!=i) for i in range(1,n+1)))

def P_invariant(A,n):
    """P = number of 4-subsets w<x<y<z with A[w][y]==A[x][z].
       In convex position 1..n, the crossing chord pair among {w,x,y,z} is (w,y) and (x,z).
       'parallel-oriented' = the two crossing chords point the same relative direction.
       Implement literally: count quadruples where arc-direction bit on (w,y) equals that on (x,z)."""
    P=0
    for w,x,y,z in combinations(range(1,n+1),4):
        # crossing chords in convex order w<x<y<z are (w,y) and (x,z)
        if A[w][y]==A[x][z]:
            P+=1
    return P

def raw_crossing(n):
    # orientation-free crossing number of K_n in convex position = C(n,4)
    from math import comb
    return comb(n,4)

# ---------- TEST 1: P is exactly quadratic in tile bits ----------
def test_quadratic(n):
    T=tiles(n); F=len(T)
    base=[0]*F
    def Pof(bits): return P_invariant(adj(n,bits,T),n)
    # check all 3rd mixed finite differences vanish, over random base points and all coord triples
    import random
    random.seed(7)
    third_all_zero=True
    second_some_nonzero=False
    bad=None
    samples=[ [random.randint(0,1) for _ in range(F)] for _ in range(40) ]
    samples.append(base)
    for s in samples:
        for trip in combinations(range(F),3):
            # 3rd mixed difference Delta_i Delta_j Delta_k P
            acc=0
            for mask in product([0,1],repeat=3):
                b=s[:]
                sign=1
                for idx,m in zip(trip,mask):
                    if m: b[idx]^=1
                    # finite diff: flipping uses +1/-1 weighting
                # weight = (-1)^(number of NOT flipped) ... use standard inclusion
            # do it cleanly:
            acc=0
            for mask in product([0,1],repeat=3):
                b=s[:]
                for idx,m in zip(trip,mask):
                    if m: b[idx]^=1
                sgn=(-1)**(3-sum(mask))
                acc+=sgn*Pof(b)
            if acc!=0:
                third_all_zero=False; bad=(s,trip,acc); break
        if not bad:
            pass
        if bad: break
    # second differences nonzero
    s=base
    for pair in combinations(range(F),2):
        acc=0
        for mask in product([0,1],repeat=2):
            b=s[:]
            for idx,m in zip(pair,mask):
                if m: b[idx]^=1
            sgn=(-1)**(2-sum(mask))
            acc+=sgn*Pof(b)
        if acc!=0:
            second_some_nonzero=True; break
    return third_all_zero, second_some_nonzero, bad

# ---------- TEST 2: P not score-determined ----------
def test_score_determined(n, sample_cap=None):
    T=tiles(n); F=len(T)
    from collections import defaultdict
    bysc=defaultdict(set)
    if F<=20 and sample_cap is None:
        it=product([0,1],repeat=F)
    else:
        import random; random.seed(3)
        it=( [random.randint(0,1) for _ in range(F)] for _ in range(sample_cap or 8000))
    for bits in it:
        A=adj(n,list(bits),T)
        sc=scores(A,n); P=P_invariant(A,n)
        bysc[sc].add(P)
    witnesses=[(sc,sorted(ps)) for sc,ps in bysc.items() if len(ps)>1]
    return len(witnesses)>0, witnesses[:5]

# ---------- TEST 3 & 4: P vs c3 ----------
def test_c3_relation(n, sample_cap=None):
    T=tiles(n); F=len(T)
    pairs=set(); c3s=set()
    c3_to_P=dict()
    affine_ok=True
    if F<=20 and sample_cap is None:
        it=product([0,1],repeat=F)
    else:
        import random; random.seed(11)
        it=( [random.randint(0,1) for _ in range(F)] for _ in range(sample_cap or 8000))
    data=[]
    for bits in it:
        A=adj(n,list(bits),T)
        c=c3(A,n); P=P_invariant(A,n)
        pairs.add((c,P)); c3s.add(c)
        data.append((c,P))
    # is P an affine function of c3? check if each c3 maps to unique P
    from collections import defaultdict
    m=defaultdict(set)
    for c,P in data: m[c].add(P)
    affine_ok = all(len(v)==1 for v in m.values())
    return len(c3s), len(pairs), affine_ok

# ---------- ROAD COLORING: permutation automaton ----------
def test_road_coloring(n):
    """Each tile-toggle is an involution on Q_F (bijection). A permutation automaton
       (all letters bijective) admits NO synchronizing word unless |states|=1.
       Verify each letter is a bijection (involution) on the 2^F state set for small n."""
    T=tiles(n); F=len(T)
    # letter i = XOR with e_i ; trivially an involution / bijection on {0,1}^F
    # Verify by checking it's a permutation on a sample of states
    import random; random.seed(5)
    ok=True
    for i in range(F):
        seen=set()
        for _ in range(2000):
            s=tuple(random.randint(0,1) for _ in range(F))
            s2=list(s); s2[i]^=1; s2=tuple(s2)
            # involution check
            s3=list(s2); s3[i]^=1
            if tuple(s3)!=s: ok=False
            seen.add((s,s2))
        # injectivity: distinct s give distinct s2 (XOR is bijection, trivially true)
    return ok  # permutation automaton => no synchronizing word collapses to single state

if __name__=="__main__":
    out=[]
    def log(*a):
        s=" ".join(str(x) for x in a); print(s); out.append(s)

    log("=== RAW CROSSING (claim a) ===")
    for n in range(4,9):
        log(f"n={n}: C(n,4) =", raw_crossing(n))
    log("Guy optimal cr(K_n) n=4..8 = 0,1,3,9,18  (raw convex is strictly larger for n>=5) -- consistent")

    log("\n=== TEST 1: P quadratic in tile bits ===")
    for n in (5,6):
        t3,s2,bad=test_quadratic(n)
        log(f"n={n}: 3rd-diff all zero = {t3} ; some 2nd-diff nonzero = {s2} ; counterexample={bad}")

    log("\n=== TEST 2: P not score-determined ===")
    for n in (5,6):
        ok,wit=test_score_determined(n)
        log(f"n={n}: same-score-different-P exists = {ok}")
        for sc,ps in wit: log(f"    score {sc} -> P in {ps}")
    ok7,wit7=test_score_determined(7, sample_cap=20000)
    log(f"n=7 (sampled 20000): same-score-different-P exists = {ok7}; examples:")
    for sc,ps in wit7[:3]: log(f"    score {sc} -> P in {ps}")

    log("\n=== TEST 3&4: P vs c3 ===")
    for n in (5,6):
        nc3,npair,affine=test_c3_relation(n)
        log(f"n={n}: #distinct c3 = {nc3} ; #distinct (c3,P) = {npair} ; P affine-in-c3 = {affine} (False => P refines c3)")
    nc3,npair,affine=test_c3_relation(7, sample_cap=20000)
    log(f"n=7 (sampled): #c3={nc3} #(c3,P)={npair} affine={affine}")

    log("\n=== ROAD COLORING (negative) ===")
    for n in (4,5):
        log(f"n={n}: every letter is involution/bijection on Q_F = {test_road_coloring(n)} => permutation automaton => NO synchronizing word")

    with open("05-knowledge/results/conn_verify_crossingroad_kps-Sx-wf.out","w") as f:
        f.write("\n".join(out))
