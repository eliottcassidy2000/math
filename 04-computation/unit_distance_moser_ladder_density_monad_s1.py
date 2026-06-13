"""
monad-explorer-2026-06-10-S1  (deep-research)
=============================================
DENSEST UNIT-DISTANCE PATCH PER MOSER-LADDER RUNG  L_t.

Background.  THM-434 proves the Moser-ladder bridge lattice
    L_t = Z[zeta6] (+) Z[zeta6].w_t,     w_t = ((2t-1) + i sqrt(4t-1))/(2t)
has EXACTLY  12 + r_E(t)  unit vectors, where r_E(t) = 6*B(t),
B(t) = sum_{d|t} chi_{-3}(d) = #ideals of norm t in Z[zeta6].
So t=3,4,9,12,16 -> 18 units;  t=13,21,28 -> 24 units;  t=49 -> 30 units; ...

THM-434's HONEST SCOPE NOTE: "This counts unit VECTORS (max possible vertex
degree), not the densest unit-distance graph the lattice supports -- a separate
(harder) question."  EVERY repo construction (u(21)=57, Engel's table) uses ONLY
the t=3 Moser lattice.  THIS SCRIPT attacks the open question:

   Does a HIGHER ladder rung (more unit vectors) support a DENSER unit-distance
   patch than t=3 at some N?  Or is t=3 special for DENSITY, not just count?

Everything EXACT: adjacency decided by |z|^2 == 1 over Q(sqrt3, sqrt(4t-1)),
Fraction coeffs in basis (1, sqrt3, sqrt d, sqrt 3d), d=4t-1.  No float decides
any adjacency.  Unit vectors are CONSTRUCTED from theory (3 families), not box-
searched, so the unit set is provably complete & exact for every t.

NEW STRUCTURAL INVARIANT computed here:  tau(t) = # unit triangles through 0
(pairs of unit vectors whose difference is also a unit vector) = local 3-clique
count = the true density-relevant rosette invariant (THM-434 gave only the count).
"""
from fractions import Fraction as F
from itertools import combinations

# ------------------------------------------------------------------ exact field
def make_field(t):
    """Return (units, mul, normsq, d) for ladder rung t. d=4t-1.
       Lattice basis (a,b,c,d_)  <->  a + b*w1 + c*w_t + d_*w1*w_t,  w1=zeta6.
       Value field basis: (1, sqrt3, sqrtD, sqrt3D), D=4t-1."""
    D = 4*t - 1
    # multiplication table on value-basis {1, sqrt3, sqrtD, sqrt3D}: (i,j)->(coeff,idx)
    MT = {
        (0,0):(1,0),(0,1):(1,1),(0,2):(1,2),(0,3):(1,3),
        (1,1):(3,0),(1,2):(1,3),(1,3):(3,2),
        (2,2):(D,0),(2,3):(D,1),(3,3):(3*D,0),
    }
    def mul(x, y):
        r = [F(0)]*4
        for i in range(4):
            if x[i]==0: continue
            for j in range(4):
                if y[j]==0: continue
                a,b = (i,j) if i<=j else (j,i)
                coeff, idx = MT[(a,b)]
                r[idx] += x[i]*y[j]*coeff
        return tuple(r)
    Z4=(F(0),)*4
    def add(x,y): return tuple(x[i]+y[i] for i in range(4))
    def smul(s,x): return tuple(s*x[i] for i in range(4))
    # complex coords (re, im) of the 4 lattice generators, in value basis
    half=F(1,2)
    twotm1=F(2*t-1)
    RE = {
        0:(F(1),F(0),F(0),F(0)),                       # 1
        1:(half,F(0),F(0),F(0)),                       # w1
        2:(twotm1/(2*t),F(0),F(0),F(0)),               # w_t
        3:(twotm1/(4*t),F(0),F(0),F(-1,1)*F(1,4*t)),   # w1 w_t : (2t-1)/4t on 1, -1/4t on sqrt3D
    }
    IM = {
        0:(F(0),F(0),F(0),F(0)),
        1:(F(0),half,F(0),F(0)),                       # w1: 1/2 on sqrt3
        2:(F(0),F(0),F(1,2*t),F(0)),                   # w_t: 1/2t on sqrtD
        3:(F(0),twotm1/(4*t),F(1,4*t),F(0)),           # w1 w_t: (2t-1)/4t on sqrt3, 1/4t on sqrtD
    }
    def coord(v):
        re=Z4; im=Z4
        for k in range(4):
            if v[k]==0: continue
            re=add(re,smul(F(v[k]),RE[k]))
            im=add(im,smul(F(v[k]),IM[k]))
        return re,im
    ONE=(F(1),F(0),F(0),F(0))
    def normsq(v):
        re,im=coord(v)
        return add(mul(re,re),mul(im,im))
    # ---- CONSTRUCT the unit vectors from theory (3 families) ----
    # triangular rosette zeta6^j
    tri = [(1,0,0,0),(0,1,0,0),(-1,1,0,0),(-1,0,0,0),(0,-1,0,0),(1,-1,0,0)]
    # w_t rosette zeta6^j * w_t
    wt  = [(0,0,1,0),(0,0,0,1),(0,0,-1,1),(0,0,-1,0),(0,0,0,-1),(0,0,1,-1)]
    # transverse: (p,q,-p,-q) for p^2+pq+q^2 = t
    trans=[]
    B=int(t**0.5)+2
    for p in range(-2*B,2*B+1):
        for q in range(-2*B,2*B+1):
            if p*p+p*q+q*q==t:
                trans.append((p,q,-p,-q))
    units = tri+wt+trans
    # sanity: every constructed vector is genuinely a unit vector
    for u in units:
        assert normsq(u)==ONE, (t,u,normsq(u))
    # no duplicates
    assert len(set(units))==len(units), t
    return units, normsq, D, len(tri), len(wt), len(trans)

# ------------------------------------------------------------------ graph tools
def U_count(points, units):
    s=set(points); e=0
    for p in points:
        for u in units:
            q=(p[0]+u[0],p[1]+u[1],p[2]+u[2],p[3]+u[3])
            if q in s: e+=1
    assert e%2==0
    return e//2

def degree(p, S, units):
    return sum(1 for u in units if (p[0]+u[0],p[1]+u[1],p[2]+u[2],p[3]+u[3]) in S)

def tau(units):
    """# unit triangles through 0 = #pairs {u,v} of units with u-v a unit vector."""
    US=set(units); c=0
    for u,v in combinations(units,2):
        d=(u[0]-v[0],u[1]-v[1],u[2]-v[2],u[3]-v[3])
        if d in US: c+=1
    return c

# ------------------------------------------------------------------ patch search
def grow_greedy(seed, target, units):
    """greedy: repeatedly add the candidate (a neighbour of current set) with the
       most neighbours already in S; tie-break toward the origin (compact)."""
    S=set([seed])
    while len(S)<target:
        cand={}
        for p in S:
            for u in units:
                q=(p[0]+u[0],p[1]+u[1],p[2]+u[2],p[3]+u[3])
                if q in S: continue
                if q in cand: continue
                cand[q]=degree(q,S,units)
        if not cand: break
        best=max(cand.items(), key=lambda kv:(kv[1], -(abs(kv[0][0])+abs(kv[0][1])+abs(kv[0][2])+abs(kv[0][3]))))
        S.add(best[0])
    return S

def local_improve(S, units, rounds=400):
    """remove the lowest-degree vertex, regrow a better one; keep if U improves.
       deterministic hill-climb over boundary swaps."""
    S=set(S); best=U_count(S,units)
    import itertools
    for _ in range(rounds):
        # candidate external points: neighbours of S not in S
        ext={}
        for p in S:
            for u in units:
                q=(p[0]+u[0],p[1]+u[1],p[2]+u[2],p[3]+u[3])
                if q not in S:
                    ext[q]=ext.get(q,0)+1
        if not ext: break
        # best external by gain
        addq=max(ext.items(), key=lambda kv:kv[1])
        # worst internal by degree
        degs={p:degree(p,S,units) for p in S}
        delp=min(degs.items(), key=lambda kv:kv[1])
        # try swap delp -> addq  (only if addq not already adjacent-counted via delp)
        if addq[0]==delp[0]: break
        newS=set(S); newS.discard(delp[0]); newS.add(addq[0])
        if len(newS)!=len(S):
            continue
        nu=U_count(newS,units)
        if nu>best:
            S=newS; best=nu
        else:
            break
    return S, best

def best_patch(units, target, seeds):
    best=-1; bestS=None
    for seed in seeds:
        S=grow_greedy(seed,target,units)
        if len(S)<target: continue
        u0=U_count(S,units)
        S2,u2=local_improve(S,units)
        u=max(u0,u2)
        if u>best:
            best=u; bestS=(S if u0>=u2 else S2)
    return best,bestS

# ------------------------------------------------------------------ main
if __name__=="__main__":
    HUB=(0,0,0,0)
    # known Engel/AMP records (t=3 Moser lattice) for calibration & comparison
    KNOWN={9:18,12:27,16:41,21:57,24:72,27:81,28:85,30:93}
    RUNGS=[3,4,9,12,13,16,21,28,49]   # valid t (d=4t-1 not 3*square); excl t=7,19
    TARGETS=[9,12,16,21,24,27,28,30]

    print("="*78)
    print("MOSER-LADDER RUNG STRUCTURE: unit-vector count and unit-triangle count tau")
    print("="*78)
    print(f"{'t':>3} {'D=4t-1':>7} {'#units':>7} {'(tri,wt,trans)':>16} {'tau(t)':>7}  note")
    rung_data={}
    for t in RUNGS:
        units,normsq,D,ntri,nwt,ntr = make_field(t)
        ta=tau(units)
        rung_data[t]=(units,ta,len(units))
        note=""
        if t==3: note="<- Moser lattice (all repo constructions)"
        if len(units)==24: note="24-unit rung (split prime =1 mod3)"
        if len(units)==30: note="30-unit RECORD rung (t=49=7^2)"
        print(f"{t:>3} {D:>7} {len(units):>7} {('('+str(ntri)+','+str(nwt)+','+str(ntr)+')'):>16} {ta:>7}  {note}")

    print()
    print("="*78)
    print("DENSEST PATCH PER RUNG (greedy multistart + hill-climb; EXACT recount)")
    print("Same search budget every rung. KNOWN = Engel/AMP best (t=3).")
    print("="*78)
    hdr="  N |"+"".join(f"{('t='+str(t)):>7}" for t in RUNGS)+" | KNOWN"
    print(hdr); print("-"*len(hdr))
    # seeds: hub plus a few short unit-vector walks (compact starts)
    for N in TARGETS:
        row={}
        for t in RUNGS:
            units=rung_data[t][0]
            seeds=[HUB]+[ (units[i]) for i in range(min(6,len(units))) ]
            b,_=best_patch(units,N,seeds)
            row[t]=b
        line=f"{N:>3} |"+"".join(f"{row[t]:>7}" for t in RUNGS)+f" | {KNOWN.get(N,'?')}"
        # mark any rung that meets/beats KNOWN
        print(line)
    print()
    print("Reading: if any t>3 column EXCEEDS the t=3 column (and KNOWN) at some N,")
    print("that rung supports a denser patch -> potential new record. If t=3 dominates,")
    print("the Moser lattice is special for DENSITY, not merely unit-vector count.")
    print("DONE.")
