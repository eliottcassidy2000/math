"""Exact controls for continuation-aware adverse budgets and shear zero-cost closure.

Reproduce: python 04-computation/cross_concepts_catalysis_20260908.py
Uses rational arithmetic and explicit exceptions (all checks survive -O).
The LRC density masks are auxiliary indicators, not canonical owner rows.
"""
from fractions import Fraction as Q
from itertools import product
from collections import Counter
from hashlib import sha256


def check(ok, label):
    if not ok:
        raise RuntimeError(label)


def add(x, y):
    return tuple(a+b for a, b in zip(x, y))


def neg(x):
    return tuple(-a for a in x)


def mean(x):
    return sum(x, Q(0))/len(x)


def center(x):
    a = mean(x)
    return tuple(t-a for t in x)


def energy(x):
    return mean(tuple(t*t for t in center(x)))


def hinge(r, x):
    return min(sum((max(Q(0), -a*(t-b)) for a,t in zip(r,x)), Q(0))/len(r)
               for b in set(x))


def transport(r, x):
    # Optional partial transport: cheapest positive-role source, most valuable
    # negative-role destination. No median, LP vertices, or hinge minimizer.
    sources = sorted([[x[i], r[i], i] for i in range(len(r)) if r[i]>0])
    sinks = sorted([[x[j], -r[j], j] for j in range(len(r)) if r[j]<0], reverse=True)
    i = j = 0
    value = Q(0)
    plan = []
    while i<len(sources) and j<len(sinks):
        low, a, u = sources[i]
        high, b, v = sinks[j]
        if high <= low:
            break
        mass = min(a,b)
        value += mass*(high-low)
        plan.append((u,v,mass))
        sources[i][1] -= mass
        sinks[j][1] -= mass
        if sources[i][1] == 0:
            i += 1
        if sinks[j][1] == 0:
            j += 1
    return value/len(r), plan


def scanner(r, x, T=Q(3)):
    P = [i for i,a in enumerate(r) if a>0]
    N = [i for i,a in enumerate(r) if a<0]
    i0, j0 = P[0], N[0]
    edges = [(i0,j) for j in N] + [(i,j0) for i in P[1:]]
    differences = {}
    for i,j in edges:
        k = tuple(-T if s==i else T if s==j else
                  2*T if r[s]>0 else -2*T if r[s]<0 else Q(0)
                  for s in range(len(r)))
        B = hinge(r, add(x,k))
        observed = len(r)*B/min(r[i],-r[j])-2*T
        check(observed == x[j]-x[i], "spanning-tree scanner edge")
        differences[i,j] = observed
    recovered = {i0: Q(0)}
    recovered.update({j: differences[i0,j] for j in N})
    recovered.update({i: recovered[j0]-differences[i,j0] for i in P[1:]})
    check(all(recovered[s] == x[s]-x[i0] for s in recovered), "scanner reconstruction")
    return len(edges)


def budget_census():
    total = probes = 0
    for p in (2,3,4):
        rows = list(product(map(Q,(-1,0,1)), repeat=p))
        local = 0
        for R in rows:
            if len(set(R))==1:
                continue
            r = center(R)
            P = [i for i,a in enumerate(r) if a>0]
            N = [i for i,a in enumerate(r) if a<0]
            for x in rows:
                B = hinge(r,x)
                Bminus = hinge(r,neg(x))
                BT, plan = transport(r,x)
                check(B==BT, "partial transport equals gauge hinge")
                zero_cone = max(x[j] for j in N)<=min(x[i] for i in P)
                check((B==0)==zero_cone, "zero cone iff")
                visible = [x[i] for i,a in enumerate(r) if a]
                check((B+Bminus==0)==(len(set(visible))==1), "symmetric nullspace")
                mad = min(mean(tuple(abs(a)*abs(t-b) for a,t in zip(r,x))) for b in set(x))
                check(B+Bminus==mad, "symmetrized median norm")
                check(B-Bminus==-mean(tuple(a*t for a,t in zip(r,x))), "antisymmetric response")
                probes += scanner(r,x)
                local += 1
        total += local
        print(f"BUDGET p={p} profile pairs={local}: transport/hinge, cone/nullspace, norm, scanner PASS")
    print(f"BUDGET total_pairs={total}; spanning_tree_scalar_queries={probes}")
    # Same present cost; a common additive continuation distinguishes them.
    R = (Q(-1),Q(0),Q(1))
    r = center(R)
    L0 = (Q(0),)*3
    L1 = R
    K = neg(R)
    check(hinge(r,L0)==hinge(r,L1)==0, "favorable current collision")
    check(hinge(r,add(L0,K))==Q(2,3) and hinge(r,add(L1,K))==0,
          "favorable common future splitter")
    # Removing invisible slots is valid only while the observer is fixed.
    invisible = (Q(0),Q(1),Q(0))
    r2 = (Q(-1),Q(1),Q(0))
    check(hinge(r,invisible)==hinge(r,neg(invisible))==0, "old observer null")
    check(hinge(r2,neg(invisible))==Q(1,3), "observer change revives null slot")
    print("CONTINUATION collision B_r(0)=B_r(r)=0; K=-r gives 2/3 versus 0 PASS")
    print("OBSERVER change r=(-1,0,1) -> (-1,1,0): old null slot becomes visible PASS")


def density_masks():
    p = 13
    R = tuple(Q(1,7) if s==0 else Q(6,91) if s in (1,12) else Q(0) for s in range(p))
    r = center(R)
    E = energy(R)
    budgets = []
    old = new = constant = 0
    profiles = []
    for mask in range(1<<12):
        L = tuple(Q(0) if s==0 or not(mask>>(s-1)&1) else Q(1,7)-R[s] for s in range(p))
        B = hinge(r,L)
        check(B==transport(r,L)[0], "density partial transport")
        a = int(bool(mask&1))+int(bool(mask&(1<<11)))
        m = (mask&((1<<11)-2)).bit_count()
        check(B == Q(min(325*m,1750-371*a+150*m),1399489), "two-count density compiler")
        C = add(R,L)
        check(energy(C)*E >= (B-E)**2, "density variance consequence")
        budgets.append(B)
        profiles.append(L)
        old += energy(L)<E
        new += B!=E
        constant += energy(C)==0
    # Moore refinement on the actual finite operation monoid mask -> mask OR bit.
    def classify(keys):
        table = {}
        return tuple(table.setdefault(key,len(table)) for key in keys)
    classes = classify(budgets)
    counts = [len(set(classes))]
    penultimate = None
    while True:
        refined = classify((classes[m], *(classes[m | (1<<s)] for s in range(12)))
                           for m in range(1<<12))
        if refined==classes:
            break
        penultimate = classes
        classes = refined
        counts.append(len(set(classes)))
    # Independent direct-context partition: no successor partition recursion.
    direct = classify(budgets)
    direct_counts = [len(set(direct))]
    contexts = sorted(range(1,1<<12),key=lambda c:(c.bit_count(),c))
    for depth in range(1,len(counts)):
        for context in contexts:
            if context.bit_count()==depth:
                direct = classify((direct[m],budgets[m | context]) for m in range(1<<12))
        direct_counts.append(len(set(direct)))
    check(direct_counts==counts, "direct-context and Moore depth classification")
    grouped = {}
    for m,c in enumerate(penultimate):
        grouped.setdefault(c,[]).append(m)
    pair = next(v[:2] for v in grouped.values() if len(v)>1)
    a6,b6 = pair
    shortest = next(c for c in contexts if budgets[a6 | c]!=budgets[b6 | c])
    check(shortest.bit_count()==len(counts)-1, "direct sharp depth witness")
    # The complement-of-one context separates every pair of distinct masks.
    full = (1<<12)-1
    check(budgets[full]==E and all(budgets[full^(1<<s)]!=E for s in range(12)),
          "analytic complete-future singleton separation")
    # Explicit same-cost state collision split by one lawful activation.
    witness = None
    for a in range(1<<12):
        for s in range(12):
            # Swapping among ten equal negative-role coordinates is currently
            # invisible; individual labelled activations need not preserve it.
            b = a ^ (1<<s)
            if budgets[a]==budgets[b]:
                for t in range(12):
                    if budgets[a | (1<<t)] != budgets[b | (1<<t)]:
                        witness = (a,b,t+1,budgets[a],budgets[a | (1<<t)],budgets[b | (1<<t)])
                        break
            if witness:
                break
        if witness:
            break
    if witness is None:
        # Equal cardinality masks, not adjacent masks, are the natural collisions.
        groups = {}
        for m,B in enumerate(budgets):
            groups.setdefault(B,[]).append(m)
        for group in groups.values():
            a = group[0]
            for b in group[1:]:
                for t in range(12):
                    if budgets[a | (1<<t)] != budgets[b | (1<<t)]:
                        witness=(a,b,t+1,budgets[a],budgets[a | (1<<t)],budgets[b | (1<<t)])
                        break
                if witness:
                    break
            if witness:
                break
    check(witness is not None, "lawful density one-step splitter exists")
    a,b,s,B,Ba,Bb = witness
    def active(mask):
        return tuple(i+1 for i in range(12) if mask>>i&1)
    semantic = "\n".join(f"{m}:{budgets[m]}:{classes[m]}" for m in range(1<<12))
    print(f"DENSITY universe=4096 labelled auxiliary activation masks; E={E}")
    print(f"DENSITY old_norm_gate={old}; exact_adverse_gate={new}; constant_targets={constant}")
    print(f"DENSITY present/future refinement class counts={counts}; final_class_sizes={dict(sorted(Counter(Counter(classes).values()).items()))}")
    print(f"DENSITY independent direct-context counts={direct_counts} PASS")
    print(f"DENSITY sharp-depth pair A={active(a6)}, B={active(b6)}, shortest context={active(shortest)}")
    print(f"DENSITY common-activation splitter A={active(a)}, B={active(b)}, activate={s}, current={B}, future=({Ba}, {Bb})")
    print(f"DENSITY semantic_sha256={sha256(semantic.encode()).hexdigest()}")


def matadd(A,B):
    return tuple(a+b for a,b in zip(A,B))


def scale(a,A):
    return tuple(a*x for x in A)


def matmul(A,B):
    return tuple(sum(A[3*i+k]*B[3*k+j] for k in range(3)) for i in range(3) for j in range(3))


def outer(q,p):
    return tuple(x*y for x in q for y in p)


ZERO = (Q(0),)*9


def is_shear(A):
    tr = A[0]+A[4]+A[8]
    minors = [A[3*i+j]*A[3*k+l]-A[3*i+l]*A[3*k+j]
              for i in range(3) for k in range(i+1,3)
              for j in range(3) for l in range(j+1,3)]
    return tr==0 and all(x==0 for x in minors) and matmul(A,A)==ZERO


def w(A):
    return (A[7]-A[5],A[2]-A[6],A[3]-A[1])


def alpha(A):
    v = w(A)
    den = sum(x*x for x in v)
    check(den>0, "nonzero vorticity")
    return Q(sum(v[i]*A[3*i+j]*v[j] for i in range(3) for j in range(3)),den)


def decomposition(M):
    S12 = (1,1,0,-1,-1,0,0,0,0)
    S23 = (0,0,0,0,1,1,0,-1,-1)
    pieces = [scale(M[0],S12),scale(M[0]+M[4],S23)]
    used = matadd(*pieces)
    for i in range(3):
        for j in range(3):
            if i!=j:
                pieces.append(tuple(M[k]-used[k] if k==3*i+j else Q(0) for k in range(9)))
    return pieces


def shear_census():
    count = 0
    for coords in product((-1,0,1),repeat=8):
        a,b,c,d,e,f,g,h = coords
        M = (a,b,c,d,e,f,g,h,-a-e)
        pieces = decomposition(M)
        check(len(pieces)==8 and all(is_shear(S) for S in pieces), "eight pure shear decomposition")
        total = ZERO
        for S in pieces:
            total = matadd(total,S)
        check(total==M, "shear reconstruction")
        count += 1
    S = (0,1,0,0,0,0,0,0,0)
    T = outer((1,1,1),(1,0,-1))
    M = matadd(S,T)
    check(is_shear(S) and is_shear(T) and not is_shear(M), "zero-cost sum leaves cone")
    check(w(M)==(1,-2,0) and alpha(M)==Q(-3,5), "two-shear stretching hostile")
    print(f"SHEAR universe=3^8={count} trace-free integer matrices; <=8 pure shear decomposition PASS")
    print(f"SHEAR hostile S={S}; T={T}; w(S+T)={w(M)}; alpha(S+T)={alpha(M)}")
    print("SHEAR consequence: d_F(S)=d_F(T)=0 but d_F(S+T)>=3*sqrt(6/7)/5>0")


if __name__ == '__main__':
    budget_census()
    density_masks()
    shear_census()
    print("PASS: exact finite controls; no canonical LRC owner or Euler trajectory realization claimed")
