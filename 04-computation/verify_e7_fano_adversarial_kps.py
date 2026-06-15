"""
verify_e7_fano_adversarial_kps.py
ADVERSARIAL re-verification of the e7-fano-7 worker claims.
We actually COMPUTE everything we can, rather than printing hardcoded numbers.

Checks:
  A. OCF: enumerate ALL tournaments n<=6, compute H by directed Ham path count AND by
     OCF (I(Omega,2)), confirm they agree, confirm 7 and 21 are absent from achievable H.
  B. Verify Omega=K3 non-realizability directly (3 odd cycles pairwise sharing a vertex,
     no independent pair) over all tournaments n<=6.
  C. Mayer / cluster check: b_2 = -|E(Omega)| ; log H vs sum b_k 2^k on sample tournaments.
  D. E7 arithmetic: rank, exponents, degrees, roots, Weyl order, 56 decompositions, 21=C(7,2).
  E. T(6) count = 56 and 35 strong + 21 non-strong, independently from enumeration.
"""

from itertools import combinations, permutations
from math import comb

# ---------------------------------------------------------------------------
# Tournament utilities. A tournament on n vertices: adjacency matrix A[i][j]=1
# iff arc i->j. For each pair exactly one direction.
# We enumerate tournaments by choosing orientation of each of C(n,2) pairs.
# ---------------------------------------------------------------------------

def tournaments(n):
    pairs = list(combinations(range(n), 2))
    m = len(pairs)
    for bits in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A

def ham_path_count(A, n):
    """Number of directed Hamiltonian paths = H(T)."""
    cnt = 0
    for perm in permutations(range(n)):
        ok = True
        for a in range(n-1):
            if not A[perm[a]][perm[a+1]]:
                ok = False
                break
        if ok:
            cnt += 1
    return cnt

def directed_odd_cycles(A, n):
    """Return list of directed odd cycles as frozensets of vertices.
    A directed cycle on vertex subset following a cyclic order.
    We enumerate all simple directed cycles of ODD length >=3.
    Represent each cycle canonically by (length, frozenset?) -- but two distinct
    cycles can share the same vertex set. The conflict graph Omega has a vertex
    PER directed odd cycle (distinct cyclic orderings count separately).
    We must enumerate distinct directed cycles (as cyclic sequences up to rotation,
    NOT reflection since direction matters)."""
    cycles = []
    # for each subset of size L (odd), enumerate cyclic orderings that are directed cycles
    for L in range(3, n+1, 2):
        for verts in combinations(range(n), L):
            # enumerate cyclic permutations of verts; fix first element to kill rotation
            rest = verts[1:]
            for p in permutations(rest):
                seq = (verts[0],) + p
                # check directed cycle seq[0]->seq[1]->...->seq[L-1]->seq[0]
                ok = True
                for a in range(L):
                    if not A[seq[a]][seq[(a+1)%L]]:
                        ok = False
                        break
                if ok:
                    cycles.append((seq, frozenset(verts)))
    return cycles

def conflict_graph(cycles):
    """Omega: vertices = directed odd cycles; edge iff share >=1 vertex."""
    nC = len(cycles)
    edges = 0
    adj = [set() for _ in range(nC)]
    for i in range(nC):
        for j in range(i+1, nC):
            if cycles[i][1] & cycles[j][1]:
                adj[i].add(j); adj[j].add(i)
                edges += 1
    return adj, edges

def independence_poly_at_2(adj, nC):
    """I(Omega,2) = sum_k alpha_k 2^k, alpha_k = # independent sets of size k.
    Compute via counting independent sets weighted by 2^size. We do this by
    DP/recursion over vertices (works fine for small Omega)."""
    # weighted independent set sum: product expansion. Use recursion with memo on remaining set.
    # For small nC just enumerate via Bron-Kerbosch-style independent set counting.
    # alpha_k counts: count independent sets of each size.
    from functools import lru_cache
    adj_bits = [0]*nC
    for i in range(nC):
        for j in adj[i]:
            adj_bits[i] |= (1<<j)
    # enumerate maximal? no -- ALL independent sets. Use recursion.
    total = [0]
    def rec(cand, weight):
        # cand = bitmask of vertices still allowed; weight = 2^(size chosen so far)
        total[0] += weight
        c = cand
        while c:
            v = (c & -c).bit_length() - 1
            c &= c-1
            # choose v: remove v and its neighbors from cand, but only vertices > v to avoid double count
            newcand = cand & ~adj_bits[v] & ~((1<<(v+1))-1)
            rec(newcand, weight*2)
    full = (1<<nC) - 1
    rec(full, 1)
    return total[0]

def H_via_OCF(A, n):
    cycles = directed_odd_cycles(A, n)
    adj, edges = conflict_graph(cycles)
    nC = len(cycles)
    return independence_poly_at_2(adj, nC), nC, edges

# ---------------------------------------------------------------------------
# CHECK A + B: OCF agreement, forbidden 7/21, Omega=K3 non-realizability
# ---------------------------------------------------------------------------
print("="*70)
print("CHECK A: OCF (I(Omega,2)) == directed-Ham-path count, and achievable H set")
print("="*70)

for n in range(2, 7):
    achievable = set()
    mismatches = 0
    k3_found = 0
    count = 0
    for A in tournaments(n):
        count += 1
        Hp = ham_path_count(A, n)
        Hocf, nC, edges = H_via_OCF(A, n)
        if Hp != Hocf:
            mismatches += 1
        achievable.add(Hp)
        # K3 detection: Omega contains an induced K3 with NO independent pair among
        # those 3 AND globally alpha_2 contribution... simplest faithful check of THM-029
        # mechanism: does there exist a tournament with EXACTLY alpha=(1, a1=3, a2=0)
        # i.e. Omega itself == K3 (3 odd cycles, all pairwise sharing, none independent)?
        if nC == 3 and edges == 3:
            k3_found += 1
    print(f"n={n}: #tournaments={count}, OCF-vs-Ham mismatches={mismatches}, "
          f"Omega==K3 cases={k3_found}")
    print(f"      achievable H (sorted) = {sorted(achievable)}")
    print(f"      7 in achievable? {7 in achievable}   21 in achievable? {21 in achievable}")

print()
print("="*70)
print("CHECK B: direct non-realizability of Omega=K3 (3 pairwise-sharing odd cycles,")
print("         no independent pair) -- THM-029 mechanism, all tournaments n<=6")
print("="*70)
for n in range(3, 7):
    found = 0
    total = 0
    for A in tournaments(n):
        total += 1
        cycles = directed_odd_cycles(A, n)
        adj, edges = conflict_graph(cycles)
        nC = len(cycles)
        # look for an induced K3 in Omega such that those 3 cycles have no
        # OTHER structure forcing alpha_1>3 or alpha_2>0... but the cleanest THM-029
        # statement: is there a tournament whose Omega is EXACTLY K3?
        if nC == 3 and edges == 3:
            found += 1
    print(f"n={n}: tournaments={total}, Omega-exactly-K3 = {found}")

print()
print("="*70)
print("CHECK C: Mayer / cluster expansion sanity")
print("  Claim: b_2 = -|E(Omega)|.  And log H ?= sum_k b_k 2^k (cumulants).")
print("="*70)
import math
# Take a few specific tournaments at n=5 and inspect.
# We compute alpha_k explicitly, then the standard Mayer cluster integrals b_k for a
# hard-core gas: log Z = sum_{l>=1} b_l z^l where Z = sum_k alpha_k z^k (alpha_0=1).
def alphas(adj, nC):
    adj_bits = [0]*nC
    for i in range(nC):
        for j in adj[i]:
            adj_bits[i] |= (1<<j)
    res = {}
    def rec(cand, size):
        res[size] = res.get(size,0)+1
        c = cand
        while c:
            v = (c & -c).bit_length()-1
            c &= c-1
            newcand = cand & ~adj_bits[v] & ~((1<<(v+1))-1)
            rec(newcand, size+1)
    rec((1<<nC)-1, 0)
    return res  # res[k] = alpha_k

def cluster_b(alpha, kmax):
    """Given Z = sum_k alpha_k z^k (alpha[0]=1), compute b_l so that
    log Z = sum_{l>=1} b_l z^l, by formal power series log."""
    # build coefficient list of Z up to kmax
    a = [alpha.get(k,0) for k in range(kmax+1)]
    a[0] = 1
    # log of power series with constant term 1
    b = [0.0]*(kmax+1)
    # Use recurrence: if L = log A, A=1+sum a_k z^k, then A' = A L', match coefficients.
    # L'_coeffs: l*b_l. A' coeff at z^{n-1} = n*a_n.
    # (A L')_m = sum_{i} a_i * ( (m-i+1)*b_{m-i+1} )? easier: standard formula.
    # b_l = a_l - (1/l) sum_{k=1}^{l-1} k b_k a_{l-k}
    for l in range(1, kmax+1):
        s = a[l]
        for k in range(1, l):
            s -= (k * b[k] * a[l-k]) / l
        b[l] = s
    return b

# pick tournaments at n=5 with varied Omega. Require alpha to reach degree>=2 so b[2] exists.
samples = []
seen_edges = set()
for A in tournaments(5):
    cycles = directed_odd_cycles(A,5)
    adj, edges = conflict_graph(cycles)
    nC = len(cycles)
    key = (nC, edges)
    if key not in seen_edges and nC>0:
        seen_edges.add(key)
        al = alphas(adj, nC)
        kmax = max(2, max(al.keys()))   # ensure we go to at least degree 2
        H = sum(al.get(k,0)*(2**k) for k in range(kmax+1))
        b = cluster_b(al, kmax)
        logH = math.log(H)
        sum_b2 = sum(b[l]*(2**l) for l in range(1, kmax+1))
        samples.append((nC, edges, al, H, b[2], -edges, logH, sum_b2))
    if len(samples) >= 10:
        break

print(f"{'nC':>3} {'|E|':>4} {'alpha_k':>22} {'H':>6} {'b2':>8} {'-|E|':>5} {'logH':>8} {'sum_b_l*2^l':>12}")
for nC, edges, al, H, b2, negE, logH, sumb in samples:
    alstr = ",".join(str(al.get(k,0)) for k in range(max(al.keys())+1))
    print(f"{nC:>3} {edges:>4} {alstr:>22} {H:>6} {b2:>8.2f} {negE:>5} {logH:>8.4f} {sumb:>9.4f}")
print()
print("Interpretation: b_2 (formal-log coeff at z^2) should equal alpha_2 - alpha_1^2/2.")
print("The CLAIM 'b_2 = -|E(Omega)|' holds ONLY if alpha_2 = C(alpha_1,2) - |E|, i.e.")
print("alpha_1^2/2 - |E| - alpha_1^2/2 = -|E| - alpha_1/2 ... check below numerically.")
for nC, edges, al, H, b2, negE, logH, sumb in samples:
    a1 = al.get(1,0); a2 = al.get(2,0)
    pred_b2 = a2 - a1*a1/2.0
    print(f"  nC={nC} a1={a1} a2={a2} |E|={edges}: b2={b2:.3f}  a2-a1^2/2={pred_b2:.3f}  "
          f"-|E|={-edges}  (a2 = C(a1,2)-|E|? {a2 == comb(a1,2)-edges})")

print()
print("="*70)
print("CHECK D: E7 arithmetic facts")
print("="*70)
def factor(n):
    f={}; d=2
    while d*d<=n:
        while n%d==0: f[d]=f.get(d,0)+1; n//=d
        d+=1
    if n>1: f[n]=f.get(n,0)+1
    return f

dim_E7=133; rank=7; roots=126; pos=63
exps=[1,5,7,9,11,13,17]; degs=[2,6,8,10,12,14,18]; h=18; weyl=2903040
print(f"dim 133 = 7*19? {133==7*19}  factor={factor(133)}")
print(f"roots 126 = 18*7? {126==18*7}; pos 63=9*7? {63==9*7}=C(9,2)? {63==comb(9,2)}")
print(f"7 in exponents? {7 in exps}; exponents sum={sum(exps)} (=63? {sum(exps)==63})")
print(f"  exponent pairs m+m'=h=18: {[(e,h-e) for e in exps]}")
print(f"  is each complement also an exponent? {all((h-e) in exps for e in exps)}")
print(f"  central exponent (self-paired, =h/2=9): {h//2} -> in exps {h//2 in exps}; "
      f"is 7 self-paired? {18-7==7} (NO, 7 pairs with 11)")
print(f"degrees = exponents+1? {degs == [e+1 for e in exps]}; degree 8=7+1 present? {8 in degs}")
print(f"|Weyl(E7)| factor = {factor(weyl)} (expect 2^10 3^4 5 7); 7^2|weyl? {weyl%49==0}")
print(f"56 = 2*28 = {2*28}? {56==2*28}; 56 = 7*8 = {7*8}? {56==7*8}")
print(f"C(7,2)=21={comb(7,2)}; C(7,3)=35={comb(7,3)}; 7+21+35={7+21+35}=63? {7+21+35==63}")
print(f"21 = 14+7 (dim G2 + vector)? {14+7==21}")
print(f"NOTE: 7+21+35=63 is the # POSITIVE ROOTS, not 56. 56 != 63.")

print()
print("="*70)
print("CHECK E: T(6) = number of 6-vertex tournament iso classes, and strong/non-strong split")
print("="*70)
# A000568: number of tournaments on n nodes: 1,1,2,4,12,56,456,...
# We verify T(6)=56 by counting iso classes via canonical form (brute over n=6 is 2^15=32768).
def canon(A, n):
    best=None
    for perm in permutations(range(n)):
        bits=0; k=0
        for i in range(n):
            for j in range(i+1,n):
                pi,pj=perm[i],perm[j]
                if A[pi][pj]: bits |= (1<<k)
                k+=1
        if best is None or bits<best: best=bits
    return best

def is_strong(A,n):
    # strongly connected: every vertex reachable from every other
    def reach(s):
        seen={s}; stack=[s]
        while stack:
            u=stack.pop()
            for v in range(n):
                if A[u][v] and v not in seen:
                    seen.add(v); stack.append(v)
        return seen
    for s in range(n):
        if len(reach(s))!=n: return False
    return True

for n in [5,6]:
    classes={}
    strong_classes=set(); all_classes=set()
    for A in tournaments(n):
        c=canon(A,n)
        all_classes.add(c)
        if c not in classes:
            classes[c]=is_strong(A,n)
    nstrong=sum(1 for v in classes.values() if v)
    print(f"n={n}: iso classes T({n}) = {len(all_classes)} ; strong = {nstrong} ; "
          f"non-strong = {len(all_classes)-nstrong}")
    if n==6:
        print(f"  Worker claim: 56 = 35 strong + 21 non-strong. "
              f"strong={nstrong}, nonstrong={len(all_classes)-nstrong}")
