#!/usr/bin/env python3
"""
bbk_investigation_s20cr.py — Why are BBK triangles ALWAYS zero?
kind-pasteur-2026-03-23-S20cr

DISCOVERY: In G_n/Z_2, no triangle has exactly 2 blue + 1 black edges.
This means: if two classes A,B are blue-connected (same SC-type) and
B,C are blue-connected, then A,C are NEVER black-connected.

Equivalently: the blue neighborhood of any node is closed under the
SC-type (if your blue-neighbor's blue-neighbor shares your SC-type,
you can't be black-connected to them).

This would follow if: blue-connected nodes that are both SC (or both NS)
can NEVER have a black edge between them — but that's the DEFINITION of
blue (same SC type) vs black (different SC type). So:
  - BBK triangle = SC-SC-NS with edges SC~SC (blue), SC~NS (black), SC~NS (???)

Wait. Let me think again:
  BBK = two blue edges + one black edge.
  Blue = same type (SC-SC or NS-NS).
  Black = different type (SC-NS).

  For a triangle (A,B,C) with edges AB=blue, BC=blue, AC=black:
    AB blue => type(A) = type(B)
    BC blue => type(B) = type(C)
    => type(A) = type(C)
    AC black => type(A) != type(C)
    CONTRADICTION!

Therefore BBK triangles are IMPOSSIBLE by the type classification.
This is a clean proof!

Similarly, BBK=0 by the same argument.
BKK can happen: A-B black, B-C black, A-C blue.
  A-B black => type(A)!=type(B)
  B-C black => type(B)!=type(C)
  => type(A) may equal type(C) (if both differ from B's type)
  A-C blue => type(A)=type(C). This is consistent!

So BKK = triangle with SC-NS-SC and the two SC connected by blue.

Let me verify this structural analysis and look at BKK triangles more closely.
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  BBK TRIANGLE IMPOSSIBILITY + BKK STRUCTURE ANALYSIS")
print("  kind-pasteur-2026-03-23-S20cr")
print("=" * 80)

# (Tournament builder from three_view_deep, abbreviated)
def tadj(n, bits):
    a = [[0]*n for _ in range(n)]; idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): a[i][j] = 1
            else: a[j][i] = 1
            idx += 1
    return a

def Hdp(a, n):
    dp = [0]*((1<<n)*n)
    for v in range(n): dp[(1<<v)*n+v] = 1
    for S in range(1, 1<<n):
        if bin(S).count('1') >= n: continue
        for v in range(n):
            if not(S&(1<<v)): continue
            val = dp[S*n+v]
            if val == 0: continue
            for u in range(n):
                if S&(1<<u): continue
                if a[v][u]: dp[(S|(1<<u))*n+u] += val
    return sum(dp[((1<<n)-1)*n+v] for v in range(n))

def canon(a, n):
    sc = [sum(a[i][j] for j in range(n)) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n): sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(set(sc))]
    if all(len(g)==1 for g in gs):
        p = [g[0] for g in gs]
        return tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
    best = None
    def gp(gs):
        if not gs: yield []; return
        for p in permutations(gs[0]):
            for r in gp(gs[1:]): yield list(p)+r
    for p in gp(gs):
        f = tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or f < best: best = f
    return best

def comp(a, n):
    return [[1-a[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

def ssq(a, n):
    return tuple(sorted(sum(a[i][j] for j in range(n)) for i in range(n)))

def c3c(a, n):
    c = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if a[i][j] and a[j][k] and a[k][i]: c += 1
                if a[i][k] and a[k][j] and a[j][i]: c += 1
    return c

def build(n):
    m = comb(n,2); total = 1<<m; t0 = time.time()
    iso = defaultdict(list)
    for b in range(total):
        a = tadj(n, b); iso[canon(a, n)].append(b)
    cl = []
    c2c = {}
    for idx, (cn, mem) in enumerate(sorted(iso.items())):
        a = tadj(n, mem[0])
        cl.append({'cid':idx,'canon':cn,'mem':mem,'adj':a,
                   'H':Hdp(a,n),'sc':canon(comp(a,n),n)==cn,
                   'score':ssq(a,n),'c3':c3c(a,n),
                   'aut':factorial(n)//len(mem),'size':len(mem)})
        c2c[cn] = idx
    for d in cl:
        cc = canon(comp(d['adj'],n),n)
        d['comp_cid'] = c2c.get(cc,-1)
    edges = set(); ecol = {}
    for d in cl:
        ci = d['cid']
        for b in d['mem']:
            for ai in range(m):
                f = b^(1<<ai); fa = tadj(n, f); fc = canon(fa, n)
                nb = c2c.get(fc)
                if nb is not None and nb != ci:
                    e = (min(ci,nb),max(ci,nb)); edges.add(e)
                    if e not in ecol:
                        ecol[e] = 'blue' if d['sc']==cl[nb]['sc'] else 'black'
    mi = {}; mid = 0
    for d in cl:
        ci = d['cid']
        if ci in mi: continue
        cp = d['comp_cid']; mi[ci] = mid
        if cp != ci and cp >= 0: mi[cp] = mid
        mid += 1
    Vm = mid
    me = set(); mec = {}
    for (a,b) in edges:
        ma, mb = mi[a], mi[b]
        if ma != mb:
            e = (min(ma,mb),max(ma,mb)); me.add(e)
            c = ecol.get((min(a,b),max(a,b)),'blue')
            if e not in mec: mec[e] = c
            elif mec[e] != c: mec[e] = 'mixed'
    mc = []; seen = set()
    for d in cl:
        mv = mi[d['cid']]
        if mv not in seen:
            seen.add(mv)
            mc.append({'mid':mv,'H':d['H'],'sc':d['sc'],'score':d['score'],
                       'c3':d['c3'],'aut':d['aut'],'size':d['size']})
    mc.sort(key=lambda x: x['mid'])
    print(f"  Built n={n}: V={Vm}, E={len(me)} [{time.time()-t0:.1f}s]")
    return Vm, mc, me, mec

# ============================================================================
# PROOF: BBK IMPOSSIBILITY
# ============================================================================

print(f"\n{'='*80}")
print("  THEOREM: BBK = 0 (two blue + one black triangle IMPOSSIBLE)")
print(f"{'='*80}")
print("""
  PROOF by type contradiction:

  In G_n/Z_2, every edge is either BLUE (same SC-type) or BLACK (different SC-type).

  Consider a triangle (A, B, C) with edges AB, BC, AC.
  Suppose AB = blue and BC = blue.
  Then: type(A) = type(B) and type(B) = type(C)
  => type(A) = type(C) [transitivity of equality]
  => AC must be BLUE (same type)

  Therefore it is impossible for AB=blue, BC=blue, AC=black.

  Equivalently: BBK = 0 for ALL n.

  This is an EXACT theorem, not a conjecture. It follows purely from
  the definition of blue/black edge coloring by SC-type.

  COROLLARY: Triangle decomposition simplifies to:
    #tri(C) = #tri(BBB) + #tri(BKK)

  since KKK = 0 (black bipartite, THM-A) and BBK = 0 (this theorem).
""")

# ============================================================================
# BKK TRIANGLE STRUCTURE
# ============================================================================

print(f"\n{'='*80}")
print("  BKK TRIANGLE STRUCTURE")
print(f"{'='*80}")
print("""
  BKK triangles have edges: AB=blue, BC=black, AC=black.
  This means: type(A)=type(C) [from AB=blue... wait, no]

  Actually:
    Edges AB=black, BC=black, AC=blue.
    AB black => type(A)!=type(B)
    BC black => type(B)!=type(C)
    AC blue => type(A)=type(C)

  So B has different type from A and C, while A,C have same type.

  If A,C are SC: B is NS. Triangle = SC-NS-SC with blue edge between two SC nodes.
  If A,C are NS: B is SC. Triangle = NS-SC-NS with blue edge between two NS nodes.

  Both are structurally possible. Let's count which types appear.
""")

for n in [3, 4, 5, 6]:
    Vm, mc, me, mec = build(n)

    is_sc = {d['mid']:d['sc'] for d in mc}
    adj = defaultdict(set)
    for e, c in mec.items():
        a, b = e
        adj[a].add(b); adj[b].add(a)

    blue_adj = defaultdict(set)
    black_adj = defaultdict(set)
    for e, c in mec.items():
        a, b = e
        if c == 'blue':
            blue_adj[a].add(b); blue_adj[b].add(a)
        elif c == 'black':
            black_adj[a].add(b); black_adj[b].add(a)

    # Count BKK triangles by type
    bkk_sc_ns_sc = 0  # SC center is NS
    bkk_ns_sc_ns = 0  # NS center is SC
    bkk_details = []

    for i in range(Vm):
        for j in range(i+1, Vm):
            for k in range(j+1, Vm):
                # Check if triangle exists
                if not (adj[i] & {j} and adj[j] & {k} and adj[i] & {k}):
                    continue

                # Count blue and black edges
                edges_ijk = [(i,j), (j,k), (i,k)]
                b_count = sum(1 for a,b in edges_ijk if b in blue_adj[a])
                k_count = sum(1 for a,b in edges_ijk if b in black_adj[a])

                if b_count == 1 and k_count == 2:  # BKK
                    # Find the blue edge
                    for a,b in edges_ijk:
                        if b in blue_adj[a]:
                            blue_pair = (a, b)
                            center = [x for x in [i,j,k] if x != a and x != b][0]
                            break

                    if is_sc[center]:
                        # Center is SC, endpoints are NS
                        bkk_ns_sc_ns += 1
                    else:
                        # Center is NS, endpoints are SC
                        bkk_sc_ns_sc += 1

                    if n <= 5:
                        bp = blue_pair
                        bkk_details.append({
                            'blue_pair': (bp[0], bp[1]),
                            'center': center,
                            'H_values': (mc[bp[0]]['H'], mc[center]['H'], mc[bp[1]]['H']),
                            'types': (is_sc[bp[0]], is_sc[center], is_sc[bp[1]])
                        })

    total_bkk = bkk_sc_ns_sc + bkk_ns_sc_ns
    print(f"\n  n={n}: BKK={total_bkk}")
    print(f"    SC-NS-SC (NS is bridge): {bkk_sc_ns_sc}")
    print(f"    NS-SC-NS (SC is bridge): {bkk_ns_sc_ns}")

    if bkk_details:
        print(f"    Details:")
        for d in bkk_details:
            print(f"      Blue edge: {d['blue_pair']}, Center: {d['center']}, "
                  f"H: {d['H_values']}, Types: {['NS' if not t else 'SC' for t in d['types']]}")

# ============================================================================
# WALK ORTHOGONALITY THEOREM
# ============================================================================

print(f"\n\n{'='*80}")
print("  WALK ORTHOGONALITY ANALYSIS")
print(f"{'='*80}")
print("""
  OBSERVATION: tr(A_B * A_K) = 0 at all tested n.

  This equals sum_{v,w} A_B[v,w] * A_K[v,w] since both matrices are symmetric.
  Since every edge is either blue or black (not both), A_B[v,w]*A_K[v,w]=0
  for all (v,w). Therefore tr(A_B * A_K) = 0 TRIVIALLY.

  More interesting: ALL odd-blue-count walk terms vanish at k=3.
  tr(BBK) = tr(BKB) = tr(KBB) = 0.

  This is equivalent to BBK triangles = 0 (for the diagonal terms)
  PLUS the statement that there are no length-3 closed walks using
  exactly 2 blue + 1 black edge in any order.

  tr(BBK) = sum_v (A_B^2 * A_K)[v,v] = sum_v sum_w A_B^2[v,w]*A_K[w,v]
  = sum_v sum_w (#blue 2-walks v->w) * (black edge w->v)

  This counts: how many ways can you go from v via 2 blue edges to w,
  then return via 1 black edge? If this is 0, it means:
  for every v with a black neighbor w, there is NO blue 2-walk from v to w.

  This follows from type analysis:
    v has a black edge to w => type(v) != type(w)
    A blue 2-walk v->u->w requires type(v)=type(u) and type(u)=type(w)
    => type(v)=type(w), contradiction.

  Therefore tr(BBK) = 0 for ALL n, and similarly tr(BKB) = tr(KBB) = 0.

  GENERAL THEOREM: For any walk of length k, if the walk uses an ODD
  number of black edges, then the diagonal contribution (closed walk
  starting and ending at the same vertex) is... NOT necessarily zero.

  Let's check k=4: BBKK, BKBK, BKKB, KBBK, KBKB, KKBB
  These all use 2 blue + 2 black. Even number of each. Should be nonzero.

  What about BBBK (3 blue, 1 black)?
""")

for n in [4, 5, 6]:
    Vm, mc, me, mec = build(n)

    is_sc = {d['mid']:d['sc'] for d in mc}
    Ab = np.zeros((Vm,Vm)); Ak = np.zeros((Vm,Vm))
    for e, c in mec.items():
        a,b = e
        if c in ('blue','mixed'): Ab[a][b]=1; Ab[b][a]=1
        if c in ('black','mixed'): Ak[a][b]=1; Ak[b][a]=1

    B = Ab; K = Ak
    B2 = B@B; K2 = K@K; BK = B@K; KB = K@B

    # k=4 decomposition: all 16 terms grouped by #blue
    terms_4 = {
        'BBBB': np.trace(B2@B2),
        'BBBK': np.trace(B2@BK),  # 3 blue, 1 black
        'BBKB': np.trace(B@BK@B) if False else np.trace(BK@B2),  # trace is cyclic: tr(ABCD)=tr(BCDA)
        'BKBB': np.trace(KB@B2),
        'KBBB': np.trace(K@B2@B),
        'BBKK': np.trace(B2@K2),
        'BKBK': np.trace(BK@BK),
        'BKKB': np.trace(BK@KB),
        'KBBK': np.trace(KB@BK),
        'KBKB': np.trace(KB@KB),
        'KKBB': np.trace(K2@B2),
        'BKKK': np.trace(BK@K2),
        'KBKK': np.trace(KB@K2),
        'KKBK': np.trace(K2@BK),
        'KKKB': np.trace(K2@KB),
        'KKKK': np.trace(K2@K2),
    }

    # Group by number of blue edges
    by_blue = defaultdict(float)
    for name, val in terms_4.items():
        nb = name.count('B')
        by_blue[nb] += val

    total_k4 = sum(terms_4.values())
    actual_k4 = np.trace((Ab+Ak)@(Ab+Ak)@(Ab+Ak)@(Ab+Ak))

    print(f"\n  n={n}: k=4 walk decomposition by #blue edges:")
    for nb in sorted(by_blue.keys()):
        print(f"    {nb} blue, {4-nb} black: {by_blue[nb]:.0f}")
    print(f"    Total: {total_k4:.0f}, actual: {actual_k4:.0f}, match: {'YES' if abs(total_k4-actual_k4)<0.5 else 'NO'}")

    # Check: are odd-blue-count terms always 0?
    odd_blue = by_blue.get(1,0) + by_blue.get(3,0)
    print(f"    Odd-blue-count terms: {odd_blue:.0f} {'(ZERO!)' if abs(odd_blue)<0.5 else ''}")

print(f"\n\n{'='*80}")
print("  MASTER THEOREM: ODD-BLUE WALK VANISHING")
print(f"{'='*80}")
print("""
  CONJECTURE (verified n=3..6, k=2,3,4):

  For closed walks of ANY length k in G_n/Z_2, the contribution from
  walks with an ODD number of blue edges is ZERO.

  Equivalently: sum over all words w in {B,K}^k with odd #B of
  tr(product of A_B and A_K in order w) = 0.

  PROOF SKETCH:
  Each vertex v has a type t(v) in {SC, NS}. A blue edge preserves type
  (t(u)=t(v)) and a black edge flips type (t(u)!=t(v)).

  A closed walk v -> w1 -> w2 -> ... -> wk -> v must return to the same
  type as the starting vertex. Each blue edge preserves type (contributes
  +0 to the type change count) and each black edge flips type (contributes
  +1). For the walk to be closed, the total number of type flips must be EVEN.

  But the number of type flips = number of black edges in the walk.
  So closed walks must use an EVEN number of black edges.

  Equivalently, closed walks use an EVEN number of blue edges (since k - #black = #blue,
  and k and #black have the same parity).

  Wait — that's not quite right. k - even = even iff k is even.
  For odd k: #black even => #blue = k - #black is odd.
  For even k: #black even => #blue = k - #black is even.

  So for ODD k, closed walks have ODD #blue (not zero contribution).
  For EVEN k, closed walks have EVEN #blue.

  At k=3 (odd): #black must be even (0 or 2), so #blue must be odd (3 or 1).
  BBB is fine (#blue=3, odd). BKK has #blue=1 (odd), also fine.
  BBK has #blue=2 (even) — but #black=1 (odd) => NOT a valid closed walk.

  This means tr(BBK) = 0 because no closed walk of length 3 can have
  exactly 1 black edge (which would be odd, violating the type constraint).

  GENERAL STATEMENT:
  tr(product of A_B^{b_i} and A_K^{k_i}) = 0 whenever the total number
  of K factors is ODD.

  This is because closed walks must use an even number of type-flipping (black) edges.
""")

print(f"\n  DONE.")
print("=" * 80)
