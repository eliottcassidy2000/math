"""
lrc14_refute_residue-orbit_M4_lean_kps-S2-wf.py

LEAN exhaustive ceiling for the M4 band-depth orbit-majority map (mod 14), n=5.
Only computes permutation-INVARIANT signature (H, #3cyc, score) -- no canon needed,
since the signature alone separates the 4 forbidden classes from each other and from
the reachable ones (each forbidden (H,#3cyc,score) is a distinct iso class on 5 verts).

Ceiling = ALL residue multisets of size 5 from {0..13} (repetition allowed: loneliness
permits residue collisions; SDR strictly stronger), sweeping ALL tie-break orders for
multisets that contain a tied pair. M4 depends only on residues mod 14, so this is the
true maximum M4 can output. Exact integer arithmetic.

Optimization: precompute, per multiset, the FORCED arcs from WIN; only branch on tied
pairs (typically 0-3 of them), enumerating 2^(#tied) orientations directly instead of
5!=120 full orderings. NOTE: not every independent 2^t orientation is realizable by a
single total tie-break order, so this is a SUPERSET of M4 (an over-approximation of the
ceiling). If even this superset misses the forbidden classes, M4 certainly does too.
We ALSO run the exact-order version (5! sweep) as a cross-check on the realized set.
"""
from math import gcd
from itertools import combinations, permutations, combinations_with_replacement, product

MOD = 14
UNITS = [a for a in range(1, MOD) if gcd(a, MOD) == 1]
def depth(r):
    r %= MOD
    return min(r, MOD - r)

WIN = [[0]*MOD for _ in range(MOD)]
for x in range(MOD):
    for y in range(MOD):
        wi = wj = 0
        for a in UNITS:
            di = depth(x*a); dj = depth(y*a)
            if di > dj: wi += 1
            elif dj > di: wj += 1
        WIN[x][y] = 1 if wi > wj else (-1 if wj > wi else 0)

m = 5
PAIRS = [(i, j) for i in range(m) for j in range(i+1, m)]
TRIPLES = list(combinations(range(m), 3))
PERMS = list(permutations(range(m)))

def H_from_arc(arc):
    # arc[(i,j)] = 1 if i->j else 0 (for i<j); reconstruct full and count Ham paths
    c = 0
    for p in PERMS:
        ok = True
        for k in range(m-1):
            a, b = p[k], p[k+1]
            if a < b:
                if arc[(a, b)] != 1: ok = False; break
            else:
                if arc[(b, a)] != 0: ok = False; break
        if ok: c += 1
    return c
def c3_from_arc(arc):
    def dir(a, b):  # 1 if a->b
        return arc[(a, b)] if a < b else (1 - arc[(b, a)])
    c = 0
    for a, b, cc in TRIPLES:
        if dir(a,b) and dir(b,cc) and dir(cc,a): c += 1
        if dir(a,cc) and dir(cc,b) and dir(b,a): c += 1
    return c
def score_from_arc(arc):
    out = [0]*m
    for (i, j) in PAIRS:
        if arc[(i, j)] == 1: out[i] += 1
        else: out[j] += 1
    return tuple(sorted(out))
def sig_from_arc(arc):
    return (H_from_arc(arc), c3_from_arc(arc), score_from_arc(arc))

FORB = {
    (9, 3, (1,1,2,3,3)),
    (13,4, (1,2,2,2,3)),
    (15,4, (1,2,2,2,3)),
    (15,5, (2,2,2,2,2)),
}

def run(exact_order):
    realized = {}      # sig -> example
    forb_hits = []
    n = 0
    for res in combinations_with_replacement(range(MOD), m):
        n += 1
        forced = {}
        tied = []
        for (i, j) in PAIRS:
            w = WIN[res[i]][res[j]]
            if w == 1: forced[(i, j)] = 1
            elif w == -1: forced[(i, j)] = 0
            else: tied.append((i, j))
        if not tied:
            sg = sig_from_arc(forced)
            if sg not in realized: realized[sg] = (res, "no-tie")
            if sg in FORB: forb_hits.append((res, "no-tie", sg))
            continue
        if exact_order:
            # exact: sweep all total orders; tie (i<j) -> 1 if order[i]<order[j]
            for order in PERMS:
                arc = dict(forced)
                for (i, j) in tied:
                    arc[(i, j)] = 1 if order[i] < order[j] else 0
                sg = sig_from_arc(arc)
                if sg not in realized: realized[sg] = (res, order)
                if sg in FORB: forb_hits.append((res, order, sg))
        else:
            # superset: independent orientations of tied pairs
            for bits in product([0, 1], repeat=len(tied)):
                arc = dict(forced)
                for (ij, b) in zip(tied, bits):
                    arc[ij] = b
                sg = sig_from_arc(arc)
                if sg not in realized: realized[sg] = (res, bits)
                if sg in FORB: forb_hits.append((res, bits, sg))
    return realized, forb_hits, n

def report(tag, realized, forb_hits, n):
    print(f"--- {tag}: enumerated {n} multisets ---")
    print(f"    signatures realized: {len(realized)}")
    for sg in sorted(realized):
        flag = "  <-- FORBIDDEN!" if sg in FORB else ""
        print(f"      {sg}{flag}")
    hit = sorted(set(s for *_, s in forb_hits))
    if forb_hits:
        print(f"    *** FORBIDDEN REALIZED: {hit}  -> REFUTED")
        for w in forb_hits[:6]:
            print("        witness:", w)
    else:
        print("    no forbidden class realized.")
    print()

if __name__ == "__main__":
    print("M4 band-depth orbit-majority, mod 14, n=5. Absolute residue ceiling.")
    print(f"UNITS={UNITS}")
    print()
    # superset first (fast, strictly larger)
    r1, f1, n1 = run(exact_order=False)
    report("SUPERSET (independent tie orientations, OVER-approx of M4)", r1, f1, n1)
    # exact tie-break orders
    r2, f2, n2 = run(exact_order=True)
    report("EXACT (realizable total tie-break orders, = true M4 ceiling)", r2, f2, n2)
    print("VERDICT:")
    if f2:
        print("  REFUTED -- forbidden class realized by exact M4 ceiling.")
    elif f1:
        print("  Forbidden appears only in the OVER-approx superset, not exact M4 ->")
        print("  M4 ceiling itself does NOT realize it (superset artifact). CONFIRMED.")
    else:
        print("  CONFIRMED -- forbidden classes unreachable even by the over-approx superset.")
