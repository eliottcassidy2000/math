"""kps-S34: THE MAGNITUDE-INDEPENDENT CRUX. For which q can 13 danger-sets cover Z/q?
danger-set of residue r at modulus q = {a : ||r a / q|| < 1/14} = {a : r a mod q in [-k,k]},
k = largest m with m/q < 1/14. Covering Z/q by 13 such sets = family blocked at q (no a/q
witness). If the max coverable q is BOUNDED by a constant B, then EVERY free modulus q>B has
a lonely witness => census q_min <= max(first_free, B) = O(log M). Compute B."""
import random
random.seed(0)

def danger_set(r, q):
    k = (q-1)//14                        # largest m with 14m < q  (m/q < 1/14)
    S = set()
    for a in range(q):
        ra = (r*a) % q
        if min(ra, q-ra) <= k and (14*min(ra,q-ra) < q):
            S.add(a)
    return S

def band(q):
    return (q-1)//14

def can_cover(q, restarts=200):
    """try to cover Z/q with 13 danger-sets (greedy + random restarts). returns (covered?, residues)."""
    universe = set(range(q))
    # precompute danger sets for all residues 1..q-1
    dsets = {r: danger_set(r, q) for r in range(1, q)}
    best_unc = q
    for _ in range(restarts):
        uncovered = set(range(1, q))     # 0 always covered by any set
        chosen = []
        # random-greedy: occasionally random first pick
        rlist = list(range(1, q))
        for step in range(13):
            if step == 0 and random.random() < 0.5:
                r = random.choice(rlist)
            else:
                # pick residue covering most uncovered
                r = max(rlist, key=lambda r: len(dsets[r] & uncovered))
            chosen.append(r)
            uncovered -= dsets[r]
            if not uncovered:
                return True, chosen
        best_unc = min(best_unc, len(uncovered))
    return False, best_unc

print("q : coverable by 13 danger-sets? (band k = floor((q-1)/14), |danger|=2k+1)")
maxcov = 0
for q in range(15, 130):
    ok, info = can_cover(q, restarts=120)
    k = band(q)
    tag = "COVER" if ok else f"gap(min_unc={info})"
    if ok: maxcov = q
    mark = "  <== coverable" if ok else ""
    # only print transitions / coverable / small gaps to keep output short
    if ok or q < 30 or (isinstance(info,int) and info <= 3):
        print(f"q={q:>3} k={k:>2} |danger|={2*k+1:>2} 13*|d|/q={13*(2*k+1)/q:>4.2f} : {tag}{mark}")
print()
print(f"MAX q coverable by 13 danger-sets (greedy search): B = {maxcov}")
print()
print("If B is a small constant => every free modulus q>B has a witness => the census")
print("q_min <= max(first_free_modulus, B) = O(log M). The compressed crux is CLOSED to")
print("`no covering above B` = a FINITE, magnitude-independent covering fact.")
