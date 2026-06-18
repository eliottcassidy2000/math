"""
Fast unbuffered driver for the M2 Weyl-phase forbidden-class refutation.
Imports the map/invariants from the main refute module, memoizes fingerprints,
and runs a BROAD adversarial search:
  (A) all primitive 5-subsets of {1..K} at optimum tau* and at all candidate times
  (B) induced 5-vertex sub-tournaments of LARGER LRC configs (n=6..13) at their
      lonely optima -- covering sets, tight AP sets, sporadic/random large sets,
      off-grid optima.
Target: realize either of the 2 n=5 forbidden classes (H=13 or H=15, 3cyc=4,
score (1,2,2,2,3)). Any single realization => REFUTED.
All decisions exact (Fraction).
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations
import importlib.util

spec = importlib.util.spec_from_file_location(
    "m2mod",
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "lrc14_refute_character-spectral_M2WeylphaseImofe_kps-S2-wf.py"))
m2 = importlib.util.module_from_spec(spec)
# prevent its __main__ from running: it's guarded, and we import as module name 'm2mod'
spec.loader.exec_module(m2)

Mgap = m2.Mgap
cand = m2.cand
m2_tournament = m2.m2_tournament
is_primitive = m2.is_primitive
FORBIDDEN_5 = m2.FORBIDDEN_5
adj_matrix = m2.adj_matrix
score_seq = m2.score_seq

# ---- memoized fingerprint keyed on frozenset of arcs (n fixed per call) ----
_fp_cache = {}
def fp_fast(n, arcs):
    key = (n, frozenset(arcs))
    v = _fp_cache.get(key)
    if v is not None: return v
    A = adj_matrix(n, arcs)
    # ham paths
    hp = 0
    rng = range(n)
    for perm in permutations(rng):
        ok = True
        for a in range(n-1):
            if not A[perm[a]][perm[a+1]]:
                ok = False; break
        if ok: hp += 1
    # 3-cycles
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # count cyclic triangles among i<j<k (either orientation)
                if A[i][j] and A[j][k] and A[k][i]: c3 += 1
                elif A[j][i] and A[k][j] and A[i][k]: c3 += 1
    sc = score_seq(A)
    v = (hp, c3, sc)
    _fp_cache[key] = v
    return v

def flush(*a):
    print(*a); sys.stdout.flush()

def scan_5subsets_at_times(speeds_sorted, times):
    """Yield fingerprints of M2 tournament on all 5-subsets of speeds at each tau in times."""
    sp = speeds_sorted
    out = []
    for sub in combinations(range(len(sp)), 5):
        subsp = [sp[i] for i in sub]
        for tau in times:
            n, arcs = m2_tournament(subsp, tau)
            fp = fp_fast(5, arcs)
            out.append((fp, subsp, tau))
    return out

def family_direct(pool, label, all_times=False):
    realized = set(); forb = []
    nconf = 0
    pool = list(pool)
    for S in combinations(pool, 5):
        if not is_primitive(S): continue
        nconf += 1
        gap, ats = Mgap(list(S))
        times = list(cand(list(S))) if all_times else ats
        sp = sorted(S)
        for tau in times:
            n, arcs = m2_tournament(sp, tau)
            fp = fp_fast(5, arcs)
            realized.add(fp)
            if fp in FORBIDDEN_5:
                forb.append((S, tau, gap, fp))
    flush(f"[{label}] primitive-configs={nconf} distinct-fps={len(realized)} forbidden={len(forb)}")
    if forb:
        for w in forb[:6]:
            flush(f"    *** REFUTED: S={w[0]} tau={w[1]} gap={w[2]} fp={w[3]}")
    return realized, forb

def family_induced(big_sets, label):
    realized = set(); forb = []
    nbig = 0
    for S in big_sets:
        S = tuple(sorted(set(S)))
        if len(S) < 5 or not is_primitive(S): continue
        nbig += 1
        gap, ats = Mgap(list(S))
        sp = sorted(S)
        for tau in ats:
            for sub in combinations(range(len(sp)), 5):
                subsp = [sp[i] for i in sub]
                n, arcs = m2_tournament(subsp, tau)
                fp = fp_fast(5, arcs)
                realized.add(fp)
                if fp in FORBIDDEN_5:
                    forb.append((S, tau, gap, subsp, fp))
    flush(f"[{label}] big-configs={nbig} distinct-n5-fps={len(realized)} forbidden={len(forb)}")
    if forb:
        for w in forb[:6]:
            flush(f"    *** REFUTED induced: bigS={w[0]} tau={w[1]} sub={w[3]} gap={w[2]} fp={w[4]}")
    return realized, forb

def main():
    allr = set(); allf = []
    flush("="*70); flush("M2 WEYL-PHASE BROAD REFUTATION DRIVER"); flush("="*70)
    flush(f"Forbidden targets: {sorted(FORBIDDEN_5)}")

    # ---- (A) direct n=5 families ----
    for K in (16, 18, 20):
        r,f = family_direct(range(1,K+1), f"direct {{1..{K}}} optimum", all_times=False)
        allr|=r; allf+=f
    for K in (16, 20):
        r,f = family_direct(range(1,K+1), f"direct {{1..{K}}} all-times", all_times=True)
        allr|=r; allf+=f

    # ---- (B) induced from larger LRC configs ----
    # B1: covering / hard families related to LRC(14)
    big = []
    big.append(tuple(range(1,12))+(13,))                  # {1..11,13} the AP-drop
    big.append(tuple(range(1,12))+(13,84))                # covering family m=1 (n=13)
    big.append(tuple(range(1,14)))                        # AP {1..13} (tight, n=13)
    big.append(tuple(range(1,12))+(13,24))                # Goddyn-Wong {1..11,13,24}
    big.append(tuple(list(range(1,12))+[13]))             # 12-set core
    for w in (14,28,42,56,84):                            # parked-runner covering cores
        big.append(tuple(list(range(1,12))+[13,w]))
    r,f = family_induced(big, "induced covering/tight n=12/13")
    allr|=r; allf+=f

    # B2: all primitive 6-subsets of {1..12} (exhaustive medium) induced->5
    six = [S for S in combinations(range(1,13),6) if is_primitive(S)]
    r,f = family_induced(six, "induced 6-subsets {1..12}")
    allr|=r; allf+=f

    # B3: all primitive 7-subsets of {1..11} induced->5
    sev = [S for S in combinations(range(1,12),7) if is_primitive(S)]
    r,f = family_induced(sev, "induced 7-subsets {1..11}")
    allr|=r; allf+=f

    # B4: random large sporadic sets (n=8..13), bigger speed range, off-grid optima
    random.seed(12345)
    rand_sets = []
    for _ in range(4000):
        n = random.randint(8,13)
        hi = random.choice([20,30,50,80,120])
        S = tuple(sorted(random.sample(range(1,hi+1), n)))
        if is_primitive(S):
            rand_sets.append(S)
    r,f = family_induced(rand_sets, "induced random sporadic n=8..13")
    allr|=r; allf+=f

    flush("")
    flush(f"GRAND TOTAL forbidden witnesses: {len(allf)}")
    flush(f"Distinct n=5 fingerprints realized across ALL families: {len(allr)}")
    flush("Realized fingerprints:")
    for fp in sorted(allr):
        flag = "  <-- FORBIDDEN REALIZED!" if fp in FORBIDDEN_5 else ""
        flush(f"   H={fp[0]:3d} c3={fp[1]} score={fp[2]}{flag}")
    missing = FORBIDDEN_5 - allr
    flush("")
    flush(f"Forbidden classes NEVER realized: {len(missing)} -> {sorted(missing)}")
    if not allf:
        flush("RESULT: claim SURVIVED this broad search (no forbidden class realized).")
    else:
        flush("RESULT: claim REFUTED (forbidden class realized; see witnesses above).")

if __name__ == "__main__":
    main()
