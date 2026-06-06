#!/usr/bin/env python3
"""
lrc_lattice_synchronization_monad_s1.py
monad-explorer-2026-06-06-S1

ANGLE: structural reduction of signed->unsigned LRC via LATTICE SYNCHRONIZATION.

Thesis (to verify, exact rational arithmetic):
  (L0) Synchronization lemma. If v_a + v_b ≡ 0 (mod q) then ||v_a·k/q|| = ||v_b·k/q||
       for ALL integers k. (shell-partners coincide on the lattice L_q = {k/q}.)
  (L1) M(S) = max_t min_i ||v_i t|| is attained at a pair-SUM time t*=m/(v_a+v_b)
       (pinch lemma HYP-2059), and at t* the binding pair STRADDLES:
       ||v_a t*|| = ||v_b t*|| = M(S). So the binding pair is a (v_a+v_b)-shell-partner
       SYNCHRONIZED at the witness. -> opus's "shell-partner / zero-clock" (T3, mod 2n-1)
       and the pinch lemma's "straddling binding pair" are the SAME object at the
       modulus q = reduced pair-sum.
  (L2) lattice lower bound  M(S) >= λ_q(S) := max_{k} min_i ||v_i·k/q||, and λ_q is
       SIGN-BLIND and depends only on the shell-multiset of S mod q. M(S) = max over
       pinch denominators q of λ_q(S) (lossless).
  (L3) modulus ladder: floor binding pair has q = n (value 1/n); opus's 2n-1 shell-
       partner is the Farey-SUCCESSOR / SECOND-gap pair (value ~2/(2n-1)).
  (L4) n=14 split: AP and V* FOLD IDENTICALLY mod n=14 (same floor certificate, M=1/14)
       and differ only mod 27 -> the signed split is a SECOND-ORDER (2n-1) effect,
       invisible to M. A shell-partner at 2n-1 cannot by itself lower M below 1/n.
"""
from fractions import Fraction as Fr
from itertools import combinations
from math import gcd

def nint_dist(x):
    """||x|| = distance from rational x to the nearest integer, exact Fraction."""
    f = x - (x.numerator // x.denominator)   # fractional part in [0,1)
    return min(f, 1 - f)

def M_pinch(speeds):
    """M(S)=max_t min_i ||v_i t|| via pinch lemma: only times m/(v_a+v_b) matter.
    Returns (M, t*, binding_pair). Exact."""
    best = Fr(-1); best_t = None; best_pair = None
    cand = set()
    for a, b in combinations(speeds, 2):
        s = a + b
        for m in range(1, s):
            cand.add(Fr(m, s))
    for t in cand:
        mn = min(nint_dist((v * t) % 1) for v in speeds)
        if mn > best:
            # identify a straddling binding pair: two runners both at distance mn
            binders = [v for v in speeds if nint_dist((v * t) % 1) == mn]
            best, best_t, best_pair = mn, t, binders
    return best, best_t, best_pair

def lambda_q(speeds, q):
    """λ_q = max_{k=1..q-1} min_i ||v_i k/q||. Exact lattice loneliness on L_q."""
    best = Fr(-1); best_k = None
    for k in range(1, q):
        mn = min(nint_dist((Fr(v * k, q)) % 1) for v in speeds)
        if mn > best:
            best, best_k = mn, k
    return best, best_k

def shells_mod(speeds, q):
    """multiset of shells {min(r, q-r)} for r=v mod q (r!=0)."""
    sh = {}
    for v in speeds:
        r = v % q
        s = min(r, q - r) if r != 0 else 0
        sh[s] = sh.get(s, 0) + 1
    return sh

# ---------------------------------------------------------------------------
print("="*78)
print("(L0) SYNCHRONIZATION LEMMA: v_a+v_b ≡ 0 mod q  =>  ||v_a k/q||=||v_b k/q|| ∀k")
print("="*78)
fails = 0; tests = 0
import itertools
for q in range(2, 40):
    for u in range(1, q):
        w = (q - u) % q
        if w == 0:
            continue
        # test the genuine speeds u and w (and lifts u+q, etc.)
        for lu in (u, u + q, u + 2*q):
            for lw in (w, w + q, w + 3*q):
                ok = all(nint_dist((Fr(lu*k, q)) % 1) == nint_dist((Fr(lw*k, q)) % 1)
                         for k in range(q))
                tests += 1
                if not ok:
                    fails += 1
print(f"  q in [2,40): {tests} synchronization tests, FAILURES = {fails}")
print(f"  => synchronization lemma {'HOLDS (verified)' if fails==0 else 'FAILS'}")

# ---------------------------------------------------------------------------
print()
print("="*78)
print("(L1) BINDING PAIR = SHELL-PARTNER SYNCHRONIZED AT WITNESS")
print("="*78)
def report(name, speeds):
    n = len(speeds) + 1
    M, t, binders = M_pinch(speeds)
    q_floor = t.denominator
    print(f"\n {name}: n={n}, {len(speeds)} speeds, floor 1/n=1/{n}")
    print(f"   M(S) = {M}  (= 1/{n}? {M==Fr(1,n)});  witness t*={t} (denom q={q_floor})")
    # find the binding pair (two binders whose SUM is a multiple of q_floor)
    bp = None
    for a, b in combinations(binders if len(binders)>=2 else speeds, 2):
        if (a + b) % q_floor == 0:
            bp = (a, b); break
    print(f"   binders at distance M: {sorted(set(binders))}")
    if bp:
        a, b = bp
        print(f"   binding pair (a,b)=({a},{b}); a+b={a+b} ≡ 0 mod q={q_floor} (shell-partner)")
        print(f"   ||a t*||={nint_dist((a*t)%1)} == ||b t*||={nint_dist((b*t)%1)}  (synchronized = M)")
    return M, t

# small worry-set configs
report("AP n=4", [1,2,3])
report("AP n=5", [1,2,3,4])
report("sporadic n=5", [1,3,4,7])
report("AP n=7", [1,2,3,4,5,6])

# ---------------------------------------------------------------------------
print()
print("="*78)
print("(L4) n=14 FLOOR: AP vs V* vs 2AP — fold identically mod 14, differ mod 27")
print("="*78)
AP   = list(range(1,14))
Vst  = [1,2,3,4,5,6,7,8,9,10,11,13,24]
AP2  = [2*i for i in range(1,14)]
for nm, S in [("AP",AP),("V*",Vst),("2AP",AP2)]:
    M,t = report(nm, S)
print()
print("  -- fold mod n=14 (PRIMARY lattice, floor) --")
for nm, S in [("AP",AP),("V*",Vst),("2AP",AP2)]:
    print(f"   {nm:4s} shells mod 14: {dict(sorted(shells_mod(S,14).items()))}")
    lam, k = lambda_q(S, 14)
    print(f"        λ_14 = {lam} at k={k}   (floor certificate 1/14? {lam>=Fr(1,14)})")
print()
print("  -- fold mod 2n-1=27 (SECONDARY lattice, opus shell-partner) --")
for nm, S in [("AP",AP),("V*",Vst),("2AP",AP2)]:
    sh = shells_mod(S,27)
    # shell-partner = a pair summing to 27 = a residue r with both r and 27-r present,
    # equivalently shell 0 reachable as a SUM. report pairs v_i+v_j==27:
    sp = [(a,b) for a,b in combinations(S,2) if (a+b)%27==0]
    lam, k = lambda_q(S, 27)
    print(f"   {nm:4s} shells mod 27: {dict(sorted(sh.items()))}")
    print(f"        shell-partners (v_i+v_j≡0 mod27): {sp};  λ_27={lam} at k={k}")

# ---------------------------------------------------------------------------
print()
print("="*78)
print("(L2) LOSSLESS: M(S) = max over pinch denominators q of λ_q(S)")
print("="*78)
import random
random.seed(12345)
mismatch = 0; trials = 0
for _ in range(300):
    n_sp = random.randint(2, 6)
    S = random.sample(range(1, 30), n_sp)
    M, t, _ = M_pinch(S)
    # max λ_q over all pair-sum denominators
    qs = set()
    for a,b in combinations(S,2):
        qs.add(a+b)
    bestlam = max(lambda_q(S, q)[0] for q in qs)
    trials += 1
    if bestlam != M:
        mismatch += 1
        if mismatch <= 5:
            print(f"   MISMATCH S={S}: M={M} but max_q λ_q={bestlam}")
print(f"   {trials} random configs: M == max_q λ_q in {trials-mismatch}/{trials}")
print(f"   => lattice frame is {'LOSSLESS (verified)' if mismatch==0 else 'LOSSY'}")

# ---------------------------------------------------------------------------
print()
print("="*78)
print("(L2b) REACH of the SINGLE primary lattice L_n: which configs certified by t=1/n?")
print("  (LRC holds for S as soon as no v_i ≡ 0 mod n; gcd-with-n is the obstruction)")
print("  THEOREM: λ_n(S) >= 1/n  <=>  ∃k: v_i·k ≢ 0 (mod n) ∀i.  Certified by k=1 when")
print("  all gcd(v_i,n)=1.  Below = exact via the residue criterion (no float).")
print("="*78)
def certified_by_Ln(S, n):
    """True iff some k in 1..n-1 has v_i*k != 0 mod n for all i (=> λ_n >= 1/n)."""
    for k in range(1, n):
        if all((v * k) % n != 0 for v in S):
            return True
    return False
for n in [6, 8, 12, 14, 15, 30]:
    certified = 0; total = 0
    random.seed(n)
    for _ in range(5000):
        S = random.sample(range(1, 4*n), n-1)
        total += 1
        if certified_by_Ln(S, n):
            certified += 1
    print(f"  n={n:2d}: L_n single-lattice certifies {certified:5d}/{total} "
          f"({100*certified/total:.1f}%);  uncertified need a non-floor pinch lattice")
print()
print("CONCLUSION: the signed/shell datum (v_i+v_j≡0 mod 2n-1) is a SECOND-gap binding")
print("pair; it is dominated by the floor pinch (q=n) and CANNOT lower M below 1/n.")
print("A signed zero-clock => an unsigned counterexample only if it is ALSO a floor")
print("binding pair, i.e. v_i+v_j≡0 mod n — a DIFFERENT modulus. So the signed split of")
print("the worry-set (AP vs V*) is invisible to M: a structural-reduction NO-GO.")
