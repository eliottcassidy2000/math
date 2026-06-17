#!/usr/bin/env python3
"""
ANGLE B: the BINDING-PAIR switch reduction for LRC(14).

CLAIM under test: M(S)=max_tau min_v ||v tau|| is ALWAYS attained at a
"crossing" tau where TWO runners are equidistant from the observer (0),
i.e. ||v_a tau|| = ||v_b tau|| = M(S) for some pair (a,b). That pair is the
BINDING PAIR. If true, LRC(14) reduces to: over O(n^2) candidate pairs,
find a crossing tau where the pair's half-gap >= 1/14 AND all other runners
are at distance >= that gap ("others-clear").

We use the EXACT M tool from the prompt (stdlib only, Fraction).
"""
from fractions import Fraction as F
from itertools import combinations

def nrm(x):
    r = x - int(x)
    r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def g(S, t):
    return min(nrm(v * t) for v in S)

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2):
            C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2):
                        C.add(F(k, d)); k += 1
    C.add(F(1, 2))
    return C

def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b:
            b = v; at = t
    return b, at

def all_optima(S):
    """All tau in candidate set achieving M(S)."""
    Mv, _ = M(S)
    return Mv, sorted(t for t in cand(S) if g(S, t) == Mv)

def binding_runners(S, t, gap):
    """Runners exactly at distance == gap from 0 at time t (the binding set)."""
    return [v for v in sorted(set(S)) if nrm(v*t) == gap]

def classify_crossing(S, t, gap):
    """
    At optimum t with value gap, identify which mechanism pins it.
    Returns dict with binding runners and the (a,b) pair structure.
    A 'pair crossing' means two runners v_a,v_b have ||v_a t||=||v_b t||=gap,
    and t = k/(v_a +- v_b). A 'half-integer' means one runner v with
    t=(2k+1)/(2v) so ||v t||=1/2 ... but that only pins if gap involves it.
    """
    bind = binding_runners(S, t, gap)
    info = {"t": t, "gap": gap, "binding_runners": bind, "type": None,
            "pairs": []}
    # check pair crossings among binding runners
    for a, b in combinations(bind, 2):
        # is t = k/(a+b) or k/(b-a)?
        for d in (a+b, abs(b-a)):
            if d > 0 and (t*d).denominator == 1:
                info["pairs"].append((a, b, d, int(t*d)))
    if len(bind) >= 2 and info["pairs"]:
        info["type"] = "binding-pair-crossing"
    elif len(bind) == 1:
        info["type"] = "single-runner-halfinteger"
    elif len(bind) >= 2:
        info["type"] = "multi-no-sum-diff"  # binding but not via sum/diff
    else:
        info["type"] = "none-binding(!)"
    return info

def others_clear(S, t, gap, pair):
    """Are all runners NOT in 'pair' at distance >= gap at time t?"""
    for v in sorted(set(S)):
        if v in pair:
            continue
        if nrm(v*t) < gap:
            return False
    return True

# ----------------------------------------------------------------------
print("="*72)
print("PART 1: For many configs, find the binding pair(s) at the optimum.")
print("="*72)

configs = {
    "{1..13} tight AP (perfect SDR)": list(range(1, 14)),
    "{1..6} small": list(range(1, 7)),
    "{1..7}": list(range(1, 8)),
    "{2,3,...,14}": list(range(2, 15)),
    "covering core {1..11,13,84}": list(range(1, 12)) + [13, 84],
    "covering {1..11,13} (no 84)": list(range(1, 12)) + [13],
    "{1,2,3,4,5}": [1,2,3,4,5],
    "{1,2,4,8,16}": [1,2,4,8,16],
    "geometric-ish {1,3,7,12}": [1,3,7,12],
    "Lonely Runner classic {1,2,3,4,5,6}": list(range(1,7)),
    "{3,5,9,11} order-3 part": [3,5,9,11],
    "single binding {5,7}": [5,7],
    "{1..13} minus 7": [v for v in range(1,14) if v != 7],
}

def summarize(name, S):
    Mv, opts = all_optima(S)
    print(f"\n### {name}")
    print(f"    S = {sorted(set(S))}")
    print(f"    M(S) = {Mv}  (= {float(Mv):.6f}),  vs 1/14 = {float(F(1,14)):.6f}")
    print(f"    LRC threshold 1/14: {'PASS (M>=1/14)' if Mv>=F(1,14) else 'FAIL (M<1/14)'}")
    print(f"    optima tau (achieving M): {opts}")
    for t in opts:
        ci = classify_crossing(S, t, Mv)
        print(f"      tau={t}: type={ci['type']}, binding runners={ci['binding_runners']}")
        for (a,b,d,k) in ci["pairs"]:
            sign = "+" if d == a+b else "-"
            oc = others_clear(S, t, Mv, {a,b})
            print(f"         pair ({a},{b}) via v_a{sign}v_b={d}, tau={k}/{d}, others-clear={oc}")
    return Mv, opts

results = {}
for name, S in configs.items():
    results[name] = summarize(name, S)

print()
print("="*72)
print("PART 2: Reformulate LRC(14) as binding-pair + others-clear.")
print("   How many candidate pair-switches must be checked? O(n^2).")
print("="*72)

def lrc_via_pairs(S, thresh=F(1,14)):
    """
    Direct reformulation: enumerate candidate pair crossings tau=k/(a+-b),
    (and half-integer single-runner candidates), and for each ask:
      gap at tau >= thresh  AND  it's a true min (others-clear automatically
      since g already = min over all). M >= thresh iff SOME candidate tau
      has g(S,tau) >= thresh.
    We separately count how many candidates come from PAIRS vs half-integers,
    to show the O(n^2) structure.
    """
    Sset = sorted(set(S)); n = len(Sset)
    pair_taus = set(); half_taus = set()
    for v in Sset:
        k = 0
        while F(2*k+1, 2*v) <= F(1,2):
            half_taus.add(F(2*k+1,2*v)); k += 1
    for i in range(n):
        for j in range(i+1, n):
            for d in (Sset[i]+Sset[j], Sset[j]-Sset[i]):
                if d > 0:
                    k = 1
                    while F(k,d) <= F(1,2):
                        pair_taus.add(F(k,d)); k += 1
    half_taus.add(F(1,2))
    # which witnesses LRC?
    best_pair = max((g(S,t) for t in pair_taus), default=F(0))
    best_half = max((g(S,t) for t in half_taus), default=F(0))
    Mv = max(best_pair, best_half)
    # is the global optimum reachable from a PAIR crossing alone?
    pair_only_M = best_pair
    return {
        "M": Mv, "M_from_pairs_only": pair_only_M,
        "M_from_halfint_only": best_half,
        "n_pair_taus": len(pair_taus), "n_half_taus": len(half_taus),
        "n_pairs_Cn2": n*(n-1)//2,
        "pair_suffices": pair_only_M == Mv,
    }

print("\nConfig | M | M(pairs only) | M(half-int only) | pairs suffice? | #pairs C(n,2) | #pair-taus | #half-taus")
for name, S in configs.items():
    r = lrc_via_pairs(S)
    print(f"  {name[:34]:34s} | M={r['M']!s:>8} | pairs={r['M_from_pairs_only']!s:>8} | "
          f"half={r['M_from_halfint_only']!s:>8} | suffice={r['pair_suffices']} | "
          f"C(n,2)={r['n_pairs_Cn2']} | #pairtaus={r['n_pair_taus']} | #halftaus={r['n_half_taus']}")

print()
print("="*72)
print("PART 3: Hard core (covering sets). Which pair binds?")
print("   Is the binding pair always {small speed, the multiple of 14}?")
print("="*72)

# Build several covering-set hard cores: include a multiple of 14.
hard_cores = {
    "{1..11,13,84}  (84=6*14)": list(range(1,12)) + [13, 84],
    "{1..13,14}     (14=1*14)": list(range(1,15)),
    "{1..13,28}     (28=2*14)": list(range(1,14)) + [28],
    "{1..11,13,14}": list(range(1,12)) + [13,14],
    "{1..13,42}     (42=3*14)": list(range(1,14)) + [42],
    "{1..9,11,13,84}": list(range(1,10)) + [11,13,84],
}

print("\nFor each: M, optima, binding runners, and whether the multiple-of-14")
print("speed participates in the binding pair.")
for name, S in hard_cores.items():
    Mv, opts = all_optima(S)
    mult14 = [v for v in sorted(set(S)) if v % 14 == 0]
    print(f"\n### {name}")
    print(f"    M(S)={Mv} ({float(Mv):.6f}), multiples of 14 in S: {mult14}")
    for t in opts:
        ci = classify_crossing(S, t, Mv)
        br = ci["binding_runners"]
        m14_in_bind = [v for v in br if v % 14 == 0]
        print(f"    tau={t}: binding={br}, type={ci['type']}, mult14 binding? {m14_in_bind or 'NO'}")
        for (a,b,d,k) in ci["pairs"]:
            sign = "+" if d==a+b else "-"
            print(f"        pair ({a},{b}) via {sign}, tau={k}/{d}")

print()
print("="*72)
print("PART 3b: Off-grid vs on-grid. The multiple-of-14 runner sits in")
print("   section 0 at EVERY a/14 (||14m * a/14|| = ||m*a|| = 0).")
print("   So on the GRID a/14 the covering core is never lonely; M is")
print("   attained OFF-grid. Verify the optimum tau is NOT of form a/14.")
print("="*72)
for name, S in hard_cores.items():
    Mv, opts = all_optima(S)
    for t in opts:
        on_grid = (t * 14).denominator == 1
        print(f"  {name[:30]:30s} tau={t!s:>8}  denom={t.denominator:<4} on-grid(=a/14)? {on_grid}")
        break  # one optimum is enough to illustrate

print()
print("="*72)
print("PART 4: Stress test the CORE CLAIM across MANY random configs:")
print("   is M ALWAYS attained at a tau with >=2 binding runners that form")
print("   a sum/diff pair? Report any counterexample.")
print("="*72)
import random
random.seed(12345)
counter_pure_single = 0
counter_no_pair = 0
tested = 0
examples_single = []
for trial in range(4000):
    n = random.randint(2, 8)
    S = random.sample(range(1, 40), n)
    Mv, opts = all_optima(S)
    if Mv == 0:
        continue
    tested += 1
    # does at least one optimum have a binding PAIR crossing?
    has_pair_opt = False
    only_single = True
    for t in opts:
        ci = classify_crossing(S, t, Mv)
        if ci["type"] == "binding-pair-crossing":
            has_pair_opt = True
            only_single = False
        elif len(ci["binding_runners"]) >= 2:
            only_single = False
    if not has_pair_opt:
        counter_no_pair += 1
        if only_single and len(examples_single) < 8:
            examples_single.append((sorted(set(S)), Mv, opts,
                [(t, binding_runners(S,t,Mv)) for t in opts]))
        counter_pure_single += 1 if only_single else 0

print(f"\nTested {tested} random configs (n=2..8).")
print(f"Configs where NO optimum is a binding-pair-crossing: {counter_no_pair}")
print(f"  ...of those, where every optimum has a SINGLE binding runner: {counter_pure_single}")
print("\nExamples where the optimum is pinned by a SINGLE runner (half-integer),")
print("NOT a pair (these are honest exceptions to a naive 'always a pair' claim):")
for S, Mv, opts, det in examples_single:
    print(f"  S={S}: M={Mv}={float(Mv):.4f}")
    for t, br in det:
        kind = "half-int (2k+1)/(2v)" if any((t*2*v).numerator%2==1 and (t*2*v).denominator==1 for v in br) else "?"
        print(f"     tau={t}, binding runner(s)={br}  [{kind}]")

print()
print("="*72)
print("PART 4b: Refined claim test. A SINGLE binding runner v at tau=(2k+1)/(2v)")
print("  means ||v tau||=1/2, the MAX possible. That only pins M when 1/2 is")
print("  the min, i.e. ALL runners are at 1/2 => M=1/2 (only S={powers giving")
print("  half}). Otherwise the min is set by a runner NOT at 1/2. Test: when")
print("  M<1/2, is the binding set ALWAYS >=2 and contain a sum/diff pair?")
print("="*72)
random.seed(999)
viol = 0; tested2 = 0; viol_ex = []
for trial in range(6000):
    n = random.randint(2, 9)
    S = random.sample(range(1, 50), n)
    Mv, opts = all_optima(S)
    if Mv == 0 or Mv == F(1,2):
        continue
    tested2 += 1
    # at least one optimum must have a sum/diff binding pair
    ok = False
    for t in opts:
        ci = classify_crossing(S, t, Mv)
        if ci["pairs"]:
            ok = True; break
    if not ok:
        viol += 1
        if len(viol_ex) < 6:
            viol_ex.append((sorted(set(S)), Mv,
                [(t, binding_runners(S,t,Mv)) for t in opts]))
print(f"\nTested {tested2} configs with 0<M<1/2.")
print(f"Configs where NO optimum has a sum/diff binding pair: {viol}")
if viol_ex:
    print("Counterexamples to 'M<1/2 => sum/diff binding pair':")
    for S, Mv, det in viol_ex:
        print(f"  S={S}: M={Mv}")
        for t, br in det:
            print(f"     tau={t}, binding={br}")
else:
    print("==> NO counterexamples: when 0<M<1/2, the optimum ALWAYS has a")
    print("    sum/diff binding pair (REFINED CLAIM holds on this sample).")
