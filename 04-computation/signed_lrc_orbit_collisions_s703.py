"""SIGNED LRC -- exhaustive small-n sign-orbit / cut-spectrum collisions.  monad-explorer-S703.

Builds on S699 (signed-LRC theory: T1 gauge-invariance, T2 sign=cut, T3 zero-clock<=>shell-partner,
T4 "sign-orbit of AP = 2,4,7,16,32,60 for n=3..8; collisions are config automorphisms").

GOAL OF THIS SCRIPT.
  (a) Re-derive T4 EXACTLY (matching S699's mod-C folded clock definition) and EXTEND the sequence.
  (b) CHARACTERIZE the collisions rigorously.  S699 called them "config automorphisms"; this script
      tests that claim against an exact automorphism check and against the UNFOLDED clock multiset.
  (c) Separate the TWO sources of an invariant being coarser than 2^{n-2}:
        (i)  UNFOLDED collisions  -> genuine additive coincidences of S (a real automorphism /
             additive-energy structure), present with or without mod C.
        (ii) FOLD-ONLY collisions -> created purely by reducing clocks mod C=2n-1 (a number-theoretic
             folding artifact, NOT an automorphism).
  (d) Exhaustive worry-set (tight, M=1/n) census at small n with signed signatures.

Definitions (match signed_lrc_cut_structure_s699c.py):
  config V = sorted distinct positive speeds (the n-1 runners; observer = speed 0).
  C = 2n-1.  A sign vector eps in {+-1}^{n-1} 2-colors the runners.
  pair (i,j) clock = eps_i V_i - eps_j V_j;  |clock| in {|Vi-Vj| (mono), Vi+Vj (bichromatic=cut)}.
  FOLDED clock (S699):  f = clock mod C; value = min(f, C-f)  in [0, C/2].
  sign-orbit = # distinct clock-multisets as eps ranges over all 2^{n-1} signs.
"""
from itertools import product, combinations
from fractions import Fraction


def gap_exact(V):
    """Exact LRC gap M(V) = sup_t min_i ||V_i t||, ||.|| = dist to nearest int.
    Candidate optimal t in (0,1): t=k/(Vi+Vj) (two runners cross) and t=(2k+1)/(2Vi) (one at 1/2)."""
    cand = set()
    for vi in V:
        for k in range(1, 2 * vi):
            cand.add(Fraction(k, 2 * vi))            # includes Vi*t = half-integer
    for i in range(len(V)):
        for j in range(len(V)):
            s = V[i] + V[j]
            for k in range(1, s):
                cand.add(Fraction(k, s))
            d = abs(V[i] - V[j])
            for k in range(1, max(d, 1)):
                cand.add(Fraction(k, d)) if d else None
    best = Fraction(0)
    for t in cand:
        if not (0 < t < 1):
            continue
        m = min(min((V[i] * t) % 1, 1 - (V[i] * t) % 1) for i in range(len(V)))
        if m > best:
            best = m
    return best


def folded_clocks(V, eps, C):
    out = []
    n1 = len(V)
    for i in range(n1):
        for j in range(i + 1, n1):
            clk = eps[i] * V[i] - eps[j] * V[j]
            f = clk % C
            out.append(min(f, C - f))
    return tuple(sorted(out))


def unfolded_clocks(V, eps):
    """|clock| with NO mod -- the raw {|Vi-Vj|, Vi+Vj} multiset."""
    out = []
    n1 = len(V)
    for i in range(n1):
        for j in range(i + 1, n1):
            out.append(abs(eps[i] * V[i] - eps[j] * V[j]))
    return tuple(sorted(out))


def cut_signature(V, A, C, fold=True):
    """A = frozenset of indices colored +; clock multiset of that cut."""
    eps = [1 if i in A else -1 for i in range(len(V))]
    return folded_clocks(V, eps, C) if fold else unfolded_clocks(V, eps)


def is_automorphism_collision(V, A1, A2):
    """Is the pair of cuts (A1,A2) explained by an affine automorphism x->a*x+b (mod nothing,
    over the integers) of the *signed configuration*?  We test the weakest reasonable notion:
    does there exist a bijection of the runner index set carrying the mono/cut edge pattern of A1
    to that of A2 while preserving the absolute speeds?  Since speeds are DISTINCT, any speed-
    preserving bijection is the identity; so a genuine 'automorphism' collision can only come from
    a value-level symmetry.  We therefore report whether the two cuts are related by GLOBAL
    structure: identical as partitions (trivially no), or related by the complement within a
    self-paired speed symmetry.  Practically: return True iff the UNFOLDED multisets also match
    (a real additive coincidence), False iff they differ (a fold-only artifact)."""
    return cut_signature(V, A1, 0, fold=False) == cut_signature(V, A2, 0, fold=False)


def analyze(V):
    n = len(V) + 1
    C = 2 * n - 1
    # enumerate cuts up to global swap: fix index 0 in '+' class
    cuts = []
    for bits in product([0, 1], repeat=len(V) - 1):
        A = frozenset([0] + [i + 1 for i, b in enumerate(bits) if b])
        cuts.append(A)
    # folded orbit
    fold_map = {}
    for A in cuts:
        sig = cut_signature(V, A, C, fold=True)
        fold_map.setdefault(sig, []).append(A)
    unfold_map = {}
    for A in cuts:
        sig = cut_signature(V, A, 0, fold=False)
        unfold_map.setdefault(sig, []).append(A)
    folded_orbit = len(fold_map)
    unfold_orbit = len(unfold_map)
    total_cuts = len(cuts)  # = 2^{n-2}
    # classify folded collisions
    fold_only = 0
    auto_real = 0
    coll_examples = []
    for sig, group in fold_map.items():
        if len(group) > 1:
            # within this folded-equal group, see how it splits under unfolded
            unf = {}
            for A in group:
                unf.setdefault(cut_signature(V, A, 0, fold=False), []).append(A)
            # pairs that are folded-equal but unfolded-DIFFERENT => fold-only
            if len(unf) > 1:
                fold_only += len(group) - len(unf)  # extra mergers due to fold
                coll_examples.append(('FOLD-ONLY', sorted(tuple(sorted(a)) for a in group)))
            else:
                auto_real += len(group) - 1
                coll_examples.append(('UNFOLDED-TOO', sorted(tuple(sorted(a)) for a in group)))
    shell = [(V[i], V[j]) for i in range(len(V)) for j in range(i + 1, len(V)) if (V[i] + V[j]) % C == 0]
    return dict(n=n, C=C, V=V, total_cuts=total_cuts, folded_orbit=folded_orbit,
               unfold_orbit=unfold_orbit, fold_collisions=total_cuts - folded_orbit,
               unfold_collisions=total_cuts - unfold_orbit, shell=shell,
               coll_examples=coll_examples)


def main():
    print("=" * 90)
    print("PART 1: AP_n sign-orbit -- reproduce + extend T4, split FOLD-ONLY vs UNFOLDED collisions")
    print("=" * 90)
    print(f"{'n':>3} {'C':>4} {'2^(n-2)':>8} {'folded':>7} {'unfold':>7} {'fold_coll':>9} {'unf_coll':>9}")
    seq = []
    for n in range(3, 17):
        V = list(range(1, n))  # AP_n = {1,...,n-1}
        r = analyze(V)
        seq.append(r['folded_orbit'])
        print(f"{n:>3} {r['C']:>4} {r['total_cuts']:>8} {r['folded_orbit']:>7} {r['unfold_orbit']:>7} "
              f"{r['fold_collisions']:>9} {r['unfold_collisions']:>9}")
    print(f"\nfolded sign-orbit sequence (n=3..10): {seq}")
    print("S699 T4 claimed: 2,4,7,16,32,60 for n=3..8")

    print("\n" + "=" * 90)
    print("PART 2: the n=5 AP collision dissected (which cuts collide, and WHY)")
    print("=" * 90)
    V = [1, 2, 3, 4]
    r = analyze(V)
    print(f"V={V} C={r['C']} folded_orbit={r['folded_orbit']} (= {r['total_cuts']} - {r['fold_collisions']})")
    for kind, group in r['coll_examples']:
        print(f"  [{kind}] colliding cuts (+ class shown): {group}")
        for A in group:
            Aset = set(A)
            eps = [1 if i in Aset else -1 for i in range(len(V))]
            print(f"      A+={str(sorted(V[i] for i in Aset)):<16} folded={folded_clocks(V,eps,r['C'])}  "
                  f"unfolded={unfolded_clocks(V,eps)}")

    print("\n" + "=" * 90)
    print("PART 3: where do collisions FIRST appear -- AP, and is it fold-only or a real automorphism?")
    print("=" * 90)
    for n in range(3, 17):
        V = list(range(1, n))
        r = analyze(V)
        comp = "PRIME" if all((2*n-1) % d for d in range(2, 2*n-1)) else "composite"
        kinds = [k for k, _ in r['coll_examples']]
        print(f"n={n}: C={2*n-1:>3} [{comp:>9}] fold_coll={r['fold_collisions']:>4} "
              f"unf_coll={r['unfold_collisions']:>4} "
              f"types={ {k:kinds.count(k) for k in set(kinds)} if kinds else 'none'}")

    print("\n" + "=" * 90)
    print("PART 4: exhaustive worry-set (tight M=1/n) census at small n + signed signatures")
    print("=" * 90)
    # enumerate all configs V (sorted distinct positive speeds) up to a speed bound; find tight ones
    for n in range(3, 8):
        C = 2 * n - 1
        target = Fraction(1, n)
        bound = 2 * n + 2  # search window for speeds
        tights = []
        for V in combinations(range(1, bound + 1), n - 1):
            V = list(V)
            if gap_exact(V) == target:
                tights.append(V)
        # dedup by gcd-primitive + report signed signature
        print(f"\n n={n} C={C}: tight configs (M=1/{n}) with speeds in [1,{bound}]:")
        seen_prim = set()
        for V in tights:
            g = 0
            for x in V:
                g = __import__('math').gcd(g, x)
            prim = tuple(x // g for x in V)
            if prim in seen_prim:
                continue
            seen_prim.add(prim)
            r = analyze(V)
            maxzero = 0
            for bits in product([0, 1], repeat=len(V) - 1):
                A = set([0] + [i + 1 for i, b in enumerate(bits) if b])
                eps = [1 if i in A else -1 for i in range(len(V))]
                z = sum(1 for x in folded_clocks(V, eps, C) if x == 0)
                maxzero = max(maxzero, z)
            tag = "SHELL!" if r['shell'] else "free"
            print(f"   V={V} prim={prim}: sign-orbit={r['folded_orbit']}/{r['total_cuts']} "
                  f"shell-partners={r['shell'] if r['shell'] else 'none'} max-zero-clocks={maxzero} [{tag}]")


if __name__ == "__main__":
    main()
