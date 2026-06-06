"""SIGNED LRC -- the MECHANISM of fold-only collisions.  monad-explorer-S703d.

S703: AP_n folded collisions exist iff C=2n-1 composite, and are all FOLD-ONLY (unfolded faithful).
Here: WHAT is a collision, structurally?

A cut chi (2-coloring, fix elt 1 in '+').  foldedclock({i,j}) = |i-j| (mono) or rho(i+j) (cut),
rho(s)=min(s,C-s).  Test two structural claims:

  CLAIM-FLIP:  every folded-collision pair of cuts is a SINGLE-FLIP (Hamming distance 1 in the
               cut hypercube) -- i.e. collisions are 'silent flips' of ONE runner.

  CLAIM-EULER: flipping runner x is silent (preserves the whole folded multiset) iff the
               multigraph G_x on the value set, with one edge {|x-y|, rho(x+y)} for each y!=x,
               is EULERIAN (every vertex even degree).  [Because flipping x toggles every incident
               edge between its mono value |x-y| and its cut value rho(x+y); the multiset is
               preserved iff the 'before' and 'after' selections induce the same multiset, which
               for a set of value-pairs happens iff each value has even total incidence = G_x even.]

Also: does the silent-flip runner x, or the divisor structure, relate to a nontrivial divisor d|C?
"""
from itertools import product, combinations
from collections import Counter


def is_prime(m):
    d = 2
    while d * d <= m:
        if m % d == 0:
            return False
        d += 1
    return m >= 2


def folded_sig(V, col, C):
    n1 = len(V)
    c = Counter()
    for i in range(n1):
        for j in range(i + 1, n1):
            if col[i] != col[j]:
                f = (V[i] + V[j]) % C
                c[min(f, C - f)] += 1
            else:
                c[abs(V[i] - V[j])] += 1
    return tuple(sorted(c.items()))


def x_is_euler_silent(V, x_idx, C):
    """Eulerian test for flipping runner x (independent of current coloring)."""
    x = V[x_idx]
    deg = Counter()
    for j, y in enumerate(V):
        if j == x_idx:
            continue
        m = abs(x - y)
        f = (x + y) % C
        c = min(f, C - f)
        deg[m] += 1
        deg[c] += 1
    return all(v % 2 == 0 for v in deg.values())


def analyze_AP(n):
    C = 2 * n - 1
    V = list(range(1, n))
    n1 = len(V)
    groups = {}
    for bits in product((0, 1), repeat=n1 - 1):
        col = (0,) + bits
        groups.setdefault(folded_sig(V, col, C), []).append(col)
    colls = [g for g in groups.values() if len(g) > 1]
    # Hamming distances within groups (up to global swap: canonicalize so col[0]=0)
    def ham(a, b):
        # min over global swap
        d1 = sum(x != y for x, y in zip(a, b))
        d2 = sum(x == y for x, y in zip(a, b))
        return min(d1, d2)
    all_flip = True
    dists = Counter()
    flip_elts = Counter()
    for g in colls:
        for a, b in combinations(g, 2):
            d = ham(a, b)
            dists[d] += 1
            if d != 1:
                all_flip = False
            else:
                # which element flipped
                for k in range(n1):
                    if (a[k] != b[k]):
                        # account for possible global-swap representation
                        pass
                diffpos = [k for k in range(n1) if a[k] != b[k]]
                if len(diffpos) == 1:
                    flip_elts[V[diffpos[0]]] += 1
    # which runners x admit a silent (Eulerian) flip
    euler_x = [V[k] for k in range(n1) if x_is_euler_silent(V, k, C)]
    return C, len(colls), all_flip, dict(dists), dict(flip_elts), euler_x


def main():
    print("=" * 80)
    print("MECHANISM of AP_n fold-collisions: single-flip? Eulerian-value-multigraph?")
    print("=" * 80)
    print(f"{'n':>3} {'C':>4} {'prime':>6} {'#collGrp':>9} {'allSingleFlip':>14} {'hamDist':>16}")
    for n in range(3, 15):
        C, ncoll, allflip, dists, felts, eulx = analyze_AP(n)
        print(f"{n:>3} {C:>4} {str(is_prime(C)):>6} {ncoll:>9} {str(allflip):>14} {str(dists):>16}")
    print("\nSilent-flip runners (Eulerian value-multigraph) and which elements flip in collisions:")
    for n in range(3, 15):
        C, ncoll, allflip, dists, felts, eulx = analyze_AP(n)
        if ncoll:
            print(f" n={n} C={C}={'prime' if is_prime(C) else 'composite'}: "
                  f"Eulerian-silent runners x = {eulx};  flipped-in-collisions = {felts}")
    print("\nDivisor check: for each composite C, do the Eulerian-silent runners x relate to d|C?")
    for n in range(3, 15):
        C = 2 * n - 1
        if is_prime(C):
            continue
        V = list(range(1, n))
        eulx = [V[k] for k in range(len(V)) if x_is_euler_silent(V, k, C)]
        divs = [d for d in range(2, C) if C % d == 0]
        print(f" n={n} C={C} divisors={divs}: Eulerian-silent x={eulx}  "
              f"(x mod small divisors: {[(x, [x % d for d in divs]) for x in eulx]})")


if __name__ == "__main__":
    main()
