"""Auditor's explicit witnesses for the corrections to THM-4807's side claims.

 (W1) The traced core alone determines T:  HP(D(T)) = HP(T) for the trace-support digraph
      D(T), and D(T) is injective on labelled tournaments (checked n <= 6).  This refutes
      "neither half suffices alone: the traced core forgets the orientation of untraced pairs".
 (W2) "Reversing it creates exactly one new cycle" is false: explicit classes where reversing
      the untraced arc u->v creates several new Hamiltonian cycles (one per Hamilton path u->v)
      and several new 3-cycles v->u->w->v.
 (W3) "The untraced arcs are exactly the places where a strong block is one vertex away from
      a factorisation": the converse fails (the strong 4-tournament; counts n <= 8).
 (W4) Case (gamma) uses at most two skipping arcs, sometimes none (i = 1 and j = k).
 (W5) 6,412,788 = 178,133 x 36 is the arc count of the strong classes at n = 9 alone.
 (W6) The n = 5 example, with all of (a)-(e) displayed.
No code from the lane is imported.
"""

import subprocess
from collections import Counter

from audit_lib import (
    all_hps, bits_iter, from_mask, from_upper_string, hp_tables, in_masks, is_strong,
    popcount, reverse_arcs, scc_ordered, traced_from_tables, traced_from_words,
)


def gentourng(n, strong=False):
    args = ["gentourng", "-q"] + (["-c"] if strong else []) + [str(n)]
    res = subprocess.run(args, capture_output=True, text=True, check=True)
    return [s.strip() for s in res.stdout.split("\n") if s.strip()]


def hp_of_digraph(n, adj):
    return set(all_hps(n, adj))


def w1():
    print("(W1) traced core D(T) versus HP(T), all labelled tournaments")
    for n in range(1, 7):
        m = n * (n - 1) // 2
        seen = {}
        bad_hp = 0
        for mask in range(1 << m):
            out = from_mask(n, mask)
            words = all_hps(n, out)
            D = traced_from_words(n, words)
            if hp_of_digraph(n, D) != set(words):
                bad_hp += 1
            key = tuple(D)
            seen.setdefault(key, []).append(mask)
        coll = sum(1 for v in seen.values() if len(v) > 1)
        print(f"    n = {n}: {1 << m} tournaments, {len(seen)} distinct trace supports D(T), "
              f"colliding supports {coll}, HP(D(T)) != HP(T) in {bad_hp} cases")


def w2():
    print("(W2) new cycles created by reversing an untraced arc u->v (classes n = 5..8)")
    shown = 0
    for n in range(5, 9):
        dist = Counter()
        for s in gentourng(n, strong=True):
            out = from_upper_string(n, s)
            inn = in_masks(n, out)
            F, B = hp_tables(n, out, inn)
            tr = traced_from_tables(n, out, F, B)
            words = None
            for u in range(n):
                for v in bits_iter(out[u] & ~tr[u]):
                    if words is None:
                        words = all_hps(n, out)
                    TD = reverse_arcs(n, out, [(u, v)])
                    # Hamiltonian cycles of T_D through v->u  <->  Hamilton words u ... v of T
                    new_ham = sum(1 for w in words if w[0] == u and w[-1] == v)
                    # brute-force check: count Hamiltonian cycles of T_D containing v->u
                    cyc = sum(1 for w in all_hps(n, TD) if w[0] == u and w[-1] == v)
                    new3 = popcount(out[u] & inn[v])
                    dist[(new_ham, new3)] += 1
                    assert new_ham == cyc
                    if shown < 2 and new_ham > 1:
                        shown += 1
                        print(f"    example n = {n}: gentourng string {s}, untraced arc {u}->{v}: reversing it "
                              f"creates {new_ham} new Hamiltonian cycles and {new3} new 3-cycles")
        print(f"    n = {n}: (new Hamiltonian cycles, new 3-cycles) over untraced arcs of strong classes: "
              f"{dict(sorted(dist.items()))}")


def w3():
    print("(W3) strong classes that have a cut vertex but no untraced arc")
    for n in range(3, 9):
        tot = 0
        cut_no_untraced = 0
        for s in gentourng(n, strong=True):
            out = from_upper_string(n, s)
            full = (1 << n) - 1
            has_cut = any(not is_strong(n, out, full & ~(1 << z)) for z in range(n))
            F, B = hp_tables(n, out)
            tr = traced_from_tables(n, out, F, B)
            untr = any(out[u] & ~tr[u] for u in range(n))
            tot += has_cut
            cut_no_untraced += (has_cut and not untr)
        print(f"    n = {n}: strong classes with a cut vertex {tot}, of which without an untraced arc {cut_no_untraced}")


def w4():
    print("(W4) case gamma at n = 3: T transitive 0->1->2 (and 0->2), T' = T with 0->2 reversed")
    out = [0b110, 0b100, 0]  # 0->1, 0->2, 1->2
    Tp = reverse_arcs(3, out, [(0, 2)])
    comps, _ = scc_ordered(3, out)
    print(f"    components {[sorted(bits_iter(c)) for c in comps]}; i = 1, j = 3 = k; word Q_2 [Q_3 ending at 2] "
          f"[Q_1 starting at 0] = (1, 2, 0): legal in T' {all((Tp[a] >> b) & 1 for a, b in [(1, 2), (2, 0)])}, "
          f"legal in T {all((out[a] >> b) & 1 for a, b in [(1, 2), (2, 0)])}; skipping arcs used: none")


def w5():
    print("(W5) arcs of strong classes")
    total = 0
    for n in range(1, 10):
        k = len(gentourng(n, strong=True)) if n > 1 else 1
        arcs = k * n * (n - 1) // 2
        total += arcs
        print(f"    n = {n}: strong classes {k}, arcs {arcs}, cumulative {total}")


def w6():
    print("(W6) the n = 5 strong class with an untraced arc")
    for s in gentourng(5, strong=True):
        out = from_upper_string(5, s)
        F, B = hp_tables(5, out)
        tr = traced_from_tables(5, out, F, B)
        untr = [(u, v) for u in range(5) for v in bits_iter(out[u] & ~tr[u])]
        if not untr:
            continue
        words = all_hps(5, out)
        (u, v), = untr
        full = 31
        cu, _ = scc_ordered(5, out, full & ~(1 << u))
        cv, _ = scc_ordered(5, out, full & ~(1 << v))
        print(f"    gentourng string {s}; out-neighbourhoods "
              f"{[sorted(bits_iter(out[i])) for i in range(5)]}; H = {len(words)}")
        print(f"    untraced arc {u}->{v}; T-u components {[sorted(bits_iter(c)) for c in cu]}; "
              f"T-v components {[sorted(bits_iter(c)) for c in cv]}")
        print(f"    Hamilton words from {u} to {v}: {[w for w in words if w[0] == u and w[-1] == v]}; "
              f"from {v} to {u}: {[w for w in words if w[0] == v and w[-1] == u]}")
        print(f"    all Hamilton words: {words}")


if __name__ == "__main__":
    w1()
    w2()
    w3()
    w4()
    w5()
    w6()
