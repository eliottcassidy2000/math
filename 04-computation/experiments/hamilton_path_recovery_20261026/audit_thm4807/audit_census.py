"""Auditor's census for THM-4807 (a), (b), (c), (e), the three decoders, the two-step
statistics, the endpoint-relation collisions and several side claims.

Universes:
  * every isomorphism class n = 1..8 (nauty gentourng; upper-triangle strings, '1' = i->j);
  * every labelled tournament n = 1..6 (bit masks).
Hamilton words are enumerated explicitly by DFS and cross-checked against subset DP.
No code from the lane is imported.

usage: python3 audit_census.py [maxclass_n=8] [maxlabelled_n=6]
"""

import subprocess
import sys
from collections import Counter, defaultdict
from itertools import combinations

from audit_lib import (
    all_hps, bits_iter, check_tournament, comp_index, decode_G1, decode_G4,
    decode_endpoint_from_words, from_mask, from_upper_string, hp_tables, in_masks,
    is_ham_word, is_legal, is_strong, path_from_tables, popcount, reverse_arcs,
    scc_ordered, traced_from_tables, traced_from_words,
)

FAIL = Counter()


def fail(tag, detail=None):
    FAIL[tag] += 1
    if FAIL[tag] <= 3:
        print("  FAILURE", tag, detail)


def gentourng(n, strong_only=False):
    args = ["gentourng", "-q"] + (["-c"] if strong_only else []) + [str(n)]
    if n == 1:
        return ["" ]
    res = subprocess.run(args, capture_output=True, text=True, check=True)
    return [line.strip() for line in res.stdout.split("\n") if line.strip()]


def analyse(n, out, stats, labelled=False):
    """Run every check on one tournament; accumulate into stats."""
    assert check_tournament(n, out)
    full = (1 << n) - 1
    inn = in_masks(n, out)
    words = all_hps(n, out)
    stats["H_odd_fail"] += (len(words) % 2 == 0)
    F, B = hp_tables(n, out, inn)
    tr = traced_from_tables(n, out, F, B)
    trw = traced_from_words(n, words)
    if tr != trw:
        fail("traced DP != traced words", (n, out))
    starts = 0
    ends = 0
    for w in words:
        starts |= 1 << w[0]
        ends |= 1 << w[-1]
    if starts != B[full] or ends != F[full]:
        fail("DP endpoints != word endpoints", (n, out))

    comps, chain = scc_ordered(n, out)
    if not chain:
        fail("tournament condensation not a chain", (n, out))
    idx = comp_index(n, comps)
    # every arc between components points forward
    for x in range(n):
        for y in bits_iter(out[x]):
            if idx[x] > idx[y]:
                fail("backward cross arc", (n, out))
    # ordered-join factorisation: each word visits components contiguously in order
    sizes = [popcount(c) for c in comps]
    bounds = [0]
    for s in sizes:
        bounds.append(bounds[-1] + s)
    for w in words:
        for t in range(len(comps)):
            seg = w[bounds[t]:bounds[t + 1]]
            if set(seg) != set(bits_iter(comps[t])):
                fail("word not an ordered join", (n, out, w))
                break
    prod = 1
    comp_traced = [0] * n  # arcs traced inside the component's own Hamilton paths
    for c in comps:
        wc = all_hps(n, out, c)
        prod *= len(wc)
        for ww in wc:
            for a, b in zip(ww, ww[1:]):
                comp_traced[a] |= 1 << b
    if prod != len(words):
        fail("H(T) != prod H(C_i)", (n, out))
    # traced arcs: inside components iff traced in component; consecutive cross all; others none
    for x in range(n):
        for y in bits_iter(out[x]):
            t = (tr[x] >> y) & 1
            if idx[x] == idx[y]:
                tc = (comp_traced[x] >> y) & 1
                if t != tc:
                    fail("intra traced mismatch", (n, out, x, y))
            elif idx[y] == idx[x] + 1:
                if not t:
                    fail("consecutive cross arc untraced", (n, out, x, y))
            else:
                if t:
                    fail("non-consecutive cross arc traced", (n, out, x, y))

    # cut positions (e, step 2)
    cuts = [m for m in range(n + 1) if len({frozenset(w[:m]) for w in words}) == 1]
    if cuts != bounds:
        fail("cut positions != component boundaries", (n, out, cuts, bounds))
    # pairwise distinct prefix windows of a rotated Hamiltonian cycle (stronger form)
    for c in comps:
        L = popcount(c)
        if L < 3:
            continue
        cyc = None
        for ww in all_hps(n, out, c):
            if (out[ww[-1]] >> ww[0]) & 1:
                cyc = ww
                break
        if cyc is None:
            fail("strong component without Hamiltonian cycle (Camion)", (n, out))
            continue
        for r in range(1, L):
            wins = {frozenset(cyc[(s + i) % L] for i in range(r)) for s in range(L)}
            if len(wins) != L:
                fail("rotation windows not pairwise distinct", (n, out, r))

    # decoders
    D = tr
    dE, iE = decode_endpoint_from_words(n, words)
    dG1, iG1 = decode_G1(n, D)
    dG4, iG4 = decode_G4(n, D)
    stats["decoded_E"] += (dE == out and iE["bad_pairs"] == 0)
    stats["decoded_G1"] += (dG1 == out and iG1["bad_pairs"] == 0)
    stats["decoded_G4"] += (dG4 == out and iG4["bad_pairs"] == 0)
    stats["decoders_agree"] += (dE == dG1 == dG4)
    for (x, y, cxy, cyx) in iG1["stats"]:
        true_xy = (out[x] >> y) & 1
        tc, fc = (cxy, cyx) if true_xy else (cyx, cxy)
        stats["twostep"][(max(cxy, cyx), min(cxy, cyx))] += 1
        stats["twostep_true_false"][(tc, fc)] += 1
        if not (tc >= 2 and fc <= 1):
            fail("two-step rule fails", (n, out, x, y, tc, fc))
    for e in iG4["eligible"]:
        stats["eligible_sep"][e] += 1
    stats["untraced_intra_block_pairs"] += len(iG1["stats"])

    # endpoint relation (labelled only)
    if labelled:
        stats["endpoint_relations"].add(frozenset((w[0], w[-1]) for w in words))
        stats["hp_sets"].add(frozenset(words))

    # (a), (b), (c) for strong T
    if not (n >= 1 and len(comps) == 1):
        return
    stats["strong"] += 1
    untraced = [(x, y) for x in range(n) for y in bits_iter(out[x]) if not (tr[x] >> y) & 1]
    stats["untraced_arcs"] += len(untraced)
    stats["strong_with_untraced"] += bool(untraced)
    has_cut = any(not is_strong(n, out, full & ~(1 << z)) for z in range(n)) if n >= 2 else False
    if has_cut and not untraced:
        stats["strong_with_cut_vertex_but_no_untraced"] += 1
    if has_cut:
        stats["strong_with_cut_vertex"] += 1
    hp_pairs = Counter((w[0], w[-1]) for w in words)
    for (u, v) in untraced:
        # (a)
        cu, _ = scc_ordered(n, out, full & ~(1 << u))
        cv, _ = scc_ordered(n, out, full & ~(1 << v))
        if len(cu) < 2 or not (cu[-1] >> v) & 1:
            fail("(a) tail half", (n, out, u, v))
        if len(cv) < 2 or not (cv[0] >> u) & 1:
            fail("(a) dual half", (n, out, u, v))
        stats["k_of_T_minus_u"][len(cu)] += 1
        stats["size_Ck"][popcount(cu[-1])] += 1
        stats["size_C1_dual"][popcount(cv[0])] += 1
        if popcount(cu[-1]) < 3 or popcount(cv[0]) < 3:
            fail("(b) |C_k|>=3 or dual", (n, out, u, v))
        # (b) explicit chord path by the recipe
        c = next(bits_iter(out[u] & cu[0]), None)
        if c is None:
            fail("(b) u has no out-neighbour in C_1", (n, out, u, v))
            continue
        Pp = []
        for t, C in enumerate(cu):
            if t == 0:
                Pp += path_from_tables(out, inn, F, B, C, start=c)
            elif t == len(cu) - 1:
                Pp += path_from_tables(out, inn, F, B, C, end=v)
            else:
                Pp += path_from_tables(out, inn, F, B, C)
        if not is_ham_word(n, out, tuple([u] + Pp)) or Pp[-1] != v:
            fail("(b) uP' not a Hamilton path from u to v", (n, out, u, v))
        if hp_pairs.get((v, u), 0):
            fail("(b) Hamilton path from v to u exists", (n, out, u, v))
        stats["num_HP_u_to_v"][hp_pairs.get((u, v), 0)] += 1
        stats["new_3cycles"][popcount(out[u] & inn[v])] += 1
        others = [a for a in untraced if a != (u, v)]
        for r in range(len(others) + 1):
            for extra in combinations(others, r):
                TD = reverse_arcs(n, out, [(u, v)] + list(extra))
                if not is_legal(TD, Pp + [u]) or not (TD[u] >> Pp[0]) & 1:
                    fail("(b) uP'u not a Hamiltonian cycle of T_D", (n, out, u, v, extra))
                if is_legal(out, Pp + [u]):
                    fail("(b) P'u legal in T", (n, out, u, v))
                stats["b_flip_sets"] += 1
        # Busch Corollary 4.8 consequence
        du, dv = popcount(out[u]), popcount(out[v])
        stats["outdeg_gap"][du - dv] += 1
        if du < dv:
            fail("d+(u) < d+(v) at an untraced arc", (n, out, u, v))
        # Busch Theorem 4.12 structure (for information)
        bjr = 0
        for z in range(n):
            if z in (u, v):
                continue
            cz, _ = scc_ordered(n, out, full & ~(1 << z))
            if len(cz) >= 4 and (cz[0] >> u) & 1 and (cz[-1] >> v) & 1:
                ok_iv = (all((out[a] >> z) & 1 for a in bits_iter(cz[1])) and
                         all((out[z] >> a) & 1 for a in bits_iter(cz[-2])))
                bjr += 1 if ok_iv else 0
        stats["bjr_separators"][bjr] += 1
        if bjr == 0:
            fail("no BJR separator z", (n, out, u, v))
    # (c) on every arc of strong T
    for u in range(n):
        cu, _ = scc_ordered(n, out, full & ~(1 << u))
        for v in bits_iter(out[u]):
            stats["c_arcs"] += 1
            pred = False
            if len(cu) >= 2 and (cu[-1] >> v) & 1:
                Ck = cu[-1]
                cond2 = (cu[-2] & ~out[u]) == 0
                cond3 = True
                rest = Ck & ~(1 << v)
                A = rest
                while A:  # all nonempty subsets A of C_k - v
                    if (F[A] & inn[u]) and (B[Ck & ~A] >> v) & 1:
                        cond3 = False
                        break
                    A = (A - 1) & rest
                pred = cond2 and cond3
            actual = not (tr[u] >> v) & 1
            if pred != actual:
                fail("(c) characterisation wrong", (n, out, u, v, pred, actual))
            else:
                stats["c_exact"] += 1


def new_stats():
    return {
        "count": 0, "strong": 0, "untraced_arcs": 0, "strong_with_untraced": 0,
        "strong_with_cut_vertex": 0, "strong_with_cut_vertex_but_no_untraced": 0,
        "decoded_E": 0, "decoded_G1": 0, "decoded_G4": 0, "decoders_agree": 0,
        "twostep": Counter(), "twostep_true_false": Counter(), "eligible_sep": Counter(),
        "untraced_intra_block_pairs": 0, "endpoint_relations": set(), "hp_sets": set(),
        "k_of_T_minus_u": Counter(), "size_Ck": Counter(), "size_C1_dual": Counter(),
        "num_HP_u_to_v": Counter(), "new_3cycles": Counter(), "b_flip_sets": 0,
        "outdeg_gap": Counter(), "bjr_separators": Counter(), "c_arcs": 0, "c_exact": 0,
        "H_odd_fail": 0,
    }


def report(label, st, labelled):
    print(f"{label}: tournaments {st['count']}, strong {st['strong']}, "
          f"untraced arcs in strong {st['untraced_arcs']} (in {st['strong_with_untraced']} strong), "
          f"untraced intra-block pairs (all) {st['untraced_intra_block_pairs']}")
    print(f"    (c) exact on {st['c_exact']}/{st['c_arcs']} arcs of strong tournaments; "
          f"(b) flip sets checked {st['b_flip_sets']}; H odd failures {st['H_odd_fail']}")
    print(f"    decoders exact: endpoint E {st['decoded_E']}, two-step G1 {st['decoded_G1']}, "
          f"deleted-reachability G4 {st['decoded_G4']}, all three agree {st['decoders_agree']}")
    print(f"    two-step (max,min): {dict(sorted(st['twostep'].items()))}")
    print(f"    two-step (true,false): {dict(sorted(st['twostep_true_false'].items()))}")
    print(f"    G4 eligible separators per pair: {dict(sorted(st['eligible_sep'].items()))}")
    if st["untraced_arcs"]:
        print(f"    k = #comp(T-u): {dict(sorted(st['k_of_T_minus_u'].items()))}; |C_k|: "
              f"{dict(sorted(st['size_Ck'].items()))}; |C'_1| of T-v: {dict(sorted(st['size_C1_dual'].items()))}")
        print(f"    #Hamilton paths u->v (= #new Hamiltonian cycles after reversing u->v): "
              f"{dict(sorted(st['num_HP_u_to_v'].items()))}")
        print(f"    #new 3-cycles v->u->w->v after reversing u->v: {dict(sorted(st['new_3cycles'].items()))}")
        print(f"    d+(u)-d+(v): {dict(sorted(st['outdeg_gap'].items()))}; "
              f"#BJR separators z (Thm 4.12 incl. (iv)): {dict(sorted(st['bjr_separators'].items()))}")
    print(f"    strong with a cut vertex {st['strong_with_cut_vertex']}, of which with no untraced arc "
          f"{st['strong_with_cut_vertex_but_no_untraced']}")
    if labelled:
        print(f"    distinct HP sets {len(st['hp_sets'])}; distinct endpoint relations "
              f"{{(first,last)}} {len(st['endpoint_relations'])}")
    sys.stdout.flush()


def main():
    maxc = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    maxl = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    for n in range(1, maxc + 1):
        st = new_stats()
        lines = gentourng(n)
        for s in lines:
            out = from_upper_string(n, s) if n > 1 else [0]
            st["count"] += 1
            analyse(n, out, st)
        # cross-check the strong count with gentourng -c
        if n > 1:
            sc = len(gentourng(n, strong_only=True))
            if sc != st["strong"]:
                fail("strong count != gentourng -c", (n, sc, st["strong"]))
        report(f"classes n = {n}", st, False)
    for n in range(1, maxl + 1):
        st = new_stats()
        m = n * (n - 1) // 2
        for mask in range(1 << m):
            out = from_mask(n, mask)
            st["count"] += 1
            analyse(n, out, st, labelled=True)
        report(f"labelled n = {n}", st, True)
    print("FAILURES:", dict(FAIL) if FAIL else "none")


if __name__ == "__main__":
    main()
