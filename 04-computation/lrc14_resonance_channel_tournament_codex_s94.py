#!/usr/bin/env python3
r"""
lrc14_resonance_channel_tournament_codex_s94.py

First scout for HYP-2867.

The low part of the Node-3 spectrum sum is

    SPEC_low(H) = sum_{1<=n<=H} 2 Re(chat(n) conj(ghat(n))).

This script treats the sign-paired channels {n,-n} as tournament vertices.
The tie Hamiltonian order is increasing n.  For i<j we orient i->j when the
pair contribution s_i+s_j is nonnegative, and j->i otherwise.  This is not
claimed as a canonical tournament invariant; it is a compact diagnostic for
whether low resonance signs are incoherent/cyclic or almost all aligned
against the proof.

It also compresses channels by residue mod 7.  Residue 0 is reported as the
apex-prime trunk.  The nonzero residues are split by the Paley cut

    chi_7 = + on QR(7)={1,2,4},  - on NQR(7)={3,5,6}.

This is the finite object suggested by HYP-+2866 and HYP-2867.
"""

import itertools
import sys
from collections import Counter, defaultdict
from math import pi

sys.path.insert(0, "04-computation")
from lrc14_spectrum_intersection_sum_kpswf12 import (  # noqa: E402
    complement,
    cover_set,
    fourier_num_of_arcs,
    intersect,
    meas,
    safe_set,
)


QR7 = {1, 2, 4}
NQR7 = {3, 5, 6}


def coeff(arcs, n):
    return fourier_num_of_arcs(arcs, n) / (2j * pi * n)


def channel_terms(P, E, H):
    gp = safe_set(P)
    covc = complement(cover_set(E))
    baseline = meas(gp) * meas(covc)
    inter = meas(intersect(gp, covc))
    spec_exact = inter - baseline
    terms = []
    for n in range(1, H + 1):
        s = 2.0 * (coeff(gp, n) * coeff(covc, n).conjugate()).real
        terms.append((n, n % 7, s))
    return gp, covc, baseline, spec_exact, terms


def additive_profile(E):
    E = list(E)
    sums = Counter(a + b for a in E for b in E)
    energy = sum(v * v for v in sums.values())
    sumset = len(sums)
    k2 = sumset / max(1, len(E))
    return energy, sumset, k2


def build_tournament(weights, order=None):
    """Return adjacency for ordered vertices.  Edge follows order iff pair sum >= 0."""
    if order is None:
        order = list(weights)
    adj = {v: set() for v in order}
    for i, a in enumerate(order):
        for b in order[i + 1 :]:
            if weights[a] + weights[b] >= 0:
                adj[a].add(b)
            else:
                adj[b].add(a)
    return adj


def score_hist(adj):
    return tuple(sorted(Counter(len(v) for v in adj.values()).items()))


def directed_triangles(adj):
    verts = list(adj)
    c = 0
    for a, b, d in itertools.combinations(verts, 3):
        edges = ((b in adj[a]), (d in adj[b]), (a in adj[d]))
        # A 3-vertex tournament is cyclic iff all three edges align around one
        # cyclic order or all three align around the reverse cyclic order.
        if all(edges) or not any(edges):
            c += 1
    return c


def scc_sizes(adj):
    verts = list(adj)
    radj = {v: set() for v in verts}
    for a in verts:
        for b in adj[a]:
            radj[b].add(a)

    seen = set()
    order = []

    def dfs(v):
        seen.add(v)
        for w in adj[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in verts:
        if v not in seen:
            dfs(v)

    seen.clear()
    sizes = []

    def rdfs(v, comp):
        seen.add(v)
        comp.append(v)
        for w in radj[v]:
            if w not in seen:
                rdfs(w, comp)

    for v in reversed(order):
        if v not in seen:
            comp = []
            rdfs(v, comp)
            sizes.append(len(comp))
    return tuple(sorted(sizes, reverse=True))


def ham_path_count(adj):
    verts = list(adj)
    n = len(verts)
    if n > 16:
        return None
    idx = {v: i for i, v in enumerate(verts)}
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            a = verts[last]
            for b in adj[a]:
                j = idx[b]
                if not (mask & (1 << j)):
                    dp[mask | (1 << j)][j] += val
    return sum(dp[-1])


def fingerprint(weights, order=None):
    adj = build_tournament(weights, order)
    return {
        "score_hist": score_hist(adj),
        "c3": directed_triangles(adj),
        "scc": scc_sizes(adj),
        "hp": ham_path_count(adj),
    }


def analyze_case(P, E, label, H=21):
    gp, covc, baseline, spec_exact, terms = channel_terms(P, E, H)
    spec_low = sum(s for _, _, s in terms)
    residue = defaultdict(float)
    for _, r, s in terms:
        residue[r] += s

    nonzero_sum = sum(residue[r] for r in QR7 | NQR7)
    paley_cut = sum(residue[r] for r in QR7) - sum(residue[r] for r in NQR7)
    qr_sum = sum(residue[r] for r in QR7)
    nqr_sum = sum(residue[r] for r in NQR7)
    trunk0 = residue[0]
    energy, sumset, k2 = additive_profile(E)

    weights_n = {n: s for n, _, s in terms}
    fp_n = fingerprint(weights_n, [n for n, _, _ in terms])
    residue_order = [1, 2, 4, 3, 5, 6, 0]
    weights_r = {r: residue[r] for r in residue_order}
    fp_r = fingerprint(weights_r, residue_order)

    baseline_f = float(baseline)
    print("=" * 96)
    print(f"{label}")
    print(f"P={P}  E={E}  H={H}")
    print(f"baseline={baseline} ({baseline_f:.8f})  SPEC_exact={spec_exact} ({float(spec_exact):+.8f})")
    print(f"SPEC_low={spec_low:+.8f}  low/baseline={spec_low / baseline_f:+.6f}")
    print(
        "residue sums mod 7: "
        + " ".join(f"r{r}={residue[r]:+.6f}" for r in [0, 1, 2, 3, 4, 5, 6])
    )
    print(
        f"QR_sum={qr_sum:+.8f}  NQR_sum={nqr_sum:+.8f}  "
        f"PaleyCut={paley_cut:+.8f}  PaleyCut/baseline={paley_cut / baseline_f:+.6f}  "
        f"nonzero/baseline={nonzero_sum / baseline_f:+.6f}  "
        f"residue0_trunk/baseline={trunk0 / baseline_f:+.6f}"
    )
    print(f"additive profile: energy={energy}  |E+E|={sumset}  K2={k2:.4f}")
    print(
        "channel tournament: "
        f"score_hist={fp_n['score_hist']} c3={fp_n['c3']} "
        f"scc={fp_n['scc']} hp={fp_n['hp'] if fp_n['hp'] is not None else 'skipped(n>16)'}"
    )
    print(
        "residue tournament: "
        f"score_hist={fp_r['score_hist']} c3={fp_r['c3']} "
        f"scc={fp_r['scc']} hp={fp_r['hp']}"
    )
    worst = sorted(terms, key=lambda row: row[2])[:5]
    best = sorted(terms, key=lambda row: row[2], reverse=True)[:5]
    print("most negative channels: " + ", ".join(f"n={n}(r{r}) {s:+.6f}" for n, r, s in worst))
    print("most positive channels: " + ", ".join(f"n={n}(r{r}) {s:+.6f}" for n, r, s in best))
    return {
        "label": label,
        "low_over_base": spec_low / baseline_f,
        "nonzero_over_base": nonzero_sum / baseline_f,
        "paley_over_base": paley_cut / baseline_f,
        "trunk0_over_base": trunk0 / baseline_f,
        "c3": fp_n["c3"],
        "scc": fp_n["scc"],
        "energy": energy,
        "k2": k2,
    }


def main():
    print("# LRC14 resonance-channel tournament scout (codex S94)")
    print("# Vertices are sign-paired low Fourier channels {n,-n}; H=21.")
    print("# Edge follows increasing n iff the two-channel signed sum is nonnegative.")
    cases = [
        ([1, 2, 3, 12, 13], list(range(8)), "k=8 consec / routed bank"),
        ([1, 2, 3, 12], list(range(9)), "k=9 consec / routed bank"),
        ([1, 3, 4, 5], list(range(9)), "exact consec R' floor row"),
        ([1, 2, 3], list(range(10)), "k=10 consec / routed bank"),
        ([5, 7, 11], list(range(10)), "coprime P independence-favourable"),
        ([1, 2, 6], [0, 4, 6, 8, 10, 12, 14, 15, 16, 17], "wide d>=2 routed row"),
    ]
    rows = [analyze_case(P, E, label) for P, E, label in cases]
    print("=" * 96)
    print("SUMMARY")
    print(
        f"{'case':<38}{'low/base':>12}{'paley/base':>13}"
        f"{'nonzero':>11}{'r0/base':>11}{'c3':>8}{'scc':>28}{'K2':>9}"
    )
    for row in rows:
        print(
            f"{row['label'][:38]:<38}{row['low_over_base']:>12.6f}"
            f"{row['paley_over_base']:>13.6f}{row['nonzero_over_base']:>11.6f}"
            f"{row['trunk0_over_base']:>11.6f}{row['c3']:>8}"
            f"{str(row['scc']):>28}{row['k2']:>9.4f}"
        )
    print("\nInterpretation:")
    print("  The low sum is negative in the AP-like rows, but the Paley cut is not the")
    print("  damage source there: it is positive in those cases.  The negative mass splits")
    print("  between the nonzero-shell mean and the residue-0 trunk.  The")
    print("  residue-compressed tournament is small enough to track exactly, while the")
    print("  21-channel tournament records whether signs are cyclic/incoherent or nearly")
    print("  transitive.  HYP-2867's next proof target is to bound the nonzero mean plus")
    print("  residue-0 trunk, using PaleyCut as a balance witness and routing coherent")
    print("  sign words to Freiman/GAP atlases.")


if __name__ == "__main__":
    main()
