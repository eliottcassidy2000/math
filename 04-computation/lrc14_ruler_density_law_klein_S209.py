#!/usr/bin/env python3
"""
lrc14_ruler_density_law_klein_S209.py

HYP-5731(d): the LIVE-MULTIPLIER DENSITY LAW on pair-sum rulers, and the
resonance ledger that governs its deficit.

Claim to test: for a pair-sum ruler q of a covering 13-set S (q > Vmax so no
zero residue), the live-multiplier density is

    LM(q)/(q-1)  ≈  (6/7)^{c(q)}    with c(q) = #±-classes of {v_l mod q}
                                     (merges v_i ≡ -v_j lift the density),

and the DEFICIT from that baseline is governed by the LOW-HEIGHT RESONANCES
mod q: integer vectors n (small support/height) with  Σ n_l v_l ≡ 0 (mod q).
 - m = 0 relations (exact integer relations of S) hit EVERY ruler;
 - m ≠ 0 relations (Σ n_l v_l = m·q) pin q — each such relation taxes ~one ruler.

If confirmed, the a-priori certificate C5 = Fourier ledger mod q + budget
pigeonhole over the 78 rulers (covering ⟹ E2/E3 budget bounded away from the
AP extremum ⟹ some ruler resonance-thin ⟹ LM > 0 ⟹ lonely).

Measures, per ruler q (on adversarial worst instances + @91 + 2 census sets):
  d(q)    = LM/(q-1) − (6/7)^{c(q)}            (signed deficit)
  R3(q)   = #height-≤3 resonances mod q:
            pair-sums v_i+v_j ≡ 0, Schur v_i+v_j ≡ v_k, triples v_i+v_j+v_k ≡ 0,
            diff-pairs v_i+v_j ≡ v_k+v_l (height 4, sampled)   [weighted count]
  corr(d, R3) across rulers within each instance.
"""

from math import gcd
import itertools

QS = list(range(2, 15))


def LM_exact(S, q):
    def r_safe(r):
        return 14 * r >= q and 14 * (q - r) >= q
    bad = bytearray(q)
    bad[0] = 1
    reps = []
    seen = set()
    for v in S:
        r = v % q
        key = min(r, (q - r) % q)
        if key in seen:
            continue
        seen.add(key)
        reps.append(r)
    if any(r == 0 for r in reps):
        return None, None
    for r in reps:
        g = gcd(r, q)
        rr, qq = r // g, q // g
        inv = pow(rr, -1, qq)
        for m in range(qq):
            s = m * g
            if not r_safe(s):
                p0 = (m * inv) % qq
                for t in range(g):
                    bad[p0 + t * qq] = 1
    return q - sum(bad), len(reps)


def resonance_score(S, q):
    """weighted low-height resonance count mod q (excluding the defining
    pair-merges, which are counted via c(q) already):
       w=1.0 : v_i + v_j + v_k ≡ 0  (mod q)
       w=1.0 : v_i + v_j ≡ v_k      (mod q)   [Schur mod q]
       w=0.5 : v_i + v_j ≡ v_k + v_l (mod q)  [E2 mod q, i<j, k<l, disjoint-ish]
       w=0.5 : 2*v_i ≡ v_j + v_k or 2*v_i ≡ -v_j - v_k (mod q)
    """
    n = len(S)
    sc = 0.0
    vs = S
    for i in range(n):
        for j in range(i + 1, n):
            sij = vs[i] + vs[j]
            for k in range(n):
                if k == i or k == j:
                    continue
                if (sij + vs[k]) % q == 0:
                    sc += 1.0
                if (sij - vs[k]) % q == 0:
                    sc += 1.0
    # E2 mod q (sampled fully for 13 speeds: C(13,2)^2/2 ~ 3k -- fine)
    pair_sums = {}
    for i in range(n):
        for j in range(i + 1, n):
            pair_sums.setdefault((vs[i] + vs[j]) % q, []).append((i, j))
    for r, plist in pair_sums.items():
        m = len(plist)
        if m > 1:
            sc += 0.5 * (m * (m - 1) // 2)
    return sc


def analyze(S, name):
    S = sorted(S)
    V = S[-1]
    rulers = sorted({a + b for i, a in enumerate(S) for b in S[i:] if a + b > V})
    rows = []
    for q in rulers:
        lm, c = LM_exact(S, q)
        if lm is None:
            continue
        base = (6 / 7) ** c
        dens = lm / (q - 1)
        d = dens - base
        r3 = resonance_score(S, q)
        rows.append((q, c, dens, base, d, r3))
    # correlation d vs r3
    import statistics
    ds = [r[4] for r in rows]
    r3s = [r[5] for r in rows]
    if len(ds) > 2 and statistics.pstdev(ds) > 0 and statistics.pstdev(r3s) > 0:
        mu_d, mu_r = statistics.mean(ds), statistics.mean(r3s)
        cov = sum((a - mu_d) * (b - mu_r) for a, b in zip(ds, r3s)) / len(ds)
        corr = cov / (statistics.pstdev(ds) * statistics.pstdev(r3s))
    else:
        corr = float('nan')
    dens_all = [r[2] for r in rows]
    base_all = [r[3] for r in rows]
    print(f"\n{name}: V={V}, supra-Vmax rulers analyzed: {len(rows)}")
    print(f"  live density: mean {statistics.mean(dens_all):.4f} vs iid-class baseline "
          f"mean {statistics.mean(base_all):.4f}; |d| mean {statistics.mean(map(abs, ds)):.4f}")
    print(f"  corr(deficit, low-height resonance score) = {corr:+.3f}")
    worst = sorted(rows, key=lambda r: r[4])[:3]
    best = sorted(rows, key=lambda r: -r[4])[:3]
    print("  3 most-deficient rulers:  " + "; ".join(
        f"q={q} c={c} dens={de:.3f} base={ba:.3f} R3={r3:.1f}" for q, c, de, ba, dd, r3 in
        [(q, c, de, ba, dd, r3) for q, c, de, ba, dd, r3 in worst]))
    print("  3 most-enriched rulers:   " + "; ".join(
        f"q={q} c={c} dens={de:.3f} base={ba:.3f} R3={r3:.1f}" for q, c, de, ba, dd, r3 in best))
    lm_max = max((r[2] * (r[0] - 1)) for r in rows)
    print(f"  max LM = {int(round(lm_max))}; ALL rulers live: {all(r[2] > 0 for r in rows)}")


def main():
    print("=" * 78)
    print("The live-multiplier density law on pair-sum rulers (klein-S209)")
    print("=" * 78)
    cases = [
        ([12, 33, 46, 47, 68, 73, 79, 81, 85, 87, 91, 112, 120], "adversarial-worst V=120 (0 certs, maxLM=30)"),
        ([31, 33, 45, 48, 73, 76, 82, 86, 98, 102, 103, 104, 120], "adversarial V=120 (0 certs)"),
        ([62, 66, 69, 102, 109, 118, 120, 126, 130, 136, 159, 185, 200], "adversarial V=200 (0 certs)"),
        ([9, 16, 24, 33, 40, 47, 54, 62, 65, 70, 77, 84, 91], "7-structured @91 (hard cluster)"),
        ([8, 9, 10, 12, 14, 189, 200, 220, 230, 252, 260, 270, 280], "mid-band-free composite V=280 (S208)"),
    ]
    for S, name in cases:
        analyze(S, name)
    print("\nInterpretation: if corr < 0 consistently (more low-height resonances mod q")
    print("=> lower live density) and density ≈ (6/7)^c everywhere with small |d|,")
    print("then the C5 Fourier-ledger certificate + resonance-budget pigeonhole over")
    print("rulers is the right a-priori shape (m≠0 relations pin ~one ruler each).")
    print("DONE.")


if __name__ == '__main__':
    main()
