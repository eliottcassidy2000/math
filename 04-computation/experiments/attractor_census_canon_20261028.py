"""opus-2026-10-28-S68: keyword census of the canon - which concept family bridges the others? (HYP-9437)

The lexicon below is FIXED and was committed with the HYP-9437 pre-registration before the first run.

Each theorem file is read: title, front matter and the first BODY_CHARS characters. It is assigned the set of families
whose regex matches, case-insensitively. Co-occurrence counts n_ij are compared with a degree-preserving bipartite
null model: curveball swaps keep, for every theorem, how many families it touches, and for every family how many
theorems touch it. This controls the length effect, since long files mention everything.

z_ij = (n_ij - mean_null) / sd_null.
lift-degree of family i = #{ j != i : z_ij >= 2 }.

Usage: python3 attractor_census_canon_20261028.py [ID_MIN] [ID_MAX] [NULL_SAMPLES]
Defaults: 4400 4839 2000 (the pre-registered recent window).
"""
import glob
import math
import os
import random
import re
import sys

ID_MIN = int(sys.argv[1]) if len(sys.argv) > 1 else 4400
ID_MAX = int(sys.argv[2]) if len(sys.argv) > 2 else 4839
NULL = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
BODY_CHARS = 6000
SEED = 20261028

FAMILIES = {
    "MODULAR": r"\bPSL|\bSL ?_?2|SL\(2|PGL|modular group|modular curve|modular form|modular function|\bFarey|Stern.?Brocot"
               r"|Calkin.?Wilf|theta group|Gamma ?\(2\)|Gamma_0|\bX_0\(|\bX\(\d|triangle group|Hauptmodul|\bcusp|newform"
               r"|Hecke|Berggren|Pythagorean|Lorentz|light.?cone|\bPell\b|Markov|continued fraction|convergent|Minkowski"
               r"|Dedekind eta|\beta-quotient|\beta\(|j-invariant|Klein quartic|icosahed|Eisenstein series",
    "QUADRATIC": r"class number|genus theory|Gaussian integer|Gaussian prime|\bZ\[i\]|Eisenstein integer|Z\[omega\]"
                 r"|Z\[w\]|quadratic form|idoneal|Heegner|norm form|imaginary quadratic|real quadratic|Hilbert class"
                 r"|sums? of two squares",
    "FINITE_GEOMETRY": r"\bPaley|quadratic residue|Legendre symbol|\bFano|Heawood|Petersen|Kneser|Steiner|Hamming"
                       r"|\bK_?5\b|K_?\{?3, ?3\}?|Kuratowski|Clebsch|strongly regular|projective plane|finite field"
                       r"|\bF_q\b|\bGF\(",
    "BINARY_ALGEBRA": r"matroid|\bF_2\b|F2-linear|binary code|parity|Walsh|Hadamard|character|square class|syndrome"
                      r"|\bXOR\b|Boolean",
    "DIVISOR": r"\bsigma\(|divisor|aliquot|F ?= ?U ?\+ ?S|perfect number|abundant|totient|Moebius|Möbius|Mobius"
               r"|Liouville|squarefree|radical",
    "COLLATZ": r"Collatz|Syracuse|3x ?\+ ?1|3n ?\+ ?1|2-adic|3-adic|log_?2 ?\(?3|parity vector|Terras|stopping time",
    "TOURNAMENT": r"tournament|Hamiltonian path|Redei|Rédei|score sequence|\bH\(T\)",
    "ZETA": r"Bernoulli|\bzeta|Apery|Apéry|L-function|L-value|Dirichlet|Kloosterman|Gauss sum|Jacobi sum|polylog"
            r"|Riemann hypothesis",
    "ABC": r"abc[- ](conjecture|triple|hit|quality)|\bABC\b|Masser|Oesterl|\brad\(|Baker|S-unit|Pillai|Mih[aă]ilescu"
           r"|Szpiro|Mason.?Stothers|Fermat.?Catalan|\bBeal\b|Catalan'?s? conjecture",
    "OBSERVER": r"observer|guarded|carrier|Myhill|Nerode|minimal automaton|deterministic quotient",
}
REGEX = {k: re.compile(v, re.I) for k, v in FAMILIES.items()}
NAMES = list(FAMILIES)


def load(id_min, id_max):
    rows = []
    for f in sorted(glob.glob(os.path.join("01-canon", "theorems", "THM-*.md"))):
        m = re.match(r"THM-(\d+)", os.path.basename(f))
        if not m:
            continue
        i = int(m.group(1))
        if not (id_min <= i <= id_max):
            continue
        if "RESERVED" in os.path.basename(f):
            continue
        txt = open(f, encoding="utf-8", errors="replace").read()
        fam = frozenset(k for k in NAMES if REGEX[k].search(txt[:BODY_CHARS]))
        rows.append((i, os.path.basename(f), fam))
    return rows


def cooc(sets):
    n = {k: 0 for k in NAMES}
    nn = {}
    for s in sets:
        for a in s:
            n[a] += 1
        ss = sorted(s)
        for x in range(len(ss)):
            for y in range(x + 1, len(ss)):
                nn[(ss[x], ss[y])] = nn.get((ss[x], ss[y]), 0) + 1
    return n, nn


def curveball(sets, rng, rounds):
    """degree-preserving bipartite randomisation (Strona et al. curveball)."""
    S = [set(s) for s in sets]
    m = len(S)
    for _ in range(rounds):
        a, b = rng.randrange(m), rng.randrange(m)
        if a == b:
            continue
        A, B = S[a], S[b]
        onlyA, onlyB = list(A - B), list(B - A)
        if not onlyA or not onlyB:
            continue
        pool = onlyA + onlyB
        rng.shuffle(pool)
        k = len(onlyA)
        common = A & B
        S[a] = common | set(pool[:k])
        S[b] = common | set(pool[k:])
    return S


def main():
    rows = load(ID_MIN, ID_MAX)
    sets = [r[2] for r in rows]
    N = len(sets)
    n, nn = cooc(sets)
    print("window THM-%d..THM-%d: %d theorem files (RESERVED stubs excluded)" % (ID_MIN, ID_MAX, N))
    print("family sizes (share of files):")
    for k in sorted(NAMES, key=lambda k: -n[k]):
        print("  %-16s %4d  (%5.1f%%)" % (k, n[k], 100.0 * n[k] / max(N, 1)))
    multi = sum(1 for s in sets if len(s) >= 2)
    print("files touching >= 2 families: %d (%.1f%%); mean families per file %.2f"
          % (multi, 100.0 * multi / max(N, 1), sum(len(s) for s in sets) / max(N, 1)))
    rng = random.Random(SEED)
    pairs = [(NAMES[x], NAMES[y]) if NAMES[x] < NAMES[y] else (NAMES[y], NAMES[x])
             for x in range(len(NAMES)) for y in range(x + 1, len(NAMES))]
    acc = {p: [] for p in pairs}
    cur = [set(s) for s in sets]
    for t in range(NULL):
        cur = curveball(cur, rng, 5 * N)
        _, nn0 = cooc(cur)
        for p in pairs:
            acc[p].append(nn0.get(p, 0))
    z = {}
    for p in pairs:
        xs = acc[p]
        mu = sum(xs) / len(xs)
        sd = math.sqrt(sum((x - mu) ** 2 for x in xs) / max(len(xs) - 1, 1))
        z[p] = (nn.get(p, 0) - mu) / sd if sd > 0 else 0.0
    print("\nall 45 pairs by z (observed n_ij, lift, z):")
    for p in sorted(pairs, key=lambda p: -z[p]):
        lift = nn.get(p, 0) * N / (n[p[0]] * n[p[1]]) if n[p[0]] and n[p[1]] else float("nan")
        print("  %-16s %-16s n=%4d  lift=%5.2f  z=%6.2f" % (p[0], p[1], nn.get(p, 0), lift, z[p]))
    deg = {k: sum(1 for p in pairs if k in p and z[p] >= 2) for k in NAMES}
    strength = {k: sum(z[p] for p in pairs if k in p and z[p] > 0) for k in NAMES}
    print("\nlift-degree (# partners with z >= 2) and positive-z strength:")
    for k in sorted(NAMES, key=lambda k: (-deg[k], -strength[k])):
        print("  %-16s degree %d   strength %6.2f" % (k, deg[k], strength[k]))
    top3 = sorted(pairs, key=lambda p: -z[p])[:3]
    core = {"MODULAR", "QUADRATIC", "FINITE_GEOMETRY"}
    print("\nPRE-REGISTERED CHECKS (HYP-9437), window THM-%d..%d:" % (ID_MIN, ID_MAX))
    p1 = n["MODULAR"] >= 0.20 * N
    p2 = deg["MODULAR"] == max(deg.values())
    p3 = any(set(p) <= core for p in top3)
    print("  P1 MODULAR touches >= 20%% of files: %s (%.1f%%)" % (p1, 100.0 * n["MODULAR"] / N))
    print("  P2 MODULAR has the largest lift-degree (ties count): %s (MODULAR %d, max %d: %s)"
          % (p2, deg["MODULAR"], max(deg.values()), [k for k in NAMES if deg[k] == max(deg.values())]))
    print("  P3 a top-3 pair lies inside {MODULAR, QUADRATIC, FINITE_GEOMETRY}: %s (top3 = %s)" % (p3, top3))


if __name__ == "__main__":
    main()
