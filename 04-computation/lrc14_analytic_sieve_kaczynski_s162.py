#!/usr/bin/env python3
"""S162: analytic sieve / Kaczynski synthesis for the LRC14 route.

This script is a synthesis aid, not a proof.  It mines the local research
workspace for analytic-number-theory proof motifs named in the prompt:
prime sums, Mobius/totient packet sums, large-sieve/circle-method lanes,
exponential sums, explicit formulas/saddle methods, smoothing choices, and the
Kaczynski boundary-function thread.  It also computes small explicit tables for
sum mu(n)/n and sum mu(n)^2/phi(n), since those weights are proof-design
signals rather than decorative analogies.

Tournament Analysis vertices are proof modules, not runners.  The pairwise
observable is a retention vector:

    LRC predicate retention,
    exact arithmetic packet retention,
    local/global decomposition strength,
    exponential-sum control,
    boundary/resonance exception handling,
    smoothing/adaptivity,
    auditability.

The switch orients A -> B when A wins more coordinates; ties follow the listed
Hamiltonian path.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
import argparse
import math
import re


REPO = Path(__file__).resolve().parents[1]
RESULT = REPO / "05-knowledge" / "results" / "lrc14_analytic_sieve_kaczynski_s162.out"

SCAN_DIRS = [
    REPO / "00-navigation",
    REPO / "01-canon" / "theorems",
    REPO / "04-computation",
    REPO / "05-knowledge" / "hypotheses",
    REPO / "07-reflections",
    REPO / "poke-forum" / "posts",
]

MOTIFS: dict[str, tuple[str, ...]] = {
    "prime_sum_goldbach": (
        "prime", "goldbach", "helfgott", "vinogradov", "hardy-littlewood",
        "lambda lambda", "von mangoldt",
    ),
    "mobius_mu_sums": (
        "mobius", "moebius", "mu(n)", "mu^2", "mu2", "mertens",
        "dirichlet convolution", "mobius inversion",
    ),
    "totient_phi_packets": (
        "totient", "phi(n)", "phi(", "euler-totient", "coprime density",
        "jordan", "primitive packet", "exact-period",
    ),
    "large_sieve_circle": (
        "large sieve", "circle method", "major arc", "minor arc",
        "quadratic sieve", "selberg", "upper-bound sieve",
    ),
    "exponential_sum_core": (
        "exponential sum", "character sum", "gauss sum", "ramanujan sum",
        "fourier", "toeplitz", "fejer", "polya-vinogradov", "weil",
    ),
    "smoothing_saddle_formula": (
        "smoothing", "smooth", "saddle", "parabolic cylinder",
        "explicit formula", "explicit formulas", "mellin", "laplace",
    ),
    "kaczynski_boundary": (
        "kaczynski", "bagemihl", "fatou", "boundary function",
        "curvilinear convergence", "ambiguous points", "nontangential",
    ),
    "lrc_source_kernel": (
        "lrc14", "lrc(14)", "source-kernel", "qdiv", "endpoint",
        "boundary-moment", "true-wide", "source-spectrum", "labelled packet",
    ),
}


@dataclass(frozen=True)
class ProofModule:
    name: str
    vector: tuple[int, ...]
    note: str


MODULES = [
    ProofModule(
        "lrc_labelled_source_kernel",
        (5, 5, 4, 4, 5, 4, 5),
        "keeps qdiv, exact M/Farey, Haar boundary, endpoint owners, C27/K33, dual exits",
    ),
    ProofModule(
        "exponential_sum_engine",
        (4, 4, 4, 5, 4, 4, 4),
        "the core analytic place where minor arcs, Fourier/Toeplitz, Ramanujan modes, and Weil/PV estimates live",
    ),
    ProofModule(
        "kaczynski_boundary_resonance",
        (4, 4, 3, 4, 5, 4, 4),
        "Fatou/Kaczynski boundary values plus Bagemihl ambiguous approach sets; LRC resonances are the exceptional arcs",
    ),
    ProofModule(
        "mobius_totient_packet_ledger",
        (3, 5, 4, 3, 3, 3, 5),
        "mu, phi, Ramanujan, and Div(D) packet laws before scalarization",
    ),
    ProofModule(
        "circle_large_sieve_decomposition",
        (3, 4, 5, 4, 3, 3, 4),
        "major/minor arcs, local products, and large-sieve upper bounds as proof decomposition",
    ),
    ProofModule(
        "smoothing_saddle_explicit_formula",
        (3, 3, 4, 4, 4, 5, 3),
        "choice of smoothing, saddle/explicit-formula kernels, and complex analytic localization",
    ),
    ProofModule(
        "upper_bound_quadratic_sieve",
        (2, 4, 4, 3, 2, 3, 3),
        "upper-bound sieve as safe majorant; useful only if residual side channels remain named",
    ),
    ProofModule(
        "raw_scalar_prime_sum",
        (1, 1, 2, 1, 1, 1, 5),
        "raw sums over primes or one numerical main term without packet labels",
    ),
]

TIE_PATH = [m.name for m in MODULES]


def normalize(text: str) -> str:
    return text.lower().replace("ö", "o").replace("μ", "mu").replace("ϕ", "phi")


def iter_markdown_and_code() -> list[Path]:
    files: list[Path] = []
    self_path = Path(__file__).resolve()
    for root in SCAN_DIRS:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if path.resolve() == self_path:
                continue
            if path.is_file() and path.suffix.lower() in {".md", ".py", ".out", ".txt"}:
                files.append(path)
    return files


def motif_census() -> tuple[Counter[str], Counter[tuple[str, str]], dict[str, list[Path]]]:
    counts: Counter[str] = Counter()
    co: Counter[tuple[str, str]] = Counter()
    examples: dict[str, list[Path]] = defaultdict(list)
    for path in iter_markdown_and_code():
        try:
            text = normalize(path.read_text(encoding="utf-8", errors="ignore"))
        except OSError:
            continue
        present: list[str] = []
        for motif, needles in MOTIFS.items():
            if any(needle in text for needle in needles):
                present.append(motif)
                counts[motif] += 1
                if len(examples[motif]) < 8:
                    examples[motif].append(path.relative_to(REPO))
        for a, b in combinations(sorted(present), 2):
            co[(a, b)] += 1
    return counts, co, examples


def mobius_phi_tables(limit: int = 20000) -> list[tuple[int, float, float, float, float]]:
    mu = [1] * (limit + 1)
    phi = list(range(limit + 1))
    is_prime = [True] * (limit + 1)
    primes: list[int] = []
    mu[0] = 0
    if limit >= 1:
        is_prime[0] = is_prime[1] = False
    for i in range(2, limit + 1):
        if is_prime[i]:
            primes.append(i)
            mu[i] = -1
            phi[i] -= phi[i] // i
        for p in primes:
            if i * p > limit:
                break
            is_prime[i * p] = False
            if i % p == 0:
                mu[i * p] = 0
                phi[i * p] = phi[i] * p
                break
            mu[i * p] = -mu[i]
            phi[i * p] = phi[i] * (p - 1)

    rows = []
    s_mu_over_n = 0.0
    s_mu2_over_phi = 0.0
    s_mu2_over_n = 0.0
    s_phi_over_n = 0.0
    checkpoints = {10, 30, 100, 300, 1000, 3000, 10000, limit}
    for n in range(1, limit + 1):
        s_mu_over_n += mu[n] / n
        if mu[n] != 0:
            s_mu2_over_phi += 1.0 / phi[n]
            s_mu2_over_n += 1.0 / n
        s_phi_over_n += phi[n] / n
        if n in checkpoints:
            rows.append((n, s_mu_over_n, s_mu2_over_phi, s_mu2_over_n, s_phi_over_n / n))
    return rows


def tournament_fingerprint() -> tuple[dict[int, int], int, list[int], list[list[str]], int, list[str]]:
    n = len(MODULES)
    scores = [0] * n
    adj = [[False] * n for _ in range(n)]
    tie_rank = {name: i for i, name in enumerate(TIE_PATH)}
    for i, a in enumerate(MODULES):
        for j, b in enumerate(MODULES):
            if i == j:
                continue
            wins = sum(x > y for x, y in zip(a.vector, b.vector))
            losses = sum(x < y for x, y in zip(a.vector, b.vector))
            if wins > losses or (wins == losses and tie_rank[a.name] < tie_rank[b.name]):
                adj[i][j] = True
        scores[i] = sum(adj[i])
    score_hist = dict(sorted(Counter(scores).items()))
    cycles = 0
    for i, j, k in combinations(range(n), 3):
        if adj[i][j] and adj[j][k] and adj[k][i]:
            cycles += 1
        if adj[i][k] and adj[k][j] and adj[j][i]:
            cycles += 1

    # Kosaraju SCC sizes.
    seen = [False] * n
    order: list[int] = []

    def dfs(v: int) -> None:
        seen[v] = True
        for w in range(n):
            if adj[v][w] and not seen[w]:
                dfs(w)
        order.append(v)

    for v in range(n):
        if not seen[v]:
            dfs(v)
    radj = [[adj[j][i] for j in range(n)] for i in range(n)]
    seen = [False] * n
    components: list[list[int]] = []

    def rdfs(v: int, comp: list[int]) -> None:
        seen[v] = True
        comp.append(v)
        for w in range(n):
            if radj[v][w] and not seen[w]:
                rdfs(w, comp)

    for v in reversed(order):
        if not seen[v]:
            comp: list[int] = []
            rdfs(v, comp)
            components.append(comp)
    components.sort(key=lambda comp: (-len(comp), min(comp)))
    sizes = [len(comp) for comp in components]
    scc_members = [[MODULES[i].name for i in comp] for comp in components]

    # Held-Karp count of Hamiltonian paths following oriented edges.
    dp: dict[tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp.get((mask, last), 0)
            if not val:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + val
    full = (1 << n) - 1
    hp = sum(dp.get((full, last), 0) for last in range(n))
    path = sorted(MODULES, key=lambda m: -scores[TIE_PATH.index(m.name)])
    return score_hist, cycles, sizes, scc_members, hp, [m.name for m in path]


def render() -> str:
    counts, co, examples = motif_census()
    mu_rows = mobius_phi_tables()
    score_hist, cycles, scc, scc_members, hp, path = tournament_fingerprint()

    lines: list[str] = []
    lines.append("LRC14 analytic sieve / Kaczynski synthesis (codex S162)")
    lines.append("=" * 63)
    lines.append("")
    lines.append("Motif census")
    lines.append("------------")
    for motif, count in counts.most_common():
        lines.append(f"{motif:30s} {count:5d}")
    lines.append("")
    lines.append("Top motif co-occurrences")
    lines.append("------------------------")
    for (a, b), count in co.most_common(20):
        lines.append(f"{a} + {b}: {count}")
    lines.append("")
    lines.append("Selected local anchors")
    lines.append("----------------------")
    for motif in sorted(MOTIFS):
        lines.append(f"{motif}:")
        for path_obj in examples.get(motif, [])[:6]:
            lines.append(f"  - {path_obj.as_posix()}")
    lines.append("")
    lines.append("Explicit Mobius/totient weight table")
    lines.append("------------------------------------")
    lines.append("N      sum_mu/n        sum_mu2/phi     sum_mu2/n      avg_phi/n")
    for n, a, b, c, d in mu_rows:
        lines.append(f"{n:5d}  {a: .8f}   {b: .8f}   {c: .8f}   {d: .8f}")
    lines.append("")
    lines.append("Readout:")
    lines.append("- sum mu(n)/n is a cancellation diagnostic; it drifts toward 0 and should not be used as a positive density.")
    lines.append("- sum mu(n)^2/phi(n) is a slow sieve-capacity ledger; it records squarefree primitive packets rather than endpoint ownership.")
    lines.append("- avg phi(n)/n stabilizes near 6/pi^2, the coprime-density floor behind Farey/exact-period packets.")
    lines.append("")
    lines.append("Tournament Analysis")
    lines.append("-------------------")
    lines.append("vertices: proof modules, not runners or primes")
    lines.append("pairwise observable: LRC/exact/local-global/exponential/boundary/smoothing/audit retention")
    lines.append(f"score_hist: {score_hist}")
    lines.append(f"directed_3cycles: {cycles}")
    lines.append(f"SCC_sizes: {scc}")
    lines.append("SCC_members:")
    for comp in scc_members:
        lines.append("  - " + ", ".join(comp))
    lines.append(f"Hamiltonian_path_count: {hp}")
    lines.append("tie Hamiltonian path:")
    lines.append("  " + " > ".join(TIE_PATH))
    lines.append("retention path:")
    lines.append("  " + " > ".join(path))
    lines.append("")
    lines.append("Module notes")
    lines.append("------------")
    for module in MODULES:
        lines.append(f"- {module.name}: vector={module.vector}; {module.note}")
    lines.append("")
    lines.append("Synthesis")
    lines.append("---------")
    lines.append("1. The Helfgott/Vinogradov template is not just 'add a variable'; it is main packet + exponential-sum minor-arc control + finite exceptional verification.")
    lines.append("2. For LRC14, qdiv/Farey/Ramanujan/totient data are the arithmetic main packet; endpoint owners and Haar boundary status are the singular-series labels.")
    lines.append("3. Kaczynski/Bagemihl boundary theory names the true-wide exceptional set: resonance approach arcs where Fatou-style decorrelated limits are insufficient.")
    lines.append("4. Large-sieve or upper-bound-sieve bounds are useful only as majorants after the residual packet labels survive; otherwise they repeat the scalar-quotient failure.")
    lines.append("5. The core missing lemma is an explicit exponential-sum/resonance estimate with adaptive smoothing: off-resonance minor arcs decay; resonant arcs reduce to finite Freiman/K33/AP-GW packets.")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="write the stored result file")
    args = parser.parse_args()
    text = render()
    print(text, end="")
    if args.write:
        RESULT.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
