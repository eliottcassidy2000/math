#!/usr/bin/env python3
"""
Lane tao_minus_sheet_ladder (wave 6, 2026-09-22).

(A) Tao, arXiv:1909.03562 ("Almost all orbits of the Collatz map attain almost
    bounded values"): fetch the abstract page and the PDF (cached in the scratch
    directory), convert with pdftotext, and count the terms that the pasted
    conversation attributes to the paper.  Only counts, line numbers and one
    short quotation are printed (the paper is not reproduced).
(B) The minus-sheet Syracuse offset: Syr_-^n(N) = 3^n 2^{-|a|} N - F_n(a) with the
    same offset polynomial F_n as the plus sheet, and the valuation-word counts
    of the two sheets coincide exactly.  Hence Syrac_-(Z/3^nZ) = -Syrac_+(Z/3^nZ).
(C) The "Delta = 4 prime ladder" 3,7,11,17, the pasted Lean finite statement, the
    prime-target claim, the PPT (77,36,85), and the inherited identity
    (2^11 - 3^7)(-17) = 2363 = 17*139 (cycle gate n_0 = bB/Delta).
(D) The 196 = 14^2 identity, the guards (G12) identity, and OEIS A110979.
(E) Basin census of 3n-1 on odd n <= 10^6 (what "attains bounded values" means
    on the minus sheet) versus 3n+1.
(F) Session-lead square-sum graph probe: components, degree <= 1 vertices,
    Hamiltonian paths for n <= 32, the two pasted paths, the Lean loopless
    witness.
(G) Verdict table.

Every claim is labelled in the accompanying note; this script only prints the
exact evidence.  Explicit raise everywhere (no bare assert).
"""
import os
import re
import subprocess
import sys
import time
from math import gcd, isqrt

from sympy import isprime, primerange

SCRATCH = "/private/tmp/claude-501/-Users-e-Documents-GitHub-math/e197ec98-d8f9-4475-947b-5af87889cf35/scratchpad/tao"
ARXIV_ID = "1909.03562"
UA = "Mozilla/5.0 (research)"


def check(cond, msg):
    if not cond:
        raise RuntimeError("CHECK FAILED: " + msg)


def hr(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


# ---------------------------------------------------------------------------
# (A) Tao's paper: fetch, text-convert, count terms
# ---------------------------------------------------------------------------
def fetch_tao_text():
    os.makedirs(SCRATCH, exist_ok=True)
    abs_html = os.path.join(SCRATCH, "abs.html")
    pdf = os.path.join(SCRATCH, "tao.pdf")
    txt = os.path.join(SCRATCH, "tao.txt")
    if not os.path.exists(abs_html):
        subprocess.run(["curl", "-s", "-A", UA, "-o", abs_html,
                        "https://arxiv.org/abs/" + ARXIV_ID], check=False)
    if not os.path.exists(pdf):
        subprocess.run(["curl", "-s", "-L", "-A", UA, "-o", pdf,
                        "https://arxiv.org/pdf/" + ARXIV_ID], check=False)
    if not os.path.exists(txt) and os.path.exists(pdf):
        subprocess.run(["pdftotext", pdf, txt], check=False)
    title = None
    if os.path.exists(abs_html):
        with open(abs_html, "r", encoding="utf-8", errors="replace") as fh:
            m = re.search(r'citation_title" content="([^"]*)"', fh.read())
            if m:
                title = m.group(1)
    lines = None
    if os.path.exists(txt):
        with open(txt, "r", encoding="utf-8", errors="replace") as fh:
            lines = fh.read().split("\n")
    return title, lines


def section_A():
    hr("(A) Tao arXiv:1909.03562 -- term census of the paper text")
    title, lines = fetch_tao_text()
    print("user agent for curl:", UA)
    print("abstract page citation_title:", title)
    if lines is None:
        print("PAPER TEXT UNAVAILABLE (curl or pdftotext failed); section A is SCOPE-only")
        return None
    print("pdftotext line count:", len(lines))
    ver = [l for l in lines[:10] if l.startswith("arXiv:")]
    print("version stamp line:", ver[0] if ver else "not found")

    def count(pattern, flags=re.I):
        rx = re.compile(pattern, flags)
        hits = [(i + 1) for i, l in enumerate(lines) if rx.search(l)]
        return hits

    # Map-variant vocabulary.  pdftotext renders 3^{n-1} as "3n-1" (with a
    # unicode minus), so every raw hit is classified by what follows it.
    raw = re.compile(r"3\s?[nNx]\s?[-−–]\s?1")
    raw_hits = []
    for i, l in enumerate(lines):
        for m in raw.finditer(l):
            tail = l[m.end():m.end() + 6]
            kind = "power 3^(n-1) (followed by Z, 2^-a or a coset)" if re.match(
                r"\s?(Z|2|−|mod|\(|/|\)|,)", tail) else "OTHER"
            raw_hits.append((i + 1, kind))
    print("raw '3n-1'-shaped hits:", len(raw_hits))
    kinds = {}
    for _, k in raw_hits:
        kinds[k] = kinds.get(k, 0) + 1
    for k in sorted(kinds):
        print("   ", k, ":", kinds[k])
    other = [ln for ln, k in raw_hits if k == "OTHER"]
    print("hits not explained as 3^(n-1) (line numbers):", other)

    terms = [
        ("3x-1 / 3n-1 map as words", r"\b3x\s?-\s?1\b|three\s?n\s?minus|3n\s?-\s?1 problem"),
        ("negative", r"\bnegative\b"),
        ("5n+1 / 5x+1", r"5\s?[nNx]\s?\+\s?1"),
        ("generalization(s)", r"generali[sz]ation"),
        ("martingale", r"martingale"),
        ("entropy decrement", r"entropy decrement"),
        ("dyadic", r"\bdyadic\b"),
        ("renewal process", r"renewal process"),
        ("two-dimensional renewal", r"two-dimensional renewal"),
        ("Syracuse random variable", r"Syracuse random variable"),
        ("first passage", r"first passage"),
        ("stabilis/stabiliz", r"stabili[sz]"),
        ("3-adic", r"3-adic"),
        ("characteristic function", r"characteristic function"),
        ("Fourier", r"Fourier"),
        ("Plancherel", r"Plancherel"),
        ("Geom(2)", r"Geom\(2\)"),
        ("Pascal", r"Pascal"),
        ("triangle", r"triangle"),
        ("logarithmic density", r"logarithmic density"),
        ("natural density", r"natural density"),
        ("Bourgain", r"Bourgain"),
        ("Korec", r"Korec"),
        ("Krasikov", r"Krasikov"),
        ("Terras", r"Terras"),
        ("entropy (any)", r"entropy"),
        ("Renyi / collision entropy", r"Renyi|Rényi|collision entropy"),
        ("Baker", r"Baker"),
        ("Lagarias", r"Lagarias"),
    ]
    print()
    print("term                          lines  first-line  ")
    for name, pat in terms:
        h = count(pat)
        print("%-30s %5d  %s" % (name, len(h), (h[0] if h else "-")))
    # One short quotation (under fifteen words) from Section 1.3, giving the
    # scope of the stabilisation-to-theorem implication.
    q = "does not use any particular properties of the Syracuse map beyond (1.19), (1.20)"
    ql = [i + 1 for i, l in enumerate(lines) if "beyond (1.19), (1.20)" in l]
    print()
    print("quotation (Section 1.3, %d words):" % len(q.split()), repr(q))
    print("quotation located at pdftotext line(s):", ql)
    check(len(ql) >= 1, "quotation must be located in the text")
    # Structural anchors (line numbers only)
    anchors = {
        "Theorem 1.3": r"^Theorem 1\.3 ",
        "Theorem 1.6 (Syracuse form)": r"^Theorem 1\.6 ",
        "Conjecture 1.5": r"^Conjecture 1\.5 ",
        "Proposition 1.9 (valuation distribution)": r"^Proposition 1\.9 ",
        "Proposition 1.11 (stabilisation of first passage)": r"^Proposition 1\.11 ",
        "Proposition 1.14 (fine-scale mixing)": r"^Proposition 1\.14 ",
        "Proposition 1.17 (decay of characteristic function)": r"^Proposition 1\.17 ",
        "Remark 1.15 (Shannon-entropy heuristic)": r"^Remark 1\.15",
        "Section 3 (reduction to stabilisation)": r"^3\. Reduction to stabilisation",
        "Section 6 (reduction to Fourier decay)": r"^6\. Reduction to Fourier decay",
        "Remark 5.1 (Korec recovery from the (1.19) argument)": r"^Remark 5\.1",
    }
    alpha = [i + 1 for i, l in enumerate(lines) if l.strip().startswith("α := 1.001")]
    print("stabilisation parameter alpha := 1.001 set at line:", alpha[0] if alpha else "not found",
          "; affine iteration formula is the paper's eq. (1.7)")
    refs = [i + 1 for i, l in enumerate(lines) if l.strip() == "References"]
    print()
    print("bibliography ('References' heading) starts at line:", refs[-1] if refs else "not found",
          "(hits at or after it are citation titles, not statements of the paper)")
    for k, pat in anchors.items():
        h = count(pat, 0)
        print("%-52s line %s" % (k, h[0] if h else "not found"))
    return lines


# ---------------------------------------------------------------------------
# (B) Minus-sheet Syracuse offset symmetry
# ---------------------------------------------------------------------------
def v2(m):
    return (m & -m).bit_length() - 1


def syr(N, b):
    m = 3 * N + b
    a = v2(m)
    return m >> a, a


def offset_F(a):
    """Tao's n-Syracuse offset F_n(a) = sum_{i=1}^n 3^{n-i} 2^{-a[i,n]} as a Fraction."""
    from fractions import Fraction
    n = len(a)
    tot = Fraction(0)
    for i in range(n):
        tot += Fraction(3 ** (n - 1 - i), 2 ** sum(a[i:]))
    return tot


def section_B():
    hr("(B) Minus-sheet Syracuse offset: Syr_b^n(N) = 3^n 2^{-|a|} N + b F_n(a)")
    from fractions import Fraction
    NMAX, nmax = 1 << 14, 6
    checked = 0
    for b in (+1, -1):
        for N in range(1, NMAX, 2):
            x, word = N, []
            for n in range(1, nmax + 1):
                x, a = syr(x, b)
                word.append(a)
                rhs = Fraction(3 ** n * N, 2 ** sum(word)) + b * offset_F(word)
                check(rhs == x, "offset identity b=%d N=%d n=%d" % (b, N, n))
                checked += 1
    print("offset identity verified for both signs, odd N < %d, n <= %d: %d instances"
          % (NMAX, nmax, checked))

    # Valuation-word counts: for every word (a_1..a_n) with |a| <= J, the number
    # of odd N in [1, 2^(J+2)) realising it is 2^(J+1-|a|) on BOTH sheets.
    J = 12
    K = J + 2
    words = {+1: {}, -1: {}}
    for b in (+1, -1):
        for N in range(1, 1 << K, 2):
            x, s, w = N, 0, []
            while True:
                x, a = syr(x, b)
                s += a
                if s > J:
                    break
                w.append(a)
                key = tuple(w)
                words[b][key] = words[b].get(key, 0) + 1
    check(words[+1] == words[-1], "word counts differ between sheets")
    bad = [(w, c) for w, c in words[+1].items() if c != 2 ** (K - 1 - sum(w))]
    check(not bad, "word count not 2^(K-1-|a|): %s" % bad[:3])
    print("valuation words with |a| <= %d: %d words, counts identical on both sheets,"
          % (J, len(words[+1])))
    print("   every word (a_1..a_n) realised by exactly 2^(%d-|a|) odd N in [1, 2^%d)"
          % (K - 1, K))
    print("   => the word law is exactly Geom(2)^n truncated at |a| <= J on both sheets")
    print("Corollary: Syrac_-(Z/3^nZ) = -Syrac_+(Z/3^nZ) as random variables mod 3^n,")
    print("   and dTV(-X, -X + Unif(3^m Z/3^n Z)) = dTV(X, X + Unif(3^m Z/3^n Z)).")
    # Explicit check of the mod-3^n negation on the finite word set, n = 3.
    n = 3
    mod = 3 ** n
    dist = {+1: {}, -1: {}}
    for w, c in words[+1].items():
        if len(w) != n:
            continue
        # Tao's Syrac uses reversed labels; mod 3^n the residue of F_n(a) is
        # F_n(a) * 2^{|a|} * inverse(2^{|a|}) -- compute directly as an integer mod 3^n.
        num = sum(3 ** (n - 1 - i) * 2 ** (sum(w[:i])) for i in range(n))  # F_n * 2^{|a|}
        inv = pow(2, -sum(w), mod)
        r = (num * inv) % mod
        for b in (+1, -1):
            dist[b][(b * r) % mod] = dist[b].get((b * r) % mod, 0) + c
    neg = {(-r) % mod: c for r, c in dist[+1].items()}
    check(neg == dist[-1], "mod 3^3 negation of the offset distribution")
    print("n = %d, |a| <= %d: the minus-sheet offset residue histogram mod %d is the"
          " negation of the plus-sheet one (checked, %d residues hit, multiples of 3 hit: %d)"
          % (n, J, mod, len(dist[+1]), sum(1 for r in dist[+1] if r % 3 == 0)))


# ---------------------------------------------------------------------------
# (C) The Delta = 4 ladder, the Lean finite statement, prime targets, PPT, 2363
# ---------------------------------------------------------------------------
def section_C():
    hr("(C) 'Delta = 4 prime ladder' 3,7,11,17 and companions")
    ladder = [3, 7, 11, 17]
    diffs = [ladder[i + 1] - ladder[i] for i in range(3)]
    print("ladder", ladder, "consecutive differences", diffs,
          "-> arithmetic progression:", len(set(diffs)) == 1)
    # Longest AP of primes with common difference 4.
    best = []
    for p in primerange(2, 10 ** 5):
        run = [p]
        while isprime(run[-1] + 4):
            run.append(run[-1] + 4)
        if len(run) > len(best):
            best = run
    print("longest run of primes with difference 4 below 10^5:", best, "length", len(best))
    # mod 3 proof: among a, a+4, a+8 the residues mod 3 are a, a+1, a+2.
    print("mod-3 residues of (a, a+4, a+8):", [(a % 3, (a + 4) % 3, (a + 8) % 3) for a in range(3)])
    shifted = [p + 4 for p in ladder]
    print("pasted Lean sets: {3,7,11,17} + 4 =", shifted, "== {7,11,15,21}:", shifted == [7, 11, 15, 21])
    print("primality of the target set:", [(q, bool(isprime(q))) for q in shifted])
    print("factorisations of the composites:", {15: "3*5", 21: "3*7"})
    # Prime targets under the odd maps T_b(n) = (3n+b)/2^v.
    for b in (+1, -1):
        first = None
        for n in range(1, 1000, 2):
            t, _ = syr(n, b)
            if t != 1 and t != n and isprime(t):
                first = (n, t)
                break
        print("b=%+d: first odd n>0 whose Syracuse target is a prime other than 1 or n:" % b, first)
        pre3 = [n for n in range(1, 1 << 16, 2) if syr(n, b)[0] == 3]
        print("b=%+d: odd n < 2^16 with Syracuse target 3:" % b, pre3, "(3*2^k = 3n%+d is impossible mod 3)" % b)
    # PPT (77, 36, 85)
    check(77 ** 2 + 36 ** 2 == 85 ** 2, "PPT")
    m, n_ = 9, 2
    print("77^2 + 36^2 = %d = 85^2; Euclid parameters (m,n) = (%d,%d): m^2-n^2=%d, 2mn=%d, m^2+n^2=%d"
          % (77 ** 2 + 36 ** 2, m, n_, m * m - n_ * n_, 2 * m * n_, m * m + n_ * n_))
    print("   (m,n) = ((7+11)/2, (11-7)/2) = (%d,%d); generic: any coprime odd p<q give leg pq"
          % ((7 + 11) // 2, (11 - 7) // 2))
    # 2363
    d = 2 ** 11 - 3 ** 7
    print("2^11 - 3^7 =", d, "; (-17)*(%d) =" % d, -17 * d, "; 17*139 =", 17 * 139)
    check(-17 * d == 2363 and 2363 == 17 * 139, "2363 identity")
    # the seven-cycle of 3n-1 with valuation word (1,1,1,2,1,1,4)
    x, word, cyc = 17, [], [17]
    for _ in range(7):
        x, a = syr(x, -1)
        word.append(a)
        cyc.append(x)
    print("3n-1 orbit of 17 (odd Syracuse form):", cyc, "valuations", word,
          "|a| =", sum(word), "L =", len(word))
    check(cyc[-1] == 17, "seven-cycle closes")
    B = 2363
    Delta = d
    print("cycle gate n_0 = b*B/Delta with b=-1, B=%d, Delta=%d: n_0 = %d; q = |Delta|/gcd(B,|Delta|) = %d"
          % (B, Delta, (-1 * B) // Delta, abs(Delta) // gcd(B, abs(Delta))))
    # Check B is indeed the cycle carry: 2^11 * 17 = 3^7 * 17 + b*B  => B = (2^11-3^7)*17/b
    print("carry check: (2^11 - 3^7)*17 =", d * 17, "= b*B with b=-1 =", -B)
    check(d * 17 == -B, "carry identity")


# ---------------------------------------------------------------------------
# (D) 196 = 14^2, (G12), OEIS A110979
# ---------------------------------------------------------------------------
def section_D():
    hr("(D) 196 = 14^2 and the prime-sum-minus-one squares")
    odd_primes = [p for p in primerange(3, 40)]
    k = 11
    s = 1 + sum(odd_primes[:k])
    print("1 + sum of first %d odd primes %s = %d = %d^2" % (k, odd_primes[:k], s, isqrt(s)))
    check(s == 196 and isqrt(s) ** 2 == s, "196")
    c = [(p - (2 * i + 1)) // 2 for i, p in enumerate(odd_primes[:k], 1)]
    print("c_i = (p_i - (2i+1))/2 (odd composites from 3 through p_i):", c, "sum", sum(c))
    print("(G12): (k+1)^2 + 2*sum c_i = %d + %d = %d" % ((k + 1) ** 2, 2 * sum(c), (k + 1) ** 2 + 2 * sum(c)))
    check((k + 1) ** 2 + 2 * sum(c) == s, "G12")
    # A110979: squares equal to the sum of the first m primes minus 1
    LIM = 20_000_000
    tot, out = 0, []
    for m, p in enumerate(primerange(2, LIM), 1):
        tot += p
        r = isqrt(tot - 1)
        if r * r == tot - 1:
            out.append((m, p, tot - 1, r))
    print("primes below %d scanned; (m, p_m, sum-1, root):" % LIM)
    for row in out:
        print("   m=%d  p_m=%d  sum-1=%d  root=%d" % row)
    cited = [1, 4, 9, 16, 196, 839056, 7796654478001]
    print("OEIS A110979 (lookup 2026-09-22) data:", cited)
    print("recomputed values:", [r[2] for r in out], "match the cited prefix:",
          [r[2] for r in out] == cited[:len(out)])
    print("196 is term number %d of A110979; the next term %d needs m = 504 primes (p_504 = 3607)"
          % (cited.index(196) + 1, 839056))


# ---------------------------------------------------------------------------
# (E) Basin census on both sheets, odd n <= 10^6
# ---------------------------------------------------------------------------
def section_E():
    hr("(E) Basin census: what 'attains bounded values' means on each sheet")
    LIM = 10 ** 6
    for b in (+1, -1):
        basin = {}
        cycle_min = {}
        counts = {}
        t0 = time.time()
        for n in range(1, LIM + 1, 2):
            x = n
            steps = 0
            while True:
                if x < n and x in basin:
                    root = basin[x]
                    break
                if b == -1 and x in (1, 5, 17):
                    root = x
                    break
                if b == +1 and x == 1:
                    root = 1
                    break
                x, _ = syr(x, b)
                steps += 1
                if steps > 5000:
                    raise RuntimeError("orbit cap exceeded at n=%d b=%d" % (n, b))
            basin[n] = root
            counts[root] = counts.get(root, 0) + 1
        total = sum(counts.values())
        print("b=%+d: odd n <= %d: basins %s (fractions %s) [%.1fs]"
              % (b, LIM, dict(sorted(counts.items())),
                 {r: round(c / total, 6) for r, c in sorted(counts.items())}, time.time() - t0))
    print("reading: on the minus sheet every tested orbit attains the bounded value 1, 5 or 17;")
    print("   'almost all orbits attain almost bounded values' is EMPIRICALLY TRUE there while")
    print("   'reaches 1' fails on a set of density about 2/3 of the odd integers tested.")


# ---------------------------------------------------------------------------
# (F) Square-sum graph probe (session lead)
# ---------------------------------------------------------------------------
def sq_adj(n):
    adj = {v: set() for v in range(1, n + 1)}
    for x in range(1, n + 1):
        for y in range(x + 1, n + 1):
            s = x + y
            r = isqrt(s)
            if r * r == s:
                adj[x].add(y)
                adj[y].add(x)
    return adj


def components(adj):
    seen, comps = set(), []
    for v in adj:
        if v in seen:
            continue
        stack, comp = [v], set()
        while stack:
            u = stack.pop()
            if u in comp:
                continue
            comp.add(u)
            stack.extend(adj[u] - comp)
        seen |= comp
        comps.append(sorted(comp))
    return comps


def ham_path(adj, budget):
    """Backtracking Hamiltonian path search with the deg'-pruning rule.
    Returns (path or None, timed_out)."""
    n = len(adj)
    verts = sorted(adj, key=lambda v: len(adj[v]))
    t0 = time.time()
    sys.setrecursionlimit(10000)

    def rec(cur, visited, path):
        if len(path) == n:
            return list(path)
        if time.time() - t0 > budget:
            raise TimeoutError
        ones = 0
        for u in adj:
            if u in visited:
                continue
            d = sum(1 for w in adj[u] if w not in visited) + (1 if u in adj[cur] else 0)
            if d == 0:
                return None
            if d == 1:
                ones += 1
                if ones > 1:
                    return None
        for u in sorted(adj[cur] - visited, key=lambda w: sum(1 for z in adj[w] if z not in visited)):
            visited.add(u)
            path.append(u)
            r = rec(u, visited, path)
            if r is not None:
                return r
            path.pop()
            visited.discard(u)
        return None

    try:
        for s in verts:
            r = rec(s, {s}, [s])
            if r is not None:
                return r, False
        return None, False
    except TimeoutError:
        return None, True


def section_F():
    hr("(F) Square-sum graph Q_n on {1..n}: session-lead probe re-verified")
    comp_counts = {}
    for n in range(1, 33):
        adj = sq_adj(n)
        comps = components(adj)
        comp_counts[n] = len(comps)
    print("component counts n=1..32:", comp_counts)
    check(all(comp_counts[n] == 3 for n in range(4, 13)), "3 components for 4<=n<=12")
    check(comp_counts[13] == 2 and all(comp_counts[n] == 1 for n in range(14, 33)), "2 at 13, 1 for 14..32")
    print("components at n=13:", components(sq_adj(13)))
    print("components at n=12:", components(sq_adj(12)))
    low = {}
    for n in range(14, 33):
        adj = sq_adj(n)
        low[n] = sorted(v for v in adj if len(adj[v]) <= 1)
    print("degree <= 1 vertices, n=14..32:", low)
    check(low[18] == [16, 17, 18] and low[19] == [16, 18] and all(low[n] == [18] for n in range(20, 31))
          and low[31] == [] and low[32] == [], "low-degree pattern")
    print("neighbours of 18 in Q_30:", sorted(sq_adj(30)[18]), "; in Q_31:", sorted(sq_adj(31)[18]))
    # Pasted paths
    p15 = [8, 1, 15, 10, 6, 3, 13, 12, 4, 5, 11, 14, 2, 7, 9]
    p23 = [18, 7, 9, 16, 20, 5, 11, 14, 2, 23, 13, 12, 4, 21, 15, 10, 6, 19, 17, 8, 1, 3, 22]
    for n, p in ((15, p15), (23, p23)):
        adj = sq_adj(n)
        ok = sorted(p) == list(range(1, n + 1)) and all(p[i + 1] in adj[p[i]] for i in range(n - 1))
        sums = [p[i] + p[i + 1] for i in range(n - 1)]
        print("pasted path n=%d valid: %s; consecutive sums %s" % (n, ok, sorted(set(sums))))
        check(ok, "pasted path n=%d" % n)
    # Hamiltonian existence n <= 32
    t0 = time.time()
    result = {}
    total_budget = 150.0
    for n in range(1, 33):
        remaining = total_budget - (time.time() - t0)
        if remaining <= 0:
            result[n] = "SKIPPED"
            continue
        adj = sq_adj(n)
        path, to = ham_path(adj, min(remaining, 40.0))
        result[n] = "TIMEOUT" if to else ("yes" if path else "no")
    print("Hamiltonian path existence n=1..32 [%.1fs]:" % (time.time() - t0))
    print("   ", {n: result[n] for n in result})
    yes = [n for n in result if result[n] == "yes"]
    no = [n for n in result if result[n] == "no"]
    print("   yes:", yes)
    print("   no :", no)
    print("   unresolved:", [n for n in result if result[n] not in ("yes", "no")])
    print("OEIS A090461 (lookup 2026-09-22, comments): every k >= 25 is in the sequence; cycles for k >= 32")
    # Lean loopless witness
    loops = [v for v in range(1, 33) if isqrt(2 * v) ** 2 == 2 * v]
    print("values v <= 32 with 2v a perfect square (Adj v v holds in the pasted Lean def):", loops)
    print("   minimal witness: Fin index 1 (value 2): 2 + 2 = 4 = 2^2, so 'loopless' is unprovable as stated")


# ---------------------------------------------------------------------------
# (G) Verdict table
# ---------------------------------------------------------------------------
def section_G():
    hr("(G) Verdict table")
    rows = [
        ("Tao's paper remarks on 3n-1 / other multipliers", "SCOPE", "0 map-variant hits in v7 text; every '3n-1' string is the power 3^(n-1)"),
        ("Tao's method is a 'dyadic martingale renewal process'", "REFUTED", "0 hits for 'martingale' and 'dyadic'; the renewal process is 2-dimensional in Z^2 (Section 7), used for Fourier decay"),
        ("Tao's method uses an 'entropy decrement'", "REFUTED", "0 hits; 'entropy' occurs only as a Shannon-entropy heuristic (Remark 1.15) and Renyi/collision entropy (Remark 6.1)"),
        ("Prefix-descent counts mod 2^J coincide on both sheets", "CITED", "minus_sheet_positive_control S2; re-verified here for valuation words |a| <= 12"),
        ("That coincidence = Tao's Proposition 1.9 input, not the stabilisation (1.20)", "CITED", "Remark 5.1 of the paper (line 1399): the (1.19) argument from Prop 1.9 + drift gives only Syr_min(N) <= N^theta, theta > 1/alpha, i.e. Korec-type"),
        ("Syrac_-(Z/3^nZ) = -Syrac_+(Z/3^nZ); Props 1.14 and 1.17 hold on the minus sheet", "PROVED", "offset identity + word-count equality (B); Prop 1.17 at -xi; dTV negation-invariant"),
        ("Full Theorem 1.3 transfers to 3n-1 (almost all orbits attain almost bounded values)", "OPEN", "Sections 3 and 5 not re-run with the sign flipped; no source states it; empirically true to 10^6"),
        ("Almost-all bounded values implies a single root", "REFUTED", "minus sheet: three cycles, basins about 0.33/0.32/0.35 on odd n <= 10^6"),
        ("Delta = 4 prime ladder 3,7,11,17 is an AP", "REFUTED", "differences 4,4,6; primes in AP with difference 4 have length <= 3 (mod 3)"),
        ("Lean delta_four_prime_step_invariance", "FINITE-EXACT", "true set statement; 15 = 3*5 and 21 = 3*7 are composite, so no prime reading"),
        ("3 is a 'first nontrivial prime target' of 3n+1 or 3n-1", "REFUTED", "3 has no Syracuse preimage on either sheet; first prime targets are 5 (plus, from 3) and 7 (minus, from 5)"),
        ("7 -> 11 'maps to' the PPT (77,36,85)", "SCOPE", "generic: leg pq for any coprime odd p<q via (m,n) = ((p+q)/2,(q-p)/2) = (9,2); no Collatz content"),
        ("(2^11-3^7)(-17) = 2363 = 17*139 and the gate n_0 = 17", "CITED", "catalan_elliptic (C6), counterexample_portrait table; re-verified"),
        ("'The difference 4 regulates the carry B_L'", "SCOPE", "no object: the carry of the seven-cycle is 2363 and involves no difference-4 structure"),
        ("1+3+5+...+37 = 196 = 14^2", "CITED", "collatz_guards_20260921_squarefree sec.7 (G12); re-verified; term 5 of OEIS A110979"),
        ("'Linear prime addition matches square area' has a trajectory consequence", "SCOPE", "no map found from A110979 to any orbit statement"),
        ("Q_n components 3 (4..12), 2 (13), 1 (14..32); deg<=1 vertex pattern", "FINITE-EXACT", "session lead probe re-verified"),
        ("Hamiltonian paths of Q_n for n <= 32", "FINITE-EXACT", "see (F); A090461 CITED for k >= 25"),
        ("Pasted Lean square_sum_graph is loopless", "REFUTED", "value 2: 2+2 = 4; also 8, 18, 32"),
    ]
    print("%-88s %-13s %s" % ("claim", "label", "evidence"))
    for c, l, e in rows:
        print("%-88s %-13s %s" % (c, l, e))


def main():
    t0 = time.time()
    section_A()
    section_B()
    section_C()
    section_D()
    section_E()
    section_F()
    section_G()
    print()
    print("total wall time %.1fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
