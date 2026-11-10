#!/usr/bin/env python3
"""AUDIT 2 (auditor's own code).  Rank law of a uniformly random alternating m x m matrix over F_2 conditioned on an odd
number of entries above the diagonal ("odd-edge model", README section 4.1b).

Closed forms used (and checked here against the brute force of x5_altrank.c for m <= 8):
  N(m, 2s)  = number of alternating matrices of rank 2s
            = prod_{i=1..s} 2^(2i-2)/(2^(2i)-1) * prod_{i=0..2s-1} (2^(m-i)-1)                       (classical)
  S(m, 2s)  = sum over alternating matrices A of rank 2s of (-1)^(number of edges of A)
            = (-1)^s 4^C(s,2) [n choose s]_4 ,  n = floor(m/2)
    (the edge-parity character is the pairing with J - I, an alternating matrix of the maximal rank 2n, so S is the
     value at the last class of the eigenmatrix of the association scheme of alternating forms, Delsarte-Goethals 1975)
  => P(rank = 2s | odd) - P(rank = 2s) = -S(m, 2s) / 2^C(m,2),   |.| <= 1.46 * 2^(-n^2)   (decreasing in m).

Usage: python3 x5_oddedge_check.py <file with the lines printed by x5_altrank> [x1 labelled raw file ...]"""
import sys, re
from fractions import Fraction

def gauss(n, k, q):
    if k < 0 or k > n: return 0
    num = den = 1
    for i in range(k): num *= q ** (n - i) - 1; den *= q ** (i + 1) - 1
    assert num % den == 0
    return num // den

def N_alt(m, s):
    x = Fraction(1)
    for i in range(2 * s): x *= 2 ** (m - i) - 1
    for i in range(1, s + 1): x *= Fraction(2 ** (2 * i - 2), 2 ** (2 * i) - 1)
    assert x.denominator == 1
    return x.numerator

def S_alt(m, s):
    n = m // 2
    return (-1) ** s * 4 ** (s * (s - 1) // 2) * gauss(n, s, 4)

def laws(m):
    tot = 2 ** (m * (m - 1) // 2); odd = {}; unc = {}
    for s in range(m // 2 + 1):
        N, S = N_alt(m, s), S_alt(m, s)
        assert (N - S) % 2 == 0
        odd[2 * s] = Fraction((N - S) // 2, tot // 2); unc[2 * s] = Fraction(N, tot)
    assert sum(odd.values()) == 1 and sum(unc.values()) == 1
    return odd, unc

# 1. brute force (lines with the same m are parts and are added)
ok = True; allcnt = {}
for line in open(sys.argv[1]):
    g = re.match(r'm (\d+)(.*)', line)
    if not g: continue
    m = int(g.group(1)); cnt = allcnt.setdefault(m, {})
    for h in re.finditer(r' r(\d+)_e(\d)_f(\d) (\d+)', g.group(2)):
        r, e, f, c = map(int, h.groups()); cnt[(r, e)] = cnt.get((r, e), 0) + c
for m, cnt in sorted(allcnt.items()):
    if sum(cnt.values()) != 2 ** (m * (m - 1) // 2): ok = False
    for s in range(m // 2 + 1):
        N, S = N_alt(m, s), S_alt(m, s)
        if cnt.get((2 * s, 0), 0) != (N + S) // 2 or cnt.get((2 * s, 1), 0) != (N - S) // 2: ok = False
    if any(r % 2 for (r, e) in cnt): ok = False
print("closed forms N(m,2s), S(m,2s) against brute force over all graphs, m =", sorted(allcnt), ":", "AGREE" if ok else "MISMATCH")

# 2. table
print("\n m | odd-edge law (rank: prob) | unconditioned | max abs difference over rank classes")
for m in range(2, 13):
    odd, unc = laws(m)
    d = max(abs(odd[r] - unc[r]) for r in odd)
    print("%2d | %s | %s | %.3e (= %s)" % (m, "  ".join("%d: %.7f" % (r, float(odd[r])) for r in sorted(odd, reverse=True) if odd[r]),
                                    "  ".join("%d: %.7f" % (r, float(unc[r])) for r in sorted(unc, reverse=True)), float(d), d))
# 3. the claim "less than 3e-5 in every rank class for m >= 8"
worst = max((max(abs(laws(m)[0][r] - laws(m)[1][r]) for r in laws(m)[0]), m) for m in range(8, 61))
print("\nmax over 8 <= m <= 60 of the largest difference: %.4e at m = %d ; < 3e-5 : %s" % (float(worst[0]), worst[1], worst[0] < Fraction(3, 10 ** 5)))
print("largest difference at m = 7: %.4e ; at m = 6: %.4e" % tuple(float(max(abs(laws(m)[0][r] - laws(m)[1][r]) for r in laws(m)[0])) for m in (7, 6)))
bound_ok = all(max(abs(laws(m)[0][r] - laws(m)[1][r]) for r in laws(m)[0]) <= Fraction(146, 100) / 2 ** ((m // 2) ** 2) for m in range(2, 61))
print("bound |difference| <= 1.46 * 2^(-floor(m/2)^2) for 2 <= m <= 60:", bound_ok)
print("exact odd-edge laws at m = 4, 5, 6:", {m: {r: str(p) for r, p in laws(m)[0].items() if p} for m in (4, 5, 6)})

# 3b. the packet's printed model values (out/oddedge_model.out), if given as --packet=<file>: compare to 7 decimals
pk = [a.split('=', 1)[1] for a in sys.argv[2:] if a.startswith('--packet=')]
sys.argv = sys.argv[:2] + [a for a in sys.argv[2:] if not a.startswith('--packet=')]
for path in pk:
    m = None; bad = 0; n = 0
    for line in open(path):
        g = re.match(r'm = +(\d+) odd-edge:(.*)', line)
        if g: m = int(g.group(1)); kind = 'odd'; rest = g.group(2)
        else:
            g = re.match(r' +uncond\.:(.*)', line)
            if not g or m is None: continue
            kind = 'unc'; rest = g.group(1)
        odd, unc = laws(m); law = odd if kind == 'odd' else unc
        for r, v in re.findall(r'rank (\d+): ([0-9.]+)', rest):
            n += 1
            if abs(float(law[int(r)]) - float(v)) > 6e-8: bad += 1
    print("\npacket file %s: %d printed probabilities (odd-edge and unconditioned, m = 2..11) compared with the closed forms to 7 decimals: %d disagreements"
          % (path.split('/')[-1], n, bad))

# 4. labelled census rank laws (from x1_cube output: labelled mode, or classes with weights m!/|Aut|; parts are added)
acc = {}
for path in sys.argv[2:]:
    for line in open(path):
        g = re.match(r'order (\d+) mode (\S+) classes (\d+) W (\d+)', line)
        if not g: continue
        t = acc.setdefault((int(g.group(1)), g.group(2)), [0, {}]); t[0] += int(g.group(4))
        for h in re.finditer(r' rank (\d+) (\d+) (\d+) (\d+) (\d+)', line):
            t[1][int(h.group(1))] = t[1].get(int(h.group(1)), 0) + int(h.group(4)) + int(h.group(5))
print()
for (m, mode), (W, d) in sorted(acc.items()):
    assert W == 2 ** (m * (m - 1) // 2), "incomplete weight"
    cen = {r: Fraction(c, W) for r, c in d.items()}
    odd, unc = laws(m)
    same = all(cen.get(r, 0) == odd[r] for r in odd)
    tv = sum(abs(cen.get(r, 0) - odd[r]) for r in odd) / 2
    print("order %d (%s): labelled rank law of beta_R = %s ; odd-edge law = %s ; EQUAL: %s ; total variation %.6f"
          % (m, mode, {r: str(p) for r, p in sorted(cen.items())}, {r: str(p) for r, p in sorted(odd.items()) if p}, same, float(tv)))
