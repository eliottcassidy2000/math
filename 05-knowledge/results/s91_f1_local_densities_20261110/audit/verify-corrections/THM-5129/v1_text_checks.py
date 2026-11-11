#!/usr/bin/env python3
"""Verifier of the corrections to THM-5129 (2026-11-11): text-level checks.

(1) frontmatter parses as YAML; line widths of the body and of the status block;
(2) the session's slope_lists_le73.txt equals the audit's s3_slope_lists.txt;
(3) |Delta(h)| against the number of integers of [h, 3h/2] (statement in (D));
(4) window counts of (D) recomputed from the lists; shares quoted in (D);
(5) r(N) on [3, 225) from the lists (statement (C)).
Reads files only; writes nothing.
"""
import io
import re
import sys
from fractions import Fraction

ROOT = "C:/tmp/wt-s91/"
THM = ROOT + ("01-canon/theorems/THM-5129-one-vertex-sources-finite-supply-complete-lists-permanent-"
              "exceptions-line-counts.md")
LISTS_A = ROOT + ("04-computation/experiments/s91_f1_local_densities_20261110/rooted_spectrum/"
                  "slope_lists_le73.txt")
LISTS_B = ROOT + ("05-knowledge/results/s91_f1_local_densities_20261110/audit/audit-5129/"
                  "s3_slope_lists.txt")


def read(path):
    with io.open(path, "r", encoding="utf-8") as f:
        return f.read()


def parse_lists(text, sep):
    out = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        head, tail = line.split(sep, 1)
        D = int(head.strip())
        d = {}
        for tok in tail.split():
            a, b = tok.split(":")
            d[int(a)] = int(b)
        out[D] = d
    return out


def main():
    ok = True
    text = read(THM)
    lines = text.split("\n")
    # (1) frontmatter
    assert lines[0] == "---"
    end = lines.index("---", 1)
    fm = "\n".join(lines[1:end])
    try:
        import yaml
        data = yaml.safe_load(fm)
        print("(1) YAML frontmatter parses; keys:", list(data.keys()))
        print("    title length", len(data["title"]), "; status length", len(data["status"]))
        print("    scripts entries:", len(data["scripts"]), "; related entries:", len(data["related"]))
        for k in ("scripts", "related"):
            for item in data[k]:
                if not isinstance(item, str):
                    print("    NOTE: %s item parsed as %s: %r" % (k, type(item).__name__, item))
    except ImportError:
        print("(1) PyYAML not available; frontmatter not parsed")
    widths = [(i + 1, len(l)) for i, l in enumerate(lines)]
    long_body = [(i, w) for (i, w) in widths if i > end + 1 and w > 118]
    status_start = next(i for i, l in enumerate(lines) if l.startswith("status:"))
    status_end = next(i for i, l in enumerate(lines) if l.startswith("session:"))
    long_status = [(i + 1, len(lines[i])) for i in range(status_start, status_end) if len(lines[i]) > 118]
    print("    body lines longer than 118:", long_body)
    print("    status lines longer than 118:", long_status)
    nonascii = sorted(set(ch for ch in text if ord(ch) > 127))
    print("    non-ASCII characters in the file:", [hex(ord(c)) for c in nonascii])
    for phrase in ("cannot work", "the reason", "level of distribution", "major-arc", "local side is closed",
                   "PRE-AUDIT", "Not independently audited", "NOT covered", "nine lattice points",
                   "singular series", "relative variance", "have the least possible exponent",
                   "least exponent at every order"):
        hits = [i + 1 for i, l in enumerate(lines) if phrase.lower() in l.lower()]
        print("    phrase %-34r lines %s" % (phrase, hits))

    # (2) lists
    A = parse_lists(read(LISTS_A), ":") if False else None
    A = {}
    for line in read(LISTS_A).splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        head, tail = line.split(": ", 1)
        A[int(head)] = {int(t.split(":")[0]): int(t.split(":")[1]) for t in tail.split()}
    B = parse_lists(read(LISTS_B), "|")
    same = (A == B)
    npairs = sum(len(v) for v in A.values())
    print("(2) offsets: session %d, audit %d; identical: %s; pairs (offset, odd part): %d"
          % (len(A), len(B), same, npairs))
    ok &= same and len(A) == 35 and npairs == 1066
    print("    attained offsets <= 73:", sorted(A))
    missing = [h for h in range(1, 74, 2) if h not in A]
    print("    odd numbers <= 73 without a list:", missing)
    ok &= (missing == [7, 21])

    # (3) |Delta(h)| against window sizes
    smaller = []
    for h in sorted(A):
        nint = (3 * h) // 2 - h + 1
        if len(A[h]) < nint:
            smaller.append(h)
    print("(3) attained h <= 73 with |Delta(h)| < number of integers of [h, 3h/2]:", smaller)
    ok &= (smaller == [3, 13, 35, 39, 43, 47, 63, 67])
    print("    |Delta| at 33, 65:", len(A[33]), len(A[65]), "; integers:", (99 // 2) - 33 + 1, (195 // 2) - 65 + 1)

    # (4) window counts
    def slopes_in(h, lo, hi):
        c = 0
        for d, j in A[h].items():
            v = d << j
            while v <= hi:
                if v >= lo:
                    c += 1
                v <<= 1
        return c
    quoted = ("1:1/1 3:1/2 5:3/3 9:3/5 11:3/6 13:4/7 15:4/8 17:4/9 19:5/10 23:8/12 25:7/13 27:5/14 29:9/15 "
              "31:9/16 33:13/17 35:5/18 37:12/19 39:7/20 41:14/21 43:9/22 45:12/23 47:5/24 49:15/25 51:14/26 "
              "53:15/27 55:10/28 57:11/29 59:19/30 61:12/31 63:5/32 65:24/33 67:9/34 69:26/35 71:23/36 73:27/37")
    mine = " ".join("%d:%d/%d" % (h, slopes_in(h, h, (3 * h) // 2), (3 * h) // 2 - h + 1) for h in sorted(A))
    print("(4) window table of (D) reproduced from the lists:", mine == quoted)
    ok &= (mine == quoted)
    # the table in the file itself
    body = " ".join(text.split())
    print("    the same table is in the theorem text:", quoted in body)
    ok &= quoted in body
    share = {h: Fraction(slopes_in(h, h, (3 * h) // 2), (3 * h) // 2 - h + 1) for h in A}
    full = [h for h in sorted(A) if share[h] == 1]
    print("    full windows:", full)
    r = [h for h in sorted(A) if 37 <= h <= 73]
    vals = [float(share[h]) for h in r]
    print("    37..73: mean %.4f, min %.3f at %d, max %.3f at %d"
          % (sum(vals) / len(vals), min(vals), r[vals.index(min(vals))], max(vals), r[vals.index(max(vals))]))
    r2 = sorted((float(share[h]), h) for h in A if 23 <= h <= 73)
    print("    poorest of 23..73:", [(h, round(s, 3)) for s, h in r2[:7]], "; attained offsets there:", len(r2))
    print("    shares at 9, 17, 19:", [round(float(share[h]), 2) for h in (9, 17, 19)])
    print("    shares at 35, 39, 47, 63:", [round(float(share[h]), 3) for h in (35, 39, 47, 63)])

    # (5) r(N) on [3, 225)
    def is_slope(D, U):
        if D not in A or U < 1:
            return False
        j = 0
        while U % 2 == 0:
            U //= 2
            j += 1
        return U in A[D] and j >= A[D][U]
    zeros = []
    rr = {}
    for N in range(3, 225, 2):
        c = sum(1 for D in range(1, N // 3 + 1, 2) if is_slope(D, (N - D) // 2))
        rr[N] = c
        if c == 0:
            zeros.append(N)
    print("(5) odd N in [3, 225) with r(N) = 0:", zeros)
    ok &= (zeros == [7, 11, 13, 21, 47, 67])
    tail = [rr[N] for N in range(129, 225, 2)]
    print("    r on [129, 225): min %d, mean %.2f" % (min(tail), sum(tail) / len(tail)))
    # U >= D in the lists (least slope of every odd part)
    viol = [(D, d, j) for D in A for d, j in A[D].items() if (d << j) < D]
    print("    least slopes below the offset (would contradict U >= D):", viol)
    ok &= not viol
    # 9-grid exponents and the bound max(0, s - 3)
    print("    least exponent of d = 9 at offsets 5, 9, 17, 33, 65:", [A[D].get(9) for D in (5, 9, 17, 33, 65)],
          "; bound:", [max(0, s - 3) for s in (2, 3, 4, 5, 6)])
    print("    dyadic slopes at 3, 5, 9, 17, 33, 65:", [(D, A[D].get(1)) for D in (3, 5, 9, 17, 33, 65)])
    print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
