"""aud_c_certs.py - audit task 2(c): rebuild the certificates of certs_plus.txt / certs_minus.txt and recount.

For every line "sign c v2 order m code Omask levels":
   start pair = tournament of order m with the upper-triangle code, out-set = bits of Omask;
   each level (delta, i, mm): a = i + delta + mm; append TT_a after the current base (base => tail), then a
   closer (tail => closer, closer beats exactly the current out-set, the rest of the base beats the closer);
   new out-set = tail vertices i and i + delta (1-based).
Checks, with the auditor's own programme aud_pairstate (subset recursion, mode 0) on ALL lines:
   order, r = sign, M = c 2^v2.
On a random sample of 150 lines per file plus the lines named in the audit brief:
   - a second, literal count (mode 1: depth-first generation of every path, no memo);
   - the closed tournaments T_1, T_2 are built and their Hamiltonian paths counted (mode 2): H = sign + 2^t c 2^v2.
Also: file statistics (line counts, '-' lines, one- and two-level lines, max order, max v2, v2 = 0 count).
"""
import random, subprocess, sys, os, time
from aud_core import adj_from_code, transitive, close, is_tournament

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(HERE, "aud_pairstate.exe")
NAMED = {"certs_plus.txt": [245, 853, 2023, 3001, 3427, 3773, 3985], "certs_minus.txt": [1, 7, 13, 65, 1039]}
SIGN = {"certs_plus.txt": 1, "certs_minus.txt": -1}


def parse(path):
    rows = []
    for l in open(path):
        if l.startswith("#") or not l.strip():
            continue
        p = l.split()
        sign, c, v, order, m, code, Om = map(int, p[:7])
        levels = [] if p[7] == "-" else [tuple(map(int, z.split(","))) for z in p[7].split(";")]
        rows.append((sign, c, v, order, m, code, Om, levels))
    return rows


def rebuild(m, code, Om, levels):
    A = adj_from_code(m, code)
    O = [i for i in range(m) if Om >> i & 1]
    for (delta, i, mm) in levels:
        a = i + delta + mm
        n = len(A)
        A = close(A, O, transitive(a))
        O = [n + i - 1, n + i - 1 + delta]
    return A, O


def rec(mode, A, O):
    n = len(A)
    return " ".join([str(mode), str(n)] + [str(sum(1 << j for j in range(n) if A[i][j])) for i in range(n)]
                    + [str(sum(1 << o for o in O))])


def run(lines):
    out = subprocess.run([EXE], input="\n".join(lines) + "\n", capture_output=True, text=True)
    if out.returncode != 0:
        print("PROGRAMME ERROR", out.stderr)
        sys.exit(1)
    return [l for l in out.stdout.split("\n") if l.strip()]


allok = True
rng = random.Random(20261111)
for fname in ("certs_plus.txt", "certs_minus.txt"):
    t0 = time.time()
    rows = parse(os.path.join(HERE, "..", fname))
    sign0 = SIGN[fname]
    # --- file statistics
    cs = [r[1] for r in rows]
    stat = dict(lines=len(rows), c_ok=(cs == list(range(1, 4098, 2))), signs_ok=all(r[0] == sign0 for r in rows),
                no_level=sum(1 for r in rows if not r[7]), one_level=sum(1 for r in rows if len(r[7]) == 1),
                two_level=sum(1 for r in rows if len(r[7]) == 2), more=sum(1 for r in rows if len(r[7]) > 2),
                max_order=max(r[3] for r in rows), max_v2=max(r[2] for r in rows),
                v2_zero=sum(1 for r in rows if r[2] == 0),
                max_m=max(r[4] for r in rows), min_order=min(r[3] for r in rows),
                levels_delta_ok=all(d in (1, 2) and i >= 1 and mm >= 0 for r in rows for (d, i, mm) in r[7]),
                start_order_le9=all(r[4] <= 9 for r in rows),
                nolevel_order_le9=all(r[3] <= 9 for r in rows if not r[7]))
    print(fname, stat)
    hist = {}
    for r in rows:
        hist[r[3]] = hist.get(r[3], 0) + 1
    print("   order histogram:", sorted(hist.items()))
    hist = {}
    for r in rows:
        if r[7]:
            hist[r[3]] = hist.get(r[3], 0) + 1
    print("   order histogram of tower lines:", sorted(hist.items()))
    # --- all lines, subset recursion
    built = [rebuild(r[4], r[5], r[6], r[7]) for r in rows]
    assert all(is_tournament(A) for A, O in built)
    outs = run([rec(0, A, O) for A, O in built])
    bad = 0
    states = []
    for r, o in zip(rows, outs):
        n, rr, M, h, E = map(int, o.split())
        states.append((n, rr, M, h, E))
        if not (n == r[3] and rr == r[0] and M == r[1] * 2 ** r[2]):
            bad += 1
            print("   MISMATCH", r, "->", o)
    ok = (bad == 0 and len(outs) == len(rows))
    allok &= ok
    print(("PASS  " if ok else "FAIL  ") + "%s: all %d lines recounted by the auditor's subset recursion, %d mismatches (%.0fs)"
          % (fname, len(rows), bad, time.time() - t0))
    # --- sample: literal enumeration + closed tournaments
    t0 = time.time()
    idx = set(rng.sample(range(len(rows)), 150))
    named = [k for k, r in enumerate(rows) if r[1] in NAMED[fname]]
    assert len(named) == len(NAMED[fname])
    idx |= set(named)
    idx = sorted(idx)
    outs = run([rec(1, *built[k]) for k in idx])
    bad = 0
    for k, o in zip(idx, outs):
        a, b = o.split("|")
        a = tuple(map(int, a.split()))
        b = tuple(map(int, b.split()))
        r = rows[k]
        if not (a[1:] == b and a[0] == r[3] and b[0] == r[0] and b[1] == r[1] * 2 ** r[2]):
            bad += 1
            print("   LITERAL MISMATCH", r, "->", o)
    ok = (bad == 0 and len(outs) == len(idx))
    allok &= ok
    print(("PASS  " if ok else "FAIL  ") + "%s: %d sampled lines (150 random + %d named) recounted by literal path generation, %d mismatches (%.0fs)"
          % (fname, len(idx), len(named), bad, time.time() - t0))
    t0 = time.time()
    recs = []
    for k in idx:
        A, O = built[k]
        for t in (1, 2):
            recs.append(rec(2, close(A, O, transitive(t)), []))
    outs = run(recs)
    bad = 0
    j = 0
    for k in idx:
        r = rows[k]
        for t in (1, 2):
            n, H = map(int, outs[j].split())
            j += 1
            if not (n == r[3] + t + 1 and H == r[0] + 2 ** t * r[1] * 2 ** r[2]):
                bad += 1
                print("   CLOSED MISMATCH", r, t, outs[j - 1])
    ok = (bad == 0)
    allok &= ok
    print(("PASS  " if ok else "FAIL  ") + "%s: closed tournaments T_1, T_2 of the %d sampled lines have sign + 2^t c 2^v2 Hamiltonian paths, %d mismatches (%.0fs)"
          % (fname, len(idx), bad, time.time() - t0))
    print("   named lines:")
    for k in named:
        r = rows[k]
        print("     c=%d v2=%d order=%d start(m=%d, code=%d, Omask=%d) levels=%s -> state (n, r, M, h, E) = %s"
              % (r[1], r[2], r[3], r[4], r[5], r[6], r[7] or "-", states[k]))
print("ALL PASS" if allok else "SOME FAILED")
