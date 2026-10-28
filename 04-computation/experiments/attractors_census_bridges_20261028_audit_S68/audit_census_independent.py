"""Independent re-implementation of the HYP-9437 keyword census (S68 audit, 2026-10-28).

Only the regex lexicon is taken from the session's script (read with ast from the committed file, so the registered
patterns are used verbatim); loading, the curveball null, the statistics and the verdicts are re-implemented here.

Rows are stored as bitmasks and every set operation is done on sorted integer lists, so a run is a deterministic
function of (corpus, lexicon, BODY, seed): it does not depend on PYTHONHASHSEED (the session's script shuffles
lists built from Python sets of strings, whose order changes with the per-process hash seed).

Usage (from the repository root):
  python audit_census_independent.py run    [--idmin 4400 --idmax 4839 --null 2000 --seed 20261028 --body 6000
                                             --lexicon registered|strict --exclude 4839 --tag NAME]
  python audit_census_independent.py terms  [--idmin ... --body ...]      per-alternative match counts
  python audit_census_independent.py spot   [--k 10 --seed 7]             random spot-check with match contexts
  python audit_census_independent.py robust                                the robustness grid reported in AUDIT.md
"""
import argparse
import ast
import glob
import math
import os
import random
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SESSION_SCRIPT = os.path.join(HERE, "..", "attractor_census_canon_20261028.py")


def registered_lexicon():
    src = open(SESSION_SCRIPT, encoding="utf-8").read()
    tree = ast.parse(src)
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", None) == "FAMILIES" for t in node.targets):
            return ast.literal_eval(node.value)
    raise RuntimeError("FAMILIES not found")


REG = registered_lexicon()
NAMES = list(REG)            # same family order as the session (pairs are reported alphabetically anyway)


def strict_lexicon(drop_carrier=False):
    """the registered lexicon with the false-positive sources found by the audit removed (exploratory only):
    'Hecke' matched 'checked'/'checker' (114 of the 320 window files), 'character' matched 'characteristic'/
    'characterisation', 'parity' matched file names and 'disparity', 'Markov' matched Markov chains/inequality,
    'Steiner' matched R. P. Steiner's 1977 Collatz cycle paper, 'Baker' matched the Heegner-Baker-Stark class-number
    citations, 'Minkowski' matched Minkowski sums/bounds."""
    L = dict(REG)
    rep = {
        "MODULAR": [("|Hecke|", r"|\bHecke\b|"),
                    ("|Markov|", r"|Markov (number|triple|equation|tree|spectrum|level)|Markoff|"),
                    ("|Minkowski|", r"|Minkowski(?! sum| bound)|")],
        "BINARY_ALGEBRA": [("|parity|", r"|\bparity\b|"), ("|character|", r"|\bcharacters?\b|")],
        "FINITE_GEOMETRY": [("|Steiner|", r"|\bSteiner (system|triple|quadruple|rule)|")],
        "ABC": [("|Baker|", r"|\bBaker'?s? (theorem|bound|theory|method|estimate)|\bBaker.W|")],
    }
    for k, subs in rep.items():
        for a, b in subs:
            assert a in L[k], (k, a)
            L[k] = L[k].replace(a, b)
    if drop_carrier:
        assert "|carrier|" in L["OBSERVER"]
        L["OBSERVER"] = L["OBSERVER"].replace("|carrier|", "|")
    return L


def compile_lex(L):
    return {k: re.compile(v, re.I) for k, v in L.items()}


def strip_front_matter(txt):
    if txt.startswith("---"):
        j = txt.find("\n---", 3)
        if j > 0:
            return txt[j + 4:]
    return txt


def load(idmin, idmax, body, exclude=(), nofm=False):
    rows = []
    for f in sorted(glob.glob(os.path.join("01-canon", "theorems", "THM-*.md"))):
        b = os.path.basename(f)
        m = re.match(r"THM-(\d+)", b)
        if not m:
            continue
        i = int(m.group(1))
        if i < idmin or i > idmax or "RESERVED" in b or i in exclude:
            continue
        path = os.path.abspath(f)
        if os.name == "nt" and len(path) >= 250:          # THM-4839's file name exceeds MAX_PATH in a long checkout
            path = "\\\\?\\" + path
        with open(path, encoding="utf-8", errors="replace") as fh:
            txt = fh.read()
        if nofm:
            txt = strip_front_matter(txt)
        rows.append((i, b, txt[:body] if body > 0 else txt))
    return rows


def incidence(rows, rx):
    X = np.zeros((len(rows), len(NAMES)), dtype=np.int64)
    for r, (_, _, txt) in enumerate(rows):
        for c, k in enumerate(NAMES):
            if rx[k].search(txt):
                X[r, c] = 1
    return X


def curveball_samples(X, nsamp, seed, trades_per_sample):
    """degree-preserving bipartite curveball chain; yields co-occurrence matrices of successive samples."""
    rng = random.Random(seed)
    rows = [sorted(np.nonzero(X[r])[0].tolist()) for r in range(X.shape[0])]
    m = len(rows)
    nf = X.shape[1]
    for _ in range(nsamp):
        for _ in range(trades_per_sample):
            a = rng.randrange(m)
            b = rng.randrange(m)
            if a == b:
                continue
            A, B = rows[a], rows[b]
            sa, sb = set(A), set(B)
            onlyA = [x for x in A if x not in sb]
            onlyB = [x for x in B if x not in sa]
            if not onlyA or not onlyB:
                continue
            pool = onlyA + onlyB
            rng.shuffle(pool)
            k = len(onlyA)
            common = [x for x in A if x in sb]
            rows[a] = sorted(common + pool[:k])
            rows[b] = sorted(common + pool[k:])
        M = np.zeros((m, nf), dtype=np.int64)
        for r, lst in enumerate(rows):
            M[r, lst] = 1
        yield M.T @ M, M


def census(X, nsamp, seed):
    N = X.shape[0]
    obs = X.T @ X
    nf = X.shape[1]
    S = np.zeros((nsamp, nf, nf))
    rs0, cs0 = X.sum(1), X.sum(0)
    for t, (C, M) in enumerate(curveball_samples(X, nsamp, seed, 5 * N)):
        S[t] = C
        if t == 0:
            assert (M.sum(1) == rs0).all() and (M.sum(0) == cs0).all()      # degree preservation
    mu = S.mean(0)
    sd = S.std(0, ddof=1)
    z = np.where(sd > 0, (obs - mu) / np.where(sd > 0, sd, 1), 0.0)
    pairs = [(i, j) for i in range(nf) for j in range(i + 1, nf)]
    return obs, mu, sd, z, pairs


def summarize(X, obs, z, pairs, out):
    N = X.shape[0]
    n = X.sum(0)
    zz = {(NAMES[i], NAMES[j]) if NAMES[i] < NAMES[j] else (NAMES[j], NAMES[i]): z[i, j] for i, j in pairs}
    order = sorted(zz, key=lambda p: -zz[p])
    deg = {k: sum(1 for p in zz if k in p and zz[p] >= 2) for k in NAMES}
    strength = {k: sum(zz[p] for p in zz if k in p and zz[p] > 0) for k in NAMES}
    mod_share = n[NAMES.index("MODULAR")] / N
    P1 = mod_share >= 0.20
    P2 = deg["MODULAR"] == max(deg.values())
    core = {"MODULAR", "QUADRATIC", "FINITE_GEOMETRY"}
    P3 = any(set(p) <= core for p in order[:3])
    res = dict(N=N, share={k: n[c] / N for c, k in enumerate(NAMES)}, z=zz, order=order, deg=deg, strength=strength,
               P1=P1, P2=P2, P3=P3)
    if out:
        out("  files %d; shares: %s" % (N, ", ".join("%s %.1f%%" % (k, 100 * res["share"][k])
                                                   for k in sorted(NAMES, key=lambda k: -res["share"][k]))))
        out("  top pairs: %s" % "; ".join("%s-%s %.2f" % (p[0], p[1], zz[p]) for p in order[:6]))
        out("  bottom pairs: %s" % "; ".join("%s-%s %.2f" % (p[0], p[1], zz[p]) for p in order[-4:]))
        out("  lift-degree: %s" % ", ".join("%s %d" % (k, deg[k]) for k in sorted(NAMES, key=lambda k: (-deg[k], -strength[k]))))
        out("  least strength: %s (%.2f)" % (min(NAMES, key=lambda k: strength[k]), min(strength.values())))
        out("  P1 %s  P2 %s  P3 %s" % (P1, P2, P3))
    return res


def lexicon(name):
    if name == "heckefix":
        L = dict(REG)
        L["MODULAR"] = L["MODULAR"].replace("|Hecke|", r"|\bHecke\b|")
        return L
    return {"registered": REG, "strict": strict_lexicon(), "strict_nocarrier": strict_lexicon(True)}[name]


def cmd_run(a):
    L = lexicon(a.lexicon)
    rows = load(a.idmin, a.idmax, a.body, set(a.exclude))
    X = incidence(rows, compile_lex(L))
    t0 = time.time()
    obs, mu, sd, z, pairs = census(X, a.null, a.seed)
    print("run %s: window %d..%d body %d lexicon %s seed %d null %d excluded %s (%.1fs)"
          % (a.tag, a.idmin, a.idmax, a.body, a.lexicon, a.seed, a.null, a.exclude, time.time() - t0))
    res = summarize(X, obs, z, pairs, print)
    print("  all pairs (n_obs, null mean, null sd, z):")
    for p in res["order"]:
        i, j = NAMES.index(p[0]), NAMES.index(p[1])
        print("    %-16s %-16s %4d %8.2f %6.2f %6.2f" % (p[0], p[1], obs[i, j], mu[i, j], sd[i, j], res["z"][p]))
    return res


def cmd_terms(a):
    rows = load(a.idmin, a.idmax, a.body, set(a.exclude))
    print("per-alternative match counts (files whose first %d chars match), window %d..%d, %d files"
          % (a.body, a.idmin, a.idmax, len(rows)))
    for k, pat in REG.items():
        alts = split_alternation(pat)
        cnt = []
        fam = re.compile(pat, re.I)
        famhits = sum(1 for r in rows if fam.search(r[2]))
        for alt in alts:
            rx = re.compile(alt, re.I)
            c = sum(1 for r in rows if rx.search(r[2]))
            only = sum(1 for r in rows if rx.search(r[2]) and
                       not any(re.search(o, r[2], re.I) for o in alts if o != alt))
            cnt.append((c, only, alt))
        cnt.sort(reverse=True)
        print("  %s: %d files; top alternatives (files, files matched ONLY by this alternative):" % (k, famhits))
        for c, only, alt in cnt[:9]:
            print("      %4d %4d  %s" % (c, only, alt))


def split_alternation(pat):
    out, depth, cur = [], 0, ""
    i = 0
    while i < len(pat):
        ch = pat[i]
        if ch == "\\":
            cur += pat[i:i + 2]
            i += 2
            continue
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        if ch == "|" and depth == 0:
            out.append(cur)
            cur = ""
        else:
            cur += ch
        i += 1
    out.append(cur)
    return out


def cmd_spot(a):
    rows = load(a.idmin, a.idmax, a.body, set(a.exclude))
    rng = random.Random(a.seed)
    pick = rng.sample(rows, a.k)
    rx = compile_lex(REG)
    for i, b, txt in sorted(pick):
        print("=" * 110)
        print("THM-%d  %s" % (i, b[:100]))
        for k in NAMES:
            ms = []
            for alt in split_alternation(REG[k]):
                for mm in re.finditer(alt, txt, re.I):
                    s0 = max(0, mm.start() - 50)
                    ctx = txt[s0:mm.end() + 50].replace("\n", " ")
                    ms.append((mm.start(), mm.group(0), ctx))
                    break
            if ms:
                ms.sort()
                print("  %s:" % k)
                for pos, g, ctx in ms[:4]:
                    print("      [%s] ...%s..." % (g, ctx))


def cmd_robust(a):
    """the grid reported in AUDIT.md (several minutes)."""
    out = open(os.path.join(HERE, "audit_census_robustness.out"), "w", encoding="utf-8")

    def say(s):
        print(s, flush=True)
        out.write(s + "\n")
        out.flush()

    def one(tag, idmin=4400, idmax=4839, body=6000, lex="registered", seed=20261028, null=2000, exclude=(4839,),
            nofm=False):
        L = lexicon(lex)
        rows = load(idmin, idmax, body, set(exclude), nofm)
        X = incidence(rows, compile_lex(L))
        obs, mu, sd, z, pairs = census(X, null, seed)
        say("-- %s: window %d..%d body %s lexicon %s seed %d null %d exclude %s front-matter-stripped %s"
            % (tag, idmin, idmax, body, lex, seed, null, list(exclude), nofm))
        return summarize(X, obs, z, pairs, say)

    t0 = time.time()
    table = []
    table.append(("registered (320 files, as pre-registered)", one("R0")))
    for s in (1, 2, 3, 4, 5, 6, 7, 8, 9):
        table.append(("seed %d" % s, one("seed%d" % s, seed=s)))
    table.append(("THM-4839 included (321 files)", one("incl4839", exclude=())))
    table.append(("BODY 3000", one("body3000", body=3000)))
    table.append(("BODY 12000", one("body12000", body=12000)))
    table.append(("BODY whole file", one("bodyall", body=0)))
    table.append(("registered, front matter stripped", one("nofm", nofm=True)))
    table.append(("strict lexicon", one("strict", lex="strict")))
    for s_ in (1, 2, 3):
        table.append(("strict lexicon, seed %d" % s_, one("strict_seed%d" % s_, lex="strict", seed=s_)))
    table.append(("strict lexicon without 'carrier'", one("strict_nocarrier", lex="strict_nocarrier")))
    table.append(("strict lexicon, front matter stripped", one("strict_nofm", lex="strict", nofm=True)))
    table.append(("strict lexicon, BODY 3000", one("strict3000", lex="strict", body=3000)))
    table.append(("strict lexicon, BODY 12000", one("strict12000", lex="strict", body=12000)))
    table.append(("full canon, registered (300 null)", one("full", idmin=1, idmax=4839, null=300)))
    table.append(("full canon, strict (300 null)", one("fullstrict", idmin=1, idmax=4839, null=300, lex="strict")))
    table.append(("only 'Hecke' fixed (registered otherwise)", one("heckefix", lex="heckefix")))
    say("")
    say("SUMMARY (z of key pairs; lift-degree leaders; least-strength family; verdicts)")
    keyp = [("FINITE_GEOMETRY", "TOURNAMENT"), ("COLLATZ", "MODULAR"), ("BINARY_ALGEBRA", "TOURNAMENT"),
            ("BINARY_ALGEBRA", "FINITE_GEOMETRY"), ("COLLATZ", "OBSERVER"), ("MODULAR", "OBSERVER"),
            ("FINITE_GEOMETRY", "MODULAR")]
    say("%-44s %s | %-34s | %-15s | P1 P2 P3 | top pair" % ("run", " ".join("%7s" % ("".join(w[0] for w in p[0].split("_")) + "-" + "".join(w[0] for w in p[1].split("_"))) for p in keyp), "lift-degree max", "least"))
    for name, r in table:
        mx = max(r["deg"].values())
        lead = ",".join(k[:6] for k in NAMES if r["deg"][k] == mx)
        least = min(NAMES, key=lambda k: r["strength"][k])
        say("%-44s %s | %-34s | %-15s | %s %s %s | %s-%s"
            % (name, " ".join("%7.2f" % r["z"][p] for p in keyp), "%d: %s" % (mx, lead), least,
               "T " if r["P1"] else "F ", "T " if r["P2"] else "F ", "T " if r["P3"] else "F ",
               r["order"][0][0], r["order"][0][1]))
    say("total %.0f s" % (time.time() - t0))
    out.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "terms", "spot", "robust"])
    ap.add_argument("--idmin", type=int, default=4400)
    ap.add_argument("--idmax", type=int, default=4839)
    ap.add_argument("--null", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20261028)
    ap.add_argument("--body", type=int, default=6000)
    ap.add_argument("--lexicon", default="registered")
    ap.add_argument("--exclude", type=int, nargs="*", default=[4839])
    ap.add_argument("--tag", default="")
    ap.add_argument("--k", type=int, default=10)
    a = ap.parse_args()
    {"run": cmd_run, "terms": cmd_terms, "spot": cmd_spot, "robust": cmd_robust}[a.cmd](a)


if __name__ == "__main__":
    main()
