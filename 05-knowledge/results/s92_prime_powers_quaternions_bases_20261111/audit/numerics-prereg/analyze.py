# Auditor numerics-prereg: compare my census json with a stored session output (cum_*.out or window_*.out).
# usage: python3 -B analyze.py my_json session_out [prereg_json]
import sys
import re
import json
import math
from aud_common import divisors, remainder, groups, L_m


def load(fn):
    J = json.load(open(fn))
    T = J["T"]
    mods = {}
    for ms, v in J["moduli"].items():
        m = int(ms)
        cls = v["cls"]
        typ = [{int(k): x for k, x in v["typ"][sh].items()} for sh in range(2)]
        diag = [{int(k): x for k, x in v["diag"][sh].items()} for sh in range(2)]
        mods[m] = (cls, typ, diag)
    return J, T, mods


def parse_out(fn):
    txt = open(fn).read().splitlines()
    info = {"T": None, "table": {}, "groups": {}, "dev": {}, "lines": txt}
    for ln in txt:
        mo = re.search(r"triples sheet 0: (\d+), sheet 1: (\d+)", ln)
        if mo:
            info["T"] = [int(mo.group(1)), int(mo.group(2))]
        mo = re.match(r"^\s+(\d+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s*$", ln)
        if mo:
            info["table"][int(mo.group(1))] = (mo.group(2), mo.group(3), mo.group(4))
        mo = re.match(r"^\s+m =\s+(\d+) sheet (\d): (.*)$", ln)
        if mo:
            m = int(mo.group(1))
            sh = int(mo.group(2))
            d = {}
            for part in mo.group(3).split(";"):
                part = part.strip()
                mm = re.match(r"g0=(\d+): ([0-9.]+) \[([0-9.]+), ([0-9.]+)\] \| (.*)$", part)
                g0 = int(mm.group(1))
                rest = [t.strip() for t in mm.group(5).split("|")]
                d[g0] = (mm.group(2), mm.group(3), mm.group(4), rest)
            info["groups"][(m, sh)] = d
        mo = re.match(r"^\s+m =\s+(\d+): ([0-9.]+) / ([0-9.]+)\s+flatness ([0-9.]+) / ([0-9.]+)", ln)
        if mo:
            info["dev"][int(mo.group(1))] = (mo.group(2), mo.group(3), mo.group(4), mo.group(5))
    return info


def agree(mine, printed):
    """mine: float; printed: string with fixed decimals.  True if printed is the rounding of mine
    (tolerance: 0.5000001 unit of the last printed digit)."""
    dec = len(printed.split(".")[1]) if "." in printed else 0
    return abs(mine - float(printed)) <= 0.5000001 * 10 ** (-dec) + 1e-12


def main():
    J, T, mods = load(sys.argv[1])
    S = parse_out(sys.argv[2])
    P = json.load(open(sys.argv[3])) if len(sys.argv) > 3 else None
    nbad = 0
    ncmp = 0
    print("my file %s  vs  session file %s" % (sys.argv[1], sys.argv[2]))
    print("triple counts: mine %s  session %s  %s" % (T, S["T"], "EQUAL" if T == S["T"] else "DIFFERENT"))
    if T != S["T"]:
        nbad += 1
    # which diagonal convention reproduces the stored remainder table?
    if S["table"]:
        for conv in (False, True):
            bad = 0
            worstdiff = 0.0
            for m, (a, b, nz) in S["table"].items():
                for sh, pr in ((0, a), (1, b)):
                    r = remainder(m, mods[m][0], mods[m][1], mods[m][2], sh, conv)
                    worstdiff = max(worstdiff, abs(r - float(pr)))
                    if not agree(r, pr):
                        bad += 1
            print("remainder table, diagonal %s the type counts: %d of %d entries disagree with the stored 5-decimal values (largest difference %.2e)" % (
                "INCLUDED in" if conv else "EXCLUDED from", bad, 2 * len(S["table"]), worstdiff))
        print("my remainders (diagonal excluded as in the statement of THM-5134 (B) / included), and noise sqrt(m/#T) per sheet:")
        for m, (a, b, nz) in S["table"].items():
            r0 = remainder(m, *mods[m], 0, False)
            r1 = remainder(m, *mods[m], 1, False)
            q0 = remainder(m, *mods[m], 0, True)
            q1 = remainder(m, *mods[m], 1, True)
            n0 = math.sqrt(m / T[0])
            n1 = math.sqrt(m / T[1])
            okn = agree(n0, nz) or agree(n1, nz) or agree(math.sqrt(m / (T[0] + T[1])), nz) or agree(math.sqrt(m / min(T)), nz)
            print("  %4d  excl %.5f / %.5f   incl %.5f / %.5f   stored %s / %s   noise sheet0 %.5f sheet1 %.5f stored %s %s" % (
                m, r0, r1, q0, q1, a, b, n0, n1, nz, "" if okn else "NOISE-MISMATCH"))
    # exactness: every triple of type d lies in a class with d | c - c_sigma: the classes c with
    # gcd(m, c - c_sigma) = g0 must have count >= ... (cannot be tested from marginals alone; see my_exact_check)
    # groups
    for conv in (False, True):
        bad = 0
        tot = 0
        worst = 0.0
        badlist = []
        for (m, sh), d in S["groups"].items():
            G = groups(m, mods[m][0], mods[m][1], mods[m][2], sh, conv)
            for g0, (mean, mn, mx, rest) in d.items():
                mine = G[g0]
                items = [(mine[0], mean, "mean"), (mine[1], mn, "min"), (mine[2], mx, "max"), (mine[4], rest[-1], "limit")]
                # first-order value with measured types: cum format rest = [phi, limit]; window format rest = [apriori, phi, limit]
                items.append((mine[3], rest[-2], "first-order"))
                for val, pr, name in items:
                    tot += 1
                    worst = max(worst, abs(val - float(pr)))
                    if not agree(val, pr):
                        bad += 1
                        badlist.append((m, sh, g0, name, val, pr))
        print("group lines (mean, min, max, first-order value, limit), diagonal %s: %d of %d printed numbers disagree (largest difference %.2e)" % (
            "included" if conv else "excluded", bad, tot, worst))
        for it in badlist[:12]:
            print("     m=%d sheet %d g0=%d %s: mine %.6f stored %s" % it)
    # a-priori column against the registered json, and the P1/P2 summary lines
    if P is not None:
        bad = 0
        tot = 0
        for (m, sh), d in S["groups"].items():
            for g0, (mean, mn, mx, rest) in d.items():
                key = "%d,%d,%d" % (sh, m, g0)
                if key in P["phi"] and len(rest) == 3:
                    tot += 1
                    if not agree(P["phi"][key], rest[0]):
                        bad += 1
                        print("   a-priori column differs from registered json at", key, P["phi"][key], rest[0])
        print("a-priori column of the stored output vs registered json: %d of %d disagree" % (bad, tot))
    return nbad


if __name__ == "__main__":
    main()
