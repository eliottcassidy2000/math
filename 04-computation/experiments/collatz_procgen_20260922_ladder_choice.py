#!/usr/bin/env python3
"""collatz_procgen_20260922_ladder_choice.py -- the choice ladder for q = 3, 5, 7.

Runs lane one's program collatz_procgen_20260922_exceptional_general.c (unchanged; compiled
into a temporary directory):   excgen q b MM mode [J r1 r2 ...]
   mode 0: no choice (the q x + 1 map);  mode 1: excursion x -> q x + 1 -> q^2 x + q + 1 allowed
   at every even x (graph E_q);  mode 2: excursion allowed iff x mod 2^J lies in S.
It prints the number of exceptional classes (no multiplicative descent within precision m)
at m = MM-8, MM-4, MM.

Checks:
  * mode-0 counts == the exact ballot numbers N_m(q) (collatz_procgen_20260922_ladder_ballot.py);
  * q=3 mode-1 and the two partial-choice sets reproduce lane one's independent programs
    (e_forward_dp_full, e_forward_dp) as recorded in collatz_procgen_20260922_choice_ladder.out.
Memory: the MM=26 runs need about 0.65 GB; runs are sequential.
"""
import importlib.util
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
SRC = os.path.join(HERE, "collatz_procgen_20260922_exceptional_general.c")
LANE1_OUT = os.path.join(ROOT, "05-knowledge", "results", "collatz_procgen_20260922_choice_ladder.out")

spec = importlib.util.spec_from_file_location("ballot", os.path.join(HERE, "collatz_procgen_20260922_ladder_ballot.py"))
ballot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ballot)

LINE = re.compile(r"q=(\d+) b=(-?\d+) mode=(\d) m=\s*(\d+) exceptional=(\d+) frac=(\S+) log2/m=(\S+)")


def run(exe, q, MM, mode, J=None, S=()):
    args = [exe, str(q), "1", str(MM), str(mode)]
    if mode == 2:
        args += [str(J)] + [str(r) for r in S]
    out = subprocess.run(args, capture_output=True, text=True, check=True).stdout
    res = {}
    for g in LINE.finditer(out):
        res[int(g.group(4))] = int(g.group(5))
    return res


def lane1_values():
    """q=3 values recorded by lane one: full E (section 2) and partial choice (section 3)."""
    E, part = {}, {}
    if not os.path.exists(LANE1_OUT):
        return E, part
    txt = open(LANE1_OUT).read()
    sec2 = txt.split("### 2.")[1].split("### 3.")[0]
    for m, bad in re.findall(r"m=\s*(\d+) mod 2\^\d+: bad=(\d+)", sec2):
        E[int(m)] = int(bad)
    sec3 = txt.split("### 3.")[1].split("### 4.")[0]
    for blk in sec3.split("--- J,S = ")[1:]:
        head, body = blk.split("\n", 1)
        key = tuple(int(t) for t in head.split())
        part[key] = {int(m): int(b) for m, b in re.findall(r"m=(\d+) bad=(\d+)", body)}
    return E, part


def main():
    tmp = tempfile.mkdtemp(prefix="ladder_choice_")
    exe = os.path.join(tmp, "excgen")
    subprocess.run(["clang", "-O3", "-o", exe, SRC, "-lm"], check=True)
    fails = 0
    print("=" * 78)
    print("LADDER-CHOICE. exceptional classes mod 2^m (lane-one program exceptional_general.c)")
    print("=" * 78)
    configs = [
        ("q=3 no choice (Collatz)", 3, 0, None, ()),
        ("q=3 E_S, S={6 mod 8}", 3, 2, 3, (6,)),
        ("q=3 E_S, S={2 mod 4}", 3, 2, 2, (2,)),
        ("q=3 full choice E", 3, 1, None, ()),
        ("q=5 no choice (5x+1)", 5, 0, None, ()),
        ("q=5 E_S, S={6 mod 8}", 5, 2, 3, (6,)),
        ("q=5 E_S, S={2 mod 4}", 5, 2, 2, (2,)),
        ("q=5 full choice E_5", 5, 1, None, ()),
        ("q=7 no choice (7x+1)", 7, 0, None, ()),
        ("q=7 full choice E_7", 7, 1, None, ()),
    ]
    table = {}
    for name, q, mode, J, S in configs:
        table[name] = run(exe, q, 24, mode, J, S)
    exact = {q: ballot.exact_counts(q, 26) for q in (3, 5, 7)}
    print(f"{'system':<26} {'m=16':>10} {'m=20':>10} {'m=24':>10} {'frac m=20':>10} {'frac m=24':>10} {'log2/m @24':>10}")
    for name, q, mode, J, S in configs:
        r = table[name]
        print(f"{name:<26} {r[16]:>10} {r[20]:>10} {r[24]:>10} {r[20]/2**20:>10.5f} {r[24]/2**24:>10.5f} "
              f"{(__import__('math').log2(r[24])/24):>10.3f}")
        if mode == 0:
            ok = all(r[m] == exact[q][m] for m in (16, 20, 24))
            fails += not ok
            print(f"{'':<26} mode-0 counts == exact ballot N_m({q}) at m=16,20,24: {'YES' if ok else 'NO'}")
    print()
    print("longer runs (MM=26; m = 18, 22, 26):")
    for name, q, mode in (("q=3 full choice E", 3, 1), ("q=5 full choice E_5", 5, 1), ("q=5 no choice (5x+1)", 5, 0)):
        r = run(exe, q, 26, mode)
        fr = "  ".join(f"m={m}: {r[m]} (frac {r[m]/2**m:.5f})" for m in (18, 22, 26))
        print(f"  {name:<24} {fr}")
        if mode == 0:
            ok = all(r[m] == exact[q][m] for m in (18, 22, 26))
            fails += not ok
            print(f"  {'':<24} == exact ballot counts: {'YES' if ok else 'NO'}")
        if q == 5 and mode == 1:
            r5_26 = r
    mu5 = 0.176025784562122695129
    print(f"  q=5 reference: mu_5 = Haar(Bad(5)) = {mu5:.10f} (no-choice limit, LADDER-C)")
    print(f"  q=5 ratio (full choice)/(no choice) at m=26: {r5_26[26]/exact[5][26]:.4f}")
    print()
    print("cross-checks against lane one's independent programs (choice_ladder.out):")
    E, part = lane1_values()
    if not E:
        print("  lane-one output not found; skipped")
    else:
        r22 = run(exe, 3, 22, 1)
        r24 = table["q=3 full choice E"]
        got = {**r22, **r24}
        pairs = [(m, got[m], E[m]) for m in sorted(got) if m in E]
        ok = all(a == b for _, a, b in pairs)
        fails += not ok
        print(f"  q=3 full E: general program vs e_forward_dp_full at m={[m for m,_,_ in pairs]}: "
              f"{'IDENTICAL' if ok else 'MISMATCH'} {[(a, b) for _, a, b in pairs]}")
        r26 = run(exe, 3, 26, 1)
        if 26 in E:
            ok = r26[26] == E[26]
            fails += not ok
            print(f"  q=3 full E at m=26: {r26[26]} vs lane one {E[26]}: {'IDENTICAL' if ok else 'MISMATCH'}")
        for key, (J, S) in {(3, 6): (3, (6,)), (2, 2): (2, (2,))}.items():
            rr = run(exe, 3, 22, 2, J, S)
            ref = part.get(key, {})
            pairs = [(m, rr[m], ref[m]) for m in (18, 22) if m in ref and m in rr]
            ok = pairs and all(a == b for _, a, b in pairs)
            fails += not ok
            print(f"  q=3 S={{{S[0]} mod {2**J}}}: general program vs e_forward_dp at m=18,22: "
                  f"{'IDENTICAL' if ok else 'MISMATCH'} {[(a, b) for _, a, b in pairs]}")
    shutil.rmtree(tmp)
    print(f"LADDER-CHOICE TOTAL: {'ALL CHECKS PASS' if not fails else 'SOME CHECK FAILED'}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
