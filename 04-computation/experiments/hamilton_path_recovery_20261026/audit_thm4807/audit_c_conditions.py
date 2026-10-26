"""Auditor's check that none of the three conditions of THM-4807 (c) is redundant.

For every arc u->v of every strong class n <= 8 (gentourng -c) compute
  cond1: T-u has k >= 2 components and v lies in the last one, C_k;
  cond2: u dominates C_(k-1);
  cond3: C_k has no partition into a nonempty path ending at an in-neighbour of u
         and a path starting at v;
and tabulate the patterns.  An arc is untraced iff cond1 & cond2 & cond3 (checked again).
No code from the lane is imported.
"""

import subprocess
from collections import Counter

from audit_lib import bits_iter, from_upper_string, hp_tables, in_masks, scc_ordered, traced_from_tables


def gentourng(n):
    res = subprocess.run(["gentourng", "-q", "-c", str(n)], capture_output=True, text=True, check=True)
    return [s.strip() for s in res.stdout.split("\n") if s.strip()]


def main():
    for n in range(3, 9):
        pat = Counter()
        bad = 0
        for s in gentourng(n):
            out = from_upper_string(n, s)
            inn = in_masks(n, out)
            F, B = hp_tables(n, out, inn)
            tr = traced_from_tables(n, out, F, B)
            full = (1 << n) - 1
            for u in range(n):
                cu, _ = scc_ordered(n, out, full & ~(1 << u))
                for v in bits_iter(out[u]):
                    c1 = len(cu) >= 2 and bool((cu[-1] >> v) & 1)
                    c2 = c3 = None
                    if c1:
                        Ck = cu[-1]
                        c2 = (cu[-2] & ~out[u]) == 0
                        c3 = True
                        rest = Ck & ~(1 << v)
                        A = rest
                        while A:
                            if (F[A] & inn[u]) and (B[Ck & ~A] >> v) & 1:
                                c3 = False
                                break
                            A = (A - 1) & rest
                    untr = not (tr[u] >> v) & 1
                    pat[(c1, c2, c3, untr)] += 1
                    if untr != bool(c1 and c2 and c3):
                        bad += 1
        print(f"n = {n}: (cond1, cond2, cond3, untraced) -> arcs: "
              f"{ {str(k): v for k, v in sorted(pat.items(), key=str)} }; mismatches {bad}")
        c1only = sum(v for k, v in pat.items() if k[0])
        print(f"    arcs with cond1: {c1only}; cond1 & cond3 but not cond2: "
              f"{sum(v for k, v in pat.items() if k[0] and k[2] and not k[1])}; cond1 & cond2 but not cond3: "
              f"{sum(v for k, v in pat.items() if k[0] and k[1] and not k[2])}")


if __name__ == "__main__":
    main()
