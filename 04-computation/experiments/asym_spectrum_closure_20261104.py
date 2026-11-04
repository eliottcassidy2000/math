"""asym_spectrum_closure_20261104.py (kind-pasteur-2026-11-04-S10) -- HYP-9462 (A).

Reads
  autspec output  ("S n=.. H=.. hc=.. aut=.. sc=.." lines from asym_spectrum_autspec_20261104 autspec NMAX HMAX)
  UH values       ("UH n=.. classes=.. H: ..." lines from uh_hpath_values_20261104)
and reports, for odd values up to LIMIT:
  * the holes of the full spectrum (closure of all strong values under products, plus 1);
  * the holes of the asymmetric spectrum (closure of strong values with |Aut| = 1, plus the UH values, which are
    asymmetric by THM-4868 (E));
  * strong values attained only with a nontrivial automorphism group;
  * the |Aut| values seen, and any violation of |Aut| | H (THM-4747).
Aut of a non-strong tournament is the product of the Aut of its strong components (they are linearly ordered), and
H is multiplicative over strong components, so these closures are exactly the spectra built from the listed strong
classes.  Completeness up to LIMIT needs every strong class with H <= LIMIT, i.e. orders up to the Busch bound.

Usage: python3 asym_spectrum_closure_20261104.py AUTSPEC_FILE UH_FILE LIMIT
"""
import re
import sys


def closure(base, lim):
    s = {1}
    frontier = [1]
    while frontier:
        nxt = []
        for a in frontier:
            for b in base:
                c = a * b
                if c <= lim and c not in s:
                    s.add(c)
                    nxt.append(c)
        frontier = nxt
    return s


def main():
    autspec, uhfile, lim = sys.argv[1], sys.argv[2], int(sys.argv[3])
    rows = [tuple(int(x) for x in re.findall(r"=(\d+)", l)) for l in open(autspec) if l.startswith("S ")]
    strong_all, strong_asym, auts, viol = set(), set(), set(), 0
    for n, H, hc, aut, sc in rows:
        if H % aut:
            viol += 1
        if H > lim:
            continue
        strong_all.add(H)
        auts.add(aut)
        if aut == 1:
            strong_asym.add(H)
    uh = set()
    for l in open(uhfile):
        if l.startswith("UH n="):
            uh.update(int(x) for x in l.split("H:")[1].split() if int(x) <= lim)
    odds = set(range(1, lim + 1, 2))
    full = closure(strong_all, lim)
    asym = closure(strong_asym | uh, lim)
    uh_only = closure(uh, lim)
    print(f"strong classes read: {len(rows)}; LIMIT = {lim}")
    print("full-spectrum holes:", sorted(odds - full))
    print("asymmetric-spectrum holes (strong asymmetric of the listed orders + UH):", sorted(odds - asym))
    print("strong values attained only with a nontrivial automorphism group:", sorted(h for h in strong_all if h not in strong_asym))
    print("|Aut| values among the listed strong classes:", sorted(auts))
    print("violations of |Aut| | H:", viol)
    print("holes of the UH-product closure alone:", len(odds - uh_only), sorted(odds - uh_only)[:40])


if __name__ == "__main__":
    main()
