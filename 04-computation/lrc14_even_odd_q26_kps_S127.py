# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont45: the EVEN/ODD structural handle at the composite modulus q=26 for Route B.
#
# CONTEXT (klein-S258 finish map, opus-S239): LRC(14) reduces to Route B spread-bulk anti-concentration --
# every spread divisor-complete (DC) family clears at some non-14 q in [15,29] (=> M>1/14, band-edge opus-S235).
# opus-S238 analyzed PRIME moduli 17,19,23 via fold-classes and found NO small-prime shortcut. THIS explores a
# DIFFERENT, COMPOSITE handle opus's prime analysis missed:
#
# AT q = 26 = 2*13: an EVEN runner v=2w gives v*p mod 26 = 2*(w*p mod 13) = ALWAYS EVEN. The danger arc is
# {0,1,25} = {0,+-1}; the only EVEN danger residue is 0. So an even runner is in danger ONLY when v*p === 0
# mod 26, i.e. (if gcd(w,13)=1) only at p=13. => at q=26, for every multiplier p != 13, ALL even runners are
# SAFE. DC families are FORCED to contain many even runners (mults of 2,4,6,8,10,12,14 are needed), so at
# q=26 the danger for p!=13 is carried ENTIRELY by the ODD runners -- a much smaller set => easier to clear.
#
# TEST: (1) how many odd runners do DC families have? (2) does q=26 clear DC families, and is it predicted by
# the odd count / the odd runners missing a fold-class? (3) combined with the parity structure, is there a
# clean sub-class guarantee?
from math import gcd
from functools import reduce
from itertools import combinations
from collections import Counter

def is_DC(v):  return all(any(x % d == 0 for x in v) for d in range(2, 15))
def prim(v):   return reduce(gcd, v) == 1
def longest_run(v):
    v = sorted(v); b = m = 1
    for i in range(1, len(v)):
        if v[i] == v[i-1] + 1: m += 1; b = max(b, m)
        else: m = 1
    return b
def clears_at(v, q):
    lo = -(-q // 14)
    return any(all(lo <= (vi * p) % q <= q - lo for vi in v) for p in range(1, q))
def clear_p_at(v, q):
    lo = -(-q // 14)
    return [p for p in range(1, q) if all(lo <= (vi * p) % q <= q - lo for vi in v)]

def main():
    # sample spread (longest-run <= 7) primitive DC families
    fams = []
    for v in combinations(range(1, 30), 13):
        v = list(v)
        if prim(v) and is_DC(v) and longest_run(v) <= 7:
            fams.append(v)
    print(f"spread primitive DC families (Vmax<=29, longest-run<=7): {len(fams)}")
    if not fams:
        print("  (none in this small range; widening)"); return

    # (1) odd-runner count distribution
    oddc = Counter(sum(1 for x in v if x % 2 == 1) for v in fams)
    print(f"(1) odd-runner-count distribution: {dict(sorted(oddc.items()))}")
    print(f"    => DC forces even runners: max odd count = {max(oddc)} of 13 (>= {13-max(oddc)} even, structural)")

    # (2) q=26 clearing, and the even/odd mechanism
    q = 26
    cl26 = sum(1 for v in fams if clears_at(v, q))
    print(f"(2) clears at q=26: {cl26}/{len(fams)} ({100*cl26/len(fams):.1f}%)")
    # verify the mechanism: at q=26, even runners are in danger ONLY at p=13
    bad_mech = 0
    for v in fams:
        for p in range(1, 26):
            if p == 13: continue
            for vi in v:
                if vi % 2 == 0 and (vi * p) % 26 in (0, 1, 25):
                    bad_mech += 1
    print(f"    MECHANISM CHECK: even runner in danger at q=26, p!=13: {bad_mech} occurrences (expect 0)")
    # for the non-clearers at 26, how many odd runners? (odd runners must cover all p!=13 fold-classes)
    nonclear26 = [v for v in fams if not clears_at(v, q)]
    if nonclear26:
        oc = Counter(sum(1 for x in v if x % 2 == 1) for v in nonclear26)
        print(f"    non-clearers at q=26: {len(nonclear26)}; their odd-count dist: {dict(sorted(oc.items()))}")

    # (3) does the WHOLE window [15,29] clear every family, and where does q=26 rank?
    window = [q for q in range(15, 30) if q % 14]
    allclear = 0; q26_helps = 0
    worst = None
    for v in fams:
        cq = [q for q in window if clears_at(v, q)]
        if cq: allclear += 1
        else:
            worst = v
        if 26 in cq: q26_helps += 1
    print(f"(3) window [15,29] non-14: clears at some q: {allclear}/{len(fams)}; q=26 in the clearing set: {q26_helps}/{len(fams)}")
    if worst: print(f"    !! family clearing NOWHERE in window: {worst}")
    else: print(f"    every spread DC family clears somewhere in [15,29] (0 exceptions) -- Route B holds on sample")

if __name__ == "__main__":
    main()
