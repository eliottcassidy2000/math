#!/usr/bin/env python3
"""Supplementary adversarial check for claim B5 (kind-pasteur-2026-06-10-S1, thread B):
"the only Paley-specific residue is the SIZE |Aut(T_p)| = p(p-1)/2".

Independent checks (exact integers, fresh code):
  1. p=7: |Aut(T_7)| by BRUTE FORCE over all 5040 permutations of S_7
     (no affine-map assumption); expect 21 = 7*6/2.
     Also enumerate the 189 directed Ham paths of T_7 directly and build the
     explicit Aut-orbit partition: expect 9 orbits, all of size exactly 21
     (LEM-003 freeness on a genuine Paley instance, beyond the n<=6 sweep).
  2. p=11: brute force over S_11 is infeasible; instead verify
     (a) every affine map x -> ax+b with a in QR(11) is an automorphism
         (gives |Aut| >= 55 = 11*10/2), and
     (b) consistency with canon: H(T_11) = 95095 (THM-212 table) and
         95095 / 55 = 1729 = r(11), the documented Paley ratio; freeness
         (LEM-003) makes 95095/|Aut| an integer, and |Aut| = 55 is the
         canon value (b is flagged as consistency, not independent proof).
"""

import itertools

# ---------------- Paley tournament builder ----------------
def paley(p):
    qr = {(x * x) % p for x in range(1, p)}
    arcs = {(u, v) for u in range(p) for v in range(p)
            if u != v and (v - u) % p in qr}
    return arcs, qr

# ---------------- p = 7 ----------------
p = 7
arcs7, qr7 = paley(p)
assert all(((u, v) in arcs7) != ((v, u) in arcs7)
           for u in range(p) for v in range(u + 1, p)), "tournament check"

aut7 = [sig for sig in itertools.permutations(range(p))
        if all((sig[u], sig[v]) in arcs7 for (u, v) in arcs7)]
print(f"p=7: QR = {sorted(qr7)}")
print(f"p=7: |Aut(T_7)| brute force over S_7 = {len(aut7)}  (claim p(p-1)/2 = 21)")

paths7 = []
def extend(seq, used):
    if len(seq) == p:
        paths7.append(tuple(seq))
        return
    for w in range(p):
        if not used[w] and (seq[-1], w) in arcs7:
            used[w] = True
            seq.append(w)
            extend(seq, used)
            seq.pop()
            used[w] = False

for s in range(p):
    used = [False] * p
    used[s] = True
    extend([s], used)
H7 = len(paths7)
pset = set(paths7)
seen = set()
orbit_sizes = []
for q in paths7:
    if q in seen:
        continue
    orb = {tuple(sig[v] for v in q) for sig in aut7}
    assert orb <= pset
    seen |= orb
    orbit_sizes.append(len(orb))
print(f"p=7: H(T_7) = {H7} (canon: 189), orbits = {len(orbit_sizes)}, "
      f"all of size 21: {all(s == 21 for s in orbit_sizes)}")
ok7 = (len(aut7) == 21 and H7 == 189 and len(orbit_sizes) == 9
       and all(s == 21 for s in orbit_sizes))
print(f"p=7 verdict: {'PASS' if ok7 else 'FAIL'}")

# ---------------- p = 11 ----------------
p = 11
arcs11, qr11 = paley(p)
aff = []
for a in sorted(qr11):
    for b in range(p):
        sig = tuple((a * x + b) % p for x in range(p))
        if all((sig[u], sig[v]) in arcs11 for (u, v) in arcs11):
            aff.append(sig)
print(f"p=11: affine maps x->ax+b, a in QR, that are automorphisms: "
      f"{len(aff)} of {len(qr11) * p}  (claim all 55 = p(p-1)/2)")
H11_canon = 95095
print(f"p=11: canon H(T_11) = {H11_canon} = 5*7*11*13*19: "
      f"{H11_canon == 5 * 7 * 11 * 13 * 19}")
print(f"p=11: H/55 = {H11_canon // 55} (= 1729, the documented r(11)): "
      f"{H11_canon % 55 == 0 and H11_canon // 55 == 1729}")
ok11 = (len(aff) == 55 and H11_canon % 55 == 0 and H11_canon // 55 == 1729)
print(f"p=11 verdict: {'PASS' if ok11 else 'FAIL'} "
      f"(lower bound |Aut|>=55 independent; |Aut|=55 exactly is canon-consistent)")

print(f"OVERALL: {'PASS' if ok7 and ok11 else 'FAIL'}")
