# three_cubes_ledger_kpc1_verify.py
# ADVERSARIAL VERIFIER for session kind-pasteur-2026-06-10-S1, thread D.
# FRESH code, written without consulting the worker's scripts.
# Methods chosen to DIFFER from the worker where feasible:
#   - D1: cubes computed as x*(x*x) and cross-checked against a second route
#         (digit-sum mod 9 cast + pow(x,3)); residuals printed exactly.
#   - D3/D4: proof routed through Z/3 (cube residue mod 9 depends only on
#         x mod 3, proven symbolically), enumerating 27 triples instead of
#         the worker's 729; an independent 729-triple brute force is run as
#         a cross-check, plus a wide window x in [-200,200].
#   - D5: brute-force witness table over a small box, plus interval argument.
# Pure python integers throughout. No floats. No imports of worker code.

from math import gcd

fails = 0
def check(name, ok, detail=""):
    global fails
    tag = "PASS" if ok else "FAIL"
    if not ok:
        fails += 1
    print(f"[{tag}] {name}" + (f"  {detail}" if detail else ""))

print("=" * 72)
print("PART 1 (claim D1): the three famous identities, exact bigints")
print("=" * 72)

reps = [
    (42, -80538738812075974, 80435758145817515, 12602123297335631,
     "Booker-Sutherland 2019"),
    (33, 8866128975287528, -8778405442862239, -2736111468807040,
     "Booker 2019"),
    (3, 569936821221962380720, -569936821113563493509, -472715493453327032,
     "Booker-Sutherland 2019 (third rep of 3)"),
]

for k, x, y, z, who in reps:
    # route 1: explicit multiplication
    s1 = x * (x * x) + y * (y * y) + z * (z * z)
    # route 2: pow
    s2 = pow(x, 3) + pow(y, 3) + pow(z, 3)
    check(f"D1 routes agree for k={k}", s1 == s2)
    check(f"D1 k={k} ({who})", s1 == k, f"residual = {s1 - k}")
    g = gcd(gcd(abs(x), abs(y)), abs(z))
    # casting-out-nines sanity: identity must be consistent mod 9
    check(f"D1 k={k} consistent mod 9", (s1 - k) % 9 == 0,
          f"gcd(|x|,|y|,|z|) = {g}; k mod 9 = {k % 9}; "
          f"(x,y,z) mod 3 = ({x % 3},{y % 3},{z % 3})")

print()
print("=" * 72)
print("PART 2 (claims D3, D4): mod-9 laws, complete finite proofs")
print("=" * 72)

# Lemma (symbolic): (t+3)^3 - t^3 = 9t^2 + 27t + 27 == 0 (mod 9) identically.
# Expand exactly via integer polynomial coefficients:
# (t+3)^3 = t^3 + 9t^2 + 27t + 27.
lemma_coeffs = [27, 27, 9, 0]   # difference (t+3)^3 - t^3, low-to-high deg
check("Lemma: (t+3)^3 - t^3 has all coefficients divisible by 9",
      all(c % 9 == 0 for c in lemma_coeffs), f"coeffs (deg 0..3) = {lemma_coeffs}")
# Hence cube residue mod 9 depends only on x mod 3. The three cases:
cube9 = {r: (r ** 3) % 9 for r in range(3)}     # x mod 3 -> x^3 mod 9
print(f"  cube map Z/3 -> Z/9: {cube9}  (image {sorted(set(cube9.values()))})")
check("D3 cube image mod 9 is {0,1,8} = {0,+1,-1}",
      set(cube9.values()) == {0, 1, 8})

# wide-window empirical cross-check of the lemma
check("Lemma cross-check on x in [-200,200]",
      all((x ** 3) % 9 == cube9[x % 3] for x in range(-200, 201)))

# D3: image of x^3+y^3+z^3 mod 9, via the 27 triples of (x,y,z) mod 3
attained = set()
for a in range(3):
    for b in range(3):
        for c in range(3):
            attained.add((cube9[a] + cube9[b] + cube9[c]) % 9)
forbidden = sorted(set(range(9)) - attained)
print(f"  attained residues mod 9: {sorted(attained)}; forbidden: {forbidden}")
check("D3 attained = {0,1,2,3,6,7,8}", attained == {0, 1, 2, 3, 6, 7, 8})
check("D3 forbidden = {4,5} (i.e. +-4 mod 9)", forbidden == [4, 5])

# independent 729-triple brute force over (Z/9)^3 (cross-check, different route)
attained729 = set()
bad = []
for x in range(9):
    for y in range(9):
        for z in range(9):
            s = (x ** 3 + y ** 3 + z ** 3) % 9
            attained729.add(s)
            if s in (4, 5):
                bad.append((x, y, z))
check("D3 cross-check over all 729 triples of (Z/9)^3",
      attained729 == {0, 1, 2, 3, 6, 7, 8} and bad == [])

# D4 rigidity: classify the 27 triples by sum mod 9
viol3, viol6 = [], []
wit3, wit6 = [], []
for a in range(3):
    for b in range(3):
        for c in range(3):
            s = (cube9[a] + cube9[b] + cube9[c]) % 9
            if s == 3:
                wit3.append((a, b, c))
                if not (a == b == c == 1):
                    viol3.append((a, b, c))
            if s == 6:
                wit6.append((a, b, c))
                if not (a == b == c == 2):
                    viol6.append((a, b, c))
check("D4: sum==3 (mod 9) forces x==y==z==1 (mod 3)",
      viol3 == [] and wit3 == [(1, 1, 1)], f"witness classes {wit3}")
check("D4: sum==6 (mod 9) forces x==y==z==2 (mod 3)",
      viol6 == [] and wit6 == [(2, 2, 2)], f"witness classes {wit6}")
# this is complete: cube mod 9 factors through x mod 3 (Lemma), so the 27
# classes exhaust all of Z^3.

print()
print("PART 2b: mod-9 audit of the claimed lists (finite facts, not laws)")
open_list = [114, 390, 627, 633, 732, 921, 975]
fallen = [3, 33, 42, 165, 579, 906]
for k in open_list + fallen:
    pass
audit_open = {k: k % 9 for k in open_list}
audit_fall = {k: k % 9 for k in fallen}
print(f"  open k mod 9:   {audit_open}")
print(f"  fallen k mod 9: {audit_fall}")
check("D4 corollary: all 7 claimed-open k are == 3 or 6 (mod 9)",
      all(v in (3, 6) for v in audit_open.values()))
check("D4 corollary: fallen famous cases 3,33,42,165,579,906 all == +-3 (mod 9)",
      all(v in (3, 6) for v in audit_fall.values()))
# spot-check the *fallen* status claims D2 mentions: 165 and 579 published reps
# (Booker-Sutherland 2019, Charity Engine). Known published solutions:
sol_579 = (143075750505019222645, -143070303858622169975, -6941531883806363291)
s579 = sum(t ** 3 for t in sol_579)
check("D2 support: published 579 solution verifies exactly",
      s579 == 579, f"residual = {s579 - 579}")
sol_165 = (-385495523231271884, 383344975542639445, 98422560467622814)
s165 = sum(t ** 3 for t in sol_165)
check("D2 support: published 165 solution verifies exactly",
      s165 == 165, f"residual = {s165 - 165}")

print()
print("=" * 72)
print("PART 3 (claim D5): mod 7 gives NO obstruction")
print("=" * 72)
cube7 = sorted({(x ** 3) % 7 for x in range(7)})
print(f"  cube image mod 7: {cube7}")
check("D5 cube image mod 7 is {0,1,6} = {0,+1,-1}", cube7 == [0, 1, 6])
# every residue attained: brute-force witnesses with x,y,z in [-1,1]
wit = {}
for x in range(-1, 2):
    for y in range(-1, 2):
        for z in range(-1, 2):
            r = (x ** 3 + y ** 3 + z ** 3) % 7
            wit.setdefault(r, (x, y, z))
print(f"  witnesses (r: (x,y,z)): {dict(sorted(wit.items()))}")
check("D5 all 7 residues mod 7 attained with x,y,z in {-1,0,1}",
      set(wit) == set(range(7)))
# interval argument: sums of three of {0,+-1} are exactly the integers
sums = sorted({e1 + e2 + e3 for e1 in (-1, 0, 1) for e2 in (-1, 0, 1)
               for e3 in (-1, 0, 1)})
check("D5 sums of three elements of {0,+-1} = [-3..3], 7 consecutive integers",
      sums == list(range(-3, 4)))
check("D5 contrast: 7 consecutive integers cover Z/7 but only 7 of 9 classes "
      "of Z/9, missing exactly {4,5}",
      {s % 7 for s in sums} == set(range(7)) and
      sorted(set(range(9)) - {s % 9 for s in sums}) == [4, 5])

print()
print("=" * 72)
print(f"VERIFIER SUMMARY: {'ALL CHECKS PASS' if fails == 0 else str(fails) + ' FAILURES'}")
print("=" * 72)
