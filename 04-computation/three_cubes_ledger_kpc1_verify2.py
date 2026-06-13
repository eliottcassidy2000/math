# =====================================================================
# ADVERSARIAL VERIFIER (independent pass 2) -- session kind-pasteur-2026-06-10-S1
# Independent check of thread D ("three-cubes ledger") claims D1, D3, D4, D5.
# FRESH code; methods deliberately different from the worker:
#   * D1: cubes computed by repeated multiplication x*x*x (no pow()), plus
#     an independent modular fingerprint mod three large primes via pow(,,p).
#   * D3/D5: residues harvested from ACTUAL integers in [-40, 40] (not from
#     residue classes 0..m-1), then residue-class exhaustion over the
#     fundamental domain {-4..4} (mod 9) / {-3..3} (mod 7) as completeness
#     cross-check, plus the shift lemma making finite exhaustion a proof.
#   * D4: full 729-triple exhaustion over (Z/9)^3 via {-4..4}^3.
# Pure python ints throughout. No floats anywhere.
# =====================================================================
import sys

fails = 0
def check(name, cond, detail=""):
    global fails
    line = ("PASS " if cond else "FAIL ") + name
    if detail:
        line += "  |  " + detail
    print(line)
    if not cond:
        fails += 1

print("=" * 70)
print("PART 1: bigint identities (D1)")
print("=" * 70)
cases = [
    (42, -80538738812075974, 80435758145817515, 12602123297335631,
     "Booker-Sutherland 2019"),
    (33, 8866128975287528, -8778405442862239, -2736111468807040,
     "Booker 2019"),
    (3, 569936821221962380720, -569936821113563493509, -472715493453327032,
     "Booker-Sutherland 2019 (3rd nontrivial rep of 3)"),
]
for k, x, y, z, who in cases:
    s = x*x*x + y*y*y + z*z*z          # repeated multiplication, no pow()
    check("D1: x^3+y^3+z^3 == %d (%s)" % (k, who), s == k,
          "sum=%d residual=%d" % (s, s - k))
    # independent fingerprint route: modular exponentiation, three primes
    for p in (10**9 + 7, 2**61 - 1, 999999937):
        sp = (pow(x % p, 3, p) + pow(y % p, 3, p) + pow(z % p, 3, p)) % p
        check("D1: k=%d fingerprint mod %d" % (k, p), sp == k % p)

print()
print("=" * 70)
print("PART 2: mod-9 forbidden values (D3)")
print("=" * 70)
# (a) cube residues mod 9 from actual integers -40..40
cubres9 = sorted({(t*t*t) % 9 for t in range(-40, 41)})
check("D3: cube residues mod 9 == {0,1,8}", cubres9 == [0, 1, 8],
      "got %s" % cubres9)
# (b) completeness lemma: cube mod 9 depends only on x mod 3
#     (x+3)^3 - x^3 = 9x^2 + 27x + 27 == 0 (mod 9) -- verify on a full period
lemma = all(((x + 3)**3 - x**3) % 9 == 0 for x in range(-9, 10))
check("D3: (x+3)^3 == x^3 (mod 9) on a full period (completeness lemma)", lemma)
# (c) attained sums: triple loop over a fundamental domain {-4..4}^3 of (Z/9)^3
attained9 = set()
for x in range(-4, 5):
    for y in range(-4, 5):
        for z in range(-4, 5):
            attained9.add((x*x*x + y*y*y + z*z*z) % 9)
check("D3: attained sums mod 9 == {0,1,2,3,6,7,8}",
      sorted(attained9) == [0, 1, 2, 3, 6, 7, 8], "got %s" % sorted(attained9))
check("D3: forbidden set == {4,5} (i.e. +-4 mod 9)",
      sorted(set(range(9)) - attained9) == [4, 5])

print()
print("=" * 70)
print("PART 3: +-3 mod 9 rigidity (D4)")
print("=" * 70)
viol3 = viol6 = n3 = n6 = 0
for x in range(-4, 5):
    for y in range(-4, 5):
        for z in range(-4, 5):
            s = (x*x*x + y*y*y + z*z*z) % 9
            if s == 3:
                n3 += 1
                if not (x % 3 == y % 3 == z % 3 == 1):
                    viol3 += 1
            elif s == 6:
                n6 += 1
                if not (x % 3 == y % 3 == z % 3 == 2):
                    viol6 += 1
check("D4: s==3 (mod 9) => x==y==z==1 (mod 3), all 729 triples",
      viol3 == 0, "triples hitting 3: %d, violations: %d" % (n3, viol3))
check("D4: s==6 (mod 9) => x==y==z==2 (mod 3), all 729 triples",
      viol6 == 0, "triples hitting 6: %d, violations: %d" % (n6, viol6))
check("D4: non-vacuous (some triples do hit 3 and 6)", n3 > 0 and n6 > 0,
      "n3=%d n6=%d" % (n3, n6))

# finite-list audit: open k < 1000 and famous fallen cases, all == +-3 mod 9
open_list = [114, 390, 627, 633, 732, 921, 975]
fallen = [3, 33, 42, 165, 579, 906]
for k in open_list + fallen:
    check("D4: k=%d == +-3 (mod 9)" % k, k % 9 in (3, 6), "k mod 9 = %d" % (k % 9))

print()
print("=" * 70)
print("PART 4: mod 7 gives NO obstruction (D5)")
print("=" * 70)
cubres7 = sorted({(t*t*t) % 7 for t in range(-40, 41)})
check("D5: cube residues mod 7 == {0,1,6}", cubres7 == [0, 1, 6],
      "got %s" % cubres7)
attained7 = set()
for x in range(-3, 4):
    for y in range(-3, 4):
        for z in range(-3, 4):
            attained7.add((x*x*x + y*y*y + z*z*z) % 7)
check("D5: attained sums mod 7 == ALL of Z/7", attained7 == set(range(7)),
      "got %s" % sorted(attained7))
# interval argument made precise: sums of three values in {-1,0,1} give
# exactly {-3..3}: 7 consecutive integers, pairwise distinct mod 7.
ivals = sorted({a + b + c for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)})
check("D5: {a+b+c : a,b,c in {0,+-1}} == {-3..3}", ivals == list(range(-3, 4)))
check("D5: {-3..3} pairwise distinct mod 7 (7 values, 7 classes)",
      len({v % 7 for v in ivals}) == 7)
# contrast (the D3 mechanism): {-3..3} mod 9 misses exactly {4,5}
check("D3/D5 contrast: {-3..3} mod 9 misses exactly {4,5}",
      sorted(set(range(9)) - {v % 9 for v in ivals}) == [4, 5])

print()
print("=" * 70)
print("VERIFIER SUMMARY: %s (%d failures)" % ("ALL CHECKS PASS" if fails == 0
      else "FAILURES PRESENT", fails))
print("=" * 70)
sys.exit(0 if fails == 0 else 1)
