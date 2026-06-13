#!/usr/bin/env python3
# three_cubes_ledger_kpc1.py -- session kind-pasteur-2026-06-10-S1, Thread D (HYP-2370)
# The three-cubes ledger: bigint verification of the 2019 Booker / Booker-Sutherland
# solutions, the mod-9 forbidden-values law (PROVED by finite exhaustion), the mod-7
# non-obstruction (PROVED), the +-3 mod 9 rigidity law behind the stubborn cases
# (PROVED by finite exhaustion), and the 2026 open-list audit.
# Pure python ints throughout. No floats. No external libraries.

import sys

PASS = "PASS"
FAIL = "FAIL"
failures = 0

def check(label, lhs, rhs):
    global failures
    ok = (lhs == rhs)
    if not ok:
        failures += 1
    if isinstance(lhs, int) and isinstance(rhs, int):
        print(f"  [{PASS if ok else FAIL}] {label}: residual = {lhs - rhs}")
    else:
        print(f"  [{PASS if ok else FAIL}] {label}: got {lhs}, expected {rhs}")
    return ok

print("=" * 78)
print("PART 1 -- BIGINT VERIFICATION OF THE 2019 SOLUTIONS (exact integer arithmetic)")
print("=" * 78)

# 42 = (-80538738812075974)^3 + 80435758145817515^3 + 12602123297335631^3
#   Booker & Sutherland, September 2019, Charity Engine (~1.3M core-hours).
x, y, z = -80538738812075974, 80435758145817515, 12602123297335631
s = x**3 + y**3 + z**3
print(f"42: x={x}, y={y}, z={z}")
print(f"    x^3+y^3+z^3 = {s}")
check("42 = x^3+y^3+z^3 (Booker-Sutherland 2019)", s, 42)

# 33 = 8866128975287528^3 + (-8778405442862239)^3 + (-2736111468807040)^3
#   Booker, March 2019 ("Cracking the problem with 33", Res. Number Theory 5:26).
x, y, z = 8866128975287528, -8778405442862239, -2736111468807040
s = x**3 + y**3 + z**3
print(f"33: x={x}, y={y}, z={z}")
print(f"    x^3+y^3+z^3 = {s}")
check("33 = x^3+y^3+z^3 (Booker 2019)", s, 33)

# 3 = 569936821221962380720^3 + (-569936821113563493509)^3 + (-472715493453327032)^3
#   Booker & Sutherland, September 2019 -- the third nontrivial representation of 3,
#   answering Mordell's 1953 question (first new rep of 3 beyond (1,1,1), (4,4,-5)).
x, y, z = 569936821221962380720, -569936821113563493509, -472715493453327032
s = x**3 + y**3 + z**3
print(f"3:  x={x}, y={y}, z={z}")
print(f"    x^3+y^3+z^3 = {s}")
check("3 = x^3+y^3+z^3 (Booker-Sutherland 2019, third rep)", s, 3)

print()
print("=" * 78)
print("PART 3a -- MOD 9 FORBIDDEN-VALUES LAW (PROVED by finite exhaustion)")
print("=" * 78)
# x^3 mod 9 depends only on x mod 9 (indeed only on x mod 3: (x+3k)^3 = x^3 + 9(...)),
# so exhausting x in 0..8 is a complete proof, not a sample.
cubes_mod9 = sorted({(x * x * x) % 9 for x in range(9)})
print(f"cube residues mod 9 = {cubes_mod9}  (i.e. {{0, +1, -1}} mod 9)")
assert cubes_mod9 == [0, 1, 8]

sums_mod9 = sorted({(a + b + c) % 9 for a in cubes_mod9 for b in cubes_mod9 for c in cubes_mod9})
missing = sorted(set(range(9)) - set(sums_mod9))
print(f"attainable residues of x^3+y^3+z^3 mod 9 = {sums_mod9}")
print(f"FORBIDDEN residues mod 9 = {missing}  (i.e. +-4 mod 9)")
check("forbidden set mod 9 == {4,5}", missing, [4, 5])
print("PROVED: x^3+y^3+z^3 is never == 4 or 5 (mod 9).  Proof: cube map factors")
print("through Z/9 (in fact (x+3)^3 == x^3 mod 9), so the 9-case + 27-triple")
print("exhaustion above is a complete case analysis, not a sampled pattern.")

print()
print("=" * 78)
print("PART 3b -- +-3 MOD 9 RIGIDITY LAW (why the stubborn k all look alike)")
print("=" * 78)
# Claim: if x^3+y^3+z^3 == 3 (mod 9) then x == y == z == 1 (mod 3);
#        if x^3+y^3+z^3 == 6 (mod 9) then x == y == z == 2 (mod 3).
# Exhaust all 729 triples (x,y,z) in (Z/9)^3 -- complete proof for the same reason.
ok3 = True
ok6 = True
for x in range(9):
    for y in range(9):
        for z in range(9):
            r = (x**3 + y**3 + z**3) % 9
            if r == 3 and not (x % 3 == y % 3 == z % 3 == 1):
                ok3 = False
            if r == 6 and not (x % 3 == y % 3 == z % 3 == 2):
                ok6 = False
print(f"  [{'PASS' if ok3 else 'FAIL'}] sum == 3 (mod 9)  ==>  x == y == z == 1 (mod 3)   [729/729 triples]")
print(f"  [{'PASS' if ok6 else 'FAIL'}] sum == 6 (mod 9)  ==>  x == y == z == 2 (mod 3)   [729/729 triples]")
if not (ok3 and ok6):
    failures += 1
print("PROVED: for k == +-3 (mod 9) every solution is 'rigid': all of x,y,z lie in")
print("a single residue class mod 3.  This kills ~26/27 of the naive search space")
print("and is the standard explanation of why the hard cases concentrate there.")

print()
print("=" * 78)
print("PART 3c -- MOD 7 GIVES NO OBSTRUCTION (PROVED)")
print("=" * 78)
cubes_mod7 = sorted({(x * x * x) % 7 for x in range(7)})
print(f"cube residues mod 7 = {cubes_mod7}  (i.e. {{0, +1, -1}} mod 7)")
assert cubes_mod7 == [0, 1, 6]
sums_mod7 = sorted({(a + b + c) % 7 for a in cubes_mod7 for b in cubes_mod7 for c in cubes_mod7})
print(f"attainable residues of x^3+y^3+z^3 mod 7 = {sums_mod7}")
check("all 7 residues attained mod 7", sums_mod7, list(range(7)))
print("PROVED, with the clean reason: cube residues mod 7 are {0,+1,-1}, so sums of")
print("three of them realize every integer in [-3,3] -- that is 7 consecutive values,")
print("hence ALL of Z/7.  (Contrast mod 9: [-3,3] is only 7 of the 9 classes, and the")
print("two missed classes are exactly +-4.)  Same {0,+-1} cube image, opposite verdict:")
print("the modulus size decides.  Mod 7 imposes nothing.")

print()
print("=" * 78)
print("PART 2 -- THE 2026 OPEN LEDGER FOR x^3+y^3+z^3 = k, k < 1000  [WEB-VERIFIED]")
print("=" * 78)
open_list = [114, 390, 627, 633, 732, 921, 975]
print(f"OPEN (no representation known), k < 1000: {open_list}")
print("smallest open k = 114")
print("Source (accessed 2026-06-10): Wikipedia 'Sums of three cubes', verbatim:")
print("  'The only remaining unsolved cases up to 1,000 are the seven numbers")
print("   114, 390, 627, 633, 732, 921, and 975, and there are no known primitive")
print("   solutions (i.e. gcd(x,y,z)=1) for 192, 375, and 600.'")
print("NOTE: 579 is NOT open -- it fell in September 2019 (Booker-Sutherland, same")
print("Charity Engine run as 42; Wikipedia verbatim: '...as well as solutions for")
print("several other previously unknown cases including n=165 and 579 for n<=1000.')")
print("The human dispatch's candidate list included 579: corrected here.")
print()
print("mod-9 audit of the open list (and the famous fallen cases):")
all_pm3 = True
for k in open_list:
    r = k % 9
    tag = "+3" if r == 3 else ("-3 (==6)" if r == 6 else str(r))
    if r not in (3, 6):
        all_pm3 = False
    print(f"  k = {k:4d}: k mod 9 = {r}  [{tag}]")
check("every open k < 1000 is == +-3 (mod 9)", all_pm3, True)
for k in (3, 33, 42, 165, 579, 906):
    print(f"  fallen k = {k:4d}: k mod 9 = {k % 9}")
print("VERIFIED(list of 7): all seven open values are == +-3 (mod 9) -- exactly the")
print("rigid class of PART 3b.  The famous fallen cases 3, 33, 42, 165, 579, 906 were")
print("also all == +-3 (mod 9).  This is a structural observation about a finite list,")
print("NOT a law about which k are hard (MISTAKE-028/036/055 discipline).")

print()
print("=" * 78)
print("HEATH-BROWN COMPLETENESS -- status: CONJECTURE")
print("=" * 78)
print("CONJECTURE (Heath-Brown 1992, 'The density of zeros of forms for which weak")
print("approximation fails', Math. Comp. 59):  every k != +-4 (mod 9) has a")
print("representation x^3+y^3+z^3 = k, and in fact infinitely many.")
print("Status 2026: open.  No representation is known for k = 114; no algorithm is")
print("even known to decide representability for a single given k.")
print("Mod 9 is the only known congruence obstruction: the form x^3+y^3+z^3 - k has")
print("zeros in every Z_p (p != 3) and in Z_3 iff k != +-4 (mod 9); failure for")
print("k == +-4 (mod 9) is a pure 3-adic (mod 9) phenomenon, cf. PART 3a.")

print()
print("=" * 78)
if failures == 0:
    print("ALL CHECKS PASS (0 failures)")
else:
    print(f"*** {failures} CHECK(S) FAILED ***")
sys.exit(0 if failures == 0 else 1)
