"""AUDIT D1 (THM-5055): arithmetic check of the case analyses, for every order up to NMAX (default 400), using
interval reasoning for the generic cases and exhaustive enumeration of the O(n) exceptional values.
  Lemma 8 (J = L(j) >= 13): cases A1/A2/A3/B1/B3 are exhaustive, their bounds hold (1 <= i <= J in A3,
    i <= J or i = J + 3 in B3, i odd), the inductive calls are legal (J - 1 >= L(i) + 8 in A3, J >= L(i) + 8 in B3,
    Lemma 7 order J + 1 >= 8 in A2), and the start orders J+4, J+2, J+1, J+6, J+1 give [J + 8, oo) via Lemma 6.
  Section 5 (n >= 13): the four ranges partition the odd k <= M(n); each exceptional k satisfies the stated bounds
    (j <= n + 1 with n >= L(j) + 8; j' <= n with n - 1 >= L(j') + 8; c in {1,3,5,7,11,13} with order >= 8) and
    X(n) is the only excluded value of range 4.
  Well-foundedness: Lemma 8 for j <= n + 1 uses the theorem only at orders <= L(n + 1) - 1 < n.
Also the four inequalities used, numerically up to 10^6."""
import sys
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 400
def M(n): return 2**n - n if n % 2 else 2**n - n - 1
def Xv(n): return 2**n - n - 2 if (n % 2 == 1 and n >= 5) else None
def inF(n, k):
    if k < 1 or k % 2 == 0 or k > M(n) or k in (3, 7, 13): return False
    return k != Xv(n)
def L(k): return k.bit_length()
bad = 0
def fail(*a):
    global bad; bad += 1; print("FAIL", *a)
# inequalities
for J in range(13, 10**6):
    if not (J - 1 >= L(J) + 8 and J >= L(J + 3) + 8): fail("ineq J", J); break
for n in range(12, 10**6):
    if not n >= L(n + 1) + 8: fail("ineq n>=L(n+1)+8", n); break
for n in range(13, 10**6):
    if not n - 1 >= L(n) + 8: fail("ineq n-1>=L(n)+8", n); break
# Lemma 6 reach: from a start order s, orders s, s+2, s+3, ... ; need [J+8, oo)
def reach_ok(s, J): return s + 2 <= J + 8 or s == J + 8
cases = {}
for J in range(13, NMAX + 1):
    # Case A: x = 2^J - j odd in [1, 2^(J-2) - 1].  Generic: x in F(J-2) (A1).
    excA = sorted(set([x for x in (3, 7, 13)] + list(range(M(J - 2) + 2, 2**(J - 2), 2)) + ([Xv(J - 2)] if Xv(J - 2) else [])))
    for x in excA:
        assert x % 2 == 1 and 1 <= x < 2**(J - 2)
        j = 2**J - x
        if inF(J - 2, x): fail("A: exceptional x in F", J, x)
        elif x in (3, 7, 13):
            cases['A2'] = cases.get('A2', 0) + 1
            if not (J + 1 >= 8 and reach_ok(J + 2, J)): fail("A2", J, x)
        else:
            cases['A3'] = cases.get('A3', 0) + 1
            i = 2**(J - 2) - x
            if not (i % 2 == 1 and 1 <= i <= J and j == 3 * 2**(J - 2) + i and J - 1 >= L(i) + 8 and i < j and reach_ok(J + 1, J)):
                fail("A3", J, x, i)
    # every odd x < 2^(J-2) not in excA must be in F(J-2): x <= M(J-2), x != X(J-2), x not in {3,7,13}
    if not (M(J - 2) < 2**(J - 2)): fail("M", J)
    cases['A1'] = cases.get('A1', 0) + (2**(J - 3) - len(excA))
    if not reach_ok(J + 4, J): fail("A1 reach", J)
    # Case B: z = y + 2 = 2^J - j + 2, odd, in (2^(J-2) + 2, 2^(J-1) + 2)
    lo, hi = 2**(J - 2) + 3, 2**(J - 1) + 1
    excB = sorted(set(z for z in list(range(M(J - 1) + 2, hi + 1, 2)) + ([Xv(J - 1)] if Xv(J - 1) else []) if lo <= z <= hi))
    for z in excB:
        j = 2**J - z + 2
        if not (2**(J - 1) < j < 3 * 2**(J - 2)): fail("B range", J, z, j)
        if inF(J - 1, z): fail("B: exceptional z in F", J, z)
        if z <= 13: fail("B: z <= 13", J, z)
        i = j - 2**(J - 1)
        ok = i % 2 == 1 and (i <= J or (z == Xv(J - 1) and i == J + 3)) and J >= L(i) + 8 and i < j and reach_ok(J + 1, J)
        cases['B3'] = cases.get('B3', 0) + 1
        if not ok: fail("B3", J, z, i)
    # generic z in [lo, hi] odd, not in excB, must be in F(J-1): need lo > 13 and z <= M(J-1) and z != X(J-1)
    nB = (hi - lo) // 2 + 1
    cases['B1'] = cases.get('B1', 0) + (nB - len(excB))
    if not (lo > 13 and reach_ok(J + 6, J)): fail("B1", J)
    # total count per J: 2^(J-3) values in case A (x odd < 2^(J-2)) and nB in case B; together 2^(J-2)
    if 2**(J - 3) + nB != 2**(J - 2): fail("count", J)
print("Lemma 8 cases for 13 <= J <=", NMAX, ":", {k: (v if v < 10**9 else '%.3g' % v) for k, v in sorted(cases.items())})
# Section 5
for n in range(13, NMAX + 1):
    # ranges as intervals of odd k: R1 = [1, M(n-1)] minus X(n-1); R2 = (M(n-1), 2^(n-1)) u {X(n-1)};
    # R3 = (2^(n-1), 3*2^(n-2) - 3]; R4 = [3*2^(n-2) - 1, M(n)]
    if not (M(n - 1) < 2**(n - 1) and 2**(n - 1) + 1 <= 3 * 2**(n - 2) - 3 and 3 * 2**(n - 2) - 1 <= M(n)): fail("ranges", n)
    # range 2 exceptional values
    for k in list(range(M(n - 1) + 2, 2**(n - 1), 2)) + ([Xv(n - 1)] if Xv(n - 1) else []):
        j = 2**(n - 1) - k
        if not (j % 2 == 1 and (j <= n - 1 or (k == Xv(n - 1) and j == n + 1)) and n >= L(j) + 8): fail("R2", n, k)
        if not (L(j) - 1 < n): fail("wellfounded", n, j)
    # range 3: x = k - 2^(n-1), 1 <= x <= 2^(n-2) - 3 odd
    for x in list(range(M(n - 2) + 2, 2**(n - 2) - 2, 2)) + ([Xv(n - 2)] if Xv(n - 2) else []):
        jp = 2**(n - 2) - x
        if not (jp % 2 == 1 and jp <= n and n - 1 >= L(jp) + 8): fail("R3b", n, x)
    for x in (3, 7, 13):
        if not (n >= 8 and x <= 2**(n - 2) - 3): fail("R3c", n, x)
    # range 4: y = k - 3*2^(n-2) + 2 in [1, M(n-2)]
    if M(n) - 3 * 2**(n - 2) + 2 != M(n - 2): fail("R4 top", n)
    for y in (3, 7, 13):
        k = 3 * 2**(n - 2) + y - 2
        if not (y - 2 in (1, 5, 11) and n - 1 >= 8 and k <= M(n) and inF(n, k)): fail("R4b", n, y)
    if Xv(n - 2) is not None:
        k = 3 * 2**(n - 2) - 2 + Xv(n - 2)
        if k != Xv(n): fail("R4c", n)
    elif Xv(n) is not None: fail("X(n) without X(n-2)", n)
print("Section 5 ranges for 13 <= n <=", NMAX, "checked")
print("failures:", bad)
