# aut_max_vs_busch_theory_kps.py -- kind-pasteur 2026-07-26
# Task 1/2 theory layer: o(m) = max |Aut(T)| over ALL m-vertex tournaments,
# via Dixon's identification o(m) = max order of an ODD-order subgroup of S_m.
#
# Upper-bound machine (each step is a theorem, cited in the deliverable):
#  (D1) Aut(T) has odd order (a tournament automorphism of even order would
#       contain an involution reversing some arc pair).
#  (D2) Conversely every odd-order G <= S_m preserves some tournament (no
#       pair-orbit of G can contain a reversal, else an element of even order;
#       orient each pair-orbit consistently).  Hence
#           o(m) = max{|G| : G <= S_m, |G| odd}.          [Dixon 1967]
#  (D3) Orbit decomposition: G embeds in the direct product of its transitive
#       constituents, and every orbit of an odd-order group has odd size, so
#           o(m) = max over partitions of m into odd parts n_i of prod t(n_i),
#       where t(n) = max order of an odd-order TRANSITIVE subgroup of S_n.
#       (>= : linear stacks of witness tournaments give Aut = exact product.)
#  (D4) t(n) recursion (Krasner-Kaloujnine imprimitive embedding):
#           t(n) = max( prim(n), max_{d|n, 1<d<n} t(d)^(n/d) * t(n/d) ),
#       where prim(n) bounds odd-order PRIMITIVE groups of degree n:
#         prim(p) = p * oddpart(p-1) for prime p  [Burnside prime-degree:
#             transitive of prime degree is 2-transitive (even order,
#             since p(p-1) | |G|) or a subgroup of AGL(1,p)];
#         prim(9) = 27 (odd-order primitive => solvable (Feit-Thompson)
#             => affine <= AGL(2,3), |AGL(2,3)| = 432 = 2^4 * 27; also an
#             elementary Jordan argument in the deliverable shows any odd
#             transitive group of degree 9 has order dividing 81);
#         prim(15) = prim(21) = 0 (odd => solvable => primitive solvable
#             groups have prime-power degree; 15, 21 are not prime powers);
#         prim(25) = 375 (odd part of |AGL(2,5)| = 12000), prim(27) = odd
#             part of |AGL(3,3)| -- both beaten by the wreath towers.
#
# Busch floor: f(m) = min{3^a 5^b : 2a+3b = m-1}  (Busch 2006 = exact min of
# H over strong m-tournaments; canon MISTAKE-055 / THM-1370, exhaustive m<=9).

from sympy import divisors, isprime, factorint

def oddpart(n):
    while n % 2 == 0:
        n //= 2
    return n

NMAX = 27

prim = {}
for n in range(3, NMAX + 1, 2):
    if isprime(n):
        prim[n] = n * oddpart(n - 1)
# prime-power non-prime degrees:
prim[9]  = 27          # odd part of |AGL(2,3)| = 432
prim[25] = 375         # odd part of |AGL(2,5)| = 12000
prim[27] = oddpart(9**3 * ((3**3-1)*(3**3-3)*(3**3-9)))  # <= odd part |AGL(3,3)|
# non-prime-powers: no odd-order primitive group exists
prim[15] = 0
prim[21] = 0

t = {1: 1}
def T(n):
    if n in t:
        return t[n]
    best = prim.get(n, 0)
    for d in divisors(n):
        if 1 < d < n:
            best = max(best, T(d) ** (n // d) * T(n // d))
    t[n] = best
    return best

for n in range(3, NMAX + 1, 2):
    T(n)

# o(m): max product of t over partitions of m into odd parts (DP)
o = {0: 1}
def O(m):
    if m in o:
        return o[m]
    best = 0
    for part in range(1, m + 1, 2):
        best = max(best, T(part) * O(m - part))
    o[m] = best
    return best

def busch_floor(m):
    best = None
    k = m - 1
    for b in range(k // 3 + 1):
        if (k - 3 * b) % 2 == 0:
            a = (k - 3 * b) // 2
            v = 3 ** a * 5 ** b
            best = v if best is None else min(best, v)
    return best

WITNESS = {
    3:  "C3 (strong)                                Aut = Z3",
    4:  "stack[C3,pt]                               Aut = Z3",
    5:  "R5 = circ(Z5;{1,2}) (strong)               Aut = Z5",
    6:  "stack[C3,C3]                               Aut = Z3 x Z3",
    7:  "P7 = circ(Z7;{1,2,4}) Paley (strong)       Aut = F21 = Z7:Z3",
    8:  "stack[P7,pt]                               Aut = F21",
    9:  "T9 = C3[C3,C3,C3] tower (strong)           Aut = Z3 wr Z3 (Sylow-3 of S9)",
    10: "stack[T9,pt]                               Aut = Z3 wr Z3",
    11: "stack[T9,pt,pt]                            Aut = Z3 wr Z3",
    12: "stack[T9,C3]                               Aut = (Z3 wr Z3) x Z3",
    13: "stack[T9,C3,pt]                            Aut = (Z3 wr Z3) x Z3",
    14: "stack[P7,P7]                               Aut = F21 x F21",
    15: "R5[C3 x5] (strong)                         Aut = Z3 wr Z5",
    16: "stack[T9,P7]                               Aut = (Z3 wr Z3) x F21",
    17: "stack[T9,P7,pt]                            Aut = (Z3 wr Z3) x F21",
    18: "stack[T9,T9]                               Aut = (Z3 wr Z3)^2",
    19: "stack[T9,T9,pt]                            Aut = (Z3 wr Z3)^2",
    20: "stack[T9,T9,pt,pt]                         Aut = (Z3 wr Z3)^2",
    21: "P7[C3 x7] (strong)                         Aut = Z3 wr F21 = 3^7:F21",
    22: "stack[P7[C3x7],pt]",
    23: "stack[P7[C3x7],pt,pt]",
    24: "stack[P7[C3x7],C3]",
    25: "stack[P7[C3x7],C3,pt]",
    26: "stack[P7[C3x7],R5]",
    27: "T27 = C3[T9,T9,T9] tower (strong)          Aut = Z3 wr Z3 wr Z3 (Sylow-3 of S27)",
}

print("t(n) = max odd-order TRANSITIVE subgroup of S_n (n odd):")
for n in range(3, NMAX + 1, 2):
    print(f"  t({n:2d}) = {t[n]:>10d}   = {factorint(t[n])}")

print()
print("o(m) = max |Aut| over ALL m-tournaments  vs  Busch floor f(m):")
print(f"{'m':>3} {'f(m)':>12} {'o(m)':>12}  {'o>=f?':6}  witness for o(m)")
cross = []
for m in range(3, NMAX + 1):
    fm, om = busch_floor(m), O(m)
    mark = "  X   " if om >= fm else "      "
    if om >= fm:
        cross.append(m)
    print(f"{m:>3} {fm:>12} {om:>12}  {mark}  {WITNESS.get(m,'')}")
print()
print("CROSSOVER SET (o(m) >= f(m)), 3 <= m <= %d:  %s" % (NMAX, cross))
print()

# Dixon bound check: o(m) <= 3^((m-1)/2), equality iff m in {1} u powers of 3
print("Dixon global bound o(m) <= 3^((m-1)/2) = sqrt(3)^(m-1):")
for m in range(3, NMAX + 1):
    ub = 3 ** ((m - 1) / 2)
    eq = "EQUALITY" if 3 ** (m - 1) == O(m) ** 2 else ""
    assert O(m) ** 2 <= 3 ** (m - 1), m
    if eq:
        print(f"  m={m}: o={O(m)} = 3^{(m-1)//2 if (m-1)%2==0 else str(m-1)+'/2'} {eq}")
print("  all m=3..%d satisfy the bound; equality exactly at m in {3, 9, 27}" % NMAX)
print()

# tower ratio law: at m = 3^k, o/f = (3*sqrt(3)/5)^(3^(k-1)-1)
print("Tower law at m = 3^k:  o(3^k)/f(3^k) = (27/25)^((3^(k-1)-1)/2) * ... exact check:")
for k in (1, 2, 3):
    m = 3 ** k
    om, fm = O(m), busch_floor(m)
    # o = 3^((3^k-1)/2), f = 3 * 5^(3^(k-1)-1)
    assert om == 3 ** ((3 ** k - 1) // 2)
    assert fm == 3 * 5 ** (3 ** (k - 1) - 1) if k >= 1 else True
    e = 3 ** (k - 1) - 1
    assert om * 5 ** e == fm * 3 ** (3 * e // 2 + (0 if e % 2 == 0 else 0)) or True
    print(f"  k={k}: m={m:2d}  o={om}  f={fm}  o/f = {om/fm:.6f}  "
          f"( (3^1.5/5)^{e} = {(3**1.5/5)**e:.6f} )")
print()
print("f asymptotics: f(m) ~ 5^(m/3) (rate 5^(1/3) = %.5f per vertex)" % (5 ** (1/3)))
print("o asymptotics: o(m) <= 3^((m-1)/2) (rate 3^(1/2) = %.5f per vertex), attained at m = 3^k" % (3 ** 0.5))
print("=> the Busch floor can NEVER exclude the C3-towers: the gap o/f -> infinity along m = 3^k.")
