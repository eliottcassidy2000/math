# A projective prime-divisor atlas, with twists and valuations retained

2026-10-28. **PROVED elementary arithmetic**, **CITED** density theorem,
and **FINITE-EXACT** controls through prime5000. No Lucas-prime frequency
prediction follows from the finite census. Same-agent audit.

## 1. The uniform family containing u -> 6/u

For an odd prime r>3, let J_6 act on P1(F_r) by u -> 6/u, exchanging
zero and infinity. Its matrix is [[0,6],[1,0]], its determinant is -6,
and its fixed points solve u^2=6. Consequently

    number of fixed points = 1+(6/r),
    J_6 belongs to PSL(2,r) iff (-6/r)=1.

The parentheses denote Legendre symbols. These are different tests. The
supplementary laws and quadratic reciprocity give the complete table:

| r mod24 | 1 | 5 | 7 | 11 | 13 | 17 | 19 | 23 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| (6/r) | + | + | - | - | - | - | + | + |
| (-6/r) | + | + | + | + | - | - | - | - |
| fixed points | 2 | 2 | 0 | 0 | 0 | 0 | 2 | 2 |

With T:u->u+1,

    <T,J_6> = PSL(2,r) if (-6/r)=1,
                PGL(2,r) otherwise.

Proof: powers of T give all upper unipotents, and conjugation by J_6 gives
all lower unipotents, with nonzero scaling of their parameter. These generate
SL2 by row elimination. Projectivizing yields PSL2; J_6 either lies there
or supplies the other determinant squareclass of PGL2. This proof uses
prime fields: powers of a single translation do not give every parameter
in an arbitrary extension field.

The same argument works with6 replaced by any nonzero c, retaining the
two characters(c/r) and(-c/r). At r=2,3 the displayed six-mirror is singular
and is excluded. This is an exact arithmetic stratification of the operation,
not a statement that primes are statistically independent of another sieve.

## 2. Lucas divisibility is a midpoint in a projective cyclic group

Use the recurrence matrix

    Q=[[1,1],[1,0]],    det Q=-1,
    Q^n=[[F_(n+1),F_n],[F_n,F_(n-1)]],    trace(Q^n)=L_n.

For each odd prime r!=5, let d_r be the order of [Q] in PGL(2,F_r).
The matrix formula proves directly

    d_r = min{n>=1 : r divides F_n}.

In odd characteristic, a nonscalar invertible two-by-two matrix represents
an involution iff its trace is zero. One direction is Cayley-Hamilton;
the other follows by subtracting the scalar-square equation from it.
Therefore

    r divides L_n
      iff d_r is even and n == d_r/2 (mod d_r).

This is a complete incidence description, with no sampled period cutoff.
Each prime either never occurs, or divides exactly one odd-multiple tower
of Lucas indices. It is an arithmetic analogue of retaining the actual
orbit/transport structure in a metagraph, rather than just the number of
states. The graph here is a cyclic group with its distinguished midpoint;
it is not a tournament.

To compute d_r efficiently, work in the quadratic algebra

    F_r[a]/(a^2-a-1),     b=1-a,     theta=a/b=-a^2.

Then d_r=ord(theta). If chi_r=(5/r), Frobenius either fixes both roots
or swaps them, giving

    d_r divides r-chi_r.

Factor r-chi_r and remove factors by powering theta. This retains the
split/inert channel rather than replacing all prime clocks by r-1.
`check.py` compares this computation with direct Fibonacci recurrences
for every odd prime r<=5000 except5.

## 3. The valuation layer is an exact strengthening of support

If d_r is even, put h=d_r/2 and e_r=v_r(L_h). Then, for every n>=1,

    v_r(L_n) = e_r+v_r(n/h)  if n/h is an odd positive integer,
               0             otherwise.

Proof: at a prime of Q(sqrt5) over r, the extension is unramified and a,b
are units. Since L_n=b^n(theta^n+1), reduce to the valuation of
theta^(h t)+1 for odd t. Set eta=-theta^h, so eta=1 mod r. The elementary
odd-prime lifting law gives v(eta^t-1)=v(eta-1)+v_r(t). To see that law,
write eta=1+r^e u; a power coprime to r keeps valuation e, while a power
r increases it by exactly one by the binomial expansion. Iterate.
The valuation of the rational integer L_n agrees with the normalized
local valuation. This proves the formula in both split and inert cases.

Thus the triple(d_r, e_r, chi_r) is a reusable prime record for exact Lucas
divisibility and valuation questions. It does not retain the entire residue
L_n mod r^j; those operations need the lifted matrix or quadratic-algebra
element as well. A prime with e_r>1 cannot be silently replaced by e_r=1.

For odd n, L_n^2+4=5F_n^2 shows that every odd prime divisor r of L_n
has(5/r)=1. For an odd prime index ell>=5, d_r must then be2ell:
the midpoint condition leaves d_r=2 or2ell, and F_2=1 excludes2. Hence

    r divides L_ell  implies  r == 1 (mod 2ell), (5/r)=1.

This elementary sieve constraint applies to all prime factors, including
factors of composite Lucas values. It does not say L_ell is prime.

## 4. An actual density theorem, with the sampling population specified

**CITED:** Moree, [On the prime density of Lucas sequences (1996),
Theorem2](https://www.numdam.org/item/JTNB_1996__8_2_449_0.pdf), gives density
2/3 among primes for the primes dividing at least one ordinary Lucas number,
crediting Lagarias for this case. By the midpoint criterion, this is the
density of primes with even d_r. It is **not** the density of indices n
for which L_n is prime. Moree's treatment uses splitting and multiplicative
order information in quadratic/Kummer fields, matching the need for a
growing order atlas rather than a fixed small-state analogy.

The historical March Fibonacci/Lucas-prime notes already noticed disjoint
factor supports and the Pell identity. Their current correction in
[abc_difference_families_20261011c.md](../abc_difference_families_20261011c.md)
explicitly prevents deriving probabilistic independence from coprimality.
The projective-order and valuation records above provide actual incidence
information that can replace that heuristic. Infinitely many prime-valued
Lucas terms are not established by this work.

## 5. The two involutions at eleven are not conjugate

The inherited eleven-word reflection becomes J_6 on the twelve matching
frames of the biplane. At r=11 it has **zero** projective fixed points.
The Lucas-divisibility half-turn [Q]^5 has trace L_5=11 and determinant-1;
its eigenlines are rational, and it has **two** projective fixed points.
Thus these two involutions are not conjugate in PGL(2,11). The original
eleven-word reflection itself fixes **three words**, in a different action.

This separates three carriers that a shared symbol "mirror" might hide:
eleven words, twelve matching frames, and the recurrence's projective line.
Their common framework is projective involutions and operation transport;
their fixed-point and determinant data determine when an identification is
actually possible. More generally a Lucas half-turn at index n has matrix
square (-1)^(n+1)I. Its twist squareclass is (-1)^(n+1), not automatically6.

## 6. Why the eleven-frame repair cannot be repeated verbatim

For every odd Lucas-prime degree p>=29, `REFLECTION_GROUPS.md` proves the
word action generated by reflection and arithmetic translation contains A_p.
Since

    |A_p|=p!/2 > p(p^2-1)=|PGL(2,p)|,

there is no faithful fractional-linear reencoding of that entire group on
P1(F_p). The eleven-point exceptional action admits exactly such a reencoding;
the obstruction beyond it is an order theorem, not unsuccessful searching.
An arbitrary action on p+1 points is a different question and is not excluded.

The robust replacement has two layers: retain the guarded word operations
for source-level questions, and use prime-specific projective clocks and
valuation lifts for divisibility. Neither layer is licensed to erase the
other's guards. Their common state is an atlas of compatible records,
not one universally sufficient small symmetry diagram.
