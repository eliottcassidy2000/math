# A primitive Fourier mode separates every twice-repeated growth chart

2026-10-04 (America/Denver).

**Status:** PROVED universal separation, compiler admissibility, count bounds,
and natural-density existence below; FINITE-EXACT for the declared period,
repetition, and replay universes. Global Collatz coverage remains OPEN.
No novelty claim is made for periodic words, affine Collatz composition,
Fourier characters, composition counting, or finite-word deviation bounds.

The finite 48-chart separation extends to **all distinct marked primitive
growth charts, provided both repetition counts are at least two**. A single
primitive cyclic Fourier character removes the remaining final-letter
ambiguity. There is also a universal atlas with a well-defined natural
density, rigorously bracketed by

    0.11862865 < delta_atlas < 0.11862913.

These decimals are rounded outward and checked against the exact rational
endpoints in the saved output. This is total atlas density among all positive integers,
**not** newly added coverage relative to the previous banks. Increasing the
period cutoff cannot make this particular atlas cover every source: the
convergent source 7 is outside it for a structural reason proved below.

## Inheritance and scope

[The rational-anchor compiler](rational_anchor_returns_20261004.md) supplies
the exact first-descent cylinders and actual valuations through the exit.
[The coverage refinement](anchor_coverage_refinement_20261004.md) compares
the first mismatch of two periodic nominal words and retains the possible
last-letter exception. That finite proof did not assert the general
twice-repeated theorem proved here.

[The two-copies investigation](fibonacci_two_copies_collatz_20260926.md),
section 4, already enumerates 130 primitive growth necklaces through odd
period eight and 5,966 through period twelve. Its repaired near miss matters:
filtering only the lexicographically first rotation had omitted favorable
phases and initially produced 120 rather than 130. We recover both corrected
counts and keep marked phases separate from necklaces.

Closest mechanism: the original-source first-descent clock. Canonical hostile:
the one-copy word `12` and the two-copy word `1` both give `3 mod16`.
The missing sidecar in an unqualified word comparison is its repetition
count and enlarged final valuation. The least-used representation here is
one **primitive cyclic Fourier mode**, used as an exact zero test.

The board is a marked valuation word, a repetition count, its Fourier
annihilator, a guarded source cylinder, a density limit, and a structural
missed source. The anchor is universal chart separation; the niche is the
density of the resulting infinite atlas; the wildcard is a faithful
sine/cosine interpretation. We prove the needed elementary periodicity
statement directly, without relying on a separately cited Fine–Wilf theorem.

## 1. One Fourier mode resolves the last-letter exception

Let w,v be nonempty words of positive integers with lengths p,q. Suppose
their valid source certificates at repetitions m,n have these two properties:

1. the exact first-descent times are pm and qn;
2. all valuations before the final step are the corresponding letters of
   `w^m` and `v^n`; the last valuation may be enlarged.

These are precisely the inherited compiler's conclusions. The separation
argument needs no other feature of its rational anchor or exit budget.

**PROVED separation theorem.** If w and v are distinct marked primitive
words, and `m,n>=2`, their positive-source certificate sets are disjoint.

Suppose a source belonged to both. Its uniquely defined exact first-descent
time would be

    T=pm=qn.

The two repeated nominal words agree in positions `0,...,T-2`. Let
`zeta=exp(2*pi*i/T)`. For the first word,

    sum_(j=0,...,T-1) (w^m)_j*zeta^j
      = (sum_(a=0,...,p-1) w_a*zeta^a)
        *(sum_(b=0,...,m-1) zeta^(pb)) = 0.

The last geometric factor vanishes because `zeta^p` has exact order
`m>=2`. The corresponding sum for `v^n` vanishes for the same reason.
Their difference has only one possible nonzero coordinate, the final one:

    [(w^m)_(T-1)-(v^n)_(T-1)]*zeta^(T-1)=0.

The phase factor is nonzero, so that final letter agrees too. Thus
`w^m=v^n` as complete marked words.

For completeness, the common circular word is invariant under shifts p
and q, hence under their integer combinations, and therefore under
`gcd(p,q)`. Its first p letters consequently have that smaller period;
primitivity forces `p=gcd(p,q)=q`. Their marked words are equal, contrary
to the hypothesis. This proves disjointness at all heights.

The source is a pair of actual certificates; the Fourier target is a
length-T nominal sequence. The map preserves the repeated word and clock.
The Fourier coefficient alone loses essentially the whole word: every
nontrivial repetition vanishes there. It becomes decisive only because
the shared actual orbit has already forced equality at T-1 positions.

The complex character is equivalently its cosine and sine components.
A fixed cosine component alone can be blind: at T=4, the vectors
`(1,1,1,1)` and `(1,1,1,2)` both have cosine coefficient zero, while their
sine coefficients are 0 and -1. This is a vector-level hostile, not a claim
that both vectors meet the twice-repeated compiler hypotheses. Conversely,
this particular proof tests **vanishing**, so the exact modulus of the one
complex coefficient also suffices. No general phase-retrieval conclusion
is being drawn.

## 2. The hypotheses cannot be dropped

The inherited compiler gives exactly the same cylinder for

    w=(1,2), m=1: C(3,4),
    v=(1),   n=2: C(3,4),

where `C(r,K)` means `r mod2^K`. The source 3 follows

    3 -> 5 -> 1, exact valuations (1,4).

Its nominal words `(1,2)` and `(1,1)` differ only at the last letter;
the actual final divisions absorb that discrepancy. The first word has
no nontrivial repetition factor, so the Fourier proof correctly fails.

Primitivity is essential too. Word `(1,1)` repeated twice and word `(1)`
repeated four times are the same marked nominal word and compile to the
same cylinder `15 mod128`. The theorem distinguishes primitive charts,
not multiple spellings of one repeated root.

There is also a structural coverage hostile. The exact first descent of
7 occurs at time four:

    7 -> 11 -> 17 -> 13 -> 5,
    valuations (1,1,2,3).

If 7 belonged to any twice-repeated chart, its primitive period would be
a proper divisor of four, hence 1 or 2. Period one would require the first
three valuations to be equal. Period two would require them to have shape
`(a,b,a)`. Both contradict `(1,1,2)`. Thus 7 lies outside the universal atlas,
although its route continues `5 -> 1`. Since 3 lies inside, 7 is the smallest
positive odd source `3 mod4` missed by this atlas. Longer period cutoffs
cannot repair this miss; a less periodic or mixed-head certificate is needed.

## 3. Every favorable growth word is admissible from repetition two

For a positive valuation word w of length p, put

    A_i=a_1+...+a_i, A=A_p, P=3^p, Q=2^A,
    S=sum_(i=0,...,p-1) 3^(p-1-i)*2^A_i.

Assume total growth `P>Q`, and write the reduced anchor as
`-h/d=-S/(P-Q)`. Each remaining valuation is at least one, so

    S/Q <= sum_(j=0,...,p-1) 3^j/2^(j+1)
         = (3/2)^p-1 < 2^p <= Q.

Since `h<=S`, this proves `h<Q^2`. Consequently every favorable word,
meaning every nonempty prefix satisfies `3^i>2^A_i`, meets the inherited
compiler's domain `Q^m>h` at every `m>=2`. There is no additional unbounded
search for its first admissible repeated chart.

Define the **universal twice-repeated atlas** as the union of all compiler
cylinders of every marked primitive favorable word at all `m>=2`, using
the inherited least sufficient exit budget. The separation theorem handles
different words; `v2(d*n+h)=Am` separates repetitions of the same word.
Thus all these cylinders are pairwise disjoint.

Marked cyclic rotations are retained as distinct charts. They can have
different favorable marked phases and represent different source families.
Selecting instead the closest-to-zero rational phase gives one canonical
chart per primitive necklace, but that is a smaller atlas and is not the
counting convention used for the density below.

## 4. Exact finite counts and a summable symbolic period tail

Let `B_p` be the largest integer A with `2^A<3^p`. There are

    M_p=sum_(A=p,...,B_p) C(A-1,p-1)=C(B_p,p)

positive marked growth words of length p. Every such word has a unique
primitive root, which is itself a growth word. Therefore the number of
primitive marked growth words is

    G_p=sum_(d|p) mu(d)*C(B_(p/d),p/d),

and the number of primitive necklaces is `G_p/p`. The script checks this
independently against literal composition enumeration and cyclic rotations.
The closest-to-zero rotation of every necklace passes the favorable-prefix
condition, while all favorable marked rotations are counted separately.

| p | B_p | All marked growth words | Primitive necklaces | Favorable primitive marked charts |
|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 1 | 1 |
| 2 | 3 | 3 | 1 | 1 |
| 3 | 4 | 4 | 1 | 2 |
| 4 | 6 | 15 | 3 | 5 |
| 5 | 7 | 21 | 4 | 11 |
| 6 | 9 | 84 | 13 | 26 |
| 7 | 11 | 330 | 47 | 84 |
| 8 | 12 | 495 | 60 | 166 |
| 9 | 14 | 2,002 | 222 | 473 |
| 10 | 15 | 3,003 | 298 | 948 |
| 11 | 17 | 12,376 | 1,125 | 2,651 |
| 12 | 19 | 50,388 | 4,191 | 8,010 |

Totals through twelve are **5,966 necklaces and 12,378 marked favorable
charts**. The count through eight is the inherited corrected 130 necklaces.

The elementary bound

    number of favorable primitive marked charts <= C(B_p,p)
      <= 2^(B_p-1) < 3^p/2

is sufficient for an all-period tail. Each chart's budget gives cylinder
mass `<P^-m`, hence total mass over `m>=2` less than `1/[P(P-1)]`.
For periods `p>D`, summing the cylinder probabilities gives

    period tail < sum_(p>D) 1/[2(3^p-1)]
                <= 3/[4(3^(D+1)-1)].                (1)

At this stage these are valid probability bounds for uniform infinite
binary addresses, equivalently Haar measure on the 2-adic integers.
They do not, by themselves, prove a natural-density assertion about the
positive integers. That separate step follows next.

## 5. Why the universal union really has natural density

The distinction is necessary. The cylinders `C(n,n+10)`, one for every
positive n, have total cylinder mass at most `1/1024`, yet their union
contains **every positive integer**, so its natural density is one.
Summing a countable list of cylinder masses is not a general natural-density
argument.

For the present atlas there is an additional dynamical tail bound.
For any exact length-k positive odd valuation word with total A,

    U^k(n)=3^k*n/2^A+c,
    0<c<=C_k=(3/2)^k-1.

Put `B_k=floor(log2(3^k))` and

    epsilon_k=1-3^k/2^(B_k+1)>0.

For every source `n>C_k/epsilon_k`, total valuation `A>=B_k+1` forces
`U^k(n)<n`. The excluded sources are finite and depend only on k, not on
the particular word or its potentially unbounded valuation sum.

Therefore the set of sources with no descent through time k has upper
natural density at most

    H_k=sum_(A=k,...,B_k) C(A-1,k-1)*2^(-A-1).        (2)

There are finitely many words in (2). A word of total A gives one exact
source cylinder modulo `2^(A+1)`, with its last-oddness bit retained; distinct
words are disjoint. This proves (2) using ordinary densities of finite unions.

The bounds H_k tend to zero along an explicit subsequence. On the product
distribution `Pr(a=j)=2^-j` for positive valuations, the generating function
at `z=3/4` is `z/(2-z)=3/5`. Thus

    H_k <= (1/2)*(3/5)^k*(4/3)^B_k.

Since `log2(3)<8/5`, proved by `3^5<2^8`, this gives

    H_(5j) < (1/2)*(65536/84375)^j -> 0.             (3)

An atlas source whose marked period is `p>D` first descends at `pm>=2p>2D`.
It is therefore in the no-descent set through
`k=5*floor(2D/5)`. Equations (2), (3) prove that the omitted-period union has
upper natural density tending to zero as D increases. No monotonicity of
the numerical sequence H_k is assumed; we use the nested no-descent sets.

For each fixed finite D, the atlas with periods at most D has a natural
density: finitely many chart families are involved, and the tail `m>M` of
each lies in the shrinking parent `d*n+h=0 mod2^(A(M+1))`. Its density is
the sum of its cylinder masses, as in the earlier fixed-atlas proof.
The full atlas contains these finite-period unions and exceeds them by an
upper-density tail tending to zero. Consequently its natural density exists
and equals the limit of their densities, hence the disjoint cylinder sum.
**Only now** may (1) serve as a sharper numerical natural-density error bound.

This argument also explains why a long first-descent delay is a much stronger
tail restriction than a large, otherwise untyped cylinder modulus. It does
not prove that every positive source descends, and descent alone does not
supply a completed route to 1.

## 6. Computed atlas interval and verification scope

The complete finite mass calculation uses all 12,378 marked charts with
`p<=12` and every `2<=m<=30`: **358,962 disjoint cylinders**. Let L denote
their exact dyadic mass, retained in the output. The remaining repetitions
of these known charts contribute less than

    E_repeat=sum_(enumerated charts) 1/[(P-1)*P^30].

All omitted periods contribute less than `3/[4(3^13-1)]`. Thus

    L < delta_atlas < L+E_repeat+3/[4(3^13-1)].

The introductory numerical interval rounds these exact endpoints outward. It is
not a Monte Carlo estimate, and it does not subtract the older bank or the
previous infinite return families. The atlas omits the single-copy charts
by definition, including any additional source cylinders unique to them.

Independent finite controls include:

- Every positive composition in the stated period/growth universe, checked
  against the binomial/Mobius formulas and rotation classes.
- Exact polynomial remainders modulo cyclotomic polynomials for every
  primitive-mode geometric factor at lengths `T=2,...,120`; no approximate
  Fourier arithmetic is used.
- Every pair of distinct primitive words over `{1,2,3}` of lengths at most
  five, compared after a common nontrivial repetition.
- The single-copy, nonprimitive, cosine-only, countable-density and source-7
  hostiles above.
- Literal first-descent replay of each chart's least positive source at
  `m=2,3`: **24,756** complete controls, with source bit length at most 58.
- Independent affine-carry and contraction controls for k=1,...,6 on 256
  successive odd sources above each proved threshold; exact density-tail
  and generating-function comparisons through k=60.

Reproduce with:

    python -X utf8 -B 04-computation/experiments/periodic_chart_separation_20261004.py
    python -X utf8 -B -O 04-computation/experiments/periodic_chart_separation_20261004.py

The [saved output](periodic_chart_separation_20261004.out) contains all count
rows, exact interval endpoints, exponent histogram, tail controls and hostile
records. All acceptance checks remain active under optimized Python. A
structured family can carry its own return certificate; a proof that every
integer has such a family representation remains a separate obligation.
