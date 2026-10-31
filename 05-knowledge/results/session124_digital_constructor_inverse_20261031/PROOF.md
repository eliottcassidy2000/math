# Exact digital normal forms and inverse conditions for marked constructors

**Status: proved identities and exact inverse algorithms, tested on actual seed
witnesses. No unbounded tournament-spectrum coverage theorem is claimed.**
This work is separate from the sealed million-prefix publication packet. It
extends the structure, not the finite coverage ceiling.

## Main results

1. The proposed digital formula is correct, including the empty/all-one case.
2. The scalar function has an exact binary inversion formula and an exact,
   nonnegative radical-inverse error controlled by runs of one-bits.
3. A unique backward `q` map yields a potential invariant and exact actual-seed
   fibers. Singleton suffixes cannot increase a fixed fiber's offset set once
   its finite seed sources have been passed.
4. A triangle is a two-bit macro in the normalized `q` address. There are
   explicit local defects from three times its two-singleton shadow.
5. Combining this with the modulo-three restriction gives a **complete finite
   inverse algorithm for every non-3-divisible target in the entire specified
   singleton/triangle language**, with arbitrary word lengths. The inverse
   uses a final-state bound, not invalid pruning inside arbitrary block words.

The fifth result is an exact decision/reconstruction normal form, not a proof
that the decision is always positive.

## 1. Conventions and the digital formula

Use state `(h,F,q,g)`, where `g=h-t`. An actual seed is a nonempty tournament A
and an external vertex x with a fixed cut. The coordinates count:

* `h=H(A)`;
* `F=H(A+x)`;
* `q`: ordered two-path covers `(P,Q)` of A, Q nonempty, P possibly empty,
  with x beating the first vertex of Q;
* `g=h-t`, with t the number of full paths of A ending at a vertex beating x.

For the actual seeds here, h is positive and odd, F is positive, q is
nonnegative, and `0<=g<=h`. The proved singleton transitions are

```
S0(h,F,q,g) = (h,F+q+g,2q,0),
S1(h,F,q,g) = (h,F,2q+h,h).
```

For a chronological word `b_0,...,b_(k-1)`, define

```
B = sum_{i=0}^{k-1} (1-b_i) 2^i.
```

Thus the **least significant address bit records the earliest singleton**, and
address bit 1 means singleton 0. Trailing singleton 1s do not change B or F.
If B>0 its bit length specifies the word after those trailing 1s are removed;
B=0 represents the empty word/all-one words as far as F is concerned.

Define

```
f(0)=f(1)=1,
f(2n+1)=f(n)              (n>=0),
f(2n)=f(n)+n+(n mod 2)    (n>=1).
```

Then, for every actual seed and every finite singleton word,

> `F_word = F0 + q B + g (B mod 2) + h (f(B)-1)`.

### Proof directly from the transitions

After i singleton steps,

```
q_i = 2^i q + h sum_{j<i} 2^(i-1-j) b_j,
g_0 = g,
g_i = h b_(i-1)           (i>=1).
```

Only zero steps increase F, by `q_i+g_i`. Summing these increments proves the
formula below and hence the displayed result. Alternatively split by the first
singleton: address `2n+1` starts with 0, giving initial increment `q+g` and
remaining q equal to 2q; address `2n` starts with 1 and leaves F unchanged while
replacing `(q,g)` by `(2q+h,h)`. These give precisely the two recurrences for f.
No profile relaxation or realizability assumption beyond the starting seed is
introduced.

## 2. Exact digit formula and radical-inverse asymptotic

Let `z_i` be the binary digits of n, and let

```
S(n) = sum_{i>j} z_i z_j 2^(i-j-1),
a(n) = sum_{i>=1} z_i z_(i-1).
```

The exact formula is

> `f(n) = 1+n-z_0-S(n)-a(n)`.

An equivalent nonnegative weighted-inversion form is

```
f(n) = 1 + sum_{j: z_j=0} ( floor(n/2^(j+1)) + z_(j+1) ),
```

where only `0<=j<bit_length(n)-1` are needed. The latter follows by unrolling
the even recurrence. Expanding its two sums gives the former. In particular
`1<=f(n)<=n+1`, and f(n) is always odd.

Define the base-two radical inverse rho_2 (not a 2-adic valuation) explicitly by

```
rho_2(n) = sum_i z_i / 2^(i+1).
```

For n>0, let its maximal one-bit runs be `[a_l,b_l]`, with length `L_l`, and put

```
R_l = sum_{j>b_l} z_j 2^(b_l-j).
```

Each `0<=R_l<1/2`, because the bit immediately above a run is zero. Then

> `f(n) = n(1-rho_2(n)) + E(n)`,
>
> `E(n) = 1-z_0 + sum_l [ R_l + 2^(-L_l)(1-R_l) ]`.

This is an equality of rational numbers, not a floating-point fit. Consequently
`E(n)>0` and `E(n)<=1-z_0+3r/4`, where r is the number of one-runs. In particular

> `f(n) = n(1-rho_2(n)) + O(log(n+1))`.

To derive the error, expand `n rho_2(n)` into terms i>j, i=j, and i<j. If w is
popcount and `T=sum_{i<j}z_i z_j 2^(i-j-1)`, the correction is
`1-z_0+w/2+T-a`. Since `w-a=r`, grouping the missing higher-bit mass separately
within each one-run gives the displayed expression. All-one n has
`E(n)=2^(-bit_length(n))`; powers of two greater than one have E=3/2.

The estimate is **not** a uniform lower bound proportional to n: along
`n=2^k-1`, one has f(n)=1. It does not by itself prove density or interval
coverage.

## 3. Exact scalar inverse, without a depth cutoff

For a fixed actual seed with q>0, every address producing target N satisfies

```
0 <= B <= floor((N-F0)/q).
```

Indeed all other terms in the digital formula are nonnegative. Testing that
finite interval in the exact formula is therefore necessary and sufficient.
Each successful B gives the actual chronological word with bits `1-z_i`.

For an actual seed with q=0, x beats no vertex: a beaten vertex v would supply
a cover with Q=(v), using a Hamiltonian path on A-v as P. Hence `(F0,g)=(h,0)`.
Initial singleton zeroes are then exact neutral loops and can be deleted.
The target must have N=hm, with `m=f(B)`. For m>1, strip trailing binary ones
from B using `f(2n+1)=f(n)`. The first nontrivial outer argument is even, so

> `m=f(n)+n+(n mod 2)` for some `(m-1)/2 <= n <= m-1`.

The upper bound follows from f(n)>=1. For the lower bound, even n has
`f(n)<=n+1`; odd n has `f(n)=f((n-1)/2)<=(n+1)/2<=n`, so in both cases the right
side is at most `2n+1`. The branches also obey

```
G(2k)   = 3k + f(k) + (k mod 2),
G(2k+1) = 2k + 2 + f(k),
```

where `G(n)=f(n)+n+(n mod 2)`.

For example no n from 6 through 12 has G(n)=13. This proves **globally** that
the one-transitive-seed singleton grammar never produces 13, irrespective of
word length. It is not an exclusion for unrestricted tournaments.

`digital.inverse_addresses` implements the exact finite inverse. It checks the
actual empty-cut condition for q=0; it does not treat a formal positive vector
as a realizable seed.

## 4. The backward carry map and a singleton potential

For fixed odd h, the last singleton is uniquely read from the parity of q:

```
q even:  predecessor q/2,       last bit 0;
q odd:   predecessor (q-h)/2,   last bit 1, provided q>=h.
```

The reduced ancestry stops at 0 (remove its zero self-loop), or an odd q<h.
Every other predecessor is strictly smaller. Define the integer potential

```
Phi_h(0)=0;
Phi_h(r)=0 for odd 0<r<h;
Phi_h(2u)=Phi_h(u)+u+h(u mod 2)  for u>0;
Phi_h(2u+h)=Phi_h(u).
```

After one singleton, `g=h(q mod 2)`. At all later singleton steps,

> `F-Phi_h(q)` is invariant.

For an initial seed `(h,F0,q0,g0)`, its first-child offsets are precisely

```
first bit 1: kappa_1 = F0-Phi_h(q0),
first bit 0: kappa_0 = F0-Phi_h(q0)+g0-h(q0 mod 2).
```

Associate each offset with its **actual first-child q node**, respectively
`2q0+h` or `2q0`. For final q, that offset is available if and only if the
unique backward ancestry of q encounters its associated node. This condition
is essential: the union of all numerical offsets is not a free bank of seeds.
For any finite actual seed bank, this gives the exact nonempty-singleton fiber

```
{ F at fixed h,q,g=h(q mod 2) }
 = Phi_h(q) + { source offsets whose source node is an ancestor of q }.
```

If L is the largest first-child source q, no new offsets enter above L. Every
q>L has a unique boundary ancestor c satisfying `L<c<=2L+h` and
`predecessor(c)<=L`. Its offset set then remains unchanged on the entire
infinite descendant cylinder rooted at c. This gives an exact finite description
of infinite singleton fibers, not a finite output census.

For the 273 actual seed states of orders at most five:

| h | Distinct global offsets | Largest source q | Eventual maximal fiber size | Eventual maximal consecutive-odd interval size |
|---:|---:|---:|---:|---:|
|1|1|63|1|1|
|3|12|93|4|2|
|5|25|123|3|2|
|9|38|159|5|3|
|11|33|189|3|2|
|13|28|195|3|2|
|15|21|225|4|2|

The exact source-ancestry boundary sets are in `fiber_boundaries.json`.
For example the two seed F values 3 and 5 at `(h,q,g)=(3,9,3)` retain their
separation 2 under every singleton suffix; the width does not multiply.
This rules out *unbounded amplification of a fixed singleton fiber* as an
interval mechanism. It does not rule out overlapping different q cylinders,
or triangle-induced offsets.

## 5. Triangle blocks as two-bit macros with local defects

For a directed triangle block C_j with x beating j local vertices, the action is

```
h' = 3h,
F' = 3F+(9-3j)q+(3-j)g + (2h if j=1 or 2, else 0),
q' = 12q+3jh,
g' = jh.
```

Fix the original seed h0. After r triangles write `h=3^r h0` and `Q=q/3^r`.
Q is integral, and

```
singleton b: Q -> 2Q+h0 b;
triangle Cj: Q -> 4Q+h0 j.
```

Thus C_j is exactly the two-bit binary macro for j in the **Q address**.
Here macro 01 means the chronological singleton shadow S0 then S1, and macro
10 means S1 then S0. On states written `(h,F,q,g)`, the exact operator defects
relative to three times that two-singleton shadow are

```
C0 - 3(S0 then S0): (0,0,0,0),
C1 - 3(S0 then S1): (0,3q-g+2h,0,-2h),
C2 - 3(S1 then S0): (0,-3q+g-4h,0,2h),
C3 - 3(S1 then S1): (0,0,0,0).
```

These are linear identities comparing actual constructor operators. They do
**not** assert that an arbitrary positive vector or a scaled intermediate
profile is itself a seed. Actual construction always starts from the retained
seed and performs the stated blocks.

Potential identities useful for carry calculations are

```
Phi_(3h)(3q) = 3 Phi_h(q),
Phi_h(4q)    = Phi_h(q)+3q+h(q mod 2),
Phi_h(4q+h)  = Phi_h(q)+q+h(q mod 2),
Phi_h(4q+2h) = Phi_h(q)+2q+2h,
Phi_h(4q+3h) = Phi_h(q).
```

If the input is canonical, `g=h e`, `e=q mod 2`, and
`kappa=F-Phi_h(q)`, then immediately after one nonconstant triangle the new
potential offsets are

```
C1: 3 kappa + 3q - he + 2h;
C2: 3 kappa - 3q + he - 4h.
```

A following singleton 1 preserves these offsets; a following singleton 0
changes them by -2h for C1 or +2h for C2. Unlike pure singleton propagation,
these offsets vary with the triangle insertion q. That is a genuine additional
resource, but their actual ancestry and macro placement cannot be discarded.

## 6. Complete non-3 inverse normal form for the whole language

Modulo three, C0 and C3 annihilate the whole state. A nonconstant triangle has
state `(h',F',q',g')=(0,(3-j)g+2h,0,jh)` modulo three. A following singleton
annihilates its h,q,g coordinates, which remain zero through any further
singletons; another triangle then annihilates F too. Three consecutive
nonconstant triangles also annihilate all coordinates. Zero remains zero.

Therefore any output not divisible by three has:

* no C0 or C3;
* at most two triangle blocks;
* two triangles adjacent, if there are two.

This argument applies to every initial integer vector, so in particular to
all actual seed states. It does not classify unrestricted tournaments.

### A safe final-state bound

Remove all trailing singleton 1s, which preserve F. If nothing remains, the
output is simply the seed F. Otherwise the last block is S0, C1 or C2. For
actual input `F>=1,h>=1,g>=0`:

```
last S0: q'=2q < 2(F+q+g)=2F';
last C1: 2F'-q' = 6F+4g+h > 0;
last C2: 4F'-q' = 12F+4g+2h > 0.
```

Consequently every nonempty reduced witness has

> `q_final < 4 F_final`.

This bounds only a **reduced final state**. It does not prune an intermediate
state with large q and therefore avoids the known arbitrary-block pruning gap.

### Necessary and sufficient finite inverse procedure

For target N not divisible by three and each actual seed:

1. Accept the empty word if its F equals N.
2. For `r=0,1,2`, enumerate integers
   `0<=Q<=floor((4N-1)/3^r)`.
3. Use the unique backward h0 carry map to ask whether seed q0 is an ancestor
   of Q, and recover the chronological binary shadow word.
4. Mark no pairs (r=0), one 01/10 pair (r=1), or two **adjacent** 01/10 pairs
   (r=2) as triangle macros; every unmarked bit is a singleton.
5. Execute those actual blocks from that actual seed and test F=N.

There is one necessary endpoint exception: at seed q0=0, reduced ancestry
forgets initial zero bits. All initial singleton zeros are neutral, but a zero
belonging to a starting C1=01 macro must be retained. Thus also allow **one**
leading zero on the reduced shadow word. More are unnecessary: the first
nonconstant macro contains a one, and preceding singleton zeros are neutral.

Every accepted word is an actual construction. Conversely every non-3 output
can be reduced as above, obeys the final q bound, has one of these macro
placements, and is therefore enumerated. This proves completeness of the
inverse decision algorithm over **all finite words in this specified language**.

`digital.inverse_non3` implements this exact procedure. It can be slow; no
polynomial-time, all-prime, density, or universal positive-answer theorem is
asserted.

## 7. The precise remaining global obligation

To prove non-3 global coverage with these seeds, one must prove that **for every
allowed N not divisible by three**, at least one actual seed, integer Q in the
proved bound, and legal zero/one/two-macro placement passes the final F equation.
The identities reduce the problem to explicit carry-and-defect arithmetic but
do not establish that universal existence statement. An eventual finite-state
carry proof would have to prove its carry-state bounds and preserve source
ancestry; neither is supplied by the finite seed bank alone.

For targets divisible by three, ordinal products can help but do not remove
all obligations: for instance `3*7=21` remains forbidden, while some larger
multiples involving 7 are attained without an H=7 factor. No unrestricted
factor closure or prime-lane conclusion is inferred here.

The concrete singleton obstruction at 103061 and its repaired triangle witness
show that a nonconstant macro genuinely adds power. The frozen-fiber theorem
explains why singleton suffix length by itself cannot widen the existing
same-coordinate intervals. These are rigorous structural reductions, not an
unbounded completeness proof or a larger finite census.

## Reproduction and verification

All paths are file-relative. Python standard library only. No `assert` is used.

```sh
PYTHONDONTWRITEBYTECODE=1 python3 check.py
PYTHONDONTWRITEBYTECODE=1 python3 -O check.py
PYTHONDONTWRITEBYTECODE=1 python3 check_macro_inverse.py
PYTHONDONTWRITEBYTECODE=1 python3 -O check_macro_inverse.py
PYTHONDONTWRITEBYTECODE=1 python3 fibers.py --output /tmp/fiber_boundaries.normal.json
PYTHONDONTWRITEBYTECODE=1 python3 -O fibers.py --output /tmp/fiber_boundaries.optimized.json
```

Normal and optimized outputs match. Retained checks cover:

* all 273 actual seed states and all singleton words of lengths 0–10:
  558,831 formula/invariant cases;
* 32,768 recurrence-versus-digit-formula cases and 4,096 exact rational error
  checks;
* 17,472 actual-state triangle macro identities;
* 1,638 full literal adjacency Hamiltonian-path DPs, without the potential or
  digital formula in the path-counting recurrence;
* 8,463 actual-seed scalar inverse cases;
* the exact source-ancestry fiber recurrence and boundary freezing;
* 51,085 non-3 actual-seed pre/post-word cases, each with pre/post lengths 0–2
  and zero, one or two adjacent nonconstant triangle macros, including the
  final-state bound and the q=0 leading-zero exception.

The input `actual_seed_roots.json` is an unchanged copy of the already
independently verified 273 seed witnesses. Seed code bits enumerate i<j in
lexicographic order, bit 1 meaning i->j; cut bit 1 means x->i. This retains
actual graph constructors rather than substituting arbitrary positive states.
`SHA256SUMS` and `verification.json` record the final checked bytes.
