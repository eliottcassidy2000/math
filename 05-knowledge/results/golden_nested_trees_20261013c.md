# Golden spectra and nested recursions: retain the root, the tree boundary and both embeddings

**2026-10-13, renewal_shapes_oct11 continuation.**
**Status:** PROVED elementary identities and the golden depth-ball character
law; inherited PROVED termination of the golden control map on Z[phi];
FINITE-EXACT nested recurrences through `2^20`, Pascal-word summaries through
level 1024, golden ray census through `2^22`; NUMERICAL ray spectra;
CONDITIONAL Q covariance model with a proved Pell conjugacy and variance law;
OPEN
transport of these mechanisms to a universal Collatz certificate.
Self-audited. Classical sequence facts are credited separately below.

Companion: [the Collatz boundary-current and conditional-pole result](collatz_amplitude_boundary_20261013c.md).
The shared question is what an aggregated statistic has forgotten about the
source that generated it. Here three forgotten coordinates are decisive:
the terminal cycle, the shape of a region in the two real embeddings, and
the location of the last generation in a tree.

## 1. Two nested recurrences, with different conserved quantities

The user's recurrence is the Conway--Hofstadter sequence

    a(1)=a(2)=1,
    a(n)=a(a(n-1))+a(n-a(n-1)).                       (1)

Hofstadter's Q is instead

    Q(1)=Q(2)=1,
    Q(n)=Q(n-Q(n-1))+Q(n-Q(n-2)).                     (2)

[Pelesko, section 1](https://cs.uwaterloo.ca/journals/JIS/VOL7/Pelesko/pel11.pdf)
records (1), its slow-growth and dyadic-anchor properties, and develops
modular-clock variants. [Kubo--Vakil](https://doi.org/10.1016/0012-365X(94)00303-Z)
describe the classical finite-set compression approach. The full latter
paper was not accessible in this session; no uninspected theorem from it
is used as a new proof dependency. The Pascal-word construction below is
also described by its authors in their
[1992 discussion](https://groups.google.com/g/rec.puzzles/c/mvuuEANp9iM).
Our explicit agreement check has the finite scope stated below.

[Pinn's primary study of Q](https://arxiv.org/pdf/chao-dyn/9803012)
investigates generations and fluctuations empirically. Its reported
patterns do not establish an all-index definition of Q. The current
calculation likewise makes no all-index claim for Q.

### The copied increment tree (PROVED, elementary rederivation)

Put `d_n=a(n)-a(n-1)`, with `d_2=0,d_3=1`. For `n>=4`, let `j=a(n-1)`.
Subtract two successive instances of (1). If `d_(n-1)=0`, the first
summands agree; if `d_(n-1)=1`, the second summands agree. Therefore

    d_n = d_(n-j)  if d_(n-1)=0,
    d_n = d_j      if d_(n-1)=1.                     (3)

Inductively all these parents lie in `[2,n-1]`. Initially `a(3)=2`;
if the previous increments are 0 or 1 then `1<=a(n-1)<=n-2`, so the
next defining arguments are smaller positive indices. Formula (3) then
shows that the new increment is again 0 or 1 and keeps this bound for the
next stage. This proves well-definedness and slow growth for (1).

Thus each increment has a terminating certificate that it is copied from
one of the two root bits. The certificate is an ordered parent history,
not just its terminal value. Losing the parent index would lose the proof
of termination. The exact audit through `2^20` finds maximum parent depth
19. This finite depth observation is not itself an all-depth theorem.

### An additive leaf charge (PROVED)

The two arguments of (1) add **exactly to n**. Define `e(n)=2a(n)-n`.
Then

    e(n)=e(j)+e(n-j), e(1)=1, e(2)=0.                (4)

In the full computation tree, e(n) counts leaves labelled 1; leaves labelled
2 carry no charge. Leaf indices sum to n, and the number of leaves is a(n).
This gives `e(n)>=0` and a literal interpretation of the excess above n/2.
The split tree preserves additive mass; its feedback rule selects which
partition occurs. An arbitrary binary splitting tree would preserve (4)
but would not automatically have Conway's quantitative balance.

For Q, with `p=n-Q(n-1), q=n-Q(n-2)` and `E_Q(n)=2Q(n)-n`, the corresponding
identity, wherever its indices are defined, is

    E_Q(n)=E_Q(p)+E_Q(q)+(p+q-n).                    (5)

The extra term is `n-Q(n-1)-Q(n-2)`. Q has a local mass defect, while
Conway's tree conserves the index. This is a concrete reason that sharing
the initial values and nested syntax does not transfer a proof.
Q's first non-slow increment is `Q(12)-Q(11)=8-6=2`.

The numbers requested by the user separate the two models clearly:

| n | a(n) | Q(n) | Conway leaf charge 2a(n)-n |
|---:|---:|---:|---:|
| 20 | 12 | 12 | 4 |
| 1024 | 512 | 557 | 0 |
| 2^20=1024^2 | 524288 | 527047 | 0 |

These are exact calculations. The equal values at 20 do not identify the
recurrences. At 1024 the Conway computation tree has no unit-labelled
leaves, whereas Q retains a nonzero defect. Pelesko's stated classical
property identifies the Conway zero-charge indices as the powers of two;
the present work does not claim that classical result as new.

## 2. A useful compression that reaches a 2^1024-length object

Define finite binary words by

    W_(k,0)=1, W_(k,k)=0,
    W_(k,j)=W_(k-1,j-1) W_(k-1,j), 0<j<k,
    P_k=W_(k,0) ... W_(k,k).

The first complete rows are

    P_1=10,
    P_2=1100,
    P_3=11101000.

Induction gives `|W_(k,j)|=binom(k,j)` and, reading 1 as +1 and 0 as -1,

    sum W_(k,j)=binom(k-1,j)-binom(k-1,j-1).

Consequently P_k has length `2^k` and total signed sum zero. Its sum just
before cell j is `binom(k-1,j-1)`. Inside a cell it changes by at most the
cell's length, so every partial sum is bounded in absolute value by

    2 binom(k,floor(k/2)) = O(2^k/sqrt(k)).          (6)

This is a proved bound for these explicitly constructed words. The
independent recurrence computation verifies that P_k is the block of
Conway increments `d_(2^k+1),...,d_(2^(k+1))` for every `k<=19`.
The familiar classical identification motivates the test; we do not pass
a finite check off as a new proof of its all-level equivalence.

To find an exact maximum without expanding the word, store a cell as

    (length, total charge, greatest prefix charge,
     first position of that greatest charge).

For concatenation UV, the new maximum is

    max(maxprefix(U), total(U)+maxprefix(V)),        (7)

with the earliest position retained on a tie. Length and total charge add.
This is an associative, lossless summary for the declared maximum query.
It is not lossless for recovering the word. The Pascal DAG needs only
quadratically many joins in k, while the row itself has `2^k` letters.

The [checker](../../04-computation/experiments/nested_trees_20261013c.py)
performs **6,293,490 exact controls** and computes these summaries through
level 1024. Writing M_k for the largest prefix charge,

    sqrt(20) M_20/2^20 = 0.4234935873...,
    sqrt(1024) M_1024/2^1024 = 0.3993652412... .

The maxima and locations are exact large integers in
[the saved output](nested_trees_20261013c.json); the displayed ratios are
floating summaries. Their proximity to a Gaussian constant is a useful
consistency check, not a newly proved limit. This answers the information
compression prompt concretely: a large object can carry an exact maximum
certificate in a small DAG when its query composes by (7).

For the Collatz discrepancy the analogous interval summary retains the
signed total, prefix extrema, the translated first moment, and terminal
phase. These summaries already compose exactly in the new census. What
is missing is an all-source structural rule that bounds them under every
arithmetic refinement, rather than only computes a finite list of them.

## 3. The golden control has an exact, nonzero depth spectrum

Let `phi^2=phi+1`, and use the existing integral map

    G(x)=x/2                         if x=0 mod2,
    G(x)=(phi x+rho(phi x))/2         otherwise,

where rho takes coefficient digits in `{0,1}^2`. In coordinates,

    G(a,b)=(a/2,b/2)                 if a,b even,
           ((b+b%2)/2,(a+b+(a+b)%2)/2) otherwise.

The [previous two-embedding proof](tree_resonance_composite_pell_20261013b.md)
shows that every source in Z[phi] reaches one of the four fixed points

    R={0,phi^2,phi^3,phi^4}.

Its
height bound takes the maximum over both branches. It is not enough to
apply only the odd-branch conjugate contraction to an even source.

**Lemma G2 (PROVED).** G is exactly four-to-one on Z[phi]. For each y its
four distinct inverse images are

    2y,
    phi^(-1)(2y-d), d in {1,phi,1+phi}.              (8)

The first has zero residue modulo2. The others have the three distinct
nonzero residues, since phi is a unit. Each replays to y exactly.

For a fixed root r, let B_K(r) consist of the sources with `G^K(x)=r`.
The sets are nested because r is fixed, and

    |B_K(r)|=4^K,
    #{x:first hit of r at depth k}=3*4^(k-1), k>=1.

There is no overlap between the four root ports. Their union over all K
is the whole integer lattice, using the inherited termination theorem.
This is a fully covered, richer-than-integer control model.

**Theorem G3 (PROVED).** If `tau_G` is strict first-hit depth, then

    (1/4^K) sum_(x in B_K(r)) (-1)^tau_G(x)
      = (3/5)(-1)^K+(2/5)4^(-K).                   (9)

It is the geometric sum `1+sum_(k=1)^K 3*4^(k-1)(-1)^k`, divided by
`4^K`. The same mean holds after joining all four root ports.
More generally, for any unit-modulus zeta,

    (1/4^K) sum zeta^tau_G(x)
      = [3 zeta/(4 zeta-1)] zeta^K
        +[(zeta-1)/(4 zeta-1)]4^(-K).               (10)

This is an **exact amplitude**, not a fitted one: parity has limiting
oscillation amplitude 3/5 in the certificate-depth ensemble. Three quarters
of the new ball lie in its outer generation. Full termination and a
persistent boundary oscillation coexist.

The [inverse-tree audit](../../04-computation/experiments/golden_tree_audit_20261013c.py)
enumerates through depth 7 at all four ports, checks every inverse edge,
strict depth, nesting and the formula. The final population has 65,536
sources. This is a direct hostile to the incoming suggestion that a
one-sided contracting walk should have no spectral lines. The repaired
claim must state its ensemble. Certificate-depth balls, integer rays and
Euclidean boxes count different populations; none can silently replace
another.

The terminal forcing is equally important. On this lattice,

    t_G(x)=-t_G(Gx)+2 1_R(x).                       (11)

Any *homogeneous* antifunction of G is identically zero on all integral
sources, since it must vanish on a fixed point and all sources reach one.
Thus omitting the four root terms would erase the very observable being
studied. This is the golden version of the character defect in (1) of the
companion note, rather than a new obstruction to G's known termination.

## 4. Integer-ray spectra and the missing shape coordinate

The [ray census](../../04-computation/experiments/golden_parity_ray_20261013c.cpp)
counts strict stopping parity for `x=n` and `x=n phi^2`, `1<=n<=2^22`.
Its integer bin endpoints are stored explicitly (their floating construction
is not an interval certification of an ideal logarithmic grid). The
independent Python tree checker reproduces every bin with upper endpoint
at most `2^18`. Each ray also obeys the exact identity

    tau_G(2x)=tau_G(x)+1

for these nonzero sources: 2x is not a fixed point and G(2x)=x.
Doubling is checked separately in C++. The four fixed roots are the
stopping set, not an implicit time cap; exceeding the guard aborts.

NUMERICAL scans exhibit substantial lines near 11.5, 24.5 and 60.5 per
octave on both rays. A *hypothetical scalar closure* with independent residue
weights 1/4 and 3/4 would have

    h(u)=-h(u-1)/4-3h(u+log_2(phi)-1)/4,
    exp(-lambda)/4+3 exp((log_2(phi)-1)lambda)/4=-1.  (12)

Both delays are negative. This prevents right-half-plane roots, but it
**does not prevent left-half-plane resonances**. A simple one-delay example
`h(u)=-h(u-1)` already has undamped half-integer frequencies.
Near half-integer f with `f log_2(phi)` close to an integer, (12) has the
same type of small-defect roots as the Collatz law.

The selected predictions and one-line fits over octaves 14 through 22 are:

| base | scalar model frequency | model decay | fitted decay on n | fitted decay on n phi^2 |
|---:|---:|---:|---:|---:|
| 11.5 | 11.474655 | 0.008853 | -0.037174 | 0.010434 |
| 24.5 | 24.513963 | 0.002679 | -0.147578 | -0.000800 |
| 60.5 | 60.502559 | 0.000090 | -0.141306 | -0.037570 |

Negative fitted decay means finite-window growth. These are individual
fits in a crowded spectrum, with large residuals; they do **not** verify
scalar closure or an asymptotic growing mode. Some frequencies are close,
but the ray-dependent amplitudes and decays forbid importing (12) as an
exact law. Output: [golden spectra](golden_spectra_20261013c.json).

The appropriate geometric sidecar comes from both embeddings. Ignoring
bounded digit translations for the leading scale model, put

    u=log_2|x|, v=log_2|x'|,
    z=(u+v)/2, w=(u-v)/2, theta=log_2(phi).

Both branches decrease z by 1. The even branch fixes w; the odd branch
increases w by theta and reverses the sign of x'. Regions with independent
embedding widths therefore have a two-variable shape equation. In its
sign-even sector, a formal mode `exp(lambda z+i kappa w)` satisfies

    exp(lambda)=-[1/4+(3/4)exp(i kappa theta)].       (13)

For kappa=0 it has undamped half-integer modes in z. This is a theorem
about the displayed leading-scale model, **not** a proved counting law on
the integer rays. The odd-sign sector changes the sign of the odd term.
Discarding w silently identifies differently shaped regions.
The exact tree law (9), in contrast, needs no closure hypothesis.

The `3+1` split here has an actual arithmetic source: 2 is inert in the
golden field, so the residue field is F4, with one zero and three nonzero
classes. This provides four inverse branches and the 3/5 amplitude in (9).
It is not an explanation of every occurrence of a triple elsewhere.

## 5. Connections, limits and next targets

- **Conway -> interval controllers.** The map sends an ordered bit word to
  `(length,total,maxprefix,first-location)`. It preserves the maximum query
  by (7), loses the full word, and needs an ordered DAG as its reconstruction
  sidecar. For Collatz, the source address, terminal character and moment
  translation must also remain. The cheap test is exact concatenation on
  the current census; the universal refinement inequality remains open.
- **Golden four-ary tree -> Collatz counting.** The exact rooted inverse
  branching proves a boundary-dominated character law. Collatz inverse
  branching is residue-dependent, and its scale population is not B_K(r).
  Its source weights and missing inverse branches are the essential debt.
  A regular-tree amplitude is a control, not evidence for its universal
  entry. Both systems now have explicit root-forcing terms to retain.
- **Pell -> golden shapes.** Multiplication by phi has conjugate magnitudes
  phi and phi^(-1): norm and area data retain z but move w. Square-triangular
  Pell units in Q(sqrt(2)) have a related reciprocal-embedding mechanism
  but a different trace and lattice. Replacing one by the other would lose
  the integral embedding and the branch translations. The existing stronger
  [Pell/zeta note](pell_zeta_scale_20261011.md) is the route for actual
  continued-fraction and odd-zeta claims; no new irrationality assertion is
  made here.
- **Geometry/Jacobian sidecar.** A norm-preserving unit quotient is not a
  Euclidean isometry, and a scalar norm does not recover two polynomial
  coordinates. Thus (13)'s shape omission is exactly the kind of lost
  information that would invalidate a transfer to unit-distance colouring
  or to the already excluded fixed genus-two carrier. These are limits on
  an analogy, not new bounds on the chromatic number or planar Jacobian.

The useful new helper questions are now specific: can the Collatz signed
cut current admit an ordered compression with a provable dyadic energy
budget; can a two-embedding golden transfer reproduce both ray spectra with
the same operator and different source ports; and which boundary functional
of a nonregular inverse tree replaces the exactly soluble coefficient
`3 zeta/(4 zeta-1)`? The amplitude numerator in the companion note answers
the last question at the level of an exact Mellin functional, leaving its
arithmetic tail as an explicit open obligation.

## 6. Incoming Q covariance work: a genuine two-step Pell mechanism

**Inheritance after the final fetch.** Commits `bba55f5d6f` and `be97fc0dad`
add the [one-sided walk analysis](meta_fibonacci_one_sided_walks_20261013.md)
and [THM-4630, silver-ratio covariance model](../../01-canon/theorems/THM-4630-hofstadter-q-jitter-is-the-silver-ratio-arma-fixed-point-rho1-equals-1-minus-sqrt2.md).
The actual Q sampling assumptions remain **OPEN**. The following is a
**PROVED statement about the explicitly stationary covariance model**, not
an unconditional theorem about Q or Collatz.

Let a parent process have positive variance V and lag-one correlation r in
[-1,1]. Assume its two half-weight source samples are uncorrelated except
for the shared adjacent-index pair. The new source has variance V/2 and
correlation c=r/2 at lag one, zero at higher lags. This covariance is
admissible: its spectral factor `1+r cos(w)` is nonnegative. Pass this
source through the stationary causal stable filter

    x_n=-x_(n-1)/2-x_(n-2)/2+s_n.

The filter roots have modulus 1/sqrt(2). Its covariance equations, with
source variance sigma^2, give

    gamma_2=-(gamma_1+gamma_0)/2,
    3 gamma_0+gamma_1=(4-2c) sigma^2,
    3 gamma_1=-gamma_0+2c sigma^2.

Consequently the model's generation update is exactly

    r_new=F(r)=(r-1)/(3-r),
    V_new/V=(3-r)/4.                                  (14)

This supplies more than a fixed-point equation. F maps [-1,1] into
[-1/2,0] and has derivative `2/(3-r)^2 <= 1/2` on that initial interval.
Every admissible starting correlation therefore converges to
`r_-=1-sqrt(2)`. If `r_+=1+sqrt(2)`, its exact projective form is

    (F(r)-r_-)/(F(r)-r_+)
      =(3-2sqrt(2)) (r-r_-)/(r-r_+).                  (15)

In particular `F'(r_-)=3-2sqrt(2)=r_-^2`, the inverse of the
square-triangular Pell unit. This coincidence has an explicit integral map:
in the coordinate s=1-r, one generation is `s -> 2(s+1)/(s+2)`, and
**two generations** are

    s -> (3s+4)/(2s+3),
    [numerator,denominator]^T -> [[3,4],[2,3]]
                                [numerator,denominator]^T.  (16)

The same matrix multiplies `X+Y sqrt(2)` by `3+2sqrt(2)` and preserves
`X^2-2Y^2`. With X=2n+1 and Y=2m, norm one is precisely
`n(n+1)/2=m^2`. Thus the map from the covariance recursion to the Pell
recursion preserves projective dynamics and a quadratic norm; it does not
make a general covariance state a square-triangular integer solution.
Its needed sidecar is the unnormalised pair (X,Y), including its norm.

The variance also has a closed form, useful as a hostile control against
merely fitting a constant exponent. Put

    A=(1+sqrt(2)-r_0)/(2sqrt(2)), B=1-A,
    v_+=(2+sqrt(2))/4, v_-=(2-sqrt(2))/4.

Then for every k>=0 in this model,

    V_k/V_0=A v_+^k+B v_-^k.                         (17)

Indeed the projective matrix for (14) is `[[1,-1],[-1,3]]`;
its denominator multiplies by `3-r_k` at each step, so V_k is that
denominator divided by 4^k. Diagonalising the matrix gives (17).
A>0 for r_0 in [-1,1]. If generation k is assigned size 2^k and V is the
variance of Q(n)/n-1/2, this model gives the exponent
`1+log_2(v_+)/2` for unnormalised rms fluctuations. Connecting that
generation observer to the actual Q sequence still needs the source
sampling and stationary-envelope assumptions; none follows from (14).

**Exact audit and first failed implication.** The incoming exact index
formulas omitted two finite-size terms. With d(n)=Q(n)/n-1/2 they are

    p/n=1/2+1/(2n)-(1-1/n)d(n-1),
    q/n=1/2+1/n-(1-2/n)d(n-2).                       (18)

n=3 is the minimal hostile to the omitted terms. The main normalized Q
identity already included its correct 3/(4n) source and is unaffected.
The companion [exact checker](../../04-computation/experiments/nested_incoming_bridge_20261013c.py)
verifies (18) and the Q identity at every source through 2^16, the
covariance equations and Pell conjugacy for 33 rational correlations
including both endpoints, and matrix/variance iterates through 1024.
It rejects r=2. Its [saved output](nested_incoming_bridge_20261013c.json)
also independently confirms Q(20)=12, Q(1024)=557 and Q(1536)=808,
where the early `Q(3*2^k)=2^(k+1)` pattern fails.

**What transfers back to the other lanes.** The golden unit has
`phi'=-phi^(-1)`; the silver covariance parameter has
`r_-=-(1+sqrt(2))^(-1)`. Their traces, integral lattices and residue
branching differ. Equations (15)-(16), not a common decimal appearance,
are the precise silver/Pell bridge. For Collatz the analogous productive
question is whether a *source-conditioned* positive covariance or energy
operator closes with an explicitly bounded remainder. Discarding the
remainder or calling deterministic source samples white would reproduce
the unproved step, not bypass it. The exact signed-cut energy in the
companion note is the consequence that such an operator must bound.

The incoming unpinned Conway continuum map eventually drifts to a flat
profile; the discrete dyadic boundary remains pinned. Our ordered Pascal
DAG retains that boundary by construction. Its maximum raw charge M_k
and the incoming maximum relative deviation A_k are different queries:
`d(n)=e(n)/(2n)`, so maximization and division by n cannot be interchanged.
The constants about 0.399 and 0.14 therefore cannot be compared as if they
measured the same quantity. The primary Pinn paper cited in section 1
likewise treats the Q fluctuation exponent empirically; the new model
explains a conditional mechanism, without upgrading those observations.

## Reproduction

```sh
python3 04-computation/experiments/nested_trees_20261013c.py --output /tmp/nested.json
python3 04-computation/experiments/nested_incoming_bridge_20261013c.py --output /tmp/incoming-bridge.json
clang++ -O3 -std=c++17 04-computation/experiments/golden_parity_ray_20261013c.cpp -o /tmp/golden-ray
/tmp/golden-ray 22 /tmp/golden22.json
python3 04-computation/experiments/golden_tree_audit_20261013c.py --ray-census /tmp/golden22.json --output /tmp/golden-tree.json
OPENBLAS_NUM_THREADS=1 python3 04-computation/experiments/golden_spectra_20261013c.py --input /tmp/golden22.json --output /tmp/golden-spectra.json
```
