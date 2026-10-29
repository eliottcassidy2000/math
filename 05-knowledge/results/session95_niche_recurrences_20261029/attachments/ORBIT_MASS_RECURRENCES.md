# Source-mass recurrences: replicas, cyclic periods, and a six-way obstruction

**PROVED (specified cyclic wreath actions). FINITE-EXACT controls.**
Source pin: `d3741fe6ed02f6072dbe6bfda2f3d04193d75ce0`.
This note derives recurrences for the *probability law of an orbit*, rather
than uniformly sampling isomorphism classes. No priority beyond the cited
repository comparison, unrestricted A000568 formula, or new graph gadget is claimed.

## Inheritance and the less-used mechanism

The [recent canonical attachment carrier](../../weighted_crt_clocks_20261028/CANONICAL.md)
provides cyclic child keys and the centralizer's paired product palette.
[Arbitrary weighted ports](../../crt_weighted_ports_20261028/WEIGHTED_PORTS.md)
separates distinct orbit mass from weighted Burnside's orbit-average mass.
Session94's cycle index counts weighted invariant attachments. Those are
inputs, not the new conclusions below.

A niche predecessor is the [old two-palette gate gluing](../../first_exit_complement_games_20261020/PALETTE-GLUE.md):
two copies sharing a *named* output have counts
`sum_c f(c)g(c)`, with the displayed 24-per-colour fibres producing 576.
[The later fibre decoder and leakage audit](../../charge_palettes_boundary_20261020g/palettes.md)
retains the externally conditioned boundary and explains why an unconditional
colour-equality gate does not follow. This is the useful mechanism to recover:
shared-label gluing is a moment/collision contraction, not multiplication of
independent marginal totals. The present work applies that contraction to
canonical orbit keys, with an explicit cyclic-period correction.

Concept board: (i) source mass versus uniform class count; (ii) replicated
shared-label gluing; (iii) cyclic word period; (iv) paired-port diagonal action;
(v) six-way alignment hidden from pairwise overlaps.

## 1. The full scalar hierarchy

Let `G_(d,0)=1` on a singleton and
`G_(d,k+1)=G_(d,k)^d semidirect C_d`. Colour its `d^k` leaves independently
with the *same* finite palette law `p=(p_a)`, restricted to its positive support.
Let `mu_k(O)` be the probability of orbit O under this labelled source and put

    S_r(k;p)=sum_O mu_k(O)^r,     r=0,1,2,... .                 (1)

Thus S0 is the number of supported orbit classes, S1=1, and S2 is the
probability that two independent samples have the same orbit. Base:
`S_r(0;p)=sum_a p_a^r`. At r=0 each positive mass contributes one.

For every integer d>=2 the exact recurrence is

    S_r(k+1;p)=sum_(ell|d) ell^(r-1)
        sum_(e|ell) mobius(e) S_(r*d*e/ell)(k;p)^(ell/e).      (2)

The factors with r=0 are rational but the final result is the ordinary integer
necklace count. For ternary tournaments the particularly short form is

    S_r(k+1)=3^(r-1) S_r(k)^3
                   -(3^(r-1)-1) S_(3r)(k).                   (3)

For nonternary d this is a theorem about the specified permutation action;
it does not identify the full automorphism group of an arbitrary tournament.

### Proof, including the non-free words

A word of d child orbit keys has cyclic orbit size ell exactly when its
minimal period is ell. Write it as d/ell repetitions of a primitive ell-word.
Every rotation has the same probability because the child laws are identical.
Its orbit mass is therefore `ell * product_(i<ell) mu(a_i)^(d/ell)`.
Summing the rth powers over distinct orbits divides the primitive-word sum
by ell. Weighted period inversion gives

    sum_(primitive ell-words) product_i q_(a_i)
      =sum_(e|ell) mobius(e) (sum_a q_a^e)^(ell/e).

This identity follows by partitioning all words by their least period and
ordinary divisor Möbius inversion; no independence beyond the stated child
product law is inserted. Substitute `q_a=mu(a)^(rd/ell)` to obtain (2).
For d=3 the only periods are 1 and 3, yielding (3).

An equivalent ternary derivation is especially informative: summing all
ordered child triples at multiplier `3^(r-1)` overcounts the constant triples.
The missing correction is precisely their `S_(3r)` mass, not a vertex count.
At r=0 this becomes `(S0^3+2S0)/3`; at r=1 it becomes normalization.

For a single requested r and ternary depth k, memoization visits at most
`(k+1)(k+2)/2` pairs `(depth, r*3^j)`. Thus it needs O(k²) arithmetic
operations and stored scalars, without building the orbit-key bank. This
charges arbitrarily large rational numerators and denominators: their bit
length can grow proportionally to `3^k`. It is not an O(k²) bit-time claim.

## 2. Shared mirrored ports: the product palette is still compulsory

Let C_k be the centralizer of digit negation in the ternary action, and let
`T_r(k;p)` be the corresponding orbit-mass power sum. The inherited canonical
key separates the zero child and the *single* paired array
`h(v)=(f_1(v),f_2(alpha(v)))`. Under iid leaf labels this array has iid palette
law `p tensor p`. The two factors have disjoint source coordinates, so

    T_r(0;p)=sum_a p_a^r,
    T_r(k;p)=T_r(k-1;p) S_r(k-1;p tensor p),
    T_r(k;p)=S_r(0;p) product_(j=0)^(k-1) S_r(j;p tensor p).  (4)

At the base, `S_t(0;p tensor p)=S_t(0;p)^2`, so (2)-(4) evaluate the hierarchy
using the initial source moments without enumerating a product palette.
Replacing the paired factor by two independent orbit quotients is invalid.

For fair binary labels, collision probabilities at depths 0,1,2,3 are:

| depth | A_k | C_k |
|---:|---:|---:|
| 0 | 1/2 | 1/2 |
| 1 | 5/16 | 1/8 |
| 2 | 2635/32768 | 23/4096 |
| 3 | 1685266004645/1125899906842624 | 26822807/17592186044416 |

In particular C2 has 192 classes, but its source collision probability is
23/4096, **not** 1/192. Uniform labelled input does not yield uniform orbits.

## 3. Exact information saving and its biased-source boundary

Let H denote Shannon entropy with any fixed log base. Because iid source
mass is constant inside each cyclic orbit, conditioning on that orbit is
uniform on its rotations. For ternary child-key law mu,

    H_(A,k+1)=3 H_(A,k)-(1-S_3(k)) log 3.                    (5)

The correction is the expected log orbit size: it vanishes on constant
triples and is log3 otherwise. For the centralizer,

    H_(C,k;p)=H_(C,k-1;p)+H_(A,k-1;p tensor p).              (6)

For fair binary labels at depth2 these give
`H_A=9-(201/64)log_2(3)` and `H_C=9-(15/16)log_2(3)` bits.
These are expected coding entropies, not fixed-width memory bounds or a
claim that transporters are free.

For two positive palette laws p,q define the mixed hierarchy

    M_(r,s)(k;p,q)=sum_O mu_p(O)^r mu_q(O)^s.

The same proof replaces r in the multiplier of (2) by r+s, and replaces
`S_(r*d*e/ell)` by `M_(r*d*e/ell,s*d*e/ell)`. In particular

    M_(r,s)(k+1)=3^(r+s-1) M_(r,s)(k)^3
                     -(3^(r+s-1)-1) M_(3r,3s)(k).          (7)

The base is `sum_a p_a^r q_a^s`; the corresponding centralizer version uses
paired laws `p tensor p` and `q tensor q` in (4).

Apply this to the inherited global-converse involution on binary attachment
orbits, which has no fixed orbit. Under Bernoulli(p), its transported source
is Bernoulli(1-p). If Q2 is the collision probability after pairing converse
orbits, then exactly

    Q2=T2(p)+M_(1,1)^C(p,1-p).                              (8)

Indeed each two-element pair contributes `(a+b)^2`. For p=1/2, the masses
are equal and Q2=2T2: the reversal quotient loses exactly one Shannon bit.
For biased sources the fixed-point-free selector still needs one *worst-case*
bit, but its expected information need not be one bit. At depth0 the loss
is merely the binary entropy h(p). This separates structural encoding from
source-weighted information.

## 4. Nonidentical child sources force a six-way alignment term

Equations (2)-(7) assume identical child laws. Here is the exact repair for
three arbitrary independent child laws mu_0,mu_1,mu_2 and nu_0,nu_1,nu_2 on
a common child-key alphabet. Let

    B_ij=sum_a mu_i(a)nu_j(a),
    J6=sum_a product_i mu_i(a)nu_i(a).

Then the probability that their cyclic triple keys agree is

    sum_(s=0)^2 product_(i=0)^2 B_(i,i+s) - 2 J6.            (9)

There are three alignment events. Each has the displayed product probability.
Any intersection of two distinct alignments forces *all six child keys* equal;
the triple intersection is that same event. Inclusion-exclusion therefore
subtracts 3J6 and adds J6. This proves (9) without an invariance assumption.
The six-way term cannot in general be reconstructed from the nine B_ij.

**Exact hostile.** Use alphabet F2^5. Let the first five laws be uniform on
`x_i=0`, each with 16 points. The sixth is uniform on
`x_1+...+x_5=b`, b=0 or1. Split them into the first three and last three laws.
Every B_ij is 1/32 in both cases. But the common six-way intersection contains
one point for b=0 and none for b=1. Hence the two output collision probabilities
are respectively `1534/16^6` and `1536/16^6`. Pairwise overlap data agree while
the assembled observable differs. The checker enumerates the two weighted
cyclic-key dictionaries independently to verify both results.

This connects back to the old named-output gates: two-copy equality uses a
pairwise Gram contraction, but composing several alignment alternatives needs
a higher shared-output tensor. A physical unit-distance gadget does not become
available merely because its formal tensor has been written down. Nor does a
moment or entropy recurrence retain named transporters or arbitrary future
vertex-specific operations.

## Verification and loss ledger

Run `python check_orbit_mass.py` and `python -O check_orbit_mass.py`.
The checker compares (2) with literal cyclic word dictionaries for d=2..9
and two-/three-colour rational laws subject to the declared 25,000-word cap;
compares (3),(4),(7),(8) with all binary attachments through depth2 for fair
and biased sources; and checks the 32-symbol six-way hostile by direct key
mass accumulation. It also checks the entropy-loss coefficients by literal orbit sizes. It performs 256 active checks; optimization does not
remove them. Depth3 table values are recurrence outputs, not a 2^27 census.

Source: labelled independent palette laws. Target: exact orbit-mass moments
and entropy, under the specified wreath or centralizer action. Map: canonical
cyclic/product-palette keys. Preserved: source probabilities, shared diagonal
pairing, equality-of-key collision observables. Lost: individual samples,
named transporters and general source correlations. Required sidecar for
heterogeneous composition: joint overlaps such as J6, or the full coupled
source tensor. No integer clock is altered or identified by this construction.
