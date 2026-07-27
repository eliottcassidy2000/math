---
source: codex-2026-07-25-lrc-orientation-current-leakage
status: >
  PROOF CANDIDATE + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Pending joint promotion with its literal-pair dependency. A literal unit-coloured
  pair A-A'=m*c_* from one exact THM-2305 handoff word either has an
  endpoint in the unrestricted exclusive-owner spectrum, or every
  nontrivial deepest cyclotomic address current vanishes separately in
  both boundary orientations. The bad branch imposes at least 336 exact
  lower-field current vanishings when b>=lambda+1, including the strict
  shallow-owner application, and at least 310 in the boundary case
  b=lambda. A projection-overlap inequality and a word-Pluecker
  alternative identify two complementary stopping boundaries. A further
  exact refinement proves that an orientation current split by the
  prime-to-thirteen endpoint residue is nonzero exactly when its fine
  cell contains an exposed endpoint. Raw currents can cancel between
  those residues; a unit-two half-period pairing is an exact hostile
  witness. The later whole-word parity/zero-divisor packet sharpens the
  final handoff to an exact total-current lcm(gcd) criterion and repairs
  separate provenance at common source/pulled-target boundaries. This
  does not exclude the orientation-zero branch or prove LRC(14).
depends_on:
  - THM-2304-deepest-boundary-cyclotomic-current-separation
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - lrc14-literal-word-quadratic-unit-pair-opus-20260725
related:
  - THM-2293-quadratic-root-energy-raises-the-ancestry-shell
  - THM-2302-same-label-expiration-dichotomy-and-pure-terminal-shell-no-go
  - THM-2303-terminal-component-phase-current-and-defect-rank
  - THM-2041-frobenius-stability-of-exact-period-projectors
  - lrc14-word-current-parity-zero-divisor-opus-20260725
---

# Two frequencies split the deepest leakage current by orientation

> **FOLLOW-UP WHOLE-WORD REPAIR.** Sections 8--9 correctly identify the
> affine tomography and provenance-level zero-divisor obstruction, but
> the exact physical handoff must use the **total** word derivative
> current after summing every pulled bank. At a common source/pulled
> boundary the two right-Leibniz pieces can cancel, and the faithful
> criterion is the row-wise `lcm(gcd(F_j,N),gcd(W_j,N))!=N` theorem in
> `07-reflections/lrc14-word-current-parity-zero-divisor-opus-20260725.md`.

## 1. The incidence problem after the literal-pair theorem

Let `E` be a full exclusive-owner step set, and let one exact terminal
word select the literal subset `E_Q`. The quadratic literal-pair theorem
gives

```text
(1_(E_Q))_hat(A)!=0,            (1_(E_Q))_hat(A')!=0,

A-A'=m c_*,                     gcd(m,91)=1.        (1)
```

Both frequencies have exact thirteen-adic grade `lambda`. Write

```text
c_*=13^c u,                     13 does not divide u,

lambda<=b<c,                                           (2)
```

and assume that every other scalar boundary bank has thirteen-depth at
most `b`. The unresolved implication is precisely

```text
literal atoms of E_Q
  -> survival of one atom in the full indicator 1_E.  (3)
```

THM-2304 separates deepest endpoint currents at one frequency. The new
move is to use both frequencies in (1) before discarding their difference.

## 2. Relative conductor basis at arbitrary grade

Let `N` be prime to thirteen and contain all prime-to-thirteen endpoint
denominators. If `b>=lambda+1`, put

```text
L=Q(zeta_(N 13^(b-lambda))),

M=Q(zeta_(N 13^(c-lambda))),

R_rel=[M:L]=13^(c-b).                              (4)
```

At the boundary `b=lambda`, the base field has no thirteen-primary
conductor. Then instead

```text
L=Q(zeta_N),

M=Q(zeta_(N 13^(c-lambda))),

R_rel=[M:L]=12*13^(c-b-1).                         (5)
```

The Euler factor `12` in (5) is load-bearing; replacing it by `13` would
give a false power-basis count. In either case, every lower-bank phase at
`A` or `A'` lies in `L`, while a deepest endpoint phase has a unique
relative-address expansion in an `L`-basis

```text
1,xi,...,xi^(R_rel-1).                             (6)
```

For each nontrivial address `1<=s<R_rel`, split the deepest endpoint current
at frequency `A'` by its geometric comb-endpoint side `epsilon`:

```text
C_(s,+),                    C_(s,-) in L.           (7)
```

Coincident scalar boundaries use THM-2304's right-sided product rule;
bank provenance is assigned before the two orientations are summed.

## 3. The two-frequency determinant

A deepest endpoint has the form

```text
x=(14r+epsilon)/(14c_*),             epsilon=+1 or -1.
```

Using (1),

```text
exp(-2*pi*i*A*x)
 =exp(-2*pi*i*A'*x)
  exp(-2*pi*i*m*c_*x)

 =exp(-2*pi*i*A'*x) theta_epsilon,                  (8)
```

where

```text
theta_+=exp(-2*pi*i*m/14),

theta_-=exp(+2*pi*i*m/14).                          (9)
```

These phases lie in `L` and do not change the relative address `s`.
Suppose both literal endpoints leak completely:

```text
(1_E)_hat(A')=0,                   (1_E)_hat(A)=0.  (10)
```

Uniqueness in the basis (6) gives, for every `1<=s<R_rel`,

```text
C_(s,+)+C_(s,-)=0,

theta_+ C_(s,+)+theta_- C_(s,-)=0.                 (11)
```

The determinant is

```text
theta_--theta_+
 =2i sin(pi*m/7)!=0,                               (12)
```

because `7` does not divide `m`. Therefore

```text
C_(s,+)=C_(s,-)=0         for every 1<=s<R_rel.    (13)
```

This proves the exact dichotomy

```text
one of A,A' survives in the unrestricted owner spectrum

or

all 2(R_rel-1) nontrivial deepest orientation currents vanish.      (14)
```

For the interior depth gap `c-b>=2`, this gives

```text
b>=lambda+1:
  2(R_rel-1)=2(13^(c-b)-1)>=336;

b=lambda:
  2(R_rel-1)=2(12*13^(c-b-1)-1)>=310.              (15)
```

Thus restricted-to-unrestricted leakage is no longer an unstructured
complex cancellation. It is a finite, labelled, orientation-resolved
top-digit boundary obstruction. The strict first-depth-one shallow-owner
application has `lambda=1<b` and therefore uses the first line of (15).
The displayed currents are distinct labelled vanishing equations; no
claim is made that they are algebraically independent constraints on the
underlying scalar-cover parameters.

## 4. Why positive gap spread does not finish (14)

THM-2273/2278-type image spreading gives positive mass in multiple
deepest-safe gaps. That is interior information. Currents (7) live on
the deepest boundary bank.

Already in the abstract scalar-boundary model, two positive intervals
lying strictly inside two different deepest-safe gaps have multi-gap mass
and every local root character, but no deepest boundary endpoint. Every
current in (13) is then zero. Hence

```text
multi-gap activation
  does not imply a nonzero orientation-resolved top-digit current.   (16)
```

The needed sidecar must force boundary occupancy or a nonuniform
top-digit derivative current, not merely positive interior measure.

## 5. A quantitative projection-overlap threshold

There is a general Hilbert-space test for forcing incidence. Let

```text
N=the unrestricted rooted gauge,

W=1_Q N,

E=||W||_2^2,                    T=||N||_2^2,

C=supp(N_hat) intersection supp(W_hat).             (17)
```

Since `<N,W>=||W||_2^2`, Parseval and Cauchy--Schwarz give

```text
sum_(h in C)|W_hat(h)|^2>=E^2/T.                    (18)
```

If no difference-`d` edge in `W_hat` touches `C`, then its entire
autocorrelation coefficient is supported on `C^c`. Another
Cauchy--Schwarz inequality yields

```text
|(|W|^2)_hat(d)|
 <=E-E^2/T
 =E(1-E/T).                                        (19)
```

Therefore a shell coefficient strictly larger than the right side of
(19) forces endpoint incidence. Current LRC arguments give only
nonvanishing, with no remotely comparable lower margin; the inherited
word/full energy ratios are on the order of `10^-7` in the strict branch
and `10^-8` in the repeated branch. Inequality (19) is exact but presently
diagnostic.

## 6. The faithful leakage object is a Pluecker current

The three terminal words partition only the return set, not the whole
circle. Add the complement word, or retain every full clock atom. Pull
each word back to its literal source subset `F_sigma`, so

```text
E=disjoint union_sigma F_sigma.                     (20)
```

At the literal pair `A,A'`, put

```text
v_sigma=(1_(F_sigma))_hat(A),

z_sigma=conjugate((1_(F_sigma))_hat(A')).           (21)
```

If neither endpoint survives, then

```text
sum_sigma v_sigma=sum_sigma z_sigma=0.              (22)
```

Linear algebra gives a second exact alternative:

```text
some v_R z_S-v_S z_R is nonzero,

or

v and z lie on one rank-one leakage eigenline.       (23)
```

The selected word has both its `v` and `z` coordinates nonzero. Hence if
all minors incident to that coordinate vanish, the two full word vectors
are proportional; otherwise (23) supplies a nonzero minor. That minor is
a genuine two-word relative-phase current.
Under source translation by `tau`, it carries exactly the affine shell
phase

```text
exp(-2*pi*i(A-A')tau).                              (24)
```

It is therefore closer to THM-2303's component holonomy than to a
tournament edge. The rank-one alternative is real: THM-2304 alone does
not eliminate it, because the conjugation in (21) compares `A` with
`-A'`, so the useful difference in (1) is not the conductor being
separated.

## 7. The missing coordinate is the prime-to-thirteen endpoint residue

There is an exact positive result behind the stopping boundary in (16).
For Sections 7--8 assume the strict branch

```text
b>=lambda+1,                    R=13^(c-b).          (24a)
```

The boundary case `b=lambda` has the different relative degree in (5) and
is not included in the pure-thirteen cyclotomic polynomial below.
Write

```text
c_*=13^c u_*,                  gcd(u_*,13)=1,

A=13^lambda a_0,               gcd(a_0,13)=1.       (25)
```

Fix one nontrivial thirteen-primary address `s`, one orientation
`epsilon`, and refine its deepest endpoints by their prime-to-thirteen
residue

```text
h in Z/u_* Z.
```

After choosing the address representative, the endpoint indices in this
fine cell have the form

```text
r=r_s+R(h+u_*v),               0<=v<13^b.           (26)
```

Let `g(h+u_*v)` be the zero-one lower-bank gate: guard and unit speeds are
safe, the selected owner is dangerous, and the other blockers are safe.
Up to a fixed nonzero endpoint scalar, the fine current is

```text
C_(s,epsilon,h)
 =kappa_h sum_v g(h+u_*v)
    exp(-2*pi*i a_0 v/13^(b-lambda)).               (27)
```

Put `n=13^(b-lambda)>=13`. Reduce the sum in (27) modulo `n`; its
coefficients are nonnegative integers. If the current vanishes, its
coefficient polynomial is divisible by

```text
Phi_n(X)
 =1+X^(n/13)+X^(2n/13)+...+X^(12n/13).              (28)
```

Consequently, for every `k<n/13`, the 13 coefficients at

```text
k, k+n/13, ..., k+12n/13
```

are equal. Across those 13 columns, the phase of the exact-depth selected
owner

```text
c_j=13^lambda u_j,                 13 does not divide u_j,
```

advances by `u_j/13`; multiplication by `u_j` permutes the thirteen-grid.
Its danger interval can contain at most two of the thirteen equally spaced
digits. At least eleven coefficients are therefore zero, and equality
forces all thirteen to be zero. We have proved the exact equivalence

```text
C_(s,epsilon,h)!=0

iff

the fine cell (s,epsilon,h) contains an exposed c_* endpoint.        (29)
```

When `u_*=1`, there is only one residue and (29) immediately upgrades any
exposed nontrivial deepest endpoint to the nonzero current needed in
(14). For arbitrary `u_*`, the raw current in (7) is a weighted sum of the
fine currents over `h`; it can cancel because it forgot this residue.

That failure is literal. Take `u_*=2`, a selected lower gate
`D_(2*13^lambda)`, and `A=13^lambda a_0` with `a_0` odd. Pair each
same-orientation, same-thirteen-address endpoint `x` with `x+1/2`. The
selected-owner gate is invariant, whereas

```text
exp(-2*pi*i A(x+1/2))=-exp(-2*pi*i Ax).             (30)
```

Thus every raw orientation current can vanish despite abundant exposed
endpoints; its two fine `h=0,1` currents are nonzero opposites. This is an
exact lower-gate hostile mechanism.

It embeds in a primitive nine-speed full scalar boundary trace. Take

```text
lambda=1,        b=3,        c=5,

(c_1,c_2,c_3)=(2*13,2*13^3,2*13^5)
              =(26,4394,742586),

H=2,

(q_1,...,q_5)=(1,1+7c_3,3,3+7c_3,4)
             =(1,5198103,3,5198105,4).              (31)
```

These nine speeds are distinct, have gcd one, and have thirteen-adic
valuation profile

```text
(0,0,0,0,0,0,1,3,5).
```

On every `c_3` endpoint, translation by `1/2` preserves orientation and
thirteen-relative address. The four even banks `H,q_5,c_1,c_2` are
invariant. Each odd pair `q,q+7c_3` swaps, since the frequency shift
`7c_3` becomes exactly a half-turn at a `c_3` endpoint. Hence the complete
lower exclusive-owner gate

```text
C_H intersect (intersection_i D_(q_i)^c)
    intersect D_(c_1) intersect D_(c_2)^c           (32)
```

is half-translation invariant on both deepest boundary orientations. At
every odd frequency, antipodal phases have opposite signs. Exact integer
enumeration finds, for each orientation,

```text
49,880 exposed endpoints,

49,536 exposed endpoints in nontrivial relative addresses,                (33)
```

yet every raw orientation/address current is zero by pairing.

This is a primitive full scalar *boundary-trace* hostile, but its guard is
even. It is not a live odd-guard fixed-section LRC row, an exact THM-2305
literal-word realization, or a counterexample to LRC(14). It nevertheless
rules out any proof based solely on endpoint occupancy, primitivity,
matchings, Frobenius--Koenig transversals, or partial-cube connectivity.
The faithful sidecar is the entire prime-to-thirteen residue vector,
equivalently a Ramanujan/DFT packet, rather than its raw sum.

A component count can now be used without ambiguity. If `L_0` is the
lower gate before intersecting the deepest safe comb, has `K` positive
interval components and mass `mu`, then either oriented deepest grid has
at least

```text
c_* mu-K                                                        (34)
```

active endpoints. The trivial relative address contains exactly
`13^b u_*` endpoints. Hence

```text
c_* mu-K>13^b u_*                                               (35)
```

forces a nontrivial exposed address. Under `u_*=1`, or after retaining the
full residue packet in (27), this produces a nonzero current by (29).
Without that sidecar, (32) alone is insufficient by the half-period
hostile.

## 8. A finite shifted-frequency bank reconstructs the lost packet

The residue loss has a physical Fourier inverse. Retain the literal pair

```text
A-A'=m c_*,                   gcd(m,91)=1,            (36)
```

and for `0<=k<u_*` define the two shifted frequencies

```text
A'_k=A'+k13^b,               A_k=A+k13^b.            (37)
```

Their difference remains exactly `m c_*`. Because `b>=lambda+1`,

```text
nu_13(A'_k)=nu_13(A_k)=lambda                          (38)
```

for every `k`; the shifts cannot cancel the unit digit at grade `lambda`.
Thus the same lower-field/deepest-field separation applies at every
frequency in this bank. Multiplication by the resulting thirteen-adic unit
induces a `k`-dependent permutation `pi_k` of the relative address classes.
It preserves nontriviality but not the literal address label.

For an endpoint in the fine cell `(s,epsilon,h)`, parameterized as in
(26), the shift ratio is

```text
exp(-2*pi*i*k13^b x)
 =eta_(s,epsilon,k) exp(-2*pi*i*k h/u_*),             (39)
```

where `eta_(s,epsilon,k)` is nonzero and independent of the remaining
index `v`. Hence, at fixed orientation, the `u_*` shifted currents are,
up to harmless diagonal factors, the complete discrete Fourier transform
of

```text
(C_(s,epsilon,h))_(h mod u_*).                        (40)
```

For each fixed `k`, comparison of the two rows `A'_k,A_k` multiplies the
two orientations by

```text
theta_+=exp(-2*pi*i*m/14),

theta_-=exp(+2*pi*i*m/14),
```

whose determinant is nonzero by (12). Therefore the combined transform

```text
fine currents indexed by (epsilon,h)

   -> currents at {A'_k,A_k:0<=k<u_*}                 (41)
```

  is, after the `pi_k` relabelling and diagonal `eta` scaling, the tensor
  product of an invertible two-by-two orientation matrix and the invertible
  Fourier matrix of `Z/u_*Z`.

Apply the relative-conductor basis separately at all frequencies in
(37). If every full Fourier coefficient in the `2u_*` bank vanished, then
every nontrivial-address raw current at each frequency would vanish; the
inverse transform (41) would then kill every fine current. By (29), no
fine cell could contain an exposed endpoint. Contrapositively:

```text
an exposed endpoint in a nontrivial deepest address

  =>

some full Fourier coefficient at A'_k or A_k is nonzero.             (42)
```

This is an exact repair of the raw-current hostile. It preserves the
literal difference `m c_*`, the grade `lambda`, relative address
nontriviality (up to `pi_k`), and boundary orientation, while allowing the
common affine frequency translate to vary. The cost is a bank of size
`2u_*`, not a single prescribed pair.

In group-algebra language, (40) is precisely the prime-to-thirteen cyclic
packet that THM-2041 decomposes into exact-period idempotents. THM-2041 can
preserve any nonzero whole packet under good-characteristic Frobenius, but
it does not produce the exposed endpoint or preserve the exact terminal
word. Here (29) supplies production and (41) supplies characteristic-zero
reconstruction. The still-missing handoff is now sharply stated:

```text
transport one nonzero member of the affine bank (37) through the same
THM-2305 terminal word, or prove that unrestricted same-grade service
already suffices without that word.                                  (43)
```

## 9. Word-provenance overlap is a two-row group-algebra zero-divisor problem

Fix one nontrivial deepest bank provenance and one source address `s`
before differentiating, and pull the same `pi_k` address permutation back
in the full and word packets. Do **not** split or renormalize the two
orientations: in each physical row, aggregate them exactly as they occur
and retain every factor `eta_(s,epsilon,k)` inside the current. The word
sequence here is the contribution of this same pulled-back provenance;
jumps from other word banks are intentionally excluded. Thus the objects
below are bank-provenance current packets, not automatically the whole
word Fourier coefficients.

Put

```text
A_cyc=C[X]/(X^u_*-1).
```

For physical row `j=0,1`, corresponding respectively to the frequency
sequences `(A'_k)_k` and `(A_k)_k`, there are unique group-algebra
elements

```text
F_j(X), W_j(X) in A_cyc                                  (44)
```

whose evaluations at `zeta_u^(-k)` are the actual unnormalized full and
exact-word bank-provenance currents in that physical row. This is simply
inverse DFT notation for the two length-`u_*` sequences; it does not assert
that the `eta`-weighted sequence is itself the DFT of a single orientation
residue vector. Since Fourier evaluation identifies the semisimple
algebra `A_cyc` with the product of its character fields,

```text
some physical index (j,k) has nonzero full and word
bank-provenance currents

iff

F_0 W_0 != 0 in A_cyc  or  F_1 W_1 != 0 in A_cyc.       (45)
```

Thus provenance-level overlap is an exact *row-wise* complementary
zero-divisor problem. On the full stratum, THM-2304's unique
deepest-address basis promotes a nonzero `F_j` evaluation to a nonzero
full Fourier coefficient. On the word side, a nonzero `W_j` evaluation
can still cancel against the intentionally omitted jumps from other word
banks. Consequently (45) isolates one necessary algebraic obstruction but
does not by itself solve the exact-member handoff (43).

It is not enough to perform the orientation reconstruction first and
then test coordinatewise products there. The invertible two-by-two matrix
in (41) preserves nonzero vectors, but not overlap of their coordinates.
Already for `u_*=1`, if

```text
M=[[1,1],[theta_+,theta_-]],

f=M^(-1)(1,0),                    w=M^(-1)(0,1),        (46)
```

then both reconstructed orientation-coordinate products are nonzero,
whereas the full physical vector `(1,0)` and word physical vector `(0,1)`
have disjoint support. The earlier one-polynomial criterion after
orientation splitting would therefore have been false.

A clean sufficient condition that is uniform in the row carrying the
full current is

```text
gcd(W_0,X^u_*-1)=gcd(W_1,X^u_*-1)=1.                   (47)
```

Then both `W_j` are units, so whichever nonzero physical row is produced
by (42) shares a shifted provenance current with the exact word. If the
producing row is known, only its unit condition is needed. More generally
it is enough to prove `F_j W_j!=0` in one exact-order idempotent component
for one physical row `j`. A separate no-cancellation or ancestry sidecar
is still required to turn that word-provenance current into a nonzero
whole word coefficient.

For algebraic coefficients, choose a finite place of good residue
characteristic `p` with `p` not dividing `u_*` and with the chosen
nonzero product projection still nonzero after reduction. THM-2041 then
preserves that whole exact-order projection under absolute Frobenius. It
does not preserve a named physical character or terminal word, and it
does not create the characteristic-zero seed in (45).

This formulation identifies the cheapest hostile: construct a legal
terminal word whose inverse-DFT physical-row current polynomial is
idempotent-complementary to the corresponding full current in both rows.
The half-period trace above shows that zero divisors are plausible in the
raw scalar algebra, but it does not realize one on a THM-2305 positive
word. The next decisive finite audit should therefore compute

```text
Res_X(W_0,X^u_*-1),              Res_X(W_1,X^u_*-1)
```

or the row-wise exact-order products `e_d F_j e_d W_j` on the canonical
word packets, with deepest provenance, owner, clock, address pullback,
and physical row still attached. A successful computation must then be
paired with a word-side cancellation audit; a nonzero product alone is
not advertised as the implication (43).

## 10. Information ledger and next target

```text
source:
  a literal unit-coloured pair on one exact handoff word, the full
  exclusive-owner scalar boundary atlas, and a strict deepest depth gap;

map:
  evaluate the full endpoint derivative current at both pair
  frequencies, retain boundary orientation, and compare the two relative
  cyclotomic coordinate vectors;

preserved:
  source owner, literal pair, exact multiplier m, thirteen address,
  boundary orientation, prime-to-thirteen endpoint residue when the fine
  current is retained, nontrivial relative-address shell under the
  tomography permutation, affine shell phase, and full-set incidence;

destroyed:
  interior component geometry and, in the raw zero-current branch, the
  prime-to-thirteen residue vector and every witness away from the
  deepest boundary bank;

strongest survivor:
  the endpoint-survival / 2(R_rel-1)-orientation-vanishing dichotomy
  (14), the fine-cell endpoint equivalence (29), and the exact affine
  tomography implication (42);

cheapest decisive test:
  prove both inverse-DFT physical-row word current polynomials are units
  as in (47),
  identify the row carrying the full current and prove its word
  polynomial is a unit, or directly prove F_j W_j!=0 in one exact-order
  component; then prove that the surviving provenance current cannot
  cancel against the remaining word-bank jumps.                     (48)
```

If that residue-faithful boundary-occupancy test succeeds, (14)
immediately lands one literal endpoint in the unrestricted owner spectrum.
No scalar profile is yet excluded, and LRC(14) remains open.

## 11. Exact reproduction

The integer-only hostile audit is

```text
04-computation/lrc14_prime_to13_current_hostile_probe.py

05-knowledge/results/lrc14_prime_to13_current_hostile_probe.out.          (49)
```

It checks distinctness, primitivity, every valuation, the half-turn gate
identity at all `4*13^5` deepest endpoints across both orientations, both
counts in (33), and the
odd-frequency phase reversal. Normal and optimized Python executions
must match the stored transcript before this packet is promoted. The
independently replayed LF-normalized hashes are

```text
script:
bfe2ec8a816528a7d3cea4a8b02280b5f7792456b79fa14b2c18839a29117c2a

output:
655b3cd9f91c8c12bbfa5a8031a867a4f805f060a4b65698f019fcdb9e7e34aa.
```
