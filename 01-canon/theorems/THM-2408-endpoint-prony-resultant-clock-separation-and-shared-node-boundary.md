---
id: THM-2408
title: "Endpoint-Prony resultant separation and the shared-node cancellation boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In a common
  coefficient Hilbert space, if U=O+H, then O and U cannot both erase H:
  their squared norms sum to at least one half of the squared norm of H,
  sharply. THM-2407 supplies this identity in one lawful
  owner/source-deletion gauge, giving an unconditional quantitative branch
  alternative. For two finite exponential sequences, cancellation is
  confined exactly to common nodes with opposite coefficient vectors. A
  node of H transverse to the O-node polynomial forces U to survive in
  every bounded consecutive window; coprime node polynomials give the r+s
  bound. Applied to common-gauge finite step amplitudes, the nodes are the
  endpoint-Prony nodes in one Fourier residue progression. A strictly
  positive rational owner can nevertheless cancel a nonnegative deletion
  packet at every charged colour when their endpoint nodes coincide. Canon
  does not yet separate the relevant endpoint nodes or identify a
  Fourier-lift index with a terminal clock. No row exclusion or LRC(14)
  conclusion is claimed.
source: codex-2026-07-26-endpoint-prony-resultant-clock-separation
depends_on:
  - THM-2286-endpoint-prony-lift-bank-and-sharp-owner-multiplier-landing
  - THM-2403-clean-toothpick-unequal-slope-target-axis-imbalance
  - THM-2407-owner-or-source-deletion-target-current-dichotomy
related:
  - THM-1710-tnc-cyclotomic-refuted-resultant-replaces
  - THM-2022-gmc2-frobenius-lowest-balanced-face
  - THM-2241-monic-transverse-response-depth-and-resultant-nonproper-quotient
  - THM-2285-centered-grid-footprint-and-generic-keller-lines
  - THM-2370-deletion-martingale-drift-conservation-and-sharp-clone-hostile
  - THM-2380-cross-word-charged-target-correlation-and-pair-twist-gate
  - THM-2383-polarized-complete-subcube-gram-tomography
  - THM-2398-prime-cyclic-rational-restoration-dichotomy
script: 04-computation/lrc14_endpoint_prony_resultant_clock_separation_thm2408.py
output: 05-knowledge/results/lrc14_endpoint_prony_resultant_clock_separation_thm2408.out
script_sha256: 20fe75c8c383880d99714806125b57db9dd1c37a39a570168e2e34f596f0d9c4
output_sha256: 342db54089f49f26e6a53358d8c4c46b1274bc3f22f52229e0ce8b19dfe81677
hash_basis: working-tree bytes (LF)
---

# THM-2408 -- endpoint-Prony resultant separation

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2403 produces a lawful fully-all-safe current `H_+` with charged
target/deep coefficients. THM-2407 places it in the exact target-equivariant
Boolean split

```text
U=O+H_+,
```

where `O` is the genuine positive source-owner packet and `U` has exactly
that stationary source factor deleted. There are two different questions:

```text
can the charged difference disappear from both physical branches?

can the restored branch itself cancel at every lift of one residue?       (1)
```

The first has an exact Hilbert-space answer: no. The second is a finite
Prony/resultant question.  Once the two branches are put in one oriented
endpoint gauge, cancellation can persist only on shared exponential nodes,
and there only by exact coefficient opposition.

The phrase **one oriented endpoint gauge** is load-bearing throughout.
THM-2407 supplies it for `O,H_+,U`; it does not supply the transverse
endpoint node required later to force `U` itself.

## 1. A charged difference survives in one of its two endpoints

Let `V` be a real or complex Hilbert space and let

```text
O,H,U in V,                         U=O+H.             (2)
```

Then the parallelogram identity at the midpoint gives

```text
||O||^2+||U||^2
 =2||O+H/2||^2+||H||^2/2
 >=||H||^2/2.                                          (3)
```

Consequently

```text
max(||O||^2,||U||^2)>=||H||^2/4.                      (4)
```

For one scalar coefficient,

```text
max(|O|,|U|)>=|H|/2.                                  (5)
```

Both constants are sharp at

```text
O=-H/2,                         U=H/2.                 (6)
```

The statement remains true after any fixed orthogonal projection `Pi`:

```text
Pi U=Pi O+Pi H

implies

||Pi O||^2+||Pi U||^2>=||Pi H||^2/2.                 (7)
```

Thus target, deep-colour, fixed-triangle, owner-word, or finite-packet
coordinates may first be frozen and only then may (3) be applied.  No
phase measurement is needed once the coefficientwise identity in (2) is
lawful.

For the eligible THM-2403 transform, equation (44b) there gives

```text
||Pi H_+||^2
 >=9 rho_R^2/250994068.                              (8)
```

THM-2407 constructs physical currents `O,U` in the same normalized
transform and common gauge and proves

```text
U=O+H_+,                                              (9)
```

Therefore one of them obeys

```text
||Pi O||^2 or ||Pi U||^2
 >=9 rho_R^2/1003976272.                             (10)
```

Likewise THM-2403's finite-colour coefficient floor (44c) and (5) give
some eligible target/deep transform coefficient in one of the two
branches with magnitude at least

```text
rho_R/456976.                                        (11)
```

The nonzero coefficient in (11) has an exact `m`-then-`X` triangle
descendant in the selected branch.  No magnitude floor for that exact
triangle is claimed: THM-2403 (44c) is a floor only before this absolutely
convergent expansion.

Thus (10)--(11) are unconditional quantitative alternatives for the
THM-2407 packets. The bank selected by these norm bounds need not be the
same bank selected by THM-2407's all-twelve-colour dichotomy. In particular:

```text
Pi O=0       implies Pi U=Pi H_+;

Pi U=0       implies Pi O=-Pi H_+.                   (12)
```

A positive owner which is genuinely target-flat in the same gauge therefore
cannot cancel `H_+`; cancellation in the restored branch instead transfers
the exact opposite charged coefficient to that owner.

## 2. Finite exponential sequences

Let `V` now be a finite-dimensional complex vector space.  Let

```text
O_n=sum_(alpha in A) o_alpha alpha^n,

H_n=sum_(beta in B) h_beta beta^n,                   (13)
```

where:

- `A,B` are finite sets of distinct nonzero complex numbers;
- every displayed coefficient vector is nonzero; and
- `n` ranges over all integers.

Put

```text
U_n=O_n+H_n.                                         (14)
```

Extend missing coefficients by zero and define the residual coefficient
at a node `z in A union B` by

```text
u_z=o_z+h_z.                                         (15)
```

Let

```text
R={z in A union B:u_z!=0},

L=|R|.                                               (16)
```

Then

```text
U_n=sum_(z in R)u_z z^n.                             (17)
```

If `L>0`, the sequence cannot vanish at `L` consecutive integers.  Indeed,
after shifting the first index to zero, `L` consecutive vector equations
have scalar coefficient matrix

```text
[z^j]_(0<=j<L,z in R),                               (18)
```

whose Vandermonde determinant is nonzero.  Tensoring this invertible matrix
with the identity on `V` forces every `u_z` to vanish, contrary to (16).
Hence:

> **Finite-node window theorem.** If `U` is not the zero sequence, every
> `L` consecutive indices contain an `n` with `U_n!=0`.  The coarser
> a priori window `|A|+|B|` is always valid.           (19)

The exact permanent-cancellation locus is

```text
U_n=0 for every n

iff

A=B and o_z=-h_z for every z in A.                   (20)
```

Equivalently, if `U` vanishes on `|A union B|` consecutive indices, then
the two effective node sets agree and all node coefficients are opposite.
This converts an infinite cancellation claim into one finite Vandermonde
test.

## 3. Shift-operator and resultant separation

Let `E` be the forward shift,

```text
(EF)_n=F_(n+1),
```

and define the monic node polynomials

```text
P_O(T)=product_(alpha in A)(T-alpha),

P_H(T)=product_(beta in B)(T-beta).                  (21)
```

Since `(E-alpha)alpha^n=0`,

```text
P_O(E)O=0.                                          (22)
```

Applying the same operator to (14) gives the exact transverse-node
projection

```text
P_O(E)U_n
 =sum_(beta in B minus A)
     h_beta P_O(beta) beta^n.                        (23)
```

Let

```text
r=|A|,

t=|B minus A|.                                      (24)
```

If `t>0`, the right side of (23) is a nonzero `t`-node sequence.  Among
any `t` consecutive values of `P_O(E)U` one is nonzero, and that value is
a linear combination of `r+1` consecutive values of `U`.  Therefore:

> **Transverse-node theorem.** If `H` has a node absent from `O`, every
> block of `r+t` consecutive indices contains an `n` with `U_n!=0`.
>                                                                  (25)

The usual resultant packages the cleanest case:

```text
Res(P_O,P_H)!=0

iff

A intersect B=empty.                                (26)
```

If (26) holds and `H` is nonzero, then `U` is nonzero and every block of

```text
r+s,                         s=|B|,                  (27)
```

consecutive indices contains a survivor.

The converse is deliberately not claimed.  A zero resultant means that
some endpoint node is shared, not that cancellation occurs.  Equations
(15) and (20) say exactly what more is required: every effective node must
be shared and every coefficient vector must be opposite.

## 4. The bound is sharp in the exponential class

Let `z_1,...,z_L` be distinct nonzero scalars and put

```text
P(T)=product_(j=1)^L(T-z_j),

c_j=1/P'(z_j),

S_n=sum_(j=1)^L c_j z_j^n.                          (28)
```

Lagrange interpolation and comparison of the coefficient of `T^(L-1)`
give

```text
S_0=...=S_(L-2)=0,

S_(L-1)=1.                                          (29)
```

Thus `L-1` consecutive zeros are possible, and the window length in (19)
cannot be reduced using only the number of distinct nodes.

## 5. Common-gauge endpoint-Prony realization

Let `o,h,u:R/Z->V` be finite step functions in one fixed oriented circle
coordinate, and assume

```text
u(x)=o(x)+h(x)                     almost everywhere. (30)
```

Fix

```text
m>=2,

1<=k<=m-1,

q_n=k+mn.                                            (31)
```

For a step function `f`, let `Delta f(x)` be its oriented jump at `x`.
Direct integration by parts gives, because `q_n!=0`,

```text
A_f(n)
 :=2 pi i q_n f_hat(q_n)
 =sum_x Delta f(x) exp(-2 pi i q_n x).              (32)
```

Group jumps with the same node

```text
z=exp(-2 pi i m x)
```

and define

```text
gamma_f(z)
 =sum_(x:exp(-2 pi i m x)=z)
    Delta f(x) exp(-2 pi i kx).                     (33)
```

After zero groups are discarded,

```text
A_f(n)=sum_z gamma_f(z)z^n.                          (34)
```

Equation (30), in the same orientation, gives

```text
gamma_u(z)=gamma_o(z)+gamma_h(z)                    (35)
```

at every node.  Sections 2--3 therefore apply verbatim to the scaled
Fourier coefficients of `o,h,u` in the residue class `k mod m`.

Let `r,s` be the effective endpoint-node counts of `o,h`.  If their node
polynomials are coprime and the `h` progression is nonzero, then every
`r+s` consecutive lifts contain

```text
q congruent k mod m,

u_hat(q)!=0.                                        (36)
```

Using the centered block

```text
n=-floor((r+s)/2),...,
  ceil((r+s)/2)-1
```

gives the explicit bound

```text
0<|q|<=m ceil((r+s)/2)-1.                           (37)
```

This is the two-packet counterpart of THM-2286's one-packet
endpoint-Prony lift.  Rational endpoints make all nodes cyclotomic and
the resultant an exact algebraic-number computation, but rationality is
not needed for the Vandermonde proof.

The word **common** in (30) cannot be weakened:

- `o,h,u` must use the same circle coordinate;
- their endpoint orientations must agree;
- the same residue progression `k+mn` must be used; and
- the identity must hold before endpoint grouping and Fourier projection.

Separate spectra in unrelated gauges do not define the nodewise sum (35).

## 6. Two sharp positivity hostiles

Positivity and rationality do not replace endpoint-node separation.
On `C_p`, with normalized Fourier transform, put

```text
epsilon=1/(2p),

H=epsilon delta_0,

O=1/p-H,

U=O+H=1/p.                                          (38)
```

Here `O` is strictly positive, all three arrays are nonnegative rational,
and, for every `k!=0`,

```text
H_hat(k)=1/(2p^2),

O_hat(k)=-1/(2p^2),

U_hat(k)=0.                                         (39)
```

At `p=13` the two nonzero values are `+1/338` and `-1/338`.

There is an endpoint-level version.  On the circle take

```text
h=1_[0,1/13),

o=2-h,

u=2.                                                (40)
```

Then `o` is strictly positive, `h` is nonnegative, and all functions have
rational values and endpoints.  The nonconstant jumps of `o` and `h`
occur at the same two endpoints with opposite signs.  Hence

```text
o_hat(q)=-h_hat(q),          u_hat(q)=0              (41)
```

for every nonzero `q`; moreover `h_hat(q)!=0` whenever `13` does not
divide `q`.  This is exactly the shared-node/opposite-coefficient locus
in (20), simultaneously on all twelve nonzero mod-thirteen residue
classes.

The hostile also shows why applying Perron--Frobenius positivity,
prime-cyclic rationality, or a scalar energy inequality separately to
the two branches cannot orient their sum.

## 7. Relation to the unrelated algebraic mechanisms

The transfer from other parts of the repository is exact but conditional.

1. **NC2/GMC(2) lowest face.** Introduce a formal branch marker

   ```text
   F(lambda)=O+lambda H.
   ```

   The coefficient of `lambda` is an exposed face and cannot cancel as a
   polynomial.  THM-2022 preserves such a face because its finite-place
   operation retains the whole face sum.  Evaluating at `lambda=1` erases
   the marker, and (38)--(41) show that cancellation can then be exact.
   LRC would need two lawful marker values or a filtration which retains
   the branch degree.

2. **Centered grid/Nullstellensatz.** Since `F` has degree one, any two
   distinct `lambda` values contain at least one nonzero evaluation when
   `H!=0`.  This is the one-variable footprint mechanism of THM-2285.
   It proves the owner-or-restored-current alternative (3)--(5), not
   survival at the prescribed physical value `lambda=1`.

3. **JC/resultants.** THM-2241 detects the precise locus on which target
   sheets meet at infinity by a resultant coefficient.  Here
   `Res(P_O,P_H)` detects the precise locus on which endpoint recurrences
   can interact.  A nonzero resultant separates the branches; a zero
   resultant records only a possible collision and still needs the
   coefficient equalities in (20).

4. **Gram/polarization.** THM-2383 reconstructs signed coefficients from
   oriented references.  Equation (3) is cheaper because the actual
   linear identity `U-O=H` already supplies the orientation.  Without
   that common-gauge identity, separate squared banks again lose the
   relative phase.

Thus the reusable algebraic mechanism is:

```text
common-gauge coefficientwise branch identity
  -> one of the two endpoints retains charged energy;

plus one transverse endpoint node
  -> the restored endpoint itself survives in a bounded lift window.    (42)
```

## 8. LRC typing and stopping boundary

Apply the theorem to THM-2407 after fixing a displayed eligible target/deep
coefficient, or the entire eligible orthogonal bank. That theorem supplies
physical step amplitudes

```text
o,h_+,u
```

satisfying the first three properties below.  The fourth is the additional
separation hypothesis needed only for the bounded-lift conclusion:

1. `h_+` is the selected THM-2403 coefficient in its actual root/target
   and endpoint orientation;
2. `o` is the preselected positive owner in that same orientation;
3. `u=o+h_+` before Poisson--Abel collapse and Fourier grouping; and
4. **additional hypothesis:** in the chosen Fourier residue progression,
   the terminal operation supplies one effective endpoint-Prony node of
   `h_+` transverse to the owner node polynomial.

THM-2407 proves 1--3 before finite target/deep projection, so Section 1
already gives the nonzero charged branch and the floors (10)--(11).
Condition 4 is the remaining hypothesis for Sections 2--5 to force a
bounded exact lift specifically in the restored branch. THM-2305 retains
the owner and word but does not separate their effective endpoint nodes;
THM-2380 names the endpoint-matched cross observable but does not construct
this resultant sidecar.

Also, THM-2286's index `n` is a Fourier-lift progression.  It becomes a
family of physical terminal clocks only if an additional covariant
clock/endpoint identification proves that equivalence.  The theorem does
not silently make that identification.

Accordingly this candidate is a rigorous algebraic separator and a sharp
no-go, not a scalar-row exclusion.  The LRC ledger remains `165`, and
LRC(14) remains open.

## 9. Exact companion

The dependency-free companion:

- exhausts Gaussian-integer controls for the sharp Hilbert constants;
- checks the THM-2403 divided constants exactly;
- exhausts small rational node banks and the consecutive-zero theorem;
- verifies the shift-operator transverse-node identity;
- verifies nonzero and zero resultants by exact rational arithmetic;
- constructs the sharp `L-1`-zero Lagrange sequence;
- checks the `C_13` strictly-positive rational cancellation hostile; and
- checks the representative `m=13`, length-five centered arithmetic used
  to illustrate the general proved frequency bound.

Run

```bash
python3 04-computation/lrc14_endpoint_prony_resultant_clock_separation_thm2408.py
python3 -O 04-computation/lrc14_endpoint_prony_resultant_clock_separation_thm2408.py
```

Both transcripts must byte-match, after LF normalization,

```text
05-knowledge/results/lrc14_endpoint_prony_resultant_clock_separation_thm2408.out
```

Every truth-bearing check raises explicitly, so optimized mode executes
the same audit.

## 10. Independent hostile audit

An independent audit rederived the Hilbert constants, vector-valued
Vandermonde and shift-operator arguments, exact window lengths, endpoint
sign, centered bound, sharp Lagrange construction, and both cancellation
hostiles.  It also enforced the distinction after (11): the displayed
magnitude belongs to the finite transform coefficient, while its exact
triangle descendant is asserted only to be nonzero.  Normal, optimized,
and stored transcripts match across sixteen LF lines with the hashes in
the frontmatter.  The full documentation check passes.
