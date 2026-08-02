---
id: THM-3133
title: "Common-simple-zero Faber exclusion and the sharp odd-bipole response boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  In every normalized nonsplit polynomial exact-square-prefix Faber bank of
  reduced index R>=7, a common simple zero of the response polynomials V and
  M is impossible.  This removes THM-2831's d_F=s_F=0 hypothesis.  Consequently
  every balanced survivor in this range has no simple response zeros; the
  THM-2827 resonance then forces R=3k+2 and genuine nonsplit degree N>=2(4k+3).
  The degree bound is abstractly sharp: an explicit odd-bipole family realizes
  pole passport (d,d) and has a first-order nonvanishing coefficient
  recurrence.  The planar Jacobian conjecture and chart entry remain open.
audit: >
  An independent derivation rebuilt the Faber coefficients from the defining
  differential equation, checked the polar, pure-q, and regular-H lanes with
  every retained lower row, and verified the balanced parity consequence.
  A separate all-degree algebra audit proved the odd-bipole ODE, leading
  coefficient, clean passport, constant recurrence, and d=11 Nielsen count;
  an independent script replayed R=7..18 and odd d through 25.  The canonical
  pole-rotation orbit is now gated rather than merely printed.  Normal,
  optimized, and stored transcripts and hashes agree.
source: root/frontier-synthesis-2026-08-02
depends_on:
  - THM-2760-exact-prefix-even-faber-flux-gcd-and-smooth-boundary-exclusion
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
  - THM-2827-uniform-pole-order-faber-nonresonance-atlas-and-double-pole-exclusion
related:
  - THM-2831-cyclic-resonance-simple-zero-and-full-lower-bank-exclusion
  - THM-2840-heptic-e3-2221-chebyshev-accessory-classification
  - THM-3123-heptic-e3-remaining-accessory-classification-and-s7-monodromy
  - MISTAKE-317
script: 04-computation/jc_common_simple_zero_odd_bipole_thm3133.py
output: 05-knowledge/results/jc_common_simple_zero_odd_bipole_thm3133.out
script_sha256: 7c2de38fba53b9942e46bf3e29daae971f3229e4cbad5b9920424e365454dd16
output_sha256: 32d8a794ce2f25ff969de77e4c3419883cd9aefb04668d7c10baba5bb89836c3
hash_basis: LF-normalized bytes
---

# THM-3133 -- common-simple-zero Faber exclusion and the sharp odd-bipole boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This theorem strengthens THM-2831 at the exact place where the planar
Jacobian response atlas still admitted formal resonances.  The new local
lemma uses all three valuation lanes and no longer specializes the source
coordinates `d_F,s_F` to zero.  It then identifies the sharp abstract
balanced-response boundary, rather than mistaking the obstruction for an
unconditional nonexistence theorem.

## 1. Common-simple-zero exclusion

Fix a normalized nonsplit polynomial exact-square-prefix Faber bank with

```text
R>=7,
Q=E_(4R-2)+sum_(j=1)^(R-2) a_j E_(4j-2).                     (1)
```

Use the THM-2827 source coordinates

```text
T_F=A_src^2/V,
d_F=C_0-B_src^2/(4V),
s_F=A_src B_src/(2V)-E_0.                                  (2)
```

The normalized flux equations are

```text
phi_Q=0,              psi_Q in C,
K_Q=T_F H_Q,          A_src K_Q=lambda M,       lambda!=0,    (3)
```

where the exact Faber extraction gives

```text
H_j in Q[T_F,s_F,T_F d_F],
H_Q=H_R+sum_(j=1)^(R-2)a_jH_j.                               (4)
```

> **Common-simple-zero lemma.**  There is no finite point `alpha` such that
>
> ```text
> ord_alpha(V)=ord_alpha(M)=1.                               (5)
> ```

### Proof

Write

```text
a=ord_alpha(A_src)>=0,       b=ord_alpha(B_src)>=0,           (6)
```

allowing `b=infinity` if `B_src=0`.  The order of `A_src` is finite because
the last equation of `(3)` has nonzero right side.  From `(2)--(5)`,

```text
ord(T_F)=2a-1,               ord(K_Q)=1-a.                    (7)
```

There are three exhaustive lanes.

#### Polar lane: `a=b=0`

Pass to the quadratic local deck `U^2=V` and put

```text
q=A_src/U,        omega=B_src/(2U),        rho=q/omega^3.     (8)
```

Then

```text
ord(s_F)=-1,               ord(rho)=1.                      (9)
```

The exact prefix formulas of THM-2760 are

```text
in(phi_j)=4 in(s_F)^(j-1) P_j(rho),
in(psi_j)=4 in(s_F)^j Q_j(rho),                             (10)
```

with `P_j(0)Q_j(0)!=0`.  Hence

```text
ord(phi_j)=-(j-1),             ord(psi_j)=-j.                 (11)
```

The top row `R` is strictly below every retained `j<=R-2` in both channels.
It therefore cannot be cancelled, contradicting `phi_Q=0` (and independently
the regularity of `psi_Q`).

#### Pure-`q` lane: `a=0`, `b>=1`

Here

```text
ord(q)=-1/2,          ord(d_F)>=0,          ord(s_F)>=0.   (12)
```

The unique pure-`q` top face cycles through the three flux channels.  Exact
coefficient extraction, with the missing `R-1` row in `(1)`, gives:

```text
R=3k+1:  phi_R has unique order -2k below every retained phi_j;
R=3k:    psi_R has unique order -2k below every retained psi_j;
R=3k+2:  ord(K_Q)=-(2k+1).                                  (13)
```

The first two cases contradict `(3)`.  In the third, `(7)` requires
`ord(K_Q)=1`, again impossible.  Notice that this removes rather than merely
catalogues the old `R=3k+2` local resonance: at a common simple zero there is
no equality face.

#### Regular-`H` lane: `a>=1`

If `b=0`, then `ord(d_F)=-1`, but

```text
ord(T_F d_F)=2a-2>=0,        ord(s_F)>=0.                  (14)
```

If `b>=1`, both `d_F` and `s_F` are regular.  In either case every
variable in the polynomial ring `(4)` is regular.  Thus

```text
ord(H_Q)>=0,                  ord(K_Q)>=2a-1>=1,               (15)
```

unless `H_Q=0`, which contradicts `(3)` even more directly.  But `(7)`
requires `ord(K_Q)=1-a<=0`.  This is the last contradiction and proves the
lemma.  QED.

The source is a common response divisor, the target is its Faber flux bank,
and the map is the three-lane valuation fan.  It preserves the unique
top-versus-lower face separation.  Valuation alone loses face coefficients,
the divisor identity `A_src K_Q=lambda M`, and global monodromy; the exact
prefix constants and the polynomial ring `(4)` are load-bearing sidecars.
This is why MISTAKE-317 remains the correct hostile warning for arbitrary
pole points even though simple common zeros now close.

## 2. Balanced consequence: the first resonant cell doubles

Use the balanced factor normal form of THM-2796:

```text
V=v S D_poles T_poles^2,
M=S E T_poles,
N=s_resp+2e=sum_j p_j,                s_resp=deg(S).           (16)
```

Every root of `S` is a common simple root of `V` and `M`.  The lemma forces

```text
s_resp=0.                                                     (17)
```

THM-2827 then says that a surviving bank must have

```text
R=3k+2,           D=4k+3,           p_j=D delta_j.            (18)
```

Since `N=2e` is even and `D` is odd, `sum_j delta_j=N/D` is even.  With
`s_resp=0`, genuine nonsplitting requires at least one pole part `p_j` to be odd,
equivalently at least one `delta_j` to be odd.  The positive integers
`delta_j` therefore contain at least two odd entries and have even sum at
least two.  Consequently

```text
N>=2D=2(4k+3).                                                (19)
```

Equality in `(19)` forces the two-pole passport

```text
(p_1,p_2)=(D,D).                                              (20)
```

At the first resonance `R=8`, `D=11`; the first possible balanced response
cell is therefore `N=22`, pole passport `(11,11)`, not the cyclic `N=11`
cell.  In particular, no degree-seven balanced response from the heptic
accessory atlases THM-2840/THM-3123 can occur in a bank with `R>=7`.  This is
a chart obstruction, not a proof of the planar Jacobian conjecture and not a
construction of chart entry.

## 3. Sharp abstract odd-bipole family

The numerical bound `(19)` is sharp within the balanced response equations.
Let

```text
d=2m+1,                  T_p(x)=x^2-1,                         (21)
```

and let `E_d` be the polynomial part at infinity of the Laurent branch
`T_p^(d/2)` with leading term `x^d`:

```text
E_d(x)=sum_(j=0)^m (-1)^j binom(d/2,j)x^(d-2j).               (22)
```

This notation does **not** treat the half-integral power as a polynomial.
Termwise binomial cancellation gives the first-order identity

```text
2T_p E_d'-2d x E_d=C_d,
C_d=-2(-1)^m binom(d/2,m) != 0.                               (23)
```

Define

```text
F=E_d^2/T_p^d,
V=T_p^(d+2),
G=E_d/T_p^(d+1).                                              (24)
```

Then `(23)` gives exactly

```text
F=VG^2,
2VG'+V'G=C_d,
V(F')^2=C_d^2 F.                                              (25)
```

The ODE also proves the required clean factor typing.  A common root of
`E_d` and `T_p`, or of `E_d` and `E_d'`, would make the left side of `(23)`
zero, contradicting `C_d!=0`.  Put

```text
W_d=E_d^2-T_p^d.                                               (26)
```

The first omitted Laurent term in `(22)` shows

```text
deg W_d=d-1,                 LC(W_d)=-C_d/(d+1).              (27)
```

Differentiating `(26)` and using `(23)` also gives

```text
T_p W_d'-2d x W_d=C_d E_d.                                  (28)
```

At a root of `W_d`, both `E_d` and `T_p` are nonzero and

```text
W_d'=C_d E_d/T_p!=0.                                         (29)
```

Thus all three fibres are disjoint and have the clean passport

```text
zero:       2^d;
pole:       (d,d);
third:      (d+1,1^(d-1)).                                   (30)
```

Because `d` is odd, `(24)` is genuinely nonsplit.  It realizes `N=2d=2D`
abstractly for every odd divisor scale and proves that `(19)` cannot be
improved from the balanced response identity alone.

The constants in `(23)` themselves form an `O(1)`-state exact sequence.  If
`C_m` denotes the constant for `d=2m+1`, then

```text
C_0=-2,
C_(m+1)=-(2m+3)/(2m+2) C_m.                                  (31)
```

This product-Gamma recurrence is the nonvanishing sidecar missed by a bare
Newton-face observer.

## 4. First Faber-resonant control at `d=11`

For `d=11`, the family is

```text
E_11=x^11-(11/2)x^9+(99/8)x^7-(231/16)x^5
     +(1155/128)x^3-(693/256)x,
C_11=693/128,                 deg W_11=10.                    (32)
```

The exact clean Nielsen census has `121=11^2` marked matchings.  Independent
rotations of the two pole cycles act transitively on them, so there is one
unmarked orbit.  A representative is transitive and its generated
permutation group has exact order

```text
40874803200=2^10*11!.                                        (33)
```

Only this finite-exact order is asserted; no group-isomorphism claim is
needed.

## 5. Reproduction and scope

The exact companion independently extracts Faber rows through `R=30`, keeps
every allowed lower row, checks all three simple-zero lanes, and verifies the
nonzero top-face recurrence.  A separate block verifies `(22)--(31)` for odd
`d<=21` and performs the `d=11` Nielsen census.  Run

```text
python3 04-computation/jc_common_simple_zero_odd_bipole_thm3133.py
python3 -O 04-computation/jc_common_simple_zero_odd_bipole_thm3133.py
```

and compare byte-for-byte with the declared output.

No claim is made that the sharp family enters a Keller chart, that every
balanced response is realized by a Faber bank, or that the planar Jacobian
conjecture follows.  The common-simple-zero lemma is an all-`R` symbolic
result; the displayed Faber and Nielsen ranges are verification sidecars.

## 6. Independent hostile audit

The independent audit rebuilt the Faber coefficients without importing the
canonical companion and checked every retained lower row for `R=7..18`.
It rederived all three valuation lanes, the balanced parity argument, the
all-degree identities `(23)` and `(28)`, squarefreeness, disjoint fibres, and
the exact leading coefficient in `(27)`.  Its separate `d=11` Nielsen replay
agreed on `121` marked matchings, one rotation orbit, and the group order.
Normal and optimized runs byte-match, and the canonical orbit assertion is
now an executable gate.  No assertion statements or floating-point literals
occur in either audit engine.

**End of proof.**
