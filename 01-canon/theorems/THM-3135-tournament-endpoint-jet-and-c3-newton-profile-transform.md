---
id: THM-3135
title: "Tournament endpoint jets and the cubic C3 Gregory-Newton profile transform"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The complete
  factorial path-cover profile of a tournament is exactly the Taylor jet at
  u=1/2 of one centered permutation polynomial; its Hamiltonian-path count is
  only the endpoint value.  For C3 substitution, the THM-3121 multivariate
  kernel is conjugate to a one-dimensional Gregory-Newton transform, giving
  the complete equal-factor profile in O(N^3) arithmetic and O(N^2) memory.
  Iteration yields exact Hamiltonian counts through 729 vertices and a free
  cyclic-wreath divisibility law for every ordered path-cover count.  No
  fixed-order scalar recurrence or arbitrary-tournament complexity collapse
  is claimed.
audit: >
  An independent engine enumerated all 1,099 labelled tournaments through
  order five, checked direct path covers, permutation cuts, the inverse jet,
  and parity, then checked 27 unequal-factor C3 substitutions, two independent
  C3 engines through 27 vertices, direct profiles through nine vertices, the
  same-H hostile pair, and cyclic-wreath divisibility.  It also audited the
  finite ranges and the general O(qM^2) complexity.  Canonical normal,
  optimized, and stored transcripts and hashes agree.
source: root/frontier-synthesis-2026-08-02
depends_on:
  - THM-3121-path-cover-walk-content-substitution-kernel
related:
  - THM-1975-the-path-cover-polynomial-is-the-refined-compositional-invariant
  - THM-2412-delta-exponential-and-central-newton-layer-split
  - THM-2438-poisson-newton-ternary-half-and-harmonic-divisor-incidence
script: 04-computation/tournament_c3_endpoint_jet_newton_thm3135.py
output: 05-knowledge/results/tournament_c3_endpoint_jet_newton_thm3135.out
script_sha256: 7aa2974d36d6431315f3afa071a7a3cb63407fee1e8c7621c53a5f8917fccb8f
output_sha256: 8ba3debc0708fa5091392c6b142662c8fc2c4c4303443cf8539fcc1cdff431c8
hash_basis: LF-normalized bytes
---

# THM-3135 -- tournament endpoint jets and the cubic `C3` Newton transform

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3121 shows that cyclic tournament substitution needs the complete
path-cover profile rather than the scalar Hamiltonian-path count.  This
theorem identifies that sidecar with an endpoint Taylor jet and conjugates
the full `C3` profile update to a Gregory--Newton transform.  The resulting
recurrence is structural and all-level; it is not a fit to the first few
integer values.

## 1. Path covers are an endpoint jet

Let `T` be a tournament on `n>=1` vertices.  Write

```text
pc_T(d)=number of spanning covers by d unordered directed paths,
F_T(d)=d! pc_T(d).                                            (1)
```

Thus `F_T(d)` counts ordered path covers.  For a permutation `pi` of the
vertices, let `b_T(pi)` be the number of consecutive pairs oriented backward
along `pi`, and put `g_T(pi)=n-1-b_T(pi)`.  Define the centered permutation
polynomial

```text
E_T(u)=sum_(pi in S_n)
       (u-1/2)^b_T(pi) (u+1/2)^g_T(pi).                       (2)
```

Then

```text
sum_(d=1)^n F_T(d) z^d = z E_T(z+1/2),                       (3)

pc_T(d)=E_T^(d-1)(1/2)/(d!(d-1)!).                           (4)
```

### Proof

Order the components of a path cover and concatenate their vertex words.
This gives a permutation together with cuts between consecutive entries.
Every backward adjacency must be cut; each forward adjacency may either be
cut or retained.  For a fixed permutation, the cut generating polynomial is

```text
z^(1+b_T(pi))(1+z)^g_T(pi).                                  (5)
```

Summing `(5)` over permutations is exactly the right side of `(3)`.  The
construction is bijective: cutting at the chosen positions recovers the
ordered directed paths.  Extracting the coefficient of `z^d` proves `(4)`.

If `A_T(b)` counts permutations with `b` backward adjacencies, `(3)` gives
the exact triangular dictionary

```text
F_T(d)=sum_(b=0)^(d-1)
       binom(n-1-b,d-1-b) A_T(b),                             (6)

A_T(b)=sum_(j=0)^b (-1)^(b-j)
       binom(n-1-j,b-j) F_T(j+1).                             (7)
```

Reversing a permutation exchanges forward and backward adjacencies.  Hence

```text
A_T(b)=A_T(n-1-b),
E_T(-u)=(-1)^(n-1)E_T(u).                                    (8)
```

The full path-cover sidecar is therefore exactly the endpoint jet of `E_T`.
The scalar `H(T)=pc_T(1)=E_T(1/2)` remembers only its value.

## 2. Gregory--Newton conjugacy for `C3`

Let `S_1,S_2,S_3` have orders `n_1,n_2,n_3`.  Put

```text
q=min_i n_i,                    M=n_1+n_2+n_3,                (9)
```

and extend each `F_i=F_(S_i)` by zero outside `1<=a<=n_i`.
For `0<=s<=q`, define

```text
B_(i,s)(k)=sum_(r=0)^(n_i-s) binom(k,r) F_i(s+r),             (10)

D_s^(m)=Delta_k^m (prod_(i=1)^3 B_(i,s)(k)) at k=0.          (11)
```

Then for every `1<=d<=M`, the complete substituted profile is

```text
F_(C3[S_1,S_2,S_3])(d)
 =sum_(j=0)^min(d,q) 2^j binom(d,j)
   sum_(t=0)^(q-j) binom(d+t-1,t) D_(t+j)^(d-j).              (12)
```

### Proof

Orient `C3` cyclically and put

```text
P=(1+x)(1+y)(1+z),       U=P-1,       p=xyz.                  (13)
```

The quotient-walk series from THM-3121 is

```text
W_(C3)=(U+2p)/(1-p).                                          (14)
```

Let `L_i(x_i^a)=F_i(a)`.  Multiplying the THM-3121 path-cover law by `d!`
gives the clean functional form

```text
F_(C3[S_1,S_2,S_3])(d)=L_1 L_2 L_3(W_(C3)^d).                (15)
```

Expand

```text
(1-p)^(-d)=sum_(t>=0) binom(d+t-1,t)p^t,
(U+2p)^d=sum_(j=0)^d 2^j binom(d,j)p^j U^(d-j),               (16)

U^m=(P-1)^m=sum_(k=0)^m(-1)^(m-k)binom(m,k)P^k.              (17)
```

For `s=t+j`, the functional applied to `p^sP^k` factors as

```text
L_1L_2L_3(p^sP^k)=prod_i B_(i,s)(k).                         (18)
```

The alternating sum in `(17)` is precisely the forward difference `(11)`.
Only `s<=q` survives, yielding the exact finite ranges in `(12)`.  QED.

This is literally the Gregory--Newton transform used in THM-2412/2438, now
acting on the endpoint-jet sidecar of a tournament substitution kernel.

## 3. Complexity and the equal-factor recurrence

Precomputing `(10)--(11)` and evaluating `(12)` takes

```text
O(q M^2) arithmetic,                 O(q M) memory.            (19)
```

For three equal `N`-vertex factors this is `O(N^3)` arithmetic and `O(N^2)`
memory.  The displayed sparse three-dimensional coefficient-powering engine
used as a cross-check has `O(N^5)` worst-case cost; this is an implementation
comparison, not a lower bound for all possible kernel algorithms.

For equal factors, write

```text
B_s(k)=sum_r binom(k,r)F(s+r),
D_s^(m)=Delta^m(B_s(k)^3) at k=0.                             (20)
```

Equations `(12)` and `(20)` are an exact all-profile recurrence.  Its `d=1`
slice recovers THM-3121's diagonal law

```text
H(C3[S,S,S])
 =3 sum_(r>=1)(F(r)^3+F(r+1)F(r)^2+F(r+1)^2F(r)).            (21)
```

## 4. Iterated `C3` values and wreath divisibility

Put

```text
T_0=K_1,                 T_(ell+1)=C3[T_ell,T_ell,T_ell].     (22)
```

The recurrence gives

```text
H(T_0)=1,
H(T_1)=3,
H(T_2)=3159,
H(T_3)=417382500592116023859.                                (23)
```

The next three exact values have respectively `98`, `404`, and `1554`
decimal digits.  Their decimal-string SHA-256 values are

```text
ell=4: 97119e814f11b6cffbbfa081da65fd7a27f77928aa9231c1100b8776af46ee19
ell=5: 006ed2c3a597eefbc50dc24aa8ea5241ca943de16f0b281c1a469f2a4c4d960c
ell=6: 6cfeb448f0af07021a3c38a7a967595bc7a085502d0271cf5bf2306c7d2a5f47. (24)
```

There is also an all-profile divisibility law.  The automorphism group of
`T_ell` contains the canonical iterated cyclic wreath subgroup

```text
G_0=1,        G_ell=G_(ell-1)^3 semidirect C3,
|G_ell|=3^((3^ell-1)/2).                                    (25)
```

This subgroup acts freely on ordered spanning path covers: an element fixing
an ordered list of vertex paths fixes every vertex.  Therefore

```text
3^((3^ell-1)/2) divides F_(T_ell)(d)       for every d>=1.   (26)
```

The factorial in `F=d!pc` is essential; `(26)` is false in general for the
unordered count, as `pc_(C3)(3)=1` already shows.  The observed valuations

```text
v_3(H(T_ell))=1,5,13,40,121,364,       ell=1,...,6           (27)
```

meet the wreath floor from `ell=3` onward, but equality beyond the verified
range is not asserted.

## 5. Hostile boundaries

Two labelled order-five tournaments, edge-bit masks `40` and `76`, both have
`H=15`, but

```text
E_40(u)=120u^4+18u^2+3,       F_40=(15,78,198,240,120),
E_76(u)=120u^4+30u^2,         F_76=(15,90,210,240,120).       (28)
```

Consequently

```text
H(C3[T_40,K_1,K_1])=3F_40(1)+F_40(2)=123,
H(C3[T_76,K_1,K_1])=3F_76(1)+F_76(2)=135.                   (29)
```

Thus the endpoint value alone provably loses information required by cyclic
substitution; the first derivative already separates the minimal hostile.

For a transitive tournament, each vertex subset has one directed increasing
path, so

```text
pc_T(d)=S(n,d),              F_T(d)=d!S(n,d),                 (30)
```

the Stirling/Fubini boundary.  In `(14)`, `U` is this separable ordered-block
part, `2p` is the extra cyclic-start multiplicity, and `(1-p)^(-d)` records
complete quotient revolutions.

The transform does not extend automatically to arbitrary quotient
tournaments, compute arbitrary base profiles, or collapse `#P`-hard
Hamiltonian-path counting.  Nor is `(12)` a fixed-order scalar recurrence in
the level `ell`: the profile length triples at each iteration.

## 6. Connection contract and verification

The source is a tournament path-cover profile; the target is the profile of
its `C3` substitution.  Maximal block-run splitting is the map.  Endpoint-jet
conjugacy preserves the entire profile, while the endpoint value destroys all
higher run counts.  The full Taylor jet and the quotient's rational walk
series restore the lost data.

The link to factorial THM-3131 is typed but not a reduction.  Both use a
prime-place observer after a Newton reframe, but here `3`-adic divisibility
comes from free symmetry orbits; the Newton summands can overlap and cancel.
The link to the Jacobian lane is likewise methodological: retaining a jet can
restore data lost by a scalar accessory, but no Keller-chart map is supplied.

The canonical companion checks the kernel and Newton engines through the
81-vertex profile, exact values through 729 vertices, wreath divisibility
through 243 vertices, all `1,099` labelled tournaments through order five,
the hostile pair, and transitive Stirling controls.  Run

```text
python3 04-computation/tournament_c3_endpoint_jet_newton_thm3135.py
python3 -O 04-computation/tournament_c3_endpoint_jet_newton_thm3135.py
```

and compare with the declared transcript.

The independent audit used a separate subset/set-partition engine and an
independent unequal-factor implementation.  All `27` unequal `C3` controls,
direct profiles for `T_1,T_2`, two profile engines through `T_3`, the inverse
jet, parity, and the symmetry scope passed.  Normal, optimized, and stored
transcripts agree byte for byte.

**End of proof.**
