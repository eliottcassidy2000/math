---
id: THM-3377
title: "Path-colour deletion compiler, positive reciprocity, and skew current"
status: >
  PROVED + VERIFIED-EXACT.  The complete symmetric subset-deletion sums of the
  THM-3166 path-colour polynomial form a bivariate multiplicative ordered-join compiler.
  Its deletion derivative is the summed marked first response, its colour-one
  face is exactly the shifted THM-3372 Hamiltonian deletion transform, and its
  negative-colour face has coefficientwise positive reciprocity.  Contracting
  marked compiler responses against A-A^T gives a closed skew ordered-join
  current.  The scalar Q and current do not close without the first-deletion
  jet: two strong order-six hostiles share both but acquire different currents
  after adjoining one sink.  No injectivity, cyclic-substitution, chronology,
  arbitrary-input complexity, LRC, GMC, or JC consequence is asserted.
source: root/repository-archaeology-second-pass-2026-08-14
depends_on:
  - THM-3166-falling-factorial-order-join-path-colour-transform
  - THM-3372-multiaffine-deletion-transform-variance-and-skew-join-current
related:
  - THM-3369-skew-deletion-response-and-ordered-join-orientation-current
script: 04-computation/tournament_path_colour_deletion_compiler_skew_current_thm3377.py
output: 05-knowledge/results/tournament_path_colour_deletion_compiler_skew_current_thm3377.out
script_sha256: e524ced0394c8796407914d283b90c87b6ec5549b4a3b1134710660c83d4a202
output_sha256: 00c2749ac7cbaad38e14f72a22642955f75e28c7e612c8f5834f5577dd6bcbbd
semantic_sha256: 4403a78a56ee1c57dcb0e64865318e4f09187b26fa6a83697f852308d8e612ba
hash_basis: LF-normalized bytes
---

# THM-3377 -- the path-colour deletion deck is a multiplicative exterior response

**PROVED + VERIFIED-EXACT.**

THM-3166 makes the complete spanning path-cover profile multiplicative after
the falling-factorial path-colour transform.  Aggregating every induced
deletion by its size produces a stronger two-variable operation interface. Its
colour-one face is THM-3372, its negative-colour face is coefficientwise
positive, and its marked derivative supplies the sidecar needed by a skew
ordered-join current.

## 1. The bivariate deletion compiler

For a tournament `T`, let

```text
Q_T(t)=sum_(d=1)^|T| p_T(d)(t)_d                       (1)
```

be the THM-3166 path-colour polynomial, where `p_T(d)` counts spanning covers
by `d` unordered directed paths and `(t)_d` is the falling factorial.  Set
`Q_empty(t)=1`.  Define

```text
F_T(a,t)=sum_(X subset V(T)) Q_(T-X)(t) a^|X|.         (2)
```

Then for all disjoint tournaments `X,Y`,

```text
F_(X triangleright Y)(a,t)=F_X(a,t)F_Y(a,t).           (3)
```

*Proof.*  Every deleted set in `X union Y` splits uniquely as `A union B`.
THM-3166 gives

```text
Q_((X-A) triangleright (Y-B))=Q_(X-A)Q_(Y-B),          (4)
```

including an empty factor under the convention above.  Substitution of `(4)`
in `(2)` factors the two finite sums and proves `(3)`. ∎

Thus the coefficient of `a^k` is the complete size-`k` deletion sum of the
path-colour polynomial.  In particular,

```text
F_T(0,t)=Q_T(t).                                       (5)
```

## 2. The deletion derivative and the Hamiltonian face

For every tournament,

```text
partial_a F_T(a,t)=sum_(v in V(T)) F_(T-v)(a,t).       (6)
```

Indeed, a fixed nonempty deleted set `Z` contributes to the right side once
for each choice `v in Z`, with remaining exponent `|Z|-1`.  This is exactly
the derivative of its term in `(2)`.

Put

```text
q_T(t)=partial_a F_T(0,t)=sum_v Q_(T-v)(t).             (7)
```

Differentiating `(3)` at `a=0` gives the dual-number or first-jet law

```text
q_(X triangleright Y)=q_X Q_Y+Q_X q_Y.                 (8)
```

Equivalently, `Q_T+epsilon q_T` is multiplicative modulo `epsilon^2`.

At the colour `t=1`, every falling factorial `(1)_d` with `d>=2` vanishes,
while `Q_T(1)=p_T(1)=H(T)`.  The empty convention agrees with `H(empty)=1`.
Consequently THM-3372's diagonal deletion transform satisfies

```text
F_T(a,1)=sum_X H(T-X)a^|X|=D_T(1+a).                   (9)
```

Thus `(2)` is not a competing Hamiltonian invariant: it is a path-cover lift
whose `t=1` face is exactly the shifted full Hamiltonian deletion deck.

## 3. Coefficientwise positive negative-colour reciprocity

For a tournament `S` of order `s` and integer `m>=1`, THM-3166 proves

```text
R_S(m):=(-1)^s Q_S(-m)
      =sum_(pi in Sym(V(S))) binom(m+b_S(pi),s)>=0,     (10)
```

where `b_S(pi)` counts backward consecutive adjacencies.  Set `R_empty(m)=1`.
Applying `(10)` to every induced deletion in `(2)` gives

```text
(-1)^n F_T(-a,-m)=sum_(X subset V(T)) R_(T-X)(m)a^|X|. (11)
```

Hence the left side has nonnegative integer coefficients in `a` for every
integer `m>=1`.  This is stronger than positivity of the undeleted value: the
whole Boolean deletion grading survives the reciprocity involution
`(a,t)->(-a,-t)` after the global parity `(-1)^n`.

## 4. Exterior current of the full compiler

Let `K_T=A_T-A_T^T`.  For a variable pair `alpha=(a,z)`, define the marked row

```text
r_T(alpha)=(F_(T-v)(a,z))_(v in V(T))                  (12)
```

and the four-variable skew current

```text
Psi_T(alpha,beta)=r_T(alpha)K_T r_T(beta)^T,            (13)
```

where `beta=(b,w)`.  It is invariant under relabelling and obeys

```text
Psi_T(beta,alpha)=-Psi_T(alpha,beta),
Psi_(T^op)=-Psi_T.                                     (14)
```

Write

```text
F_Z^alpha=F_Z(a,z),       dot F_Z^alpha=partial_a F_Z(a,z),
F_Z^beta =F_Z(b,w),       dot F_Z^beta =partial_b F_Z(b,w).
```

Then

```text
Psi_(X triangleright Y)(alpha,beta)
 =F_Y^alpha F_Y^beta Psi_X(alpha,beta)
  +F_X^alpha F_X^beta Psi_Y(alpha,beta)
  +(F_Y^alpha dot F_X^alpha)(F_X^beta dot F_Y^beta)
  -(F_X^alpha dot F_Y^alpha)(F_Y^beta dot F_X^beta).   (15)
```

*Proof.*  If `x in X`, then `(X triangleright Y)-x=(X-x) triangleright Y`,
so its marked response is `F_(X-x)^alpha F_Y^alpha`; the analogous response
for `y in Y` is `F_X^alpha F_(Y-y)^alpha`.  The two diagonal blocks of

```text
K_(X triangleright Y)=[[K_X,+J],[-J,K_Y]]              (16)
```

give the first two terms of `(15)`.  Summing each all-ones cross block and
using `(6)` gives the exterior pair in the last line. ∎

For the repeated ordered join `T^r=T triangleright ... triangleright T`, the
two cross terms cancel because their factor derivatives are proportional.
Induction in `(3)` and `(15)` therefore gives the compact all-`r` law

```text
F_(T^r)(alpha)=F_T(alpha)^r,
Psi_(T^r)(alpha,beta)
 =r F_T(alpha)^(r-1)F_T(beta)^(r-1)Psi_T(alpha,beta).   (15a)
```

Thus the deletion compiler and its order current retain the same exact
polynomial--exponential repeated-join structure as the undeleted THM-3166
response.

The scalar slice

```text
psi_T(z,w)=Psi_T((0,z),(0,w))
            =sum_(u,v)K_(u,v)Q_(T-u)(z)Q_(T-v)(w)      (17)
```

therefore satisfies

```text
psi_(X triangleright Y)(z,w)
 =Q_Y(z)Q_Y(w)psi_X(z,w)+Q_X(z)Q_X(w)psi_Y(z,w)
  +[Q_Yq_X](z)[Q_Xq_Y](w)-[Q_Xq_Y](z)[Q_Yq_X](w).     (18)
```

Unlike the scalar pair `(Q,psi)`, the full pair `(F,Psi)` needs no external
first-deletion sidecar: equation `(6)` recovers it by differentiation.

## 5. Sharp order witness and the load-bearing first jet

For `K1 triangleright C3` and its reverse, the commutative compiler is the
same:

```text
F(a,t)=t^4+2t^2+a(4t^3+2t)+6a^2t^2+4a^3t+a^4.         (19)
```

But their scalar currents are opposite:

```text
psi_(K1 triangleright C3)(z,w)=6zw(w^2-z^2),
psi_(C3 triangleright K1)(z,w)=6zw(z^2-w^2).           (20)
```

Thus `psi(1,2)=+36` versus `-36` restores the source/sink bit lost by `(19)`.
At `t=1`, equations `(19)--(20)` specialize to THM-3372 after the shift
`a=y-1`.

The first deletion jet is genuinely necessary for scalar composition.  In
the pair-bit convention of the companion, the strong nonisomorphic order-six
masks `16` and `83` have the common path-colour polynomial and current

```text
Q(t)=t^6+8t^4+8t^2,                  psi(z,w)=0,        (21)
```

but

```text
q_16(t)=6t^5+24t^3+8t,
q_83(t)=6t^5+24t^3+4t.                                (22)
```

After adjoining one sink on the right, their currents differ already in the
coefficient of `z^2w^4`, which is respectively `-128` and `-160`.  Therefore
neither `q` nor the joined scalar current is determined by `(Q,psi)` alone.
The `a`-coefficient of `F` supplies exactly the missing coordinate.

The typed contract is

```text
source:  complete path-cover profile on every induced deletion + arc signs
map:     falling-factorial colour transform -> deletion compiler -> r K r^T
target:  exact ordered-join composition with one restored order channel
kept:    deletion size, every path-cover depth, negative-colour positivity
lost:    deleted-set labels after diagonalization, full SCC order, chronology
sidecar: none for (F,Psi); q=[a]F is mandatory for the scalar (Q,psi) slice
tests:   K1 triangleright C3 reversal; strong masks 16/83 joined with K1. (23)
```

No injectivity, reconstruction, cyclic-substitution, arbitrary-input path
counting, LRC, Gaussian-moment, Keller, or chronological claim follows.

## 6. Exact audit

The standard-library-only companion imports the independently audited
THM-3166 path-cover engine and the audited THM-3372 deletion engine.  It checks:

- all `1,099` labelled tournaments of orders one through five;
- `5,405` marked-deletion cells in `(6)` and the specialization `(9)`;
- `169,330` induced-deletion negative-colour cells against the permutation
  backward-adjacency formula `(10)`;
- `2,197` relabellings, converse sign and variable skewness;
- nonzero-current labelled counts `0,0,0,16,320` by order;
- all `2,553` ordered labelled factor pairs of total order at most six for
  `(3)` and the full four-variable law `(15)`; and
- the literal order-four reversal and order-six/order-seven sidecar hostiles.

Normal and optimized modes byte-match the stored eleven-line transcript.
Reproduce with

```text
python3 04-computation/tournament_path_colour_deletion_compiler_skew_current_thm3377.py
python3 -O 04-computation/tournament_path_colour_deletion_compiler_skew_current_thm3377.py
```

All hashes are pinned in the frontmatter.

**QED.**
