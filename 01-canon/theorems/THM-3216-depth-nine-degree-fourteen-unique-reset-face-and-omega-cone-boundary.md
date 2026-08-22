---
id: THM-3216
title: "Depth-nine degree-fourteen unique-reset face and omega cone boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For the support-(1,3), bank-I2 product-Gamma selector, the complete
  depth-at-most-nine degree-fourteen feasibility cone is the singleton
  delta_(1,3,3,4,5,6,7,8).  An explicit positive 18-row combination of
  lawful coarsening-upset responses vanishes at that reset and is strictly
  negative on the other 3,128 states.  Hence the same singleton statement
  holds at every horizon D>=14, including the all-degree intersection.
  Multiset complement preserves the full exterior response through the
  symmetric-function involution omega, but it does not preserve the selector
  order cone: failure occurs at degree four and mixed signs occur at degree
  six.
source: root/multiscale-newton-flag/product-gamma-width3/2026-08-02
audit: >
  The exact companion pins promoted THM-3209, reconstructs all 3,129 legal
  states without a discovery cache, verifies every upset and degree
  normalizer, clears the 18 positive rational multipliers primitively, and
  checks the full depthwise negative/zero census and exact nearest/farthest
  gaps.  It separately verifies all 631 depth-nine/depth-seven complement
  identities and derives the degree-four and degree-six omega images by an
  exact power-sum/monomial change of basis.  Normal, optimized, and stored
  replay agree byte-for-byte.  Two independent hostile audits checked the
  cone implication, exact finite evidence, complement typing, and omega
  boundary.
depends_on:
  - THM-3127-partition-refinement-strassen-upset-dual-and-filter-response
  - THM-3160-complete-pluecker-pole-holotopy-and-selector-projection-no-go
  - THM-3184-depth-seven-degree-fourteen-farkas-death
  - THM-3209-depth-eight-complete-quotient-reset-and-negative-singleton-tangent
script: 04-computation/gmc_depth_nine_unique_reset_face_omega_boundary_thm3216.py
output: 05-knowledge/results/gmc_depth_nine_unique_reset_face_omega_boundary_thm3216.out
script_sha256: 3951b7e9a5e08199fd03ab2b3827dcbeb8c1039790acd4d218713c78d44ab9cf
output_sha256: d303a7fb18b07f2091e1fe5c705b15e62ecd2f9049fe6f29a003834288a24636
hash_basis: LF-normalized bytes
---

# THM-3216 -- depth-nine degree-fourteen unique-reset face and omega cone boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3209 found an all-degree reset at depth eight and proved only that its
four one-extra-pole tangent rays point out of the selector cone.  The other
627 legal depth-nine states were not controlled.  Here an exact global
supporting cocircuit controls all of them at once.  The result is stronger
than tangent isolation: adding the complete ninth layer does not enlarge the
feasible set at degree fourteen at all.

The second half identifies the tempting but false shortcut suggested by pole
complement.  Complement is an exact symmetry of the full exterior endpoint
plane from THM-3160.  It is not a symmetry of the fine-to-coarse order cone
used by the selector.  This is the precise holotopy boundary: the flat
Pluecker object transports, while its positive same-degree projection does
not.

## 1. The finite selector cone

Use the support `(1,3)` and signed bank `I2` of THM-3209.  Its reduced-pole
multiset and distinguished quotient alphabet are

```text
P=(8,7,6,5,5,4,4,3,3,2,2,2,1,1,1,1),
Q=(1,3,3,4,5,6,7,8).                                    (1)
```

Let `S_<=9` be the nonempty submultisets of `P` of cardinality at most nine.
The exact layer sizes, from depth one through nine, are

```text
(8,33,93,200,348,507,631,678,631),       |S_<=9|=3129.   (2)
```

For a partition `mu` of `N`, retain THM-3209's response

```text
G_N^sigma(mu)
 =Phi^sigma(h_N)m_mu[Q-sigma]
  -Phi^sigma(m_mu)h_N[Q-sigma].                           (3)
```

For a nonempty proper coarsening upset `U` in the partition lattice of `N`,
put

```text
r_(N,U)(sigma)=sum_(mu in U) G_N^sigma(mu).               (4)
```

By THM-3127, a probability law `lambda` belongs to `C_D^(<=9)` exactly when

```text
sum_sigma lambda_sigma r_(N,U)(sigma)>=0                  (5)
```

for every lawful `U` and every `5<=N<=D`.

## 2. The exact eighteen-row dual

Write `1^k` for `k` trailing ones and `<g_1,...,g_s>` for the coarsening
upset with the displayed minimal generators.  The certificate uses the
following rows.  The final column is the exact positive factor `a_i` in
the factorized multiplier `theta_i=a_i/M_N`.

| `i` | `N` | `|U_i|` | minimal generators | `a_i` |
|---:|---:|---:|---|---:|
| 1 | 5 | 1 | `(5)` | `5/72755559` |
| 2 | 14 | 130 | `(3,2,1^9);(2,2,2,1^8)` | `1577093/47102691` |
| 3 | 8 | 21 | `(2,1^6)` | `53234/56508147` |
| 4 | 10 | 40 | `(3,1^7);(2,2,1^6)` | `151685/66318847` |
| 5 | 12 | 76 | `(2,1^10)` | `30203115/76749488` |
| 6 | 12 | 74 | `(2,2,1^8)` | `57683/3902927` |
| 7 | 14 | 128 | `(4,1^10);(3,2,1^9);(2^6,1,1)` | `283819/63667394` |
| 8 | 14 | 121 | `(7,1^7);(3,3,2,1^6);(2^4,1^6)` | `2861/13045678` |
| 9 | 14 | 132 | `(2,2,1^10)` | `10395173/58298789` |
| 10 | 14 | 129 | `(5,1^9);(3,3,1^8);(2,2,2,1^8)` | `1951835/98897507` |
| 11 | 14 | 127 | `(5,1^9);(4,2,1^8);(3,3,1^8);(2^4,1^6)` | `3979/1773976` |
| 12 | 10 | 41 | `(2,1^8)` | `4381327/51259840` |
| 13 | 11 | 55 | `(2,1^9)` | `6670053/32445565` |
| 14 | 9 | 29 | `(2,1^7)` | `1319165/92555776` |
| 15 | 13 | 97 | `(3,1^10);(2^4,1^5)` | `54743/8093902` |
| 16 | 11 | 54 | `(3,1^8);(2,2,1^7)` | `198899/26117585` |
| 17 | 14 | 113 | `(7,1^7);(6,2,1^6);(5,3,1^6);(5,2,2,1^5);(4,4,1^6);(4,3,2,1^5);(3,3,3,1^5);(3,3,2,2,2,1,1);(2^6,1,1)` | `2057/73733672` |
| 18 | 13 | 98 | `(3,1^10);(2,2,2,1^7)` | `1313225/43330928` |

The degree normalizers are

```text
M_5 = 1300713120,
M_8 = 9558239814808320,
M_9 = 1141773245941342464,
M_10= 98177407096144199040,
M_11= 6402822693369065077104,
M_12= 623312546828577253020720,
M_13= 39394150378793250375693600,
M_14= 2528571670236939601303479360.                         (6)
```

Thus every `theta_i` is a displayed positive rational number.  Let `L` be
the least common denominator of the eighteen `theta_i`, divide the integer
tuple `(L theta_i)` by its gcd, and call the resulting primitive positive
tuple `(c_i)`.  The exact companion prints all eighteen integers and their
SHA-256 digest.  This factorized presentation is the same primitive dual,
without hiding the natural degree scaling inside two-hundred-digit entries.

Define

```text
H(sigma)=sum_(i=1)^18 c_i r_(N_i,U_i)(sigma).              (7)
```

Direct integer evaluation on the complete universe `(2)` gives

```text
H(Q)=0,
H(sigma)<0                    for every sigma!=Q.           (8)
```

The negative/zero census by depth is

```text
(8/0,33/0,93/0,200/0,348/0,507/0,631/0,677/1,631/0).      (9)
```

There is no floating-point inference in `(8)`.  Every response coordinate,
every multiplier after clearing, and every comparison is integral.  The
primitive multiplier tuple has digest

```text
22b881266348913671f4993625bfa9a6ea7fcf7b49760ac71e4028292f9ab1de. (9a)
```

The closest negative state is

```text
(3,3,4,4,5,5,6,7,8),

H=-136976983387068608966758345593427925570435791167435166854712683791392270004153630838413137075479245536925214723447440160337957901738538321634033707008893414500105328960783600207156143126472074967127158706197965444153154158192334977812841408. (9b)
```

The farthest is the singleton `(8)`; its exact coordinate and the full
eighteen-entry primitive tuple are printed by the companion.  The giant
integers are a denominator-clearing artifact; the factorized table and `(6)`
are the readable exact certificate.

## 3. The unique reset theorem

**Theorem.** For every finite `D>=14`, and also for the all-degree
intersection,

```text
C_D^(<=9)={delta_Q}.                                      (10)
```

Indeed, let `lambda` be feasible.  Each row in `(7)` is lawful and each
`c_i` is positive, so `(5)` gives

```text
E_lambda H>=0.                                             (11)
```

But `(8)` gives `H<=0` pointwise.  Therefore equality holds in `(11)`, and
the strict part of `(8)` forces the support of `lambda` to lie in `{Q}`.
Since `lambda` is a probability law, `lambda=delta_Q`.  Conversely,
THM-3209 proves that the quotient alphabet at `Q` is empty, so every
positive-degree response at `Q` vanishes.  Hence `delta_Q` is feasible for
every horizon.  This proves `(10)`.

Two boundaries matter.

1. The unique law is a **weak boundary law**, not a strict selector: all its
   response rows vanish.
2. The theorem is global only on the finite bank `S_<=9`.  It does not
   exclude resurrection at depth ten or construct a physical stochastic
   pole carrier.

Promoted THM-3219 supplies a clean adjacent collar.  For all 63 nonempty
physical completions `Q+tau` through full depth sixteen, its degree-five
principal-upset coordinate is exactly

```text
-1440 sum_(r in tau) r^5<0.                                (11a)
```

Inside `S_<=9`, this accounts only for the four one-pole completions of `Q`;
it is the local mechanism behind row 1 of the dual.  The other seventeen
rows in `(7)` are what stitch the thousands of states outside that principal
filter into the global singleton face.  Neither result substitutes for the
other.

The normalized discovery LP had margin about `3.321763938907998e-6`.
Independent denominator rounding first became exact at cap `10^8`, while a
naive common-grid rounding failed through `5*10^9`.  This records numerical
conditioning, not a minimal-denominator theorem; `(6)--(8)` are the exact
certificate.

## 4. Complement is exact on the exterior plane

The full pole complement of the reset is

```text
R=P-Q=(1,1,1,2,2,2,4,5).                                 (12)
```

If `sigma` has depth nine, write `sigma=P-tau`; then `tau` has depth seven
and

```text
Q-sigma=-(R-tau).                                          (13)
```

For every signed-bank atom `S_i`, put `check S_i=P-S_i` in the **virtual
alphabet group** and define the complementary endpoint

```text
check Phi^tau(f)=sum_i epsilon_i f[check S_i-tau].          (14)
```

For a homogeneous symmetric function `f` of degree `N`, the antipode
identity is

```text
f[-X]=(-1)^N omega(f)[X].                                  (15)
```

Equations `(13)--(15)` give

```text
Phi^(P-tau)(f)=(-1)^N check Phi^tau(omega(f)),
f[Q-(P-tau)] =(-1)^N omega(f)[R-tau].                      (16)
```

Let

```text
E_sigma(f)=(Phi^sigma(f),f[Q-sigma]).                      (17)
```

Then `G_N^sigma(mu)=det(E_sigma(h_N),E_sigma(m_mu))`.
Because both degree-`N` columns in `(17)` acquire the same sign in `(16)`,
it cancels in the determinant, and `omega(h_N)=e_N` yields the exact
complete-to-elementary reciprocal formula

```text
G_N^(P-tau)(mu)
 =det(check E_tau(e_N),check E_tau(omega(m_mu))).           (18)
```

Thus physical prefix complement, together with this virtual endpoint
reindexing, is a genuine involutive transport of the full exterior response.
The `check S_i` need not themselves be physical pole submultisets.  This is
compatible with THM-3160's flat Pluecker pole holotopy.
At the reset, the second coordinate in `(17)` vanishes on every
positive-degree function, so the evaluation plane drops to rank at most one
and all same-degree wedges vanish.  The dual `(7)` proves that this rank-drop
state is the unique point of the depth-nine probability simplex compatible
with the degree-fourteen order inequalities.

## 5. Omega destroys the selector order cone

Formula `(18)` does **not** transport the dual certificate to depth seven,
because the selector cone is phrased in sums of monomials over coarsening
upsets, and `omega` does not preserve that cone.

The smallest clean failure is degree four.  The principal lawful upset
generated by `(2,2)` is `{(4),(2,2)}`.  Since

```text
omega(m_4+m_(2,2))=m_(2,2),                               (19)
```

its image omits the mandatory coarsest partition `(4)`.  It lies outside
both signs of the cone generated by lawful upset indicators.

At degree six the failure is already sign-indefinite.  If `U` is the
principal upset generated by `(3,1,1,1)`, exact change of basis gives

```text
omega(sum_(mu in U)m_mu)
 =-m_6+m_(4,1,1)-m_(3,3)+m_(3,1,1,1).                    (20)
```

The exact companion derives `(19),(20)` from
`p_mu=sum_lambda A_(mu,lambda)m_lambda` and
`omega(p_mu)=(-1)^(N-length(mu))p_mu`; they are not fitted identities.

Consequently the exterior complement `(18)` preserves the source and target
evaluation plane and its determinant, but destroys the order predicate
needed by `(5)`.  This is also why THM-3214's reciprocal two-jet locality
cannot simply dualize the depth-nine cocircuit into a depth-seven positive
certificate: the missing object is a cone-compatible sidecar, not another
jet coefficient.

## 6. Scope and exact evidence

This theorem proves an exact finite selector-polytope statement.  It does not
prove strict Hasse positivity, a nonlinear or Markov factorization, a
positive physical carrier, arbitrary depth, arbitrary support, or a new
Gaussian-moment conclusion.  The omega calculation is a positive theorem on
the exterior response and a no-go on its coarsening-order projection; neither
direction may be silently strengthened.

Run

```text
python 04-computation/gmc_depth_nine_unique_reset_face_omega_boundary_thm3216.py
python -O 04-computation/gmc_depth_nine_unique_reset_face_omega_boundary_thm3216.py
```

and compare LF-normalized bytes with the declared output.  The companion
uses exact integers and `Fraction`s only after importing and hash-pinning the
promoted THM-3209 companion/output.  It uses no discovery pickle, floating
point, randomness, imported executable, or assertion-sensitive test.

Independent audits separately reconstructed the 3,129-state census, every
lawful upset and its size, positivity and primitive digest of the eighteen
multipliers, and the supporting-face implication.  They rederived the
virtual complement signs and both omega cone hostiles, and confirmed the
weak-boundary, depth-ten, and no-physical-carrier scope.  The immutable owner
replays supply the exact 3,129-state sign matrix; normal and optimized output
byte-match the stored transcript and declared LF hashes.

QED.
