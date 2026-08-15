---
id: THM-3407
title: "Hadamard-core multi-toggle response, plaquette shells, and trade distance"
status: >
  PROVED (exact determinant algebra and self-contained trade-floor proof) +
  VERIFIED-EXACT (three-path Paley controls and hostile triples) + THM-3394
  DEPENDENT FINITE COROLLARY + INDEPENDENTLY AUDITED.  The Hadamard-trade and switching comparisons
  are CITED classical context; no novelty or priority claim is made.
source: hadamard-multitoggle-2026-08-15
audit: determinant-sign/type audit, direct/event/two-support replay, Boolean-Mobius audit, boundary/equality audit, and normal/-O byte-identical controls; independent nonsymmetric H8 and Paley H4/H12 proof/replay audit clean
depends_on:
  - THM-3403-hadamard-core-maxdet-smith-and-circuit-descent
  - THM-3394-twelve-formerly-missing-hadamard-orders-through-2000
dependency_scope: THM-3403 supplies the universal core inverse and binary maxdet equality criterion; THM-3394 is used only for the pinned order-668 specialization in Section 7
related:
  - THM-3393-hadamard-order-668-explicit-certificate
  - THM-447-skew-sylvester-doubling
  - THM-451-skew-tower-hadamard-chirality
script: 04-computation/hadamard_core_multitoggle_response_thm3407.py
output: 05-knowledge/results/hadamard_core_multitoggle_response_thm3407.out
script_sha256: 83f38703b126cdc9bf9358bce20803ed890a1f19d685c7644e909c1a8df45e90
output_sha256: 46b6b6176a23ebb0554c0359315f5584227576c753b5e3e1d1f6cc8ce1f1d31f
semantic_sha256: 6333b6b8292916e6917b9d04c82bb7b4bf86d6fe2c9dfda45acfce5e3eb98a09
hash_basis: working-tree bytes with LF line endings
---

# THM-3407 -- Hadamard-core multi-toggle response, plaquette shells, and trade distance

**PROVED + VERIFIED-EXACT + THM-3394 DEPENDENT FINITE COROLLARY +
INDEPENDENTLY AUDITED.**

Let `m>=1`, put `N=4m` and `v=N-1`, and normalize a real Hadamard
matrix as

~~~text
    [ 1   1^T ]
H = [           ],              B=(J-K)/2 in {0,1}^{v by v}.       (1)
    [ 1    K  ]
~~~

THM-3403 proves

~~~text
B^(-1)=-K^T/(2m),              |det B|=2m^(2m),                  (2)
~~~

and proves that this is the global binary determinant maximum, with equality
exactly when the bordered sign matrix is Hadamard.  No symmetry of `K` is
assumed anywhere below.

## 1. Two exact sparse-response compilers

Let `T` be a set of `t` distinct positions in the core.  If `M` is its
zero-one mask, put

~~~text
E=K circle M,                  B_T=B+E.                           (3)
~~~

Here `circle` is entrywise product.  The sign is load-bearing: toggling
`B_ij` adds `K_ij`, and flips the corresponding entry of `K`.

Let `R,C` be the distinct rows and columns met by `T`, in any fixed orders,
and let

~~~text
L=K[R,C],                      D=E[R,C].                          (4)
~~~

If `r=|R|` and `c=|C|`, then the normalized signed determinant response is

~~~text
rho(T):=det(B_T)/det(B)
      =det(I_r-DL^T/(2m))
      =det(I_c-L^TD/(2m)).                                      (5)
~~~

Thus the cheaper of `r` and `c` gives a support-compressed exact compiler,
including repeated toggle rows or columns.

There is also an event compiler.  Index the events as
`e_a=(i_a,j_a)`, put `delta_a=K_(i_a,j_a)`, and define

~~~text
Q_ab=delta_a K_(i_b,j_a),                  Q_aa=1.                (6)
~~~

Then

~~~text
rho(T)=det(I_t-Q/(2m))=(2m)^(-t)det(2mI_t-Q).                    (7)
~~~

Reordering the events conjugates `Q` by a permutation.  Transposition and
moving `diag(delta)` give determinant-equivalent gauges, but (6) is the
orientation convention used in every cycle formula below.

### Proof of (5)--(7)

Let `P_R,P_C` be the row- and column-inclusion matrices.  Then

~~~text
B_T=B+P_R D P_C^T,
P_C^T B^(-1)P_R=-L^T/(2m).
~~~

The generalized determinant lemma and Sylvester's identity give both forms
in (5).  Alternatively write

~~~text
E=sum_a delta_a e_(i_a)e_(j_a)^T.
~~~

The event update matrix has entries

~~~text
delta_a (B^(-1))_(j_a,i_b)
       =-delta_a K_(i_b,j_a)/(2m),
~~~

which is exactly `-Q/(2m)` and proves (7).  This also audits the transpose in
(2); replacing `K^T` by `K` would generally be wrong.

## 2. Boolean Mobius coefficients are signed core minors

For `S subseteq {1,...,t}`, in the inherited event order, put

~~~text
mu_S=(product_(a in S) delta_a)
     det K[(i_a)_(a in S),(j_a)_(a in S)],       mu_empty=1.      (8)
~~~

Principal-minor expansion of (7) gives

~~~text
rho(T)=sum_(S subseteq [t]) (-1/(2m))^|S| mu_S.                  (9)
~~~

More strongly, give every candidate event an indeterminate `z_a`.  The
complete multilinear determinant polynomial is

~~~text
F(z)=det(B+sum_a z_a delta_a e_(i_a)e_(j_a)^T)/det(B)
    =sum_(S subseteq [t]) (-1/(2m))^|S| mu_S product_(a in S)z_a. (10)
~~~

Hence the coefficient in (10) is exactly the Boolean Mobius interaction of
`S`, not a Taylor approximation.  If `S` repeats a row or column then its
minor vanishes.  Therefore only matchings in the event row-column bipartite
graph can have nonzero interaction.  Some matchings can still vanish: their
signed minor is the exact remaining sidecar.

The map here is

~~~text
labelled toggle set T  ->  multilinear response polynomial F.
~~~

It preserves every toggled determinant and forgets the toggled matrix itself.
The packet `(mu_S)_S` is response-complete; the oriented event matrix `Q` is a
convenient sufficient lift, but principal minors are not claimed to reconstruct
`Q` entrywise.

## 3. Complete two-toggle shells and multiplicities

For distinct events `e_1=(i,j)` and `e_2=(k,l)`, define the plaquette sign

~~~text
chi=Q_12 Q_21=K_ij K_kl K_kj K_il in {+1,-1}.                   (11)
~~~

If the events share a row or column then `chi=+1`.  Otherwise it is the sign
of their `2 by 2` plaquette.  Taking the determinant in (7) gives

~~~text
rho_2=((2m-1)^2-chi)/(2m)^2.                                  (12)
~~~

Thus there are exactly two determinant magnitudes:

~~~text
chi=+1:
  rho=(m-1)/m,
  |det B_T|=2(m-1)m^(2m-1),
  loss from max=2m^(2m-1);

chi=-1:
  rho=(2m^2-2m+1)/(2m^2),
  |det B_T|=(2m^2-2m+1)m^(2m-2),
  loss from max=(2m-1)m^(2m-2).                               (13)
~~~

The negative-plaquette shell is higher by exactly `m^(2m-2)`.  At `m=1`
the positive shell is singular and the negative shell has determinant one.
No nonempty two-toggle set is a second maximizer.

After subtracting the two one-bit tariffs, the exact second interaction is

~~~text
rho({1,2})-rho({1})-rho({2})+rho(empty)
       =(1-chi)/(2m)^2.                                      (14)
~~~

It is zero on positive plaquettes and `2/(2m)^2` on negative ones.

The shell multiplicities depend only on `m`.  For any two core rows, their
entrywise product has `2m-1` plus entries and `2m` minus entries: their full
Hadamard inner product is zero and the common border entry is plus.  Put
`C_v=binom(v,2)`.  Among all `binom(v^2,2)` unordered event pairs,

~~~text
n_low =2 C_v [v+(2m-1)^2],
n_high=2 C_v [2m(2m-1)].                                    (15)
~~~

Indeed the first `2v C_v` pairs share a row or column.  For each unordered
row pair, two matchings on a column pair contribute, and equal versus opposite
row-product signs give `(2m-1)^2` versus `2m(2m-1)` column pairs.

If `m` is odd, THM-3403's mod-two circuit has a determinant-level shadow:
the high determinant in (13) is odd, while the low determinant is divisible
by four.  One toggle repairs the unique mod-two corank.  A second toggle
preserves that repair exactly on a negative plaquette and restores singularity
modulo two on a positive plaquette.

## 4. Orientation first appears at three events

For three events write

~~~text
chi_ab=Q_ab Q_ba,
gamma_abc=Q_ab Q_bc Q_ca,
gamma_acb=Q_ac Q_cb Q_ba.
~~~

Then

~~~text
mu_123=det Q
      =1-chi_12-chi_13-chi_23+gamma_123+gamma_132,             (16)

rho_3=1-3/(2m)
      +[(1-chi_12)+(1-chi_13)+(1-chi_23)]/(2m)^2
      -mu_123/(2m)^3.                                         (17)
~~~

If the three events use distinct rows and columns, `mu_123` is a signed
`3 by 3` minor of `K`, so it lies in `{0,+4,-4}`.  This follows because a
`3 by 3` sign determinant is divisible by four after subtracting one row
from the other two, while Hadamard's bound is less than six.  A repeated row
or column makes the minor zero.

The orientation debt is exactly one bit.  Indeed

~~~text
gamma_123 gamma_132=chi_12 chi_13 chi_23.                    (18)
~~~

Let `s` be the number of negative `chi` values.  Equations (16) and (18)
give the complete abstract palette

| `s` | possible `mu_123` | extra datum needed beyond labelled `chi` |
|---:|:---:|:---|
| 0 | `0,-4` | `gamma_123` |
| 1 | `0` | none |
| 2 | `0,+4` | `gamma_123` |
| 3 | `+4` | none |

Thus labelled pair data determines the triple response when its plaquette
product is negative (`s` odd).  When that product is positive (`s` even),
the directed 3-cycle gain selects one of two algebraic shells.  The exact
companion exhausts all `2^6=64` oriented off-diagonal sign assignments and
checks this palette; realizability of every abstract state inside a specified
Hadamard core is not asserted.

Label pair data in the order `(12,13,23)`.  The independently constructed
Paley-I controls contain the following exact hostile pairs; coordinates are
zero-based positions in `K`:

| `N` | first triple / second triple | common pair data | `mu` | `rho` | `|det|` |
|---:|:---|:---:|:---:|:---:|:---:|
| 8 | `((0,0),(1,2),(2,1))` / `((0,0),(1,1),(2,3))` | `(-1,+1,-1)` | `4 / 0` | `7/16 / 1/2` | `14 / 16` |
| 12 | `((0,2),(1,1),(2,0))` / `((0,1),(1,0),(2,3))` | `(-1,-1,+1)` | `4 / 0` | `16/27 / 11/18` | `864 / 891` |
| 20 | `((0,2),(1,1),(2,0))` / `((0,0),(1,1),(2,3))` | `(-1,-1,+1)` | `4 / 0` | `92/125 / 37/50` | `14,375,000 / 14,453,125` |

Thus even the complete **labelled** pair-plaquette packet does not determine
the three-toggle determinant.  The first destroyed coordinate is oriented
cycle gain, equivalently the signed `3 by 3` minor.  For general `t`, the
determinant of `Q` is a signed cycle-cover sum: two-cycles are plaquettes,
while directed cycles of length at least three supply the higher orientation
sidecar.

## 5. Rank-one rectangles, the singular wall, and equality trades

Suppose every entry of a full `a by b` rectangle `R by C` is toggled.  Then
`D=L` in (5), so

~~~text
rho=det(I_a-LL^T/(2m)).                                      (19)
~~~

If `L=uw^T` has real rank one, with sign vectors `u,w`, then

~~~text
rho=1-ab/(2m)=1-2ab/N.                                      (20)
~~~

THM-3403's global binary maxdet bound gives `|rho|<=1`; consequently every
rank-one sign rectangle in a Hadamard core has

~~~text
ab<=N.                                                       (21)
~~~

The three exact regimes are

~~~text
0<ab<N:   strict determinant loss;
ab=N/2:   singular toggled core;
ab=N:     rho=-1 and another Hadamard maximizer.              (22)
~~~

A `t`-toggle set in one row or one column is the special case

~~~text
rho=1-t/(2m).                                                (23)
~~~

It is singular at `t=2m`; at `t=v` its signed response is
`-(1-1/(2m))`, the same magnitude as a one-toggle neighbor.

Every pair of lower Hadamard rows differs in `N/2=2m` core columns.  Toggling
those columns in both rows exchanges the two rows.  The changed block is a
rank-one `2 by 2m` rectangle of area `N`, proving a universal `N`-entry trade
with `rho=-1`; columns give the dual trade.

If `N>4` and four Hadamard rows form a closed quadruple, choose the normalized
border row outside it.  Orthogonality partitions the columns of the quadruple
into four sign fields, modulo simultaneous negation, each of size `m`.  Choose
a field not containing the border column.  Its `4 by m` core block has rank
one and area `N`; flipping it gives `rho=-1`.  This is the determinant-response
form of closed-quadruple switching.

## 6. Sharp equality floor and labelled Hamming distance

For every nonempty core toggle set,

~~~text
|rho(T)|<=1,
|rho(T)|=1  iff  the bordered flipped sign matrix is Hadamard. (24)
~~~

Moreover equality is impossible when `|T|<N`.  This trade floor has the
following self-contained proof.  Let two real Hadamard matrices `H,H'` differ
on exactly `t` entries.  Let `r` be the number of changed rows and choose a
changed row `i` having the minimum number `a` of changed entries.  Let `S` be
its changed columns and define `x` by

~~~text
x_c=H_ic for c in S,              x_c=0 otherwise,
y=Hx.                                                          (25)
~~~

Then `||y||^2=N||x||^2=Na` and `|y_j|<=a`.  If row `j` is unchanged,
orthogonality of that row with row `i` both before and after the flip gives
`y_j=0`.  Hence `support(y)<=r`, and

~~~text
Na=||y||^2<=r a^2,              ra>=N,              t>=ra>=N. (26)
~~~

The row-swap trades above attain equality.  Therefore, inside the labelled
cube `{0,1}^{v^2}`, the set of binary core maximizers descended with the fixed
all-plus border has exact minimum Hamming distance

~~~text
d_min=N=4m.                                                   (27)
~~~

Equation (24) follows directly from THM-3403's border determinant identity
and equality case.  Equations (26)--(27) prove the sharp distance and do not
depend on the external trade citation in Section 9.

## 7. Odd full circuits forbid closed quadruples but not higher shells

For lower core rows define the pair-product signature

~~~text
x_(rs)(j)=K_rj K_sj.                                         (28)
~~~

The event response samples these signatures:

~~~text
Q_ab=x_(i_a,i_b)(j_a).                                       (29)
~~~

A closed quadruple is exactly a collision of two disjoint unordered
pair-product signatures (including the all-plus border row when appropriate).

When `m` is odd, THM-3403 proves that the lower binary rows form one full
binary circuit.  Therefore

~~~text
S subseteq {lower rows}  ->  product_(i in S) K_i             (30)
~~~

is injective modulo complementation: two products agree exactly when their
index sets agree or are complements.  If `m>1`, then `v>=11`, while the
symmetric difference of two row pairs has size at most four.  Hence every
unordered pair-product signature among all `N` normalized rows is distinct.
There is no closed row quadruple; the column statement is identical.  The
excluded `m=1` case is sharp at `H_4`.

The binary circuit does **not** determine the sampled signs in (28), directed
cycle gains in (16), or higher real minors in (8).  It therefore forbids the
closed-quadruple `4 by m` trade for odd `m>1`, but neither classifies higher
determinant shells nor removes the universal `2 by 2m` row/column-swap trades.
The order-8/12/20 hostile triples in Section 4 are the cheapest decisive test
for the missing coordinate.

## 8. Pinned order-668 specialization

This section, and only this section, uses THM-3394.  Its first explicit
Hadamard certificate has `N=668`, `m=167`, normalized text SHA-256

~~~text
73f1de1539849e1dc7e6085cc69c563fd2965c44970263e8203384bd1a46aa63.
~~~

Any two lower rows have core-product counts `333` plus and `334` minus.
Sections 1--6 specialize to

~~~text
max                 =2*167^334,
one toggle          =333*167^333,             rho=333/334,
two-toggle low      =332*167^333,             rho=166/167,
two-toggle high     =55445*167^332,           rho=55445/55778,
high-low difference =167^332.                                  (31)
~~~

The exact shell multiplicities are

~~~text
low =49,555,629,432,
high=49,407,259,284,
sum =98,962,888,716=binom(667^2,2).                            (32)
~~~

The actual difference set of any two lower rows has `334` core columns, so
their `2 by 334` row swap is a `668`-toggle equality trade with `rho=-1`.
The companion hash-pins the THM-3394 theorem, source, data, and output and
checks the first output record; it does not reconstruct a dense `667 by 667`
binary determinant.

## 9. Exact finite companion

The standard-library companion constructs Paley-I matrices independently at
orders `4,8,12,20`.  Its explicit universe and controls are:

- all `binom(v^2,2)` event pairs at all four orders are checked by the event
  and both support compilers; direct Bareiss determinants check every pair
  through order `12` and all four structural representatives at order `20`;
- every one of the `512` masks in a fixed `3 by 3` core window is checked by
  direct determinant, event response, and both support responses, and its
  entire Boolean Mobius transform is compared with (8)--(10);
- all `64` oriented `3 by 3` event sign matrices with diagonal one are
  exhausted against the parity/orientation palette following (18);
- all same-row and same-column prefix sizes `0<=t<=v`, universal row and
  column swaps, and the explicit order-8 closed-quadruple field are checked;
- the three hostile pairs in Section 4 are frozen, including their literal
  determinants; and
- the order-668 dependency bytes, normalized matrix hash, formulas, shell
  counts, and row-swap dimensions are pinned.

The order-4 full-core window has exactly three equality masks of size four and
two of size six; its smallest is the sharp `N=4` floor.  The fixed windows at
orders `8,12,20` contain no equality mask.  These are positive and boundary
controls, not an enumeration of all Hadamard trades at those orders.

Reproduce from the repository root with

~~~bash
python3 04-computation/hadamard_core_multitoggle_response_thm3407.py
python3 -O 04-computation/hadamard_core_multitoggle_response_thm3407.py
~~~

The two outputs are byte-identical to the stored output.  The companion uses
no floats, randomness, solver, network, subprocess, dynamic import, `assert`,
or repository-file writes.

## 10. Classical context, synthesis boundary, and non-consequences

The following are citations for context, not proof dependencies and not
priority claims.

- Ó Catháin and Wanless,
  [*Trades in complex Hadamard matrices*](https://arxiv.org/abs/1502.02353),
  Theorem 4, prove that every trade in a real Hadamard matrix of order `N`
  has at least `N` entries.  Their introduction also records the universal
  `2 by N/2` row trades and `4 by N/4` closed-quadruple trades.  Equation (26)
  gives a self-contained proof of the real trade floor needed here.
- Orrick,
  [*Switching operations for Hadamard matrices*](https://arxiv.org/abs/math/0507515),
  is the classical switching and closed-quadruple context.
- Brent and Osborn,
  [*On minors of maximal determinant matrices*](https://arxiv.org/abs/1208.3819),
  and the Jacobi/Szollosi complementary-minor route are adjacent determinant-
  minor literature.  They are not used to derive (5)--(18).

The support/event compilers, Boolean response packaging, two-shell accounting,
and oriented triple controls are presented here as an elementary repository
synthesis of standard determinant tools with THM-3403.  No assertion is made
that this packaging is new in the literature.

This theorem classifies the complete two-toggle determinant shell and gives a
response-complete invariant at arbitrary toggle rank.  It does **not**
classify Hadamard equivalence classes, all equality trades, or all higher
determinant shells.  It proves no new Hadamard order, no Hadamard conjecture,
and no LRC(14), JC(2), Crouzeix, tournament, or coding-theory consequence.
