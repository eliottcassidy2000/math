# Product-Gamma pole x partition bicomplex: first active Hasse obstruction

Status: **FINITE-EXACT / hostile boundary and algebraic mechanism, not a proved dependency.**

The first failed implication is

```text
scalar pole-flag positivity of the row rational function
    ==> virtual-prefixed type-current refinement.                    (1)
```

At the active top prefix `j=d`, `(1)` already fails on the principal
coarsening upset `{(5)}`.  This note records both the exact commuting square
which made `(1)` plausible and the smallest exact witness which disproves it.

## 1. Inheritance and concept board

The closest proved mechanism is
`THM-3119-factorial-normalized-labelled-deletion-and-young-carrier-order.md`:
after factorial normalization, same-label deletion is stochastic and preserves
Young refinement.  The incoming
`THM-3120-row-pole-prefix-newton-flag-positivity.md` gives a different positive
object, namely a reduced numerator expanded in a descending pole-prefix Newton
flag.  The corrected near miss is
`THM-3122-labelled-deletion-positive-kernel-ghost-and-no-upward-induction.md`:
downward commutation cannot recover a transverse kernel component.  The
least-used sidecar is the exact upset dual in
`THM-3127-partition-refinement-strassen-upset-dual-and-filter-response.md`.

The live board is therefore:

1. factorial-normalized labelled deletion `A`;
2. pole stripping as plethystic virtual-letter subtraction;
3. the signed residual-alphabet functional `Phi`;
4. the distinguished dominant alphabet `Q`;
5. zero-mass Young-type currents `G_mu`;
6. coarsening-upset inequalities; and
7. a missing transverse Young selector.

## 2. The exact commuting square

Put `w_mu=product_i mu_i!` and

```text
Mtilde_N[X]=sum_(mu partition N) m_mu[X]/w_mu e_mu.             (2)
```

THM-3119's same-label deletion identity is

```text
A Mtilde_(N+1)[X]=p_1[X] Mtilde_N[X].                          (3)
```

Multiplication of a row generating function by `(1-Mt)` is not, in general,
deletion of an actually present root.  It is plethystic subtraction of a
virtual letter:

```text
P_M f[X]=f[X-M].                                               (4)
```

Equations `(3)--(4)` give the strict mixed square

```text
A P_M Mtilde_(N+1)[X]
  =(p_1[X]-M)Mtilde_N[X-M]
  =P_M A Mtilde_(N+1)[X].                                     (5)
```

It holds term by term before summing either signed THM-3110 bank, hence for
every pole prefix.  The exact companion checks `735` termwise squares, `882`
row-strip identities, and `30` signed-bank squares on supports `(1,2)`,
`(1,3)`, and `(2,3)`, in degrees through ten.

Literal deletion is sharply false.  At support `(1,2)` the top reduced pole
is `M=5`, but it is absent from `21/24` first-bank atoms and `21/25`
second-bank atoms.  Thus `(5)` is an algebraic plethystic commutation theorem,
not a common physical-root coupling.

## 3. The nontrivial prefixed current

For the reduced pole list `r_1>=...>=r_E`, set

```text
R_j={r_1,...,r_j},
Phi^j(f)=sum_S epsilon_S f[S-R_j],
Q^j=Q-R_j,                                                     (6)

G^j_(N,mu)=Phi^j(h_N)m_mu[Q^j]-Phi^j(m_mu)h_N[Q^j].            (7)
```

The positive scalar/base factors from THM-3115 are immaterial at one fixed
support.  Identity `h_N=sum_mu m_mu` gives

```text
sum_mu G^j_(N,mu)=0.                                          (8)
```

At `j=0`, `(7)` is the pointwise THM-3115 current.  A nonnegative
fine-to-coarse Hasse boundary would imply every coarsening-upset mass is
nonnegative by the exact Strassen/upset dual.  This is the cheapest decisive
test of whether the pole flag lifts from the scalar row to Young type.

There is a useful exact simplification.  In either signed bank, every row
vector has total signed multiplicity zero.  Power sums are additive over the
row alphabets, while both common-root deletion and subtraction of `R_j`
contribute a bank-independent term killed by total signed mass zero.  Hence

```text
Phi^j(p_N)=0                                                   (8a)
```

for every `N`, support, bank, and prefix.  The principal coarsest-upset test
therefore reduces to

```text
G^j_(N,(N))=Phi^j(h_N)p_N[Q^j].                               (8b)
```

Since `p_N[Q^j]=p_N[Q]-sum_(ell<=j)r_ell^N`, every prefix with
`Phi^j(h_N)>0` has the necessary **power budget**

```text
sum_(ell<=j) r_ell^N <= p_N[Q].                               (8c)
```

Equivalently, define

```text
j*_N=max{j: sum_(ell<=j)r_ell^N <= p_N[Q]}.                   (8d)
```

Any Hasse lift of a scalar-positive active prefix must have `j<=j*_N`.
This is a cheap symbolic selector boundary, not an empirical flow statistic.

The companion uses the complete THM-3120 universe

```text
1<=a<=10,        a<b<=min(3a+4,21),                           (9)
```

containing `115` supports and `230` bank cases.  Their reduced numerator
degrees range from `2` to `68`; the `8,241` active prefixes give `41,205`
exact currents for `5<=N<=9`.  All `1,150` `j=0` controls pass.

## 4. First active failure

The lexicographically first active obstruction is

```text
N=5,       (a,b)=(1,3),       bank I2,       j=d=5.            (10)
```

Here

```text
reduced poles = (8,7,6,5,5,4,4,3,3,2,2,2,1,1,1,1),
P(t) = 1440-16440t+51264t^2-37176t^3+84240t^4-325248t^5,
R_5 = (8,7,6,5,5).                                            (11)
```

The dominant alphabet is

```text
Q=(1,3,3,4,5,6,7,8),       Q^5=(1,3,3,4)-5                  (12)
```

Numerically, `p_5[Q]=61987`; the first four pole powers sum to `60476`,
whereas the first five sum to `63601`.  Thus `j*_5=4<d=5`, and the active
top prefix crosses the power budget by exactly `1614`.

Here `Q^5` is interpreted as a virtual alphabet.  Exact evaluation gives

```text
Phi^5(h_5)=1440,    h_5[Q^5]=-358,
Phi^5(p_5)=0,       p_5[Q^5]=-1614.                           (13)
```

Therefore the coefficient of the unique coarsest type is

```text
G^5_(5,(5))
 =Phi^5(h_5)p_5[Q^5]-Phi^5(p_5)h_5[Q^5]
 =-2324160.                                                    (14)
```

The full vector is

```text
(5)             -2324160
(4,1)             544320
(3,2)            2237760
(3,1,1)          -656640
(2,2,1)          -915840
(2,1,1,1)         972000
(1,1,1,1,1)       142560.                                    (15)
```

Since `{(5)}` is already a coarsening upset, `(14)` is a complete symbolic
dual obstruction.  The independent rational max-flow view agrees: it routes
`1,572,480` of demand `3,896,640`, leaving deficit `2,324,160`, exactly the
unpayable coarsest debt.

The full THM-3119 labelled-deletion kernel cannot supply the missing selector.
In the ordered type basis

```text
(5),(4,1),(3,2),(3,1,1),(2,2,1),(2,1,1,1),(1^5),
```

the kernel of `A:P_5->P_4` has the exact basis

```text
K1=(-1,5,-2,-8,6,0,0),
K2=( 1,-5,0,10,0,-10,4),                                    (15a)
```

where `K2` is the THM-3122 derangement ghost.  The exact deletion matrix has
rank five, so these two directions exhaust its kernel.  In raw `Lbar`
coordinates, factorial conjugacy sends them to `WK1,WK2`, with

```text
w=(120,24,12,6,4,2,1),
WK1=(-120,120,-24,-48,24,0,0),
WK2=( 120,-120,0,60,0,-20,4).                                (15b)
```

Both weighted directions have zero mass on the upset `{(5),(4,1)}`, while
the current `(15)` has mass

```text
-2324160+544320=-1779840.                                    (15c)
```

Consequently every raw current with the same deleted image retains this
negative cut.  There is no Hasse-positive preimage anywhere in the complete
factorial deletion-kernel coset, not merely no repair by one derangement
direction.  A successful selector must be genuinely transverse to deletion.
The gauge-correct fixed-fibre theorem is packaged in
`THM-3128-active-pole-prefix-labelled-deletion-no-positive-preimage.md`.

The invariant-cut mechanism is universal.  For every `N>=3`, the row of raw
block deletion `B_N` indexed by target `(N-1)` is exactly the indicator of

```text
U_N={(N),(N-1,1)}.                                            (15d)
```

Indeed, only lowering `(N)` or deleting the singleton from `(N-1,1)` can
produce the one-block target.  Hence `G(U_N)=(B_NG)_(N-1)` is constant on an
entire deletion fibre.  Whenever this top-two upset mass is negative, no
point in that fibre can be a nonnegative Hasse boundary.  The value
`-1779840` above is the first active instance of this reusable obstruction.

The boundary is not merely too many prefixes.  At `(1,2)` the active flag
has `d=2` and its prefixes `j<=2` pass; the first failure there is the
inactive `j=4`.  In `(10)`, by contrast, `j=5` is the active top flag itself.

Nor can one repair the lift by reordering the poles.  At the same support
`(1,3)` and bank `I2`, test every possible one-letter prefix `(r)`.  The
upset `P_N\{(1^N)}` has mass equal to the negative of the singleton
coefficient.  Its first exact failures are

```text
first pole r       first N       - upset mass
     1                9          230041249606656
     2                9           86086680574464
     3                8             164016803808
     4                8             842945556672
     5                8            4210741742976
     6                8           12707727217920
     7                8           24539565851424
     8                8           33677263528128.              (15e)
```

The rational-flow deficit equals the displayed number in every row.  Any
nonempty ordering must choose one of these eight poles first, so no pole
reordering can make a nonempty prefix flag refinement-positive through
degree nine.  In particular, replacing the duplicate `5` by `4` repairs the
degree-five power budget but fails on a higher upset.

## 5. What the obstruction says

The precise connection audit is

```text
source:    positive scalar pole-prefix coefficient Phi^j(h_N)
target:    positive refinement of the zero-mass Young current G^j
map:       virtual subtraction in both Phi and the dominant alphabet Q
preserves: exact degree, zero mass, factorial deletion commutation
destroys:  signs of p_N[Q-R_j] and general coarsening-upset masses
sidecar:   a transverse Young selector controlling Phi^j(m_mu)/m_mu[Q^j]
test:      the principal coarsest upset, i.e. G^j_(N,(N)).      (16)
```

Thus exact commutation plus scalar pole-flag positivity still lacks a
transverse Young selector.  THM-3122 explains abstractly why a deletion
kernel permits this loss; the upset-dual theorem turns the loss into the
one-line witness `(14)`.  The most promising repair is not another larger
flow census, but a pole/filter condition which controls the relative
symmetric-function terms in `(7)`--`(8b)`, or a selector which splits off the
offending virtual negative letter.

The failure of a Hasse certificate does **not** prove the represented
Young-gap operator is non-PSD; the refinement cone is only a sufficient
operator-positive cone.  It also does not weaken THM-3120's scalar row
certificate.

## 6. Reproduction and scope

Run

```text
python 04-computation/gmc_pole_partition_plethystic_commutation_scout.py
python -O 04-computation/gmc_pole_partition_plethystic_commutation_scout.py
python 04-computation/gmc_pole_prefix_hasse_current_scout.py
python -O 04-computation/gmc_pole_prefix_hasse_current_scout.py
```

and compare with the corresponding stored outputs in `05-knowledge/results`.
Only integers and exact rational numbers are used.

The LF-normalized SHA256 hashes are

```text
plethystic script  a096d8feb27a8efeb3d448415ab28270d04b68568bcfaf6c5a27f399ec03c445
plethystic output  5f0100d62a8f83325f22e1667a7e78a73c4d18a3524a0c48ec00ad1fadf45265
Hasse script       151edb9b8cee4807d3f08dc17af32e45420021ba30dfd116c04d9fcaf8bbd5b7
Hasse output       0b4c68b8ba9ed4d32ad51dc6573e10fdd56cea61e1344218879b0c644d9eb856
```

This reflection proves the algebraic square and gives a finite exact hostile
to one proposed lift.  It proves no arbitrary-support pole flag, all-degree
nonrow positivity, product-Gamma conjecture, Gaussian Moment Conjecture,
NC2, LRC(14), JC(2), or DC(2).
