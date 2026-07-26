---
id: THM-2370
title: "Deletion-martingale drift conservation and the sharp clone hostile"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. For any
  nested sequence of nonnegative THM-2365 tables, the noncirculant
  component lost between the first and last tables is the sum of the
  projected deletion layers. If delta is the positive difference of
  their drift norms, one layer has drift energy at least delta^2/n^2
  and the sum of all layer energies is at least delta^2/n. A lawful
  physical deletion layer therefore lands its own nonzero target, with
  a finite-colour amplitude floor delta/(n sqrt(26208)) and an exact
  91-unit deep multiplier. The n^-2 and n^-1 Hilbert constants are sharp
  even for commuting Boolean rational-interval masks whose cumulative
  products are nested, with disjoint physical deletion supports and
  individual mass 1-1/(91n): cloning the exact 90/91 THM-2367
  circulant-restoration mask gives equality. On the full
  Hilbert-valued Boolean clone cube, the normalized Walsh expansion
  has only its constant and n equal singleton terms; every mixed ANOVA
  tensor and every disjoint pair Dirichlet defect vanishes. Global
  Boolean complementation preserves the complete squared Dirichlet
  bank while swapping the target-null terminal C with the drifting
  terminal T. Thus a complete squared interaction spectrum still loses
  the absolute terminal orientation. The scalar table of squared drift
  instead has explicit nonzero degree-two terms; it still does not
  restore phase. The landed
  target belongs to the deleted packet, not automatically to the
  canonical full-owner word. A complementary sharp anchored-slice lemma
  shows that a shifted all-safe residual which vanishes on its original
  scalar-cover slice is either identically zero or has drift at least
  1/169 of its squared norm; thus a nonzero circulant ghost residual is
  impossible. Guard/unit failure routing and canonical current transfer
  remain open. No scalar-row exclusion, ledger decrement, or LRC(14)
  consequence follows.
source: codex-2026-07-25-deletion-martingale-drift
depends_on:
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
related:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2367-septimal-root-averaging-graft-and-cover-alignment
  - THM-2374-binary-allocation-complete-subcube-dirichlet-spectrum
script: 04-computation/lrc14_deletion_martingale_drift_thm2370.py
output: 05-knowledge/results/lrc14_deletion_martingale_drift_thm2370.out
script_sha256: 34d60f9111963e1400ec77a0842806f6ad415430d5be00c50429da4b5a41f75d
output_sha256: 82ed7199feff21923d79e9afdaef72e40d65cd61ab5cb4ac324ed052c0648fd4
hash_basis: working-tree bytes (LF)
---

# THM-2370 -- drift cannot disappear without entering a deletion layer

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2367 gives both a locally drifting lawful graft and a large Boolean
mask which makes the final overlap tensor exactly circulant. The natural
question is whether the drift can disappear invisibly while owner factors
are inserted.

It cannot. The missing drift must enter at least one deleted layer:

```text
initial packet = final packet disjoint-union deletion layers

noncirculant part of initial packet
  = noncirculant part of final packet
    + sum of noncirculant deletion parts.                    (1)
```

The resulting quantitative bound is optimal, including in the positive
Boolean interval category. Its limitation is equally exact: the surviving
target is a target of the **deleted packet**, which need not be a canonical
full-owner word.

## 1. The orthogonal deletion identity

Let

```text
H_0,H_1,...,H_n : F_13^3 -> R_(>=0)
```

be tables with

```text
H_j<=H_(j-1),

L_j=H_(j-1)-H_j>=0,                    j=1,...,n.   (2)
```

Assume every table has the THM-2365 diagonal zero

```text
H_j(t,s,t)=0
```

up to the null strict-open endpoints. The same is then true for every
`L_j`.

Use normalized counting measure and let `P` be THM-2365's target-action
projection

```text
(P H)(r,s,t)
 =13^(-2)sum_(u,v in F_13)H(r+v,s+u,t+v).          (3)
```

Put

```text
Q=I-P,

D(H)=||QH||_2^2
    =13^(-3)sum_(r,s,t)|H(r,s,t)-(PH)(r,s,t)|^2.   (4)
```

This is exactly the nonzero-target drift energy of THM-2365. Telescoping
(2) and applying `Q` gives

```text
QH_0-QH_n=sum_(j=1)^n QL_j.                        (5)
```

Define

```text
delta
 =(
    sqrt(D(H_0))-sqrt(D(H_n))
  )_+.                                              (6)
```

The reverse triangle inequality, the triangle inequality, and
Cauchy--Schwarz give

```text
delta
 <=||QH_0-QH_n||_2
 <=sum_j||QL_j||_2,                                (7)

max_j D(L_j)>=delta^2/n^2,                         (8)

sum_j D(L_j)>=delta^2/n.                           (9)
```

In particular, if `H_n` is circulant and `H_0` is not, then at least one
deletion layer is noncirculant. No monotonicity of `D(H_j)` is assumed or
needed.

The constants in (8)--(9) are the Hilbert-space optimum. Simultaneous
equality is possible when

```text
QH_n=0,

QL_1=...=QL_n=QH_0/n.                              (10)
```

Section 5 realizes (10) with Boolean interval masks whose cumulative
products are nested.

## 2. A lawful deletion lands its own target

Suppose now that the tables arise before integration. Let
`B_(s,t)(x)>=0` be a lawful THM-2365 base packet, and let

```text
M_1^(s,t),...,M_n^(s,t) in {0,1}
```

be lawfully co-shifted Boolean factors. Insert them successively:

```text
H_j(r,s,t)
 =integral_T
   B_(s,t)(x) Delta_r(x)
   product_(i=1)^j M_i^(s,t)(x) dx.                (11)
```

Then the deleted layer has the positive physical factorization

```text
L_j(r,s,t)
 =integral_T
   B_(s,t)(x) Delta_r(x)
   product_(i<j)M_i^(s,t)(x)
   (1-M_j^(s,t)(x)) dx.                            (12)
```

If (12) retains the same target quotient, diagonal zero, and
Poisson--Abel typing as THM-2365, that theorem applies to `L_j`
itself. In particular, deleting the deepest safe factor which enforces
the diagonal zero is outside this clause unless another retained factor
restores that zero. Write `B_j(a,b,h)` for the normalized finite transform
of a lawful layer. For the layer selected by (8),

```text
sum_(
  a!=0,
  (b,a+h)!=(0,0)
 ) |B_j(a,b,h)|^2
 >=D(L_j)/13
 >=delta^2/(13n^2).                                (13)
```

There are exactly

```text
(13-1)(13^2-1)=2016
```

eligible deep/target colours. Hence some coefficient satisfies

```text
|B_j(a,b,h)|
 >=delta/(n sqrt(13*2016))
 =delta/(n sqrt(26208)).                           (14)
```

Summing the THM-2365 energy bound over all layers also gives

```text
sum_j sum_(
  a!=0,
  (b,a+h)!=(0,0)
 ) |B_j(a,b,h)|^2
 >=delta^2/(13n).                                  (15)
```

The ordinary absolute `m`-then-`X` collapse of THM-2365 applies to a
nonzero coefficient in (13). It produces exact integers `m,X` and a
nonzero fixed-triangle target fibre for the deleted packet, with

```text
m=a mod 13,

gcd(m,91)=1,

q=(b,a+h)!=0.                                      (16)
```

Equation (16) does not say that the same triangle survives after the
deleted factor is restored.

## 3. A zero base slice cannot support a hidden circulant ghost

There is one especially useful deletion boundary on which much more can
be said. Let `Z:F_13^3->C` be any table satisfying

```text
Z(r,0,0)=0                         for every r.     (16a)
```

It need not be nonnegative for the Hilbert-space statement. Then

```text
D(Z)>=||Z||_2^2/13^2.                              (16b)
```

### Proof

Put `G=PZ`, so `G(r,s,t)=g(r-t)` for some `g`. Let `A` be the orthogonal
projection which zeros the thirteen-cell slice

```text
S={(r,0,0):r in F_13}.
```

Every difference class meets `S` once and has `13^2` cells. Consequently

```text
||AG||_2^2=(1-13^(-2))||G||_2^2.                  (16c)
```

Since `Z=AZ` and `G=PZ`,

```text
||G||_2^2
 =<Z,G>
 =<Z,AG>
 <=||Z||_2||AG||_2.
```

Thus

```text
||G||_2^2
 <=(1-13^(-2))||Z||_2^2.
```

Orthogonality of `P` now gives

```text
D(Z)
 =||Z||_2^2-||G||_2^2
 >=13^(-2)||Z||_2^2,
```

which proves (16b).

The constant is sharp. Choose any nonzero profile `g>=0` with `g(0)=0`,
put

```text
C(r,s,t)=g(r-t),

Z=AC.                                               (16d)
```

Then `Z` is nonnegative, has the THM-2365 diagonal zero, and

```text
PZ=(1-13^(-2))C,

D(Z)=||Z||_2^2/13^2.                               (16e)
```

Combining (16b) with THM-2365 gives, for a lawful physical `Z`,

```text
sum_(
  a!=0,
  (b,a+h)!=(0,0)
 )|B_Z(a,b,h)|^2
 >=||Z||_2^2/13^3.                                 (16f)
```

Thus every nonzero lawful anchored-zero packet lands its own nonzero
target and an exact `91`-unit deep triangle.

This applies directly to the all-safe shifted residual. Partition a
lawfully shifted present packet by its complete blocker word, retaining
the deepest safe factor which supplies the diagonal zero. The original
scalar cover makes the empty-word table satisfy (16a). Therefore:

```text
empty word stays absent on every target shift;

or

the reappearing empty word has target drift
  at least ||Z||_2^2/169.                           (16g)
```

There is no third possibility in which a positive shifted all-safe ghost
is nonzero and circulant. This disposes of that residual as a possible
sink for invisible drift. Its landed current is still a derived empty-word
current, not a canonical positive owner word.

## 4. The exact Boolean cancellation seed

The sharpness construction begins with the exact THM-2367 interval
hostile, reproduced independently here.

Put

```text
D=16562=2*7^2*13^2
```

and let `W` be one on every half-open grid cell

```text
[j/D,(j+1)/D)
```

except the ten ranges

```text
[16555,16562), [0,13),
[1625,1651),   [2457,2463),
[3263,3289),   [4907,4927),
[7449,7455),   [9087,9113),
[10725,10751), [12363,12389).                     (17)
```

The first two ranges form one circular interval. Exactly `182` cells are
removed, so

```text
mu(W)=90/91.                                       (18)
```

For `r,t in F_13`, define

```text
T(r,t)
 =integral_T
   d(13x-r/13)g(13x-t/13)g(7x+t/13) dx,           (19)

C(r,t)
 =integral_T
   W(x)d(13x-r/13)g(13x-t/13)g(7x+t/13) dx.       (20)
```

Embed these as three-coordinate tables independent of `s`. Exact cell
counting gives

```text
T>=C>=0,

T(t,t)=C(t,t)=0,                                   (21)

C(r,t)=
  0          if r-t=0,
  11/169     if r-t=+1 or -1,
  11/91      otherwise.                            (22)
```

Thus `C` is circulant. The unmasked table is not, and its entire
noncirculant part lies in the deleted tensor:

```text
D(T)=D(T-C)
    =852/11589168409,

D(C)=0,

Q(T-C)=QT.                                         (23)
```

The grid in (17) resolves every strict comb boundary, so the midpoint
count is an exact rational-interval computation, not a sampling
approximation.

## 5. Clone every physical cell: equality for all n

Fix any positive integer `n`. Subdivide every grid cell in Section 4
into `n` equal clone subcells. For `j=1,...,n`, let `U_j` be the union
of the `j`-th clone subcell inside each of the `182` excluded cells, and
put

```text
M_j=1-1_(U_j).                                      (24)
```

The masks are Boolean, rational-interval, target-independent, commuting,
and their deleted physical supports are pairwise disjoint. Each one has

```text
mu(M_j)=1-1/(91n),                                 (25)
```

while

```text
product_(j=1)^n M_j=W,

mu(product_j M_j)=90/91.                           (26)
```

Let `H_j` be (19) after inserting `M_1...M_j`. Every comb factor is
constant almost everywhere on an original grid cell, so the clone
decomposition is exact:

```text
H_j
 =C+(n-j)(T-C)/n,                                  (27)

L_j=H_(j-1)-H_j=(T-C)/n.                           (28)
```

Equations (23), (27), and (28) yield

```text
D(H_j)=((n-j)/n)^2 D(T),                           (29)

D(L_j)=D(T)/n^2               for every j,         (30)

sum_jD(L_j)=D(T)/n.                                (31)
```

Thus both (8) and (9) are sharp simultaneously for every `n`. This
remains true despite:

```text
nonnegativity,
Boolean masks,
nested cumulative packets,
commuting target-independent factors,
pairwise-disjoint physical deletion supports,
individual mask mass tending to one.
```

The reason disjoint physical supports do not improve `n^-2` is now
visible: integration sends all `n` clone layers to the **same** projected
table `QT/n`. Their target-drift vectors are perfectly coherent rather
than orthogonal.

Unequal clone widths make both inequalities strict. Omitting the last
clone leaves exactly `D(T)/n^2` terminal drift. These are exact hostile
controls in the companion.

## 5a. The full clone cube has no mixed ANOVA

The nested order is not hiding a higher Boolean interaction. Keep every
subset of the clone masks. For

```text
x=(x_1,...,x_n) in F_2^n,

H_x=the packet after applying M_j exactly when x_j=1.
```

The deleted clone supports are disjoint, so (27) strengthens to the
whole-cube identity

```text
H_x=C+(1-|x|/n)(T-C).                              (31a)
```

Fix the character convention

```text
chi_A(x)=(-1)^(sum_(j in A)x_j),

Hhat(A)=2^(-n)sum_x H_x chi_A(x).                  (31b)
```

Directly expanding (31a) gives

```text
Hhat(empty)=(T+C)/2,

Hhat({j})=(T-C)/(2n)                 for every j,

Hhat(A)=0                            for |A|>=2.    (31c)
```

The sign in the singleton line is positive because

```text
2^(-n)sum_x x_j(-1)^x_j=-1/2.
```

After applying the noncirculant projection `Q=I-P`, (23) turns (31c)
into

```text
QHhat(empty)=QT/2,

QHhat({j})=QT/(2n),

QHhat(A)=0                            for |A|>=2.   (31d)
```

Thus all singleton deletion vectors are not merely coherent: they are
equal and collinear.

For completeness, let `U subset {1,...,n}` and average the squared
target difference over the whole coordinate subcube:

```text
mathcal D_U
 =E_(x in F_2^n,h in F_2^U)
    ||QH_(x+h)-QH_x||_2^2.
```

The elementary `0/2` Walsh multiplier and (31d) give

```text
mathcal D_U
 =2 sum_(A:A intersect U!=empty)||QHhat(A)||_2^2
 =|U|D(T)/(2n^2).                                  (31e)
```

Consequently, for every disjoint `U,V`,

```text
mathcal D_U+mathcal D_V-mathcal D_(U union V)=0.
```

Every mixed rectangle defect and every Walsh component of degree at
least two vanishes even though

```text
D(H_0)=D(T)>0,                 D(H_(1,...,1))=D(C)=0.
```

The **complete squared** subcube bank is also insufficient to identify
which endpoint is terminal. Global Boolean complementation gives the
equally physical cube

```text
Htilde_x=H_(1-x)=C+(|x|/n)(T-C).
```

It has exactly the same values of every `mathcal D_U`, equivalently the
same **nonconstant** squared Walsh spectrum; in this example the constant
coefficient also agrees. Yet

```text
Htilde_(1,...,1)=T
```

drifts. The forgotten coordinate is the applied/unapplied orientation,
or equivalently the constant and signed singleton data relative to an
absolute terminal reference. A mixed-ANOVA or rectangle certificate
cannot supply the owner/current transfer missing in Section 7.

This is a statement about the Hilbert-valued packet/current cube
`x -> QH_x`. It must not be conflated with the scalar quadratic table

```text
e(x)=D(H_x)=D(T)(1-|x|/n)^2.
```

That scalar table has normalized Walsh coefficients

```text
ehat(empty)=D(T)(n+1)/(4n),

ehat({j})=D(T)/(2n),

ehat({i,j})=D(T)/(2n^2),

ehat(A)=0                              for |A|>=3,
```

and therefore has positive scalar pair defects

```text
J_(U,V)=|U||V|D(T)^2/(2n^4)
```

for disjoint `U,V`. Squaring before taking the Boolean transform creates
an interaction that was absent from the current itself. It does not
restore the terminal phase or turn the scalar spectrum into a canonical
owner service.

## 6. What the martingale does and does not decompose

At the overlap-table and finite-transform levels, (5) is an exact linear
decomposition. There is also a useful fixed-frequency form. Keep the same
bare right endpoint from the isolated packet `H_0`; the retained left
packet and all deletion left packets partition the original left packet.
With that common right endpoint, every fixed `(X,m)` marked current
decomposes linearly at the indicator boundary. Equivalently, one may first
partition and then Poisson-smooth each whole layer. Separately smoothing
the Boolean factors at finite radius does not preserve their pointwise
partition identity.

This does **not** make a deletion current a subcurrent of the canonical
full-owner current. The latter uses the fully masked bare endpoint.
Replacing the common right endpoint by that endpoint changes the
off-diagonal Fourier products and introduces cross-layer terms.

The Boolean owner typing gives a second obstruction. At the unshifted
cell, scalar cover may route a blocker deletion into canonical
THM-2305 owner words. After arbitrary lawful target translations, the
all-safe blocker cell can reappear: scalar-cover routing is not currently
known to commute with the target action. Therefore the theorem proves

```text
full packet drifts
  or
some deleted packet drifts,                                (32)
```

not

```text
some canonical full-owner word drifts.
```

For the usual factor order:

- deleting a guard-safe or `q_i`-safe factor lands in a guard/unit danger
  packet, generally outside the scalar safe core;
- deleting an owner-danger factor after the safe core is installed gives
  an owner-complement layer, but its shifted fibres still require a
  cover-covariant word decomposition;
- deleting a nonowner-safe blocker produces a simultaneous blocker
  danger/collision stratum, which can be refined only while the remaining
  blocker bits are retained.

Ordering therefore changes the physical meaning of the landed layer even
though (8)--(9) are order-independent.

The shifted all-safe residual in the preceding paragraph no longer needs
to be declared circulant: Section 3 shows that if it reappears, it already
lands a target.

## 7. The precise next bridge

The exact missing service is one of:

1. a target-shift-covariant deletion/owner-cover identity;
2. a proof that every guard/unit failure packet is circulant or transfers
   to a marked blocker word; or
3. a same-clock phase mechanism transferring the deleted-packet triangle
   into the canonical full-owner current.

Section 5a rules out replacing these services by a higher mixed Boolean
interaction or by the complete **squared** clone spectrum. A subcube route
would still need signed singleton orientation and an absolute terminal
charged reference.

Without such a service, the theorem does not exclude a scalar row. The
ledger remains `165`, and LRC(14) remains open.

## 8. Exact companion

The dependency-free companion uses only integer bitsets and `Fraction`
arithmetic to:

- reconstruct all `16,562` rational cells and the exact `182`-cell
  Boolean deletion;
- verify `T>=C`, both diagonal zeros, all `169` entries of the circulant
  profile (22), and the exact segment/matrix hashes;
- compute the orthogonal projection and the three drift energies in
  (23);
- verify (27)--(31) for `n=1,2,3,5,7,13`;
- check the common clone-mask mass at each tested depth;
- compute every Walsh coefficient and complete-subcube energy in
  (31a)--(31e), verify all higher components and disjoint pair defects
  vanish, check the globally complemented terminal hostile, and verify
  separately the nonzero degree-two spectrum of `x -> D(H_x)`;
- verify one unequal-clone strictness control and the omitted-final-mask
  hostile;
- enumerate the sharp anchored-zero example (16d)--(16e) over all `2,197`
  cells; and
- check the exact `2016` and `26208` coefficient counts in (13)--(14).

Run

```bash
python3 04-computation/lrc14_deletion_martingale_drift_thm2370.py
python3 -O 04-computation/lrc14_deletion_martingale_drift_thm2370.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_deletion_martingale_drift_thm2370.out
```

byte-for-byte after LF normalization. Every executable check raises
explicitly under optimized Python.

Independent audits rederived the Hilbert projection inequalities and both
equality conditions; the `1/13`, `2016`, and `26208` target ledger; the
anchored-slice principal angle and its nonnegative diagonal-zero equality
case; the lawful deletion/current typing and common-right-endpoint boundary;
and the complete Boolean clone construction. A separate implementation
reconstructed all rational cells and every displayed drift value. Normal
and optimized transcripts, stored output, LF-normalized Git hashes, and
documentation routing pass. QED.
