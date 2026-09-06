# Normalized reciprocal jets: exact tree updates and target-dependent precision

**PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.** The starting
arbitrary-multiplicity inverse theorem is the already **PROVED** incoming
[THM-4443](../../01-canon/theorems/THM-4443-arbitrary-jet-precision-and-dyadic-unit-boundary.md),
with proof in
[the higher-jet unit-boundary report](hermite-higher-jet-unit-boundary-overnight-hexagon-sep05.md).
The new contribution is a precise integral tree-shell update, a target-dependent
residue budget, an explicit all-depth compression for unequal multiplicities,
and an unbounded terminal-pair cancellation family with a necessary target
firewall. It does not re-claim the incoming dyadic three-node discovery.

## Inheritance and scope

The closest proved mechanism is THM-4443's attainment of every inverse
denominator in an actual top-row cardinal coefficient. The canonical hostile
is its dyadic isometric three-jet pair. The corrected near miss is recorded
in MISTAKE-547: the complete uniform two-jet terminal cancellation budget
cannot be transferred to a changed observer. The least-used sidecar is
the *degree of the reciprocal Taylor coefficient*, which assigns its own
cost under common dilation. Relevant exact-statement searches found no
earlier normalized shell-jet or unbounded pair-budget result in the local
canon/results; no external priority claim is made.

The concept board is: metric denominator powers; reciprocal power moments;
integral truncated products; shell depth attenuation; all-depth precision
envelopes; and local-versus-global cancellation. These are one connected
calculation, not independent experimental lanes.

Fix a prime p, distinct integer nodes x_i, positive multiplicities m_i,
and M=sum m_i. Observe Hasse orders `0,...,m_i-1` at each node on degree<M
polynomials. The target is the largest p-Smith exponent L, equivalently
the sharp extra coefficient-recovery precision. Ordinary derivatives,
partial local jets, changing coefficient modules, and the full Smith
partition are not substituted for this target.

## 1. The normalized integral state

For n>=2 set

```
h_ij=v_p(x_i-x_j),
D_i=sum_(j!=i) m_j h_ij,
f_i=max_(j!=i) h_ij,
B_i(T)=product_(j!=i)(1+p^(f_i)T/(x_i-x_j))^(-m_j)
       modulo T^(m_i)
      =sum_(l=0)^(m_i-1) b_(i,l) T^l.                 (1)
```

Every normalized reciprocal `p^f_i/(x_i-x_j)` is p-adically integral.
The negative-binomial coefficients are integers, so all b_(i,l) lie in
Z_p and b_(i,0)=1. For a singleton use f_i=D_i=0 and B_i=1.

Let Q_i be the monic integer product over all other nodes, as in THM-4443.
Then

```
Q_i(x_i+p^f_i T)^(-1)=Q_i(x_i)^(-1) B_i(T),
a_(i,l)=Q_i(x_i)^(-1) p^(-l f_i) b_(i,l).
```

Consequently the exact loss is

```
L=max_(i,0<=l<m_i) [D_i+l f_i-v_p(b_(i,l))].          (2)
```

Here v_p(0)=infinity. This is an equality because the inherited inverse
theorem proves attainment, not merely a denominator upper bound. The
valuation tree supplies D_i and f_i. The normalized product jet supplies
the missing cancellations, in precisely the observed Taylor orders.

Uniform affine changes `x_i -> u x_i+a`, with u a p-adic unit and integer
images, send b_(i,l) to u^(-l)b_(i,l); valuations and the target are
unchanged. Common multiplication by p^e leaves B_i itself unchanged.

## 2. Fixed-anchor tree recursion and power-moment form

Fix the leaf coordinate x_i. Group other leaves into shells of common
depth h, or partition each shell into its sibling blocks in the distance
tree. For such a block C define the unit-coordinate jet

```
J_(i,C)(T)=product_(j in C)(1+p^h T/(x_i-x_j))^(-m_j)
           modulo T^(m_i).
```

Disjoint blocks merge by multiplication followed by truncation. Their
contribution to B_i is `J_(i,C)(p^(f_i-h) T)`. Thus a coefficient of
degree r from depth h acquires the factor `p^(r(f_i-h))`. Multiplying all
sibling blocks along the root-to-leaf path gives exactly (1); grouping,
associativity and permutation of blocks change no coefficient.

For a three-jet state and an incoming block, write the jets as
`1+b1 T+b2 T^2` and `1+s T+t T^2`. The exact update is

```
b1' = b1+s,
b2' = b2+b1*s+t.                                    (3)
```

The cross term is necessary. Adjoining one node of multiplicity m with
normalized reciprocal alpha gives

```
b1' = b1-m*alpha,
b2' = b2-m*alpha*b1+binom(m+1,2)*alpha^2.
```

If the new node increases f_i, first replace the old jet by
`B_i(p^(f_new-f_old) T)`. Independently update the metric constant D_i by
`m*v_p(x_i-x_new)`. These formulas include unequal multiplicities.

Equivalently retain the normalized weighted power moments

```
P_(i,s)=sum_(j!=i) m_j (p^f_i/(x_i-x_j))^s.
```

Differentiating the product in (1) yields the exact recurrence

```
r*b_(i,r)=sum_(s=1)^r (-1)^s P_(i,s) b_(i,r-s).      (4)
```

In particular b1=-P1 and b2=(P1^2+P2)/2. Thus moments through the observed
order suffice, but their residues must retain cancellation. Valuations
alone fail already at terminal pair `{0,2}`: with uniform multiplicity
three and outsider 1 versus 3, the x=0 normalized moments have
`(v2(P1),v2(P2))=(0,0)` in both cases, while v2(b2) is 4 versus 2.

Equation (4) must not be executed by modular division at unchanged
precision when p divides r. Induction shows that `r! b_r` is an integer
polynomial in P1,...,Pr: multiply the recurrence by (r-1)! and use
`(r-1)!/(r-s)!` as the integer coefficient of each earlier factorial
polynomial. Hence a sufficient common precision for computing a jet
through order m_i-1 modulo p^K is to know these moments modulo
`p^(K+v_p((m_i-1)!))`. This sufficient padding is not claimed optimal.
The supplied modular compiler uses integral negative-binomial convolution,
so it incurs no division by p in the first place.

The fixed anchor is a real condition. Truncated moments at one centre
cannot be freely recentered with the same order budget: translation mixes
higher reciprocal moments, with depth attenuation controlling the required
extra terms. The present theorem gives fixed-leaf shell updates, not a
centre-free constant-size message per tree vertex or a new complexity bound.

## 3. A finite residue budget for the current target

The metric alone gives the lower bound `D=max_i D_i` through order zero.
For each coefficient put

```
K_(i,l)=max(0,D_i+l f_i-D).                           (5)
```

If K_(i,l)=0, that coefficient cannot improve D. Otherwise it is enough
to know b_(i,l) modulo p^K_(i,l). A nonzero residue determines its exact
valuation below K_(i,l); a zero residue means its candidate in (2) is
at most D. Taking the maximum of D and the candidates with nonzero
residue therefore recovers L exactly. One may replace D by any other
already certified lower bound on L to reduce the required precision.

For a full leaf jet, using
`K_i=max_l K_(i,l)=max(0,D_i+(m_i-1)f_i-D)` is sufficient. In the shell
compiler, every shell coefficient of degree r with
`r(f_i-h)>=K_i` vanishes modulo p^K_i and may be omitted. Shallower shells
therefore require fewer degrees and fewer p-adic digits. This is the
specific sense in which the sidecar is adapted to the target rather than
retaining complete unit coordinates.

This pruning is for the current depth only. It must not delete orders
permanently when the goal includes future dilations; Section 4 gives the
separate exact all-depth object.

## 4. The all-depth envelope has at most max m_i slopes

For n>=2, under `x_i -> p^e x_i`, e a nonnegative integer, the normalized
b_(i,l) are unchanged, D_i increases by `(M-m_i)e`, and f_i increases by e.
For a singleton the convention instead keeps D_i=f_i=0, B_i=1; all losses
are zero and the profile is `{0:0}`. Set

```
sigma_d=max_(i,l: M-m_i+l=d)
          [D_i+l f_i-v_p(b_(i,l))],                 (6)
```

omitting zero coefficients; an empty maximum is minus infinity. Then

```
L(e)=max_d [sigma_d+d e].                            (7)
```

All slopes lie in the integer interval `[M-max_i m_i,M-1]`, so there are
at most max_i m_i of them, regardless of the number of nodes. Equal
slopes are merged by their maximal intercept. Further lines may be
discarded if they never change the maximum at any integer e>=0.

This gives an exact order-aware finite envelope, including heterogeneous
observers. It is not a claim that one fixed residue modulus computes the
whole future envelope, nor that the packet is globally information-minimal
among arbitrary encodings of node sets. Truncation to the observed jet
order is exact; information above that order never enters the target.

The positive controls recover the incoming dyadic uniform three-jet twins:

```
(0,1,2): sigma={6:3,7:4,8:1},
(1,0,3): sigma={6:3,7:4,8:3}.
```

The order-zero line is dominated, and the remaining lines give exactly
THM-4443's Gamma-dependent law. The heterogeneous `(2,2,1)` twins give
`{3:2,4:1}` and `{3:2,4:0}`. The earlier ternary complete two-jet twins
both give `{6:18,7:22}`, retaining their equal largest loss despite
different intermediate Smith factors.

A decisive pruning hostile is uniform `(0,1,2)`: order two cannot beat
the metric baseline at e=0, but at e=4 its line gives 33, while the
lower-order lines give only 32. The current-depth residue budget must
not be mistaken for a permanent deletion from the all-depth profile.

## 5. Terminal pairs have no node-count-independent three-jet budget

For a terminal dyadic pair `{0,2}`, uniform three jets, write the
normalized reciprocals at an endpoint as `alpha_j=2/(x-y_j)`. The
order-two coefficient is

```
b2=(3/2) [3 (sum alpha_j)^2+sum alpha_j^2].           (8)
```

The incoming three-node law implies `min_pair v2(b2)<=5` in this
normalization. This statement is false with arbitrary outside nodes.
Already at

```
X={0,2,-7,-5,1,7,9}
```

the pair is terminal, all five outsiders have depth zero from it, and
both endpoint coefficients equal `1452032/33075`, with valuation 11.
Exact reciprocal inversion gives the whole observer loss 31; the local
cancellation statistic is not the global largest exponent.

There is an unbounded family. For any integer h>=2 put

```
r=2^h-1, K=h+2,
X_h={0,2} union {1+2^K*j: 0<=j<r},
all node multiplicities=3.                          (9)
```

These are distinct integer nodes. The pair has depth one, and every
outsider is odd, so the pair remains a complete terminal cluster.
At endpoint zero the internal reciprocal is -1 and all external
reciprocals are -2 modulo 2^(K+1). At endpoint two they are +1 and +2,
respectively. Thus each endpoint has the same model value

```
(3/2)[3(1+2r)^2+(1+4r)] =6(r+1)(3r+1).             (10)
```

To justify the division in this congruence, each external reciprocal
changes by a multiple of 2^(K+1); the total sum is odd, so its square
changes by a multiple of 2^(K+2). Each external square also changes
by a multiple of 2^(K+2), since the reciprocal is even. Dividing the
bracket by two therefore preserves congruence modulo 2^(K+1)=2^(h+3).
For h>=2, `v2(r+1)=h` and `v2(3r+1)=1`. Hence

```
v2(b2 at 0)=v2(b2 at 2)=h+2.                        (11)
```

This proves unbounded simultaneous terminal-pair cancellation as node
count grows, without a Hensel argument or extrapolation from a head.
It does not claim unbounded cancellation at one fixed node count or
one fixed full metric tree.

The target firewall is part of the result. At the displayed scale the
pair has D_i=3 and f_i=1. Its normalized first coefficient is odd, so
its actual largest local candidate is exactly
`max(3,4,5-(h+2))=4`; its uncorrected order-two cost is only 5. Every
outsider sees r-1 other outsiders at depth at least K, giving
`D_outside>=3(r-1)K>5`. Thus the badly cancelling pair is already
excluded by the metric baseline at this scale. The family refutes a
universal *local cancellation budget*, not a new metric-only theorem
for the global largest exponent. No all-depth pair-domination assertion
is made from this displayed-scale argument.

## 6. Limits, exact tests and connection contract

Even within geometric branch products, the first-order jet does not
determine the next coefficient. At p=2 and anchor zero, depth-zero blocks with
nodes `{-11,-231}` and `{-15,-35}`, each multiplicity three, have jets

```
1-(2/7)T+(947/17787)T^2,
1-(2/7)T+(179/3675)T^2.
```

Their first-order states agree, but their second-order states differ.
This is a formal-state obstruction, not a claim that the two complete
observers necessarily have different largest Smith exponents. The
moment-valuation hostile in Section 2 separately shows why phases must
be combined before their valuations are taken.

The source is the multiplicity-labelled node set; the map retains its
valuation tree and fixed-anchor reciprocal shell jets. The tree preserves
metric denominator powers but loses the unit cancellations. Truncated
integral products restore precisely the observed Taylor coefficients;
the target-dependent residue budget discards only coefficients proved
unable to change the current loss. The all-depth map instead retains
degree-labelled intercepts before taking an upper envelope. The lost
information includes full node positions and the full Smith partition.
The canonical twins, the permanent-pruning hostile, and the terminal-pair
family are three distinct decisive tests of those losses.

The source `overnight6_20260906_jets_sidecar.py` imports no earlier producer
mathematics. It compares shell products with independently expanded Q_i
and reciprocal inversion, exact power-moment recurrence, modular shell
convolution, node insertion, and four common dilation scales. Nine small
observers also have literal integer Hasse Smith checks; the 21-column
seven-node hostile uses independent reciprocal inversion instead of a
large full-matrix computation. The family congruences are checked through
h=12, with additional literal rational checks through h=5. The universal
proof is (9)-(11), not these finite controls.

```
python -B 04-computation/overnight6_20260906_jets_sidecar.py
python -B -O 04-computation/overnight6_20260906_jets_sidecar.py
```

Normal and optimized outputs are byte-identical LF: **596 explicit gates
PASS**. The [independent audit](overnight6_20260906_jets_independent_audit.md)
accepts the full proof and reconstructs 34 literal Hasse inverses using only
standard-library rational Gaussian elimination, with **866 explicit gates**
and matching normal/optimized outputs. It separately verifies the seven-node
global loss31, family h=2 global loss33, pair-local loss4, and all profile
constants. Its singleton and prime-scope wording clarifications are
incorporated above. Root filed the audited artifacts in the sixth checkpoint.

```
source 5e6cf9e03659d0b74cb69e059c211f4026590c942ec95a4d803c75c14c2e3f4b
output d8c60d7b394e0b35bef50ccaaee8f0bcecac808e3d5f468a07afc36ec101ab87
```

Still open: smaller centre-independent tree messages, a classification
of necessary unit residues in higher-node observers, and whether useful
global cancellation bounds survive after metric-baseline pruning. The
proved incoming inverse formula and the exact shell algorithm make these
precise questions without promising a new metric-only classification.

**Filing:** root integrated this audited report after `f5f0f7f75`;
portable reproduction paths are shown above. The exact producer and
transcript bytes remain pinned by the sixth manifest.
