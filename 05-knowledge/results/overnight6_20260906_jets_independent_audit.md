# Independent audit of normalized reciprocal jets and terminal-pair cancellation

**Status: INDEPENDENT PROOF AUDIT + EXACT CONTROLS PASS.** This audits
[the primary report](overnight6_20260906_jets_sidecar.md)
and its frozen producer. It is an independent proof audit, not a separate
claim of discovery. No repository or Git files were changed.

The accepted scope is arbitrary positive multiplicities at distinct integer
nodes, complete Hasse jets on degree below their total multiplicity, and a
specified prime. The consumer is the largest Smith exponent at that prime.
The full Smith partition, ordinary derivative observations and partial jets
are different targets.

## 1. Inheritance and exact normalization

The starting equality is already proved in
[THM-4443, arbitrary-jet precision and dyadic unit boundary](../../01-canon/theorems/THM-4443-arbitrary-jet-precision-and-dyadic-unit-boundary.md).
Its cardinal inverse gives both the upper denominator bound and attainment:
the top coefficient of the column for jet order `r` at node `i` is exactly
`a_(i,m_i-r-1)`. Thus no necessity/sufficiency reversal is hidden in the
primary use of the largest inverse denominator.

[THM-4439, complete two-jet terminal-cluster precision](../../01-canon/theorems/THM-4439-all-node-twojet-metric-precision-by-terminal-clusters.md)
is not extended to higher or unequal multiplicities. The current correction
is [MISTAKE-547](../../01-canon/MISTAKES.md): changing the observer changes
the cancellation budget. The primary correctly preserves this boundary.

For `q_i=Q_i(x_i)`, `D_i=v_p(q_i)`, and
`f_i=max_(j!=i)v_p(x_i-x_j)`, direct substitution gives

```text
Q_i(x_i+p^f_i T)^-1=q_i^-1 B_i(T),
a_(i,l)=q_i^-1 p^(-l f_i)b_(i,l).
```

All normalized reciprocals are integral in `Z_p`; the negative-binomial
coefficients are ordinary integers. Hence all `b_(i,l)` are integral and
`b_(i,0)=1`. The exact attained loss is therefore

```text
L=max_(i,l<m_i) (D_i+l f_i-v_p(b_(i,l))).
```

Zero coefficients contribute minus infinity. Singleton data give `B=1`,
loss zero and no exceptional inverse denominator.

## 2. Recursion, precision and the dilation consumer

Partitioning the product into fixed-leaf shells is exact. A shell at depth
`h` has unit reciprocals `p^h/(x_i-x_j)`; replacing its variable by
`p^(f_i-h)T` multiplies degree `r` by `p^(r(f_i-h))`. Truncating every
product at degree `m_i-1` commutes with later multiplication because no
factor has negative powers. The three-jet cross term `b1*s` is required.
Adjoining a deeper neighbor first rescales the existing jet, then multiplies
the new negative-binomial factor. The report keeps the anchor fixed.

The logarithmic derivative independently gives

```text
r b_r=sum_(s=1)^r (-1)^s P_s b_(r-s),
P_s=sum_j m_j alpha_j^s.
```

Multiplying by `(r-1)!` proves inductively that `r!b_r` is an integer
polynomial in `P_1,...,P_r`. Consequently the stated common padding
`v_p((m_i-1)!)` is sufficient when moments are used to compute the whole
jet modulo `p^K`. It is not claimed optimal. The primary explicitly avoids
division by a nonunit at unchanged precision; its actual modular compiler
uses the integral binomial factors.

Let `D=max D_i`. For a coefficient with uncorrected cost `c=D_i+l f_i`,
the residue modulo `p^max(0,c-D)` is sufficient: a nonzero residue gives its
valuation exactly, while a zero residue proves its candidate is at most
the already attained baseline. A coefficient can be pruned only for this
target and scale. The shell attenuation test is safe since all remaining
factors have integral coefficients.

For `n>=2`, dilation by `p^e` increases `D_i` by `(M-m_i)e` and `f_i`
by `e`, leaving the normalized jet fixed. Thus the individual slope is
`M-m_i+l`. These slopes lie in `[M-max m_i,M-1]`, so merging equal slopes
leaves at most `max m_i` lines. Taking their maximum is the exact all-depth
loss. The singleton has the separate profile `{0:0}`. There is no claim
that the residues sufficient at depth zero suffice for the whole envelope.

I verified the decisive pruning hostile independently from literal matrix
inverses. Uniform `(0,1,2)` has profile `{6:3,7:4,8:1}`. The slope-eight
branch cannot beat the metric baseline at depth zero, but at `e=4` it gives
33, while the lower-order lines give only 32. This proves why the current
target and the future-scale target need different retained information.

## 3. The unbounded local family and its global boundary

For uniform three jets and terminal pair `{0,2}`, put
`alpha_y=2/(x-y)`. Direct expansion gives

```text
b2=(3/2)(3(sum alpha_y)^2+sum alpha_y^2).
```

For `h>=2`, set `r=2^h-1`, `K=h+2` and put the `r` outsiders at
`1+2^K j`, `0<=j<r`. They are distinct odd integers, so `{0,2}` remains
a complete terminal cluster of depth one.

At zero the reciprocals are `-1` and `-2` modulo `2^(K+1)`; at two they
are `+1` and `+2`. Their common model coefficient is
`6(r+1)(3r+1)`. Each outsider reciprocal changes by a multiple of
`2^(K+1)`. Squaring the odd total sum introduces an extra factor two;
each outsider square has at least that extra factor as well. After division
by two the error remains divisible by `2^(K+1)=2^(h+3)`.

Since `v_2(r+1)=h` and `v_2(3r+1)=1`, the model has valuation `h+2`,
strictly below the modulus exponent. Both actual coefficients therefore
have valuation exactly `h+2`. This is an all-height proof of unbounded
simultaneous local cancellation as the node count grows. It asserts
neither fixed-node-count unboundedness nor a fixed full metric tree.

The pair has `D=3`, `f=1`, and its normalized first coefficient is odd.
Its actual local maximum is exactly

```text
max(3,4,5-(h+2))=4.
```

Every outsider has order-zero metric candidate at least
`3(r-1)(h+2)>5`. Thus this terminal pair never controls the global loss
at the displayed scale. The primary correctly limits this conclusion to
that scale. The family refutes a uniform local cancellation budget, not
a new global metric-only theorem.

The separate seven-node witness also passes: at
`(0,2,-7,-5,1,7,9)`, both pair coefficients are
`1452032/33075`, of valuation 11. The pair-local loss is 4, whereas a
literal 21-by-21 Hasse inverse has global loss 31. For the unbounded
family's `h=2` member, a literal 15-by-15 inverse has global loss 33.

## 4. Independent exact verification and scope clarifications

The [independent companion](../../04-computation/overnight6_20260906_jets_independent_audit.py)
uses only the Python standard library. It constructs each literal Hasse
matrix and performs Fraction Gauss-Jordan inversion. It imports no primary
producer and uses neither SymPy nor a Smith-form algorithm. Every inverse
entry contributes to the directly computed largest denominator; extracting
the top inverse row separately checks the claimed cardinal attainment.

There are 34 full inverse controls: eight observers at four scales each,
plus the seven-node witness and the `h=2` family member. These include
uniform and heterogeneous dyadic twins, signed nodes, multiplicity four,
primes three and five, and a singleton. The same literal data check all
shell products, moment recurrences, adaptive budgets and grouped slopes.
The family congruence is also checked by exact rational sums for `h=2..7`.

Further hostiles distinguish the information losses:

- At pair `{0,2}`, outsiders 1 versus 3 give the same valuations `(0,0)`
  of `P1,P2`, but `b2=48` versus `44/3`, of valuations 4 versus 2.
- The constant and first-order branch jets can agree while the second
  coefficients differ; this is an exact-state obstruction, not by itself
  a claim about unequal global precision.
- Knowing a numerator modulo two does not determine its half modulo two:
  numerators zero and two are the elementary nonunit-division hostile.

Normal and optimized audit outputs are identical. Both pass **866 explicit
gates**. Reproduce with

```text
python -B 04-computation/overnight6_20260906_jets_independent_audit.py
python -B -O 04-computation/overnight6_20260906_jets_independent_audit.py
```

Audit hashes:

```text
source d452112ff59443d6b7475435c22e9d458b6ba1188ab0d8d2f95786523329b287
output 83ec0c72c24e28f1b962c9327a06e88f382b1e94ab10956f1fc27b5ee285af16
semantic 56a1bba0e33516656c58c01fdec9212b8cbb92045ea24edb9eb8d7efa7feb406
```

The audited frozen primary producer and output are

```text
source 5e6cf9e03659d0b74cb69e059c211f4026590c942ec95a4d803c75c14c2e3f4b
output d8c60d7b394e0b35bef50ccaaee8f0bcecac808e3d5f468a07afc36ec101ab87
```

Two wording clarifications were requested and incorporated during review: restrict the
`f_i -> f_i+e` statement to `n>=2`, retaining singleton profile `{0:0}`;
and name `p=2` for the formal depth-zero branch witness with nodes
`{-11,-231}` and `{-15,-35}`. Neither changes the producer or any numerical
result. I reread both corrections; no mathematical correction remains outstanding.

**Filing:** root integrated this audited report after `f5f0f7f75`;
portable reproduction paths are shown above. The exact producer and
transcript bytes remain pinned by the sixth manifest.
