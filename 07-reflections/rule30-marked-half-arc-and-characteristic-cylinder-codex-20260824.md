# Rule 30 marked half-arc and characteristic-cylinder note

**Status:** UNNUMBERED RESEARCH NOTE.  The identities and Haar comparison
below have complete elementary proofs in this note.  The single-seed census
is `FINITE-EXACT` only.  Nothing here proves a Rule 30 prize, and this note is
not a canon promotion.

## Inheritance pass and live board

The closest proved mechanism is the terminal-current identity of
`THM-3471-rule30-motzkin-strip-circuit-and-innovation-carry-spectrum`, together
with the cyclic arc operator of
`THM-3481-rule30-cyclic-arc-norm-rank-and-marked-innovation-spectrum` and the
physical left-front coordinate of
`THM-4047-rule30-left-front-affine-monodromy-clock`.

The canonical hostile is that an unpointed phase law, full Walsh support, or
Haar measure on the innovation cube does not select the physical phase
origin.  The corrected near miss is to retain the physical endpoint before
integrating the current.  The least-used sidecar selected here is the
**nearest backward physical zero address**, rather than another phase average
or fixed-depth column.

The live board is:

1. the terminal phase line `T_k`;
2. its transition current `Q_k`;
3. the physical light-cone endpoint;
4. the nearest-zero address `z_k`;
5. the right-moving characteristic cylinder; and
6. fixed-length cyclic arc rank versus data-dependent stopping arcs.

All sums are in `F_2`.  For the single-seed Rule 30 evolution, use

```text
a_(t+1)(x)=a_t(x-1) + (a_t(x) or a_t(x+1)),
a_0(0)=1,  a_0(x)=0 for x!=0.
```

For `k>=2` and `-k<=h<=0`, retain the canonical definitions

```text
T_k(h)=a_(k+h)(h),
Q_k(h)=T_k(h+1)+T_k(h),
c_k=T_k(0).
```

## Exact half-arc rebase

Put

```text
m_k=floor(k/2)+1.
```

Then, for every integer `k>=2`,

```text
T_k(-m_k)=0,
c_k=xor_(h=-m_k)^(-1) Q_k(h).                       (1)
```

This shortens THM-3471's universal length-`k` terminal arc to the last
`floor(k/2)+1` current values.

The geometric boundary is sharp.  If `k=2r`, then

```text
T_k(-r)=1,
c_k=1+xor_(h=-r)^(-1) Q_k(h).                       (2)
```

### Proof

At phase `h=-a`, the terminal point is the physical cell

```text
T_k(-a)=a_(k-a)(-a).
```

It lies strictly outside the left light cone of the seed precisely when
`a>k-a`, equivalently `2a>k`.  Thus it is zero at `a=m_k`.  Telescoping the
transition identity from `-m_k` to `0` proves (1).  When `k=2r`, the point
`(t,x)=(r,-r)` is instead the extreme-left ray of the seed cone.  That ray is
identically one because its Rule 30 predecessor neighborhood is always
`001`.  This proves (2).

Applying THM-3481 at the new fixed length `m_k` gives a useful operator
corollary: the complete cyclic `m_k`-arc profile is lossless exactly when
`m_k` is odd, hence when

```text
k mod 4 is 0 or 1.                                  (3)
```

In particular, depths `k=0 mod 4` acquire a lossless shortened arc profile
even though the old length-`k` profile was even.  This is an operator
statement about the complete phase profile; the single marked scalar in (1)
still does not reconstruct that profile.

## Nearest zero and the characteristic-cylinder theorem

Define the physical stopping address

```text
z_k=min{a in {1,...,m_k}: T_k(-a)=0}.                (4)
```

It exists by (1), satisfies `z_k<=m_k`, and gives

```text
c_k=xor_(h=-z_k)^(-1)Q_k(h).                         (5)
```

More structure is available.  For every `1<=r<=m_k`, the following are
equivalent:

```text
z_k>r;
T_k(-1)=T_k(-2)=...=T_k(-r)=1;
T_k(-1)=1 and Q_k(-2)=...=Q_k(-r)=0;
(a_(k-r)(-r),...,a_(k-r)(r-2))=(1,0,...,0).          (6)
```

The current list in (6) is empty at `r=1`.  At the final transition,

```text
c_k=Q_k(-1)+1_[z_k>1].                              (7)
```

Thus the tail of the stopping address is exactly a marked zero-run cylinder
of the current and, one characteristic step earlier, exactly the spatial word
`1 0^(2r-2)`.  The address is not merely a convenient statistic: it is the
location of a physical cylinder occurrence.

### Proof of (6)

The first three clauses follow immediately from the definition of `z_k` and
`Q_k(h)=T_k(h+1)+T_k(h)`.  For the last clause, use the elementary
right-characteristic lemma:

> If `F` is Rule 30 and `x` is any row, then
> `(F^j x)(j)=1` for every `0<=j<r` if and only if
> `(x(0),...,x(2r-2))=(1,0,...,0)`.

The lemma is inductive.  The first characteristic value forces `x(0)=1`.
Given the forced prefix `1 0^(2j-2)`, the next characteristic update is one
if and only if its two right parents vanish.  The bit `x(2j-1)` reaches the
nearer parent on its unobstructed extreme-left ray.  If that bit is zero,
`x(2j)` reaches the farther parent on its unobstructed extreme-left ray.
Thus both parents vanish exactly when
`x(2j-1)=x(2j)=0`.  This adds two forced zeros at each step.  Conversely,
the word `1 0^(2r-2)` has a one on its extreme-right characteristic for `r`
steps.  Translate the lemma to the base point
`(t,x)=(k-r,-r)` to obtain the last clause of (6).  Equation (7) is the
transition identity at `h=-1`.

## Exact Haar law and the exact transfer obstruction

The cylinder theorem explains the geometric-looking stopping-radius data
without licensing a single-seed density claim.  Start Rule 30 from a
Bernoulli-`1/2` iid row.  For fixed `k` and `1<=r<=k`, let `Z_k` be the first
zero on the same finite backward terminal line, with `Z_k=infinity` if the
line has no zero.  Uniform Bernoulli measure is preserved by the
left-permutive Rule 30 map.  At time `k-r`, the word in (6) is therefore one
specified word among `2^(2r-1)` equally likely words.  Consequently

```text
P_Haar(Z_k>r)=2^(-(2r-1)),                            (8)
P_Haar(Z_k=1)=1/2,
P_Haar(Z_k=r)=3*2^(-(2r-1))             (2<=r<=k).   (9)
```

Appending the center is one more characteristic step.  Hence

```text
P_Haar(c_k=1 | Z_k=1)=3/4,
P_Haar(c_k=1 | Z_k=r)=1/4                 (r>=2),    (10)
P_Haar(c_k=1,Z_k=1)=3/8,
P_Haar(c_k=1,Z_k=r)=3*2^(-(2r+1))         (r>=2).    (11)
```

These identities recover Haar fairness of the center after summing over the
radius classes.  They do **not** transfer to the single-seed temporal orbit.
The two probability spaces are different:

```text
Haar statement:       random row, fixed marked cylinder;
Rule 30 prize target: one deterministic seed, time average at a physical mark.
```

For the fixed seed, (6) turns a proposed transfer into the precise missing
theorem: prove temporal equidistribution of the orbit through the widening
fixed-window cylinders `1 0^(2r-2)`.  Neither odd support, maximal ANF degree,
full Walsh support of one current readout, nor phase-Haar iid supplies joint
cylinder frequencies at the physical mark.  This is the minimal hostile to
an iid shortcut.

## Finite exact probe, hostile, and positive signal

The companion enumerates the exact single-seed universe

```text
2<=k<=2^18-1=262143.
```

It finds

```text
z : count
1 : 131120
2 :  98428
3 :  24449
4 :   6074
5 :   1570
6 :    383
7 :     90
8 :     20
9 :      8
```

The exact first witness to `z_k=9` is

```text
(k,z_k,c_k)=(79883,9,0).
```

So the cheap hostile refutes the tempting universal bound `z_k<=8`.  No
radius above nine occurs in this finite universe, but boundedness and
unboundedness are both open.

The positive signal is sharper than the raw histogram: its main ratios are
close to the exact Haar laws (9)--(11), and the center split is close to the
Haar conditional `3/4` versus `1/4`.  The stored output records exact
cross-multiplied discrepancies, so this is not being reported as equality,
convergence, or a prize-density theorem.

## Loss ledger and next target

| map | preserves | destroys / does not supply | required sidecar or debt |
|---|---|---|---|
| `T_k -> Q_k` | every transition and holonomy | additive owner bit | one physical zero endpoint |
| length-`k` arc -> length-`m_k` arc | marked center | older physical transitions | seed light-cone geometry |
| half arc -> `(z_k,Q_k[-z_k,-1])` | center and terminal one-run | current before the stopping address | the stopping address `z_k` |
| fixed-length arc -> variable nearest arc | marked evaluation | linear cyclic-operator structure after conditioning | nonlinear cylinder class |
| Haar cylinder law -> single-seed time average | the formal target frequency | physical orbit and temporal correlations | a cylinder-discrepancy theorem |

The next precise object is

```text
D_r(N)=#{k: r<=k<=N and
              (a_(k-r)(-r),...,a_(k-r)(r-2))=1 0^(2r-2)}
       -(N-r+1)/2^(2r-1).                            (12)
```

For fixed `r`, bounding `D_r(N)` is a lawful marked-cylinder discrepancy
problem.  To reach center balance one also needs the signed joint version
with the final current bit `Q_k(-1)`, because (7) shows that a marginal law
for `z_k` alone is insufficient.  The nested Haar cylinders have conditional
factor `1/4`; this gives a natural filtration/martingale comparison, but no
deterministic martingale estimate is currently proved.

There is a restrained cross-thread analogy: `z_k` plays the same *type* of
role as an endpoint-owner address in LRC or a retained filtration level in
factorial work—it must be kept before scalar integration.  No map to an LRC
loneliness predicate or factorial constant-term predicate is established, so
no theorem transfers across those subjects.

## Verification universe and reproduction

The companion uses two independent complete Rule 30 implementations through
time `512`, with `131327` direct transition-current checks, `511` full-arc
checks, `511` shortened-half-arc checks, `256` sharp even-boundary checks, and
`66047` current-zero-run / earlier-cylinder equivalence checks.  It also
exhausts all `43690` dependency words through radius eight and runs the packed
single-seed census through `262143` with an explicit search cap of `64` and
zero unresolved rows.  Ordinary and optimized Python must match the stored
output byte for byte.

```bash
python3 04-computation/rule30_marked_half_arc_zero_radius_probe_20260824.py
python3 -O 04-computation/rule30_marked_half_arc_zero_radius_probe_20260824.py
```

Artifacts:

- `04-computation/rule30_marked_half_arc_zero_radius_probe_20260824.py`
- `05-knowledge/results/rule30_marked_half_arc_zero_radius_probe_20260824.out`
