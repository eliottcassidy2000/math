---
id: THM-2120
title: "Thirteen-division phase rigidity empties the generic rank-eight terminal-blocker branch"
status: >
  PROVED. In THM-2116's rank-eight branch, assume the unique 13-content
  blocker is the terminal c_*=13u, with g and u independent, g nonzero modulo
  13, and seven residual terminals transverse to g modulo 13. An almost-
  everywhere terminal cover of the guard complement is impossible. After the
  surjective change Y=13X, two residual singleton windows outside the blocker
  band contradict THM-2116's thirteen-point incidence ledger. Avoiding every
  such pair forces an exact character dichotomy: either all seven residuals
  lie on the u-line, or one lies off it and the other six equal +/-u. The
  second case has danger-union measure at most 3/7<5/7; the first is excluded
  by settled LRC(9) on the one-dimensional terminal fiber. Thus the two
  extremal colored-toothpick patterns of THM-2116 cannot persist almost
  everywhere. This closes only the unique independent terminal-blocker branch;
  guard blockers, multiple 13-blockers, dependent g,u, nontransverse residues,
  and ranks nine through twelve remain.
source: codex-2026-07-22-LRC-thirteen-division-phase-rigidity
depends_on:
  - THM-2116
related:
  - THM-2098
  - THM-2114
  - THM-2115
external: Settled Lonely Runner Conjecture for at most eight integer speeds (LRC(9) in the total-runner convention).
---

# THM-2120 -- thirteen-division phase rigidity

Retain the notation and hypotheses of THM-2116. Thus `K` is a connected
rank-two torus with character lattice `Gamma`,

```text
c_*=13u,                    g,u Q-independent,         (1)
```

`g mod 13` is nonzero, and `c_1,...,c_7` are the residual terminal
characters. Put

```text
epsilon=1/14.                                         (2)
```

Assume, toward a contradiction, that the terminal-danger bands cover

```text
C={X in K:||g.X||>1/7}                                (3)
```

up to a null set `E`. The strict-band convention is used throughout;
replacing a strict inequality by its closed endpoint changes only a null set.

## 1. Division by thirteen turns orbit phases into ordinary characters

Apply the surjective isogeny

```text
Y=13X.                                                 (4)
```

For every residual character and for the blocker,

```text
13(c_i.X)=c_i.Y,                 c_*.X=u.Y             (5)
```

in `R/Z`. Every `Y` has a guard-safe thirteenth root. Indeed, start with any
`X_0` satisfying `13X_0=Y` and add `z in K[13]`. Since `g mod 13` is nonzero,
the values of `g.z` run through the full grid

```text
{j/13:j in F_13}.                                     (6)
```

One point of every translate of this grid lies within `1/26` of `1/2`, and
therefore has circle norm at least `6/13>1/7`.

Fix distinct residual labels `i,j`. If the open set

```text
O_ij={Y: ||c_i.Y||<epsilon,
            ||c_j.Y||<epsilon,
            ||u.Y||>epsilon}                          (7)
```

were nonempty, it would have positive Haar measure. Take a guard-safe lift of
one of its points. Strictness and a local inverse branch of (4) give a
nonempty open set

```text
P_ij={X:13X in O_ij and ||g.X||>1/7}.                 (8)
```

Let `v` be THM-2116's guard-kernel direction and discard

```text
E_v=union_(r in F_13) (E-rv/13).                      (9)
```

This finite union is null, so choose `X in P_ij\E_v`. Equation (5) says the
blocker is safe and residuals `i,j` both satisfy the strict singleton phase
criterion. Every point of the orbit `X+rv/13` lies outside `E`, so its seven
residual danger sets cover all thirteen labels. But two singleton sets leave
total incidence at most

```text
1+1+5*2=12,                                           (10)
```

which cannot cover thirteen points. This contradiction proves the global
phase containment

```text
{Y:||c_i.Y||<epsilon and ||c_j.Y||<epsilon}
       subset {Y:||u.Y||<=epsilon}                    (11)
```

for every residual pair. This is stronger than the local diamond obstruction:
the full torus, rather than one principal lift, remembers finite kernels.

## 2. Exact two-character rigidity

We use the following character lemma.

> **Lemma.** Let `a,b,u` be nonzero characters of a connected two-torus and
> let `0<epsilon<1/3`. If
> ```text
> {||a.Y||<epsilon and ||b.Y||<epsilon}
>       subset {||u.Y||<=epsilon},                    (12)
> ```
> then:
> 1. if `a,b` are Q-independent, `u` is one of `+/-a,+/-b`;
> 2. if `a,b` are Q-dependent, then `a,b,u` are Q-proportional.

**Independent case.** The torus map

```text
phi=(a,b):K -> T^2                                    (13)
```

is surjective with finite kernel `L`. On `L`, (12) says `||u.Y||<=epsilon`.
If `u(L)` were a nontrivial finite subgroup of the circle, it would contain a
point of distance at least `1/3` from zero. Hence `u` kills `L`.

The character exact sequence for (13), equivalently Smith normal form, now
gives integers `m,n` with

```text
u=ma+nb.                                               (14)
```

Surjectivity lets `(a.Y,b.Y)` range over a small real square about `(0,0)`.
If `S=|m|+|n|>1`, choose `0<t<epsilon` with

```text
epsilon<tS<1/2                                        (15)
```

and take the two coordinates to be `t sign(m),t sign(n)`, using zero when a
coefficient vanishes. Then (12) fails. Thus `|m|+|n|<=1`; since `u` is
nonzero, `u=+/-a` or `u=+/-b`.

**Dependent case.** Write `a=r alpha`, `b=s alpha` for the primitive generator
`alpha` of their saturated character line and nonzero integers `r,s`. If the
restriction of `u` to the connected circle `ker(alpha)` were nonzero, it would
be surjective. A point with `a=b=0` and `u=1/2` would contradict (12). Hence
`u` kills `ker(alpha)`, so it is an integer multiple of `alpha`. This proves
the lemma.

The integer conclusion in (14) is the load-bearing improvement over a purely
Lie-algebraic drift. Locally, writing `u=Aa+Bb` gives only the absolute-diamond
condition `|A|+|B|<=1`. Testing every finite sheet promotes the coefficients
to integers and collapses that diamond to its four vertices. In the dependent
case, common proportionality is necessary but not sufficient; the surviving
one-dimensional interval inclusion must still be retained.

## 3. Seven residuals collapse to two one-dimensional ledgers

Apply the lemma to (11) for every residual pair. If every `c_i` is rationally
proportional to `u`, call this the **all-line case**. Otherwise choose `c_1`
off the `u`-line. For every `j>=2`, the characters `c_1,c_j` cannot be
dependent, since the dependent clause would put `u` on their common line. The
independent clause then says

```text
c_j=+u or c_j=-u.                                     (16)
```

Thus there are only two possibilities:

```text
all seven c_i lie in Q u;                             (17a)
exactly one c_i lies off Q u and the other six are +/-u. (17b)
```

The colored steps cannot remain freely arranged: global phase avoidance
collapses them to one character line.

## 4. Both ledgers are impossible

In case (17b), the blocker `13u`, the six characters `+/-u`, and the one
exceptional character contribute at most three distinct danger bands. Each
has Haar measure `1/7`, so their union has measure at most

```text
3/7 < 5/7=measure(C).                                 (18)
```

They cannot cover (3), even almost everywhere.

In case (17a), let `alpha` be the primitive generator of the saturated line
`Gamma intersect Q u`. The eight terminal characters `c_*,c_1,...,c_7` are
nonzero integer multiples

```text
13k alpha,m_1 alpha,...,m_7 alpha.                    (19)
```

After discarding repetitions and signs, the settled LRC for at most eight
integer speeds supplies `t in T` with

```text
||d t||>=1/9>1/14              for every speed d.     (20)
```

Because `g` and `alpha` are Q-independent, `(g,alpha):K->T^2` is surjective.
Choose `X` with `g.X=1/2` and `alpha.X=t`. Equation (20) makes all eight
terminals strictly safe and the guard is strictly safe. The inequalities have
margin, so a neighborhood of `X` is uncovered, contradicting the almost-
everywhere cover.

Both cases are impossible. Therefore THM-2116's unique independent terminal
`13`-blocker branch is empty. QED.

## 5. Frontier effect and Tournament Analysis

The exact scope is important. This theorem does not handle a blocker in the
guard, two or more terminal blockers, `g,u` dependent, residuals that fail
mod-13 transversality, or terminal ranks above eight. Those are now the live
prime-invoice branches.

The challenged assumption is that the seven colored toothpicks should be
classified as a finite graph on `F_13`. That graph is only one fiber. The map
`Y=13X` assembles every fiber simultaneously, and the common-kernel test in
the character lemma detects information no single chord graph retains.

Candidate tournament vertices were residual labels, character rays, finite
kernel sheets, singleton windows, and coverage obligations. Orienting residual
labels by projective slope gives a tie Hamiltonian path, but score histograms,
cycles, SCCs, edge flips, and path counts all forget whether `u` descends
integrally through a two-character map. The faithful carrier is

```text
(character lattice; finite kernels of (c_i,c_j); u-line;
  thirteen-division lift; null-set sidecar).          (21)
```

It preserves the LRC predicate through (4)--(11), while the finite orbit alone
destroys the global kernel constraint that forces (17).
