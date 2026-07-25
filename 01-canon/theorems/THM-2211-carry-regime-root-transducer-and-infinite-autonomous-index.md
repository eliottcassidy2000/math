---
id: THM-2211
title: "Carry-regime root transducer and infinite autonomous index"
status: >
  PROVED + HOSTILE-AUDITED. A scalar guard, five unit masks, and three
  divided blockers admit an exact finite root transducer once the next
  nine floor-carry digits are supplied as edge labels. For fixed coefficient
  residues, at most 39^9 states retain every labelled sheet owner and the
  signed anti-defect law. The carry input cannot be folded into any finite
  autonomous exact state: already one unit mask has infinite zero-root
  continuation index. The separation can be completed by four legal unit
  chords so that one continuation covers the guard-safe sheets and the
  other leaves one sheet uncovered. Thus exact all-depth scalar propagation
  must retain unbounded winding/current or an unbounded carry stream. This
  is not a no-go for a finite recognizer of only the final zero-Haar
  predicate under restricted continuations.
source: codex-2026-07-24-scalar-carry-transducer
depends_on: []
related:
  - THM-2174-endpoint-phase-scale-obstruction
  - THM-2197-scalar-chord-coverage-has-a-boolean-deficiency-quotient
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier
  - THM-2204-scalar-depth-223-thirteen-lift-capacity-law
  - THM-2205-scalar-depth-113-exact-lift-capacity-exclusion
  - THM-2218-labelled-guard-hole-fourier-and-signed-lift-energy
---

# THM-2211 -- the carry-regime root transducer

Let

```text
W=(H,q_1,...,q_5,s_1,s_2,s_3)                       (1)
```

be positive integer coefficients. The actual deep blockers are `13s_j`.
Assume

```text
13 does not divide H q_1 q_2 q_3 q_4 q_5.             (U)
```

The divided blockers `s_j` may have arbitrary thirteen-adic depth. Fix a
generic phase `y`, away from every band endpoint. For each `w in W`, retain

```text
r_w=w mod 13,
a_w(y)=floor(wy) mod 13,                             (2)
```

and the three-valued regime of `theta_w={wy}`. For terminals and blockers,

```text
L=(0,1/14),      M=(1/14,13/14),      R=(13/14,1),  (3)
```

while for the guard the two cuts are `1/7` and `6/7`.

## 1. Exact finite transition with carry input

For the rooted sheets

```text
x_k=(y+k)/13,                 k in F_13,
m_w(k)=a_w(y)+r_w k mod 13,                           (4)
```

one has the exact identity

```text
{w x_k}=(m_w(k)+theta_w)/13.                          (5)
```

Indeed, write `wy=n+theta_w` and `w=13d+r_w`; the omitted terms in (5)
are integers. Consequently the next regime is determined by the current
finite state and the root digit.

For a terminal or blocker, the transition table is

```text
             current regime
m             L       M       R
0             L       L       M
1,...,11      M       M       M
12            M       R       R.                    (6)
```

For the guard it is

```text
             current regime
m             L       M       R
0             L       L       L
1             L       L       M
2,...,10      M       M       M
11            M       R       R
12            R       R       R.                    (7)
```

The next floor coordinate

```text
c_w(k)=floor(w x_k) mod 13                            (8)
```

is the carry digit not determined by (2)--(3). Supply the nine-component
vector `c(k)` as the edge label. Then (6)--(8) give the complete next state.
For fixed coefficient residues `r_w`, every coefficient has only
`13*3=39` pairs `(a_w,regime)`, so the state set has size at most

```text
39^9.                                                (9)
```

This is an exact transducer, not a probabilistic abstraction. For every
root it outputs whether the guard is safe, which unit masks own that sheet,
and whether each divided blocker is dangerous at the parent. Hence it
retains the labelled owner counts needed for endpoint deletion in THM-2197
and the full Hasse-jet incidence represented in THM-2201.

There is also an exact fibre acceptance rule. If some `s_j` is dangerous
at `y`, then

```text
{13s_j x_k}={s_j y}
```

for every `k`, so that actual blocker covers the whole root fibre. If all
three divided blockers are safe, the scalar cover holds on the fibre
exactly when every guard-safe root has at least one unit owner.

Put

```text
R=1_(C_H)-sum_(i=1)^5 1_(D_(q_i)),
C_H={x:||Hx||>1/7}.                                  (10)
```

The tables give the signed transfer invariant

```text
sum_(k=0)^12 R(x_k)=-R(y).                            (11)
```

On an accepted blocker-safe state this forces the exact words

```text
R(y)=0  ->  0^13,
R(y)=1  ->  one labelled -1 at kappa and twelve 0's. (12)
```

Thus the finite transducer contains THM-2198's one-anti-defect digit and
updates it without losing ownership.

## 2. Infinite autonomous continuation index

The carry labels in (8) cannot be removed. Set

```text
q_j=5+650*13^j,             y=1/50,
y_r=y/13^r.                                           (13)
```

All `q_j` are positive thirteen-units and

```text
q_j mod 13=5,
q_j y=1/10+13^(j+1).                                 (14)
```

Their initial rooted masks therefore coincide.

For `j<ell`, put `h=j+1` and

```text
q_ell-q_j
 =50*13^h u,             u=13^(ell-j)-1,
13 does not divide u.                                 (15)
```

Along the common zero-root path,

```text
(q_ell-q_j)y_r=13^(h-r)u.                            (16)
```

For every `r<h`, the fractional phases agree and the floor coordinates
agree modulo thirteen, so the complete labelled root masks agree. At
`r=h`, the fractional phases still agree but the floor coordinates differ
by `u mod 13`. Since both coefficient residues are five, their nonempty
one- or two-sheet masks differ by the nonzero translation

```text
-5^(-1)u=-8u mod 13.                                 (17)
```

No nonempty proper subset of `F_13` of size one or two is invariant under a
nonzero translation.

Suppose an autonomous finite state exactly propagated every positive unit
mask under the zero-root transition. Infinitely many pairs `(q_j,y)` would
give two equal initial states. Determinism would make all their future
states and outputs equal, contradicting (16)--(17) at their first split.
Thus the exact zero-root continuation quotient already has infinite
Myhill--Nerode index for one labelled terminal.

## 3. The separation changes a legal five-mask cover

This is not merely a difference between decorative sheet labels. At the
first split `h=j+1`,

```text
q_j y_h=1+epsilon,          q_ell y_h=13^(ell-j)+epsilon,
0<epsilon<1/14.                                      (18)
```

Their rooted masks are the singleton sheets `5` and `0`, respectively.
Put `D=50*13^h` and choose

```text
H_0=D/2+2.
```

This is an odd thirteen-unit, `D/7<H_0<6D/7`, and its nine guard-safe root
sheets are

```text
B={1,2,3,4,5,8,9,10,11}.                            (19)
```

Pair the eight sheets of `B\{5}` and realize the four pairs by four
positive thirteen-unit terminal coefficients, chosen mutually distinct
and distinct from `H_0` and the singleton coefficient.

The realization is exact. For a desired rooted chord `{u,v}`, choose

```text
r=-(v-u)^(-1) mod 13,
a=-ru mod 13.                                       (20)
```

Choose a residue `rho` with

```text
D/14<rho<13D/14,          rho=r mod 13,
```

and choose `A=a mod 13`. Then

```text
q=AD+rho
```

has middle regime at `y_h`, residue `r`, floor digit `a`, and rooted danger
chord exactly `{u,v}`. Adding distinct multiples of `13D` preserves all
these data while making the four coefficients positive and distinct.
Finally choose three distinct divided blockers `s_j` in `(D/7,6D/7)`;
they are safe at `y_h=1/D`, so none of the actual blockers `13s_j` covers
the missed sheet.

The singleton at sheet `5` completes a cover of `B`. Replacing it by the
singleton at sheet `0` leaves sheet `5` uncovered: if `0 in B`, it was
already paired; if `0 notin B`, it covers no required sheet. Hence any
continuation quotient compositional under legal four-mask contexts must
distinguish the two infinite-index states.

For a frozen hostile control, take

```text
h=1, D=650, y_h=1/650, H_0=327,
unit singletons: 655 -> {5},       8455 -> {0},
four chords:
  701 -> {1,2},  2001 -> {3,4},
 5251 -> {8,9},  6551 -> {10,11},
divided blockers: (326,328,329).                     (21)
```

Every displayed inequality is strict. The first five-mask row covers
exactly the guard-safe block `B`; replacing `655` by `8455` leaves sheet
`5` uncovered.

## 4. Scope and next target

Exact all-depth scalar propagation must therefore retain either:

1. infinitely many autonomous states, such as exact scale/current or full
   coefficient winding; or
2. the unbounded carry-digit input stream of (8), on top of the finite
   regime state.

This does not contradict THM-2174's one/two-state fixed-core quotient,
which remembers only the positive-/zero-Haar observable on a restricted
ray. Nor does it prove that no finite state can recognize only the final
zero-safe predicate under every restricted legal five-mask continuation.
It proves the stronger local fact actually needed for exact labelled root
propagation and compositional scalar coverage.

The remaining scalar task is therefore not another smaller static chord
quotient. THM-2207 eliminates the final depth-three word, and
THM-2213/2215/2219 close every depth-four profile. THM-2224/2226 further
reduce the ledger to `458`, THM-2227 excludes six special residues, and
THM-2229's positive-set Bellman leaves exactly `240`. Those accepted carry
words still require bounded relations, signed guard-hole correlation,
endpoint current, or owner/incidence refinements beyond the arbitrary-
coupling relaxations.

An independent hostile audit checked every regime-table boundary, the
carry/state separation, the first-split arithmetic, the nonzero sheet
translation, the four-chord continuation construction, and the scope
distinction from a bare zero-Haar recognizer. QED.
