---
id: THM-3291
title: "Two-block TV homogenization rigidity"
status: >
  PROVED + VERIFIED-EXACT; INDEPENDENT IMMUTABLE AUDIT PENDING.
  The total variation distance between two two-trial binomials has the exact
  closed form `TV(Bin(2,x),Bin(2,y)) = |x-y|(1+|x+y-1|)`.  For two
  single-observation blocks the homogenization bound
  `delta_N <= delta_1+delta_2-delta_1 delta_2` follows from exactly two
  ingredients: a box constraint `|x+y-1| <= 1-(delta_1+delta_2)/2` and AM-GM.
  Its equality locus is classified and is RIGID: apart from the trivial locus
  `delta_1=delta_2=0`, equality holds precisely when the two block gaps have
  the same sign, are EQUAL, and every block is pinned to a common face of the
  unit square (`min(p_i,q_i)=0` for both `i`, or `max(p_i,q_i)=1` for both).
  There is no equality with exactly one vanishing block gap.  The equivalent
  reading is that the affinity `1-delta` is supermultiplicative across blocks.
audit: >
  The exact rational companion verifies the closed form on a complete 25x25
  grid, sweeps the COMPLETE 17^4 = 83,521-point rational grid for violations
  (none), matches the nontrivial equality set against the predicted rigid face
  in both directions (62 = 62), confirms there is no mixed equality, checks
  each of the three proof ingredients separately on its own complete grid, and
  carries a hostile control showing that doubling the multiplicative
  correction breaks the bound at 11,520 points.  Normal and `-O` replay are
  byte-identical.  A concurrent session (boxeph, 2026-08-03) reached the same
  closed form and equality classification by a different route (sign-resolved
  symbolic certificates plus a ragged-denominator hostile scan); its 351-point
  equality set is exactly this file's 62 nontrivial points together with the
  degenerate faces, so the two classifications cross-confirm.  Independent
  immutable audit is pending.
source: death-star-gvc3-counterexample-2026-08-03
depends_on: []
related:
  - THM-3290-archimedes-flatness-and-the-gmc3-gvc3-counterexample-family
  - 05-knowledge/results/tv-fusion-homogenization-lemma-boxeph.md
external:
  - "The general block inequality `delta_N <= delta_I + delta_J - delta_I
    delta_J` for the homogenization map is CITED, not proved here: A.
    Kontorovich, *TV homogenization inequalities*, arXiv:2601.04079v3
    (v1 2026-01-07, v3 2026-02-25), math.PR.  That paper's point is that
    homogenization -- sending each Bernoulli parameter to the cumulative mean
    -- is NOT a data-processing map, so this bound is not the usual coupling
    bound.  What is proved below is the two-single-observation-block case and,
    beyond the cited statement, the exact equality classification."
script: 04-computation/two_block_tv_homogenization_rigidity_thm3291.py
output: 05-knowledge/results/two_block_tv_homogenization_rigidity_thm3291.out
script_sha256: 6a7ff7e89a714907cd08abb2b68cab6e2e7e412b5b899c1a475636d8bdd7321a
output_sha256: f19438cb639eb288ea75acbb18025dedf0fbc806d493fdc3595fc3f1085d2ac7
hash_basis: LF-normalized bytes
---

# THM-3291 -- two-block TV homogenization rigidity

**PROVED + VERIFIED-EXACT; INDEPENDENT IMMUTABLE AUDIT PENDING.**

## 1. The closed form

For `x,y in [0,1]`,

```text
TV(Bin(2,x), Bin(2,y)) = |x-y| * (1 + |x+y-1|).                    (1)
```

*Proof.*  The three signed pmf differences are `(y-x)(2-x-y)`,
`2(x-y)(1-x-y)` and `(x-y)(x+y)`.  Since `x,y in [0,1]` the factors `2-x-y`
and `x+y` are nonnegative, so half the sum of absolute values is

```text
(1/2)|x-y| [ (2-x-y) + 2|1-x-y| + (x+y) ] = |x-y|(1+|1-x-y|).      (2)
```

QED

## 2. The two-block bound

Let block `i` carry the single Bernoulli parameters `p_i` and `q_i`, put
`delta_i = |p_i - q_i|`, and let

```text
x=(p_1+p_2)/2,   y=(q_1+q_2)/2,   delta_N = TV(Bin(2,x),Bin(2,y)).  (3)
```

**Theorem.** `delta_N <= delta_1 + delta_2 - delta_1 delta_2`.

*Proof.*  Write `e_i = p_i - q_i`, so `x-y=(e_1+e_2)/2` and by `(1)`

```text
delta_N = (|e_1+e_2|/2)(1+s),      s := |x+y-1| in [0,1].          (4)
```

**Opposite signs.**  Then `|e_1+e_2| <= max(delta_1,delta_2)`, and since
`(1+s)/2 <= 1`,

```text
delta_N <= max(delta_1,delta_2)
        <= delta_1+delta_2-delta_1 delta_2,                        (5)
```

the last step because the difference is `min(delta)*(1-max(delta)) >= 0`.

**Same signs.**  Then `|e_1+e_2| = delta_1+delta_2`.  Because
`p_i,q_i in [0,1]` and `|p_i-q_i|=delta_i`, the sum `p_i+q_i` is confined to
`[delta_i, 2-delta_i]`; averaging the two blocks gives the **box constraint**

```text
s = |x+y-1| <= 1 - (delta_1+delta_2)/2.                            (6)
```

Substituting `(6)` into `(4)` and then applying **AM-GM**,

```text
delta_N <= ((delta_1+delta_2)/2)(2-(delta_1+delta_2)/2)
        = (delta_1+delta_2) - (delta_1+delta_2)^2/4
        <= (delta_1+delta_2) - delta_1 delta_2.                    (7)
```

QED

So the bound is exactly a box constraint plus AM-GM; nothing about binomials
beyond `(1)` is used.

## 3. The equality locus is rigid

**Theorem.**  Assume `delta_1 delta_2 != 0`.  Equality holds in the two-block
bound **iff** all three of the following hold:

```text
(R1)  e_1 and e_2 have the same sign;
(R2)  delta_1 = delta_2;
(R3)  min(p_i,q_i)=0 for both i,  or  max(p_i,q_i)=1 for both i.   (8)
```

Moreover there is **no** equality when exactly one of `delta_1,delta_2` is
zero.

*Proof.*  In the opposite-sign branch the first inequality of `(5)` is strict
unless one gap vanishes, which is excluded.  In the same-sign branch equality
forces equality in both `(6)` and the AM-GM step of `(7)`.  AM-GM is tight iff
`delta_1=delta_2`, giving `(R2)`.  The box constraint `(6)` is tight iff each
`p_i+q_i` attains the same endpoint of `[delta_i,2-delta_i]`, i.e. iff for
every `i` the smaller of `p_i,q_i` is `0`, or for every `i` the larger is `1`;
that is `(R3)`.  For the last claim, if `delta_2=0` the bound reads
`delta_N <= delta_1` and `(4)` gives `delta_N=(delta_1/2)(1+s)`, so equality
needs `s=1`, i.e. `x+y in {0,2}`, i.e. all four parameters equal `0` or all
equal `1` -- which forces `delta_1=0`.  QED

The face `(8)` has codimension two inside the same-sign chamber (one equation
`delta_1=delta_2`, one face condition), so equality is genuinely exceptional.
The complete `17^4` rational grid finds exactly 62 nontrivial equality points
and exactly 62 predicted face points, and they coincide.

## 4. Reading, and the honest contrast with THM-3290

The bound is equivalent to supermultiplicativity of the **affinity**:

```text
1 - delta_N >= (1 - delta_1)(1 - delta_2).                         (9)
```

That form is what makes it iterate across a block decomposition.

It is worth recording next to THM-3290 precisely because the two results point
in opposite directions, and the reason is a type distinction rather than a
technique:

- `TV` is a **positive** functional.  Homogenizing -- averaging the parameters
  over a block -- can only contract it, and `(9)` says the contraction is
  exactly supermultiplicative.  Information is degraded, never annihilated.
- The Gaussian moment `L` of THM-3290 is a **signed** functional.  Averaging
  over the sphere can annihilate every moment of `P` while leaving `x^2 P^m`
  visible; there the average manufactures a vanishing that no single slice has.

This is a contrast, **not a bridge**: there is no map here from one setting to
the other, no shared preserved predicate, and no transfer is claimed.  What
transfers is only the warning that "orbit average" behaves categorically
differently on signed and positive functionals.

## 5. Scope

Proved: `(1)`, the two-block bound for two single-observation blocks, and the
equality classification `(8)`.  **Not proved:** the general block inequality
`delta_N <= delta_I + delta_J - delta_I delta_J` for arbitrary blocks, which is
the cited external theorem; any statement about Poisson binomials of larger
blocks; any LRC, GMC, or AMM 12592 consequence.

Run

```text
python 04-computation/two_block_tv_homogenization_rigidity_thm3291.py
python -O 04-computation/two_block_tv_homogenization_rigidity_thm3291.py
```

and compare LF-normalized bytes with the declared output.  Exact rational
arithmetic on complete grids; no floating point, random sampling, imported
executable, or assertion-sensitive test.

**QED.**
