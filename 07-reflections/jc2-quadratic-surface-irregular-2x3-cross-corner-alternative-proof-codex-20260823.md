# THM-3786 cross-corner alternative proof and hostile replay

**Status: INDEPENDENT SHORTER ALTERNATIVE PROOF + VERIFIED-EXACT HOSTILE
REPLAY.**  The canonical source theorem is already proved and independently
hostile-audited as
[THM-3786](../01-canon/theorems/THM-3786-quadratic-etale-surface-irregular-two-by-three-support-no-go.md).
This note gives a second proof of its irregular-support conclusion.  It uses
only the proved homogeneous anatomy and lower-support no-gos of
[THM-3783](../01-canon/theorems/THM-3783-quadratic-tower-etale-surface-maximal-polynomial-observable.md),
then replaces the canonical three-arm hub peel by a forced-cross-corner
argument.  It is not a new theorem or an arbitrary Darboux-pair no-go: the
common-step `d=e=f` lane survives.

## Exact rederivation

Write the exact active supports, after subtracting output constants, as

```text
A={a,a+d},                 B={b,b+e,b+e+f},
d,e,f>0.
```

For homogeneous components `x^u p(s),x^v q(s)`, THM-3783 gives

```text
J=2x^(u+v+3)(u p q'-v p' q).                         (A1)
```

If `u,v` are nonzero and have opposite signs, a vanishing bracket would give
in `k(s)`

```text
u q'/q = v p'/p.                                     (A2)
```

At any irreducible divisor the two sides say `u ord(q)=v ord(p)`.  The
multiplicities are nonnegative and the weights have opposite signs.  Hence
both profiles would have to be constant; a negative-weight legal profile is
divisible by `s^2-1`, contradiction.  If `u=0!=v`, vanishing in (A1) forces
the weight-zero profile to be constant.  Thus for genuinely active profiles
a commuting endpoint pair has either two weights of the same strict sign or
two zero weights.

A homogeneous scalar bracket has `u+v=-3`.  Two negative profiles make (A1)
vanish at `s=1`; a zero/negative pair does too.  In the remaining orientation
`u>0>v`, nonzero value at `s=1` forces the negative profile to have a simple
zero there.  Since it is divisible by
`(s^2-1)^ceil(-v/4)`, this leaves only `(u,v)=(1,-4)` and its swap.  There the
profile equation is

```text
2(pq'+4p'q)=constant,
```

whose leading coefficient contains `deg(q)+4deg(p)>0`, in positive degree.
Therefore a singleton contribution can never be the scalar bracket.

The lowest and highest pair-sum buckets are singletons, so their brackets
vanish.  The low endpoint weights cannot both be nonnegative, because no pair
sum could be `-3`.  Combining this with the exact commutation alternatives
above forces

```text
a<0,                         b<0.                    (A3)
```

The high endpoint weights cannot both be nonpositive.  If they were, all six
component brackets would vanish at `s=1`: negative profiles contain
`s^2-1`, negative/zero brackets still retain the negative profile as a
factor, and zero/zero brackets vanish identically.  The high endpoint
commutation alternatives then force

```text
a+d>0,                       b+e+f>0.                (A4)
```

Consequently the two cross-corner brackets `(A_low,B_high)` and
`(A_high,B_low)` are nonzero by (A2).  Neither can occupy a singleton bucket,
even when that bucket has pair sum `-3`, because the homogeneous scalar case
was just excluded.

Relative to `a+b`, the complete labeled pair-sum list is

```text
0:(0,0), e:(0,1), e+f:(0,2),
d:(1,0), d+e:(1,1), d+e+f:(1,2).                    (A5)
```

Positivity of the gaps makes this list exhaustive.  The corner at `e+f`
collides only if `d=f` (with `d+e`) or `d=e+f` (with `d`).  The corner at `d`
collides only if `d=e` (with `e`) or `d=e+f` (with `e+f`).  Requiring both
collisions therefore gives exactly

```text
d=e=f                            or d=e+f.           (A6)
```

There are no triple collisions.  On `d=e+f`, the increments `e` and `d+e`
are singletons, so

```text
{A_low,B_middle}=0,              {A_high,B_middle}=0. (A7)
```

The middle weight is positive, negative, or zero.  In the first case the
first bracket joins opposite signs; in the second the second does; in the
third either equation forces its nonconstant zero-weight profile to be
constant.  All cases contradict (A7).  Transposing (A5) gives the literal
output-swap `3 by 2` proof with no change except the sign of the total bracket.

Thus every exact irregular `2 by 3` or output-swapped `3 by 2` support cell is
empty.  The sole combinatorial survivor is `d=e=f`.

## Hostile boundaries and controls

- If subtracting an output constant deletes a weight-zero component, the
  support shrinks to `1 by 3` or at most `2 by 2`; THM-3783 already excludes
  those shapes.  Otherwise every retained weight-zero profile is
  nonconstant, exactly as required in (A2) and (A7).
- A scalar at an endpoint or any other singleton is a homogeneous scalar
  bracket and is excluded before the collision census; it is not silently
  treated as a zero bucket.
- The common-step controls
  `A={2-t,2}, B={-5,t-5,2t-5}` for every `t>=3` cross both endpoint signs and
  put pair sum `-3` in a double bucket.  They confirm that collision logic
  alone must not claim the common-step lane is empty.  They are not Darboux
  constructions.
- The independent replay
  `04-computation/jc2_quadratic_surface_irregular_2x3_cross_corner_audit_20260823.py`
  executes 2,397,670 optimization-safe checks, including 373,248 gap triples,
  the transposed/output-swap census, constant and nonconstant zero-weight
  controls, 25,527 endpoint-sign controls, 1,223,168 hub-sign controls, and 70
  unaligned common-step positive controls.  Every substantive test uses an
  explicit runtime gate, so optimization cannot erase the audit.

## Reproduction and evidence

From the repository root, the following PowerShell block runs the normal and
optimized replays and compares both with the LF frozen transcript:

```powershell
$script = '04-computation/jc2_quadratic_surface_irregular_2x3_cross_corner_audit_20260823.py'
$frozen = '05-knowledge/results/jc2_quadratic_surface_irregular_2x3_cross_corner_audit_20260823.out'
$normal = (& python -B $script) -join "`n"
$optimized = (& python -B -O $script) -join "`n"
$stored = (Get-Content -Raw -LiteralPath $frozen).Replace("`r`n", "`n").TrimEnd("`n")
if (($normal -cne $optimized) -or ($normal -cne $stored)) { throw 'THM-3786 cross-corner replay mismatch' }
$normal
```

The command returns `PASS`, and all three line-normalized outputs agree.
Raw LF working-tree byte SHA-256 hashes are

```text
3a25d7788793714b4365471c615de7056138cc11b02d4b225abe897712a75245
  jc2_quadratic_surface_irregular_2x3_cross_corner_audit_20260823.py
17e9e2b4f3a3c60bf2a281bfb46c8adf674f58b405f2cb5fea4bc16f87adadd2
  jc2_quadratic_surface_irregular_2x3_cross_corner_audit_20260823.out
```

## Exact scope and common-step boundary

For exact active supports

```text
supp(F)={a,a+d},       supp(G)={b,b+e,b+e+f},       d,e,f>0,
```

on the THM-3783 quadratic etale surface, this alternative proof establishes
the THM-3786 implication `{F,G}=1 => d=e=f`; output swap gives the literal
`3 by 2` form.  If target translation deletes a weight-zero component, the
support falls into THM-3783's homogeneous-output or complete `2 by 2` no-go.

The common-step condition is necessary here, not sufficient.  The 70
positive controls show that the collision shape genuinely survives the
argument; they are not Darboux pairs.  Apart from THM-3783's named aligned
orientation, common-step `2 by 3`, larger supports, arbitrary Darboux pairs,
`JC(2)`, and `DC(2)` remain open.

