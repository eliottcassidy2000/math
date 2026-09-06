# THM-4448 independent audit: general shore attachment

**Verdict: PASS.**  The attachment theorem is correctly typed for arbitrary
three distinct positive tails, including multiples of three and tails that
coincide with speeds in the scaled body.  The exact event atlases, strict
endpoint behavior, filtered and unfiltered pair-gap constants, and both
parametric hostile families also pass an implementation independent of the
primary artifact.  This proves no additional row exclusion and does not prove
`LRC(14)`.

## Proof audit

Let the ten-body be `B union hU`, with `|U|=r` and `|B|=10-r`, where
`1<=r<=9`.  The physical base `3B union T` has **at most** `13-r` distinct
moving speeds: a tail divisible by three may coincide with a member of `3B`,
but that only decreases the count.  Thus cited `LRCUpTo13` supplies clearance
at least `1/(14-r)`.  Relative to the target `1/14`, the available margin is
exactly

```text
1/(14-r)-1/14 = r/[14(14-r)].
```

Put `y_0=3x_0`.  Over any proper closed quotient arc there is one unique local
inverse branch `sigma` through `x_0`.  Moving quotient distance `s` costs at
most `b s` for a body speed `3b`, while `sigma` moves distance `s/3` and costs
at most `t s/3` for a tail.  This gives

```text
rho_r=min(r/[14(14-r)max(B)], 3r/[14(14-r)max(T)]).
```

No inverse of `t mod 3` occurs here.  Intrinsically,
`K_t(y)={x:3x=y, ||tx||<1/14}`; when `3|t`, this set can be empty or the whole
three-point fibre.  The selected branch remains safe throughout the protected
arc, which is exactly what is needed to show the arc lies outside `F_T`.

The packet danger set satisfies `D_(hU)={y:hy in D_U}`.  Each open component
of length `ell` pulls back to components of length `ell/h`.  If the protected
closed arc were wholly dangerous, it would lie in one such open component.
The hypothesis `2h rho>=delta(U)` rules this out even at equality: a closed arc
of the same length cannot be contained in an open arc without two positive
endpoint gaps.  A proper open subset of the circle may have a component of
length exactly one if its complement is a singleton; using `ell<=1` handles
that boundary directly, while lower-dimensional LRC supplies a positive safe
arc and hence `ell<1` in the intended range.

## Exact independent checks

The companion script imports no THM-4448 primary code.  It independently:

- reconstructs the complete failure cells and masses `6/77` for `(1,5,11)`
  and `11/140` for `(2,11,20)` from literal three-sheet masks;
- retains the equality-safe singleton components of `(1,13)` at `1/14` and
  `13/14`, obtaining strict gap `1/7` and the deliberately incorrect
  closed-merge control `15/91`;
- verifies full-failure-overlap controls for both tails;
- enumerates all `19,314` coprime `p<q`, `p+q<=356`, with unique maximum
  `15/98` at `(1,14)`;
- separately applies THM-3818's primitive-sum filter (every prime divisor is
  `2 mod 3`, exponent at most two), obtaining `5,855` ratios and unique maximum
  `29/196` at `(1,28)`;
- checks all nine LRC margin identities and 45 exact pullback scale controls;
- checks the overlapping nonunit-tail control `T=(3,6,9)` on its local inverse
  branch; and
- checks exact representative values from both hostile progressions.

The number `46,837` has the distinct THM-3818 scope of positive-scale divisor
seam triples `(p,q,s)` with `s>=2`; it is not a ratio count.  Together with
the `5,855` scale-one triples it gives `52,692`.

For `N=53+2310k` and `N=121+210k`, the step is divisible by
`2*3*5*7`, the base is coprime to that product, and the relevant `N y_*`
residue is frozen.  Hence for every `k>=0`, `N` is coprime to `1,...,8`, the
reduced distinguished pair is `(1,4)`, and its cross height is `N`.  The base
inequalities

```text
6/(7*53)<3/154,    6/(7*121)<1/140
```

strictly place the marked components inside the claimed failure cells and
only strengthen with `k`.  Both families satisfy the cofinal cone, so they are
hostile only to prescribed-component selection.

## Reproduction

From the repository root:

```powershell
python -B 04-computation/lrc14_general_shore_attachment_thm4448_independent.py
python -B -O 04-computation/lrc14_general_shore_attachment_thm4448_independent.py
```

Both commands exit zero and are line-identical to
`05-knowledge/results/lrc14_general_shore_attachment_thm4448_independent.out`.
Hashes are SHA-256 of raw LF repository bytes:

```text
04-computation/lrc14_general_shore_attachment_thm4448_independent.py
  989ed36ca551e0e31d29071579e13e4a5ba696d2ea3776f8a61d4683f374fb26
05-knowledge/results/lrc14_general_shore_attachment_thm4448_independent.out
  b02e4827a56d453df05836bc8d6ed13eacfb70253d0ffa742d5d96bc0aa5f916
```
