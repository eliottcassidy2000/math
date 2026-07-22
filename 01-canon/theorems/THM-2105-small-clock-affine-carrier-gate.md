---
id: THM-2105
title: "Small-clock affine carrier gate on the guard half-fiber"
status: >
  PROVED. Let p,q be independent integer characters, let the guard be p,
  and write all transverse terminals as c_i=a_i p+n_i q. If their strict
  radius-1/14 danger bands cover the guard-safe half-fiber p.X=1/2, then for
  every 2<=m<=7 the affine congruence classes
  n_i*j+m*a_i=0 mod 2m cover Z/(2m). For each odd prime ell in {3,5,7},
  this forces either one universal carrier with a_i even and 2ell|n_i, or
  two half-carriers with ell|n_i, n_i odd, and opposite parities of a_i.
  The modulus-four case has a separately classified singleton exception.
  This is an exact necessary gate for the all-transverse fiber, not a proof
  that the complete higher-rank LRC row has an escape.
source: codex-2026-07-22-LRC-small-clock-affine-carriers
depends_on:
  - THM-2097
related:
  - THM-645
  - THM-2069
  - THM-2072
  - THM-2104
---

# THM-2105 -- small-clock affine carrier gate

Let `p,q:T^2->T` be independent integer characters and write

```text
g=p,
c_i=a_i p+n_i q,              a_i in Z, n_i in Z\{0}.  (1)
```

Assume the terminal strict-danger bands cover the guard-safe half-fiber:

```text
for every beta in T, some i satisfies
||a_i/2+n_i beta||<1/14.                                (2)
```

Condition (2) is necessary if the mixed safe cell for the guard `p` and the
listed transverse terminals is empty, because every point with `p.X=1/2` is
strictly guard-safe and independence of `p,q` realizes every `q.X=beta`.

## 1. Exact clocks through denominator fourteen

For every integer `m` with `2<=m<=7`, define

```text
S_i(m)={j in Z/(2m): n_i j+m a_i=0 mod 2m}.            (3)
```

Then (2) forces the exact finite cover

```text
Z/(2m)=union_i S_i(m).                                  (4)
```

### Proof

In (2), sample

```text
beta=j/(2m),                 j in Z/(2m).               (5)
```

The `i`-th terminal phase is

```text
a_i/2+n_i beta=(m a_i+n_i j)/(2m).                     (6)
```

It lies on the `1/(2m)` grid. Since `m<=7`, every nonzero grid residue has
circle distance at least

```text
1/(2m)>=1/14.                                           (7)
```

The inequality in (2) is strict, so (6) is dangerous exactly when its
numerator vanishes modulo `2m`. This is precisely membership in (3), and
sampling every `j` proves (4). Notice that equality in (7) at `m=7` is safe;
the strict endpoint convention is load-bearing. QED.

The result is slightly stronger than a prime test: the composite clocks
`m=4` and `m=6` are also mandatory. Formula (3), rather than a lossy count of
prime divisors, is the exact carrier.

## 2. Complete odd-prime classification

Fix an odd prime

```text
ell in {3,5,7}.                                         (8)
```

For one label `i`, put `d_i=gcd(n_i,2ell)`. The elementary linear-congruence
law says

```text
S_i(ell) is empty if d_i does not divide ell*a_i,
and otherwise S_i(ell) is one coset of size d_i.         (9)
```

Since `d_i` divides `2ell`, there are only four cases.

1. **Universal carrier:** if `d_i=2ell`, then `S_i(ell)` is the whole clock
   exactly when `a_i` is even, and is empty when `a_i` is odd.
2. **Half-carrier:** if `d_i=ell`, then `n_i=ell u_i` with `u_i` odd, and

   ```text
   S_i(ell)={j:j=a_i mod 2}.                            (10)
   ```

3. **Noncarrier:** if `d_i=1` or `2`, every nonempty `S_i(ell)` is contained
   in the two clock points `{0,ell}`.

For completeness, (10) follows after dividing the congruence by `ell`:
`u_i j=-a_i mod 2`, and signs agree modulo two. If `d_i=1`, the unique
solution is `0` for even `a_i` and `ell` for odd `a_i`; if `d_i=2`,
solvability forces `a_i` even and the solutions are `{0,ell}`.

Remove the universal-carrier case. The noncarriers touch no clock point
outside `{0,ell}`. Because `ell>=3`, both the even and odd parity halves have
points outside this pair. Cover (4) therefore requires one half-carrier of
each parity. We have proved the exact dichotomy

```text
either some i has a_i even and 2ell|n_i,
or there exist i,k with
  ell|n_i, ell|n_k, n_i,n_k odd, and a_i!=a_k mod 2.    (11)
```

This holds separately at `ell=3,5,7`.

## 3. The modulus-four exception, classified

At `ell=2`, the clock is `Z/4` and the same congruence is

```text
n_i j+2a_i=0 mod 4.                                    (12)
```

Its solution sets are exactly

```text
4|n_i, a_i even:                 {0,1,2,3};
4|n_i, a_i odd:                  empty;
n_i=2 mod 4, a_i even:           {0,2};
n_i=2 mod 4, a_i odd:            {1,3};
n_i odd, a_i even:               {0};
n_i odd, a_i odd:                {2}.                  (13)
```

Consequently a cover of the four-clock requires either a universal first
line, or an odd-parity half-carrier `{1,3}` together with a cover of
`{0,2}`. The latter cover may be one even-parity half-carrier `{0,2}`, or two
odd-`n_i` singleton labels of opposite `a_i` parity. This singleton option is
why the clean odd-prime dichotomy (11) must not be asserted at two.

## 4. Immediate family closures

Three useful contrapositives of (11) are exact.

- If, for one `ell in {3,5,7}`, no quotient coefficient is divisible by
  `ell`, then the half-fiber contains a sampled mixed-safe point.
- If `ell|n_i` occurs only on even `n_i` and no such label has both even
  `a_i` and `2ell|n_i`, the half-fiber again has a sampled escape.
- If there is no universal `2ell`-carrier, then the row must expose a named
  pair of odd quotient coefficients divisible by `ell`, with opposite lift
  parity. It is not enough merely to record that two coefficients are
  divisible by `ell`.

In particular every purely dyadic quotient family fails the `ell=3` clock
immediately. This supplies a shorter finite-clock explanation of the dyadic
subfamily already closed more generally by THM-2104's continuous valuation
clock. THM-2104 remains stronger on constant positive valuation layers where
a raw small clock may be swallowed by a universal carrier.

## 5. Code-wheel interpretation and limits

For fixed `m`, homogenize (3) as the zero-coordinate condition

```text
<(n_i,m a_i),(j,1)>=0 mod 2m.                           (14)
```

Thus (4) is an **affine slice** of the low-weight evaluation-code wheel from
THM-2069. The last coordinate is frozen to one, so projectivizing or deleting
the lift parity destroys the predicate. At prime clocks the universal and
half-carriers are precisely the rank-zero and parity-coset degeneracies of
this affine slice.

This theorem does not by itself close arbitrary coefficient rows. A single
label with even `a_i` and `n_i` divisible by every `2m`, `2<=m<=7`, swallows
all six sampled clocks. This is the same sensor-blindness phenomenon as
THM-2072's common-multiple no-go for fixed owner-clock banks. The correct next
step is adaptive: use a modulus selected from the actual quotient profile,
or combine the named carrier labels from (11) with collision/tree and
deletion persistence. Adding finitely many fixed clocks without such a
sidecar cannot be a universal proof.

The theorem also concerns the complete band list on the half-fiber. Extra
guard-proportional terminals have `n_i=0`; they may be dangerous at
`p.X=1/2` and must be retained separately. Therefore (4) must not be imported
unchanged into THM-2098's mixed vertical/transverse branch.

## 6. Assumption challenge and Tournament Analysis

The challenged assumption is that the useful small-prime datum is only the
valuation multiset. Equation (11) depends on the **incidence between an
`ell`-carrier and its lift parity**. Forgetting that incidence turns two
carriers on the same half-clock into a false certificate.

Candidate tournament vertices were terminals, clock residues, carrier
labels, parity halves, congruence obligations, and valuation walls. The
faithful finite object is the bipartite incidence graph between terminal
labels and residues of `Z/(2m)`. A useful scheduler may orient two carriers
by which covers a residue first, with numerical residue order as a tie
Hamiltonian path, but scores, cycles, SCCs, and edge flips do not determine
whether the union in (4) is complete. The quotient preserved by (11) keeps
the two parity halves and universal-carrier flag; it destroys the composite
clock incidences at `m=4,6`, which remain mandatory sidecars.

This is a proved modular gate and a new finite filter. It does not show that
the surviving carrier configurations exist geometrically, propagate through
peeling, or prove LRC(14).
