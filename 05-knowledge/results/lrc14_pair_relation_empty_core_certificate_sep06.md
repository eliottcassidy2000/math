# The actual-zero-coordinate addendum removes the support restriction

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**

Let `w=(a,b,c)` be primitive, distinct, positive, odd, nonzero modulo three,
and sorted. If there is any nonzero integer relation `v.w=0` of `l1` norm
at most twenty, then the complete THM-4414 network projections satisfy

```text
min_i E_i(w)<=6/77,
```

with equality only at `(1,5,11)`. If its complete live support is empty,
multi-ray, or on a primitive ray of norm other than four, then

```text
E_i(w)<6/77 for every i.
```

This combines the independently audited
[full-support theorem](lrc14_empty_core_certificate_sep06.md) with the
support-two result proved here. A relation can be divided by its content,
which decreases its norm, so primitivity need not be assumed in the
headline. A relation with only one nonzero coordinate is impossible because
all speeds are positive.

Consequently every hypothetical failure of the universal network target
must satisfy

```text
min_{0!=v in Z^3, v.w=0} ||v||_1 >=22.                 (1)
```

The lower bound is twenty-two, not twenty-one, because all three speeds are
odd and every integer relation has even coefficient sum modulo two. This
is a necessary condition on a failure of a sufficient local certificate;
it does not assert a counterexample or supply entry, synchronization, or
LRC(14).

## Inheritance and exact interval retained

The raw-carrier, error-cube, affine integer defect, and owner-residue arguments
are inherited from Sections 1-2 of the full-support note and its cited
THM-4386/4398 sources. An actual zero coefficient is distinguished from a
nonzero coefficient divisible by three. The new operation keeps the former
and uses the now-free error coordinate to get an explicit rectangle width.
The full-support source and its independent audit remain frozen.

The inherited hostile `(1,5,7)` shows why the norm-four live ray must remain
an exception to the every-projection assertion. The map here is the same
complete affine carrier interval; it preserves every owner-live raw point,
and discards individual roof values only after the common cap `3/(7c)`.

## 1. The eleven possible support-two coefficient patterns

After a coordinate permutation and overall sign, a primitive support-two
relation is `v=(p,-q,0)`, with `1<=p<q` and `gcd(p,q)=1`.
The corresponding speeds have the form

```text
w=(qs,ps,t), s,t>0.
```

The coefficients `p,q` are odd and units modulo three, because the speeds
are odd and units modulo three. Equal coefficients would be `p=q=1` and
would force equal speeds, excluded by hypothesis. Under `p+q<=20`, the
complete list is

```text
(1,5), (1,7), (1,11), (1,13), (1,17), (1,19),
(5,7), (5,11), (5,13), (7,11), (7,13).                 (2)
```

The zero coordinate makes this the one-zero-mod-three case of the owner
dichotomy. The complete allowed defect list is

```text
D={delta in Z: |delta|<(3/14)(p+q), 3 does not divide delta}.
```

For each such defect exactly one scalar multiplier class modulo three is
owner-live on the affine carrier line `C_delta+Z(p,-q,0)`.

## 2. Exact rectangle width and count envelope

Write `r=3/14`. The error slice has the description

```text
p e_1-q e_2=delta, |e_i|<=r,
e_2 in [max(-r,(-pr-delta)/q), min(r,(pr-delta)/q)],
e_3 in [-r,r].
```

Let `J_delta` be the length of the displayed `e_2` interval. Equivalently,

```text
J_delta=min(2pr,r(p+q)-|delta|)/q.                     (3)
```

On the carrier line, a scalar coordinate is

```text
(w cross e)_1/p = s e_3-(t/p)e_2.
```

The two error coordinates vary independently over a rectangle, so its
exact scalar width is

```text
T_delta=2rs+(t/p)J_delta.
```

With `c=max(qs,ps,t)`, we have `s<=c/q` and `t<=c`. Hence

```text
T_delta <=c(2r/q+J_delta/p).
```

The open-interval one-residue count from the full-support theorem gives

```text
N <A_(p,q)c+B_(p,q),
A_(p,q)=(1/3) sum_(delta in D) (2r/q+J_delta/p),
B_(p,q)=|D|.                                         (4)
```

Every projection is at most `3N/(7c)`. Thus `A<2/11` again gives all three
strict certificates for `c>=B/(2/11-A)`. This derivation does not divide by
the zero coordinate and needs no polytope extremum assumption.

## 3. Complete finite head

All eleven rational slopes in (4) are strictly below `2/11`:

| Pair | `A` | `B` | Cutoff `B/(2/11-A)` |
|---|---:|---:|---:|
| `(1,5)` | `2/21` | 2 | `231/10` |
| `(1,7)` | `4/49` | 2 | `539/27` |
| `(1,11)` | `8/77` | 4 | `154/3` |
| `(1,13)` | `8/91` | 4 | `2002/47` |
| `(1,17)` | `8/119` | 4 | `2618/75` |
| `(1,19)` | `34/399` | 6 | `13167/212` |
| `(5,7)` | `6/49` | 4 | `539/8` |
| `(5,11)` | `2/21` | 4 | `231/5` |
| `(5,13)` | `116/1365` | 4 | `30030/727` |
| `(7,11)` | `50/539` | 4 | `539/12` |
| `(7,13)` | `68/637` | 6 | `21021/263` |

The largest cutoff is below eighty. The exact verifier generates every
odd ternary-unit integer pair `(s,t)` below each cutoff with `gcd(s,t)=1`,
forms the sorted distinct triple `(ps,qs,t)`, and retains all relation
presentations. This is complete because every pair with reduced ratio
`p:q` has that form.

There are `209` pair-pattern/head memberships and `182` unique triples:
`21` empty, `16` norm-four, `113` other one-ray, and `32` multi-ray supports.
All exact selected projections satisfy the target, with the sole equality
`(1,5,11)`. All three are strict whenever the support is outside the
norm-four class. The infinite tail is strict by (4), proving the headline
after adjoining the independently audited full-support result.

The combined coefficient count is `73+11=84`: seventy-three full-support
patterns, including the inherited norm-four exception, and eleven
support-two patterns. This is a coefficient presentation atlas, not a
partition of speed triples.

## Reproduction and exact dependency scope

```bash
python3 -B 04-computation/lrc14_pair_relation_empty_core_certificate_sep06.py
python3 -B -O 04-computation/lrc14_pair_relation_empty_core_certificate_sep06.py
```

The [producer](../../04-computation/lrc14_pair_relation_empty_core_certificate_sep06.py)
imports the frozen full-support producer's raw-carrier and projection routines
as declared dependencies. Every finite head is also independently compared
with the literal six-sheet contact engine from the incoming one-ray verifier.
It checks every affine defect and strict layer count. The
[output](lrc14_pair_relation_empty_core_certificate_sep06.out) records all
rational interval lengths, constants, and proof counts after `3,144`
explicit checks.

The [independent referee](../../04-computation/lrc14_pair_slice_audit_empty_core_three_ray_sep06.py)
recomputes the widths by both polygon clipping and a closed formula,
enumerates the entire `2,910`-row eligible universe through height `79`,
selects the heads by direct gcd-pair tests, and rebuilds carriers by a full
integer box. It reproduces all `209` memberships, `182` unique heads,
projections and physical masses. Producer and referee normal/optimized runs
each have byte-identical output.

```bash
python3 -B 04-computation/lrc14_pair_slice_audit_empty_core_three_ray_sep06.py
python3 -B -O 04-computation/lrc14_pair_slice_audit_empty_core_three_ray_sep06.py
```

Frozen raw-byte SHA256:

```text
source 78cef0ea9036493a7e85c79394ff9e9eb1547822b79d99c75819859b2073289e
output d75aa683996ecdd5d96ec4c833eb1309f74a0f17180c6c302f663cdb5bef1927
independent source 86ae1d3094b545737b7dea498dcd9969744d8dcfb3b90729eba79b457d136da1
independent output a74639810c314759d13512a36e15658a9a5485200314ac73d44f0c6d17246b75
shared semantic head digest e57cd661d136ac45511b795207fcf3ab2dc1df5a7f16f28487cb56407eba097d
```

These finite heads are part of the proof. No statement bounds the shortest
relation by twenty for arbitrary eligible triples. Longer-relation control,
the universal network ceiling, entry, synchronization, and LRC(14) remain
open.
