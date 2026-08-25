# LRC14 `d=2,3,4` pack-safe/spoiled-phase intersection audit

**Status:** `FINITE-EXACT` on the three finite exception-profile banks declared
by the THM-4041, THM-4032, and THM-4030 audits.  This is **not** an exhaustive
test of the physical THM-3818 survivor branch: none of the low control rows in
these banks passes the inherited THM-3818 filters.

## Exact typed sets and endpoint conventions

Work on `R/Z`, represented linearly by `[0,1]` with `0~1`, and write
`||z||=min({z},1-{z})`.  For a divided pack `H`, divisor `d`, and exception
tuple `E`, define

```
G(H)       = {y : ||h y|| >= 1/14 for every h in H},
D(w,j)     = {y : ||w(y+j)/d|| < 1/14},
Sigma_d(E) = intersection over j=0,...,d-1 of union over w in E of D(w,j),
R_d(H,E)   = G(H) \ Sigma_d(E).
```

Thus `G(H)` is closed, every danger interval `D(w,j)` and `Sigma_d(E)` are
open, and `R_d(H,E)` is closed.  Every wall and open cell is classified with
exact `Fraction` arithmetic; interval records have the form
`(left,left_closed,right,right_closed)`.

The lift implication is exact.  If `y` belongs to `R_d(H,E)`, some label `j`
survives.  For `x=(y+j)/d`,

```
||d h x|| = ||h(y+j)|| = ||h y|| >= 1/14     (h in H),
||w x|| >= 1/14                               (w in E).
```

Hence the full row `dH union E` has a lonely phase.  This is the predicate the
intersection audit tests, not merely a scalar certificate count.

## Universes

The named rows are retained as positive/hostile controls; they are not called
THM-3818 survivors.

| `d` | pack `H` | exact declared profile bank | profiles | certificate-positive |
|---:|---|---|---:|---:|
| 2 | `{1,...,10,12}` | unordered distinct positive odd pairs `alpha<beta<=79` | 780 | 759 |
| 3 | `{1,...,10}` | unordered triples from `{1<=w<=23 : 3 does not divide w}` | 560 | 312 |
| 4 | `{1,...,10}` | `(2r,a,b)`, `r in {1,3,5,7,9,11}`, odd `1<=a<b<=19` | 270 | 67 |

For `d=2`, positivity is independently checked against
`alpha/gcd(alpha,beta)+beta/gcd(alpha,beta)>7`.  For `d=3,4`, every positivity
decision is checked both by direct open wall cells and by the canonical affine
defect certificate functions from THM-4032/4030.

## Result

No certificate-positive profile contains the pack-safe set in its spoiled set.
In fact the exact minimum surviving measure is positive in every bank:

| `d` | positive profiles | `G(H) subset Sigma` | exact minimum `measure(R)` | unique minimizer |
|---:|---:|---:|---:|---|
| 2 | 759 | 0 | `145/1764` | `(11,13)` |
| 3 | 312 | 0 | `629/8820` | `(2,11,13)` |
| 4 | 67 | 0 | `839/8820` | `(2,11,13)` |

A cheaper exact certificate uses only two closed-safe endpoints per bank.  In
each row below, every positive profile leaves at least one of the two points
outside `Sigma`:

| `d` | two-point subset of `G(H)` | profiles spoiling the first point | profiles spoiling both |
|---:|---|---:|---:|
| 2 | `{1/14,5/56}` | 30 | 0 |
| 3 | `{1/14,13/140}` | 4: `(1,13,14)`, `(1,14,16)`, `(2,13,14)`, `(2,14,16)` | 0 |
| 4 | `{15/98,1/14}` | 1: `(14,1,13)` | 0 |

For the required canonical hostile controls the spoiled selector phase and a
full-row survivor are both reproduced:

| `d` | exceptions | spoiled safe `y` | surviving safe `y` | full phase `x` | full-row clearance |
|---:|---|---:|---:|---:|---:|
| 2 | `(1,11)` | `1/11` | `229/280` | `229/560` | `5/56` |
| 3 | `(1,5,11)` | `2/11` | `3/14` | `1/14` | `1/14` |
| 4 | `(2,9,11)` | `1/11` | `9/11` | `21/22` | `1/11` |

These controls demonstrate why certificate positivity alone cannot close the
row: the certificate spoils at least one pack-safe phase, but not every
pack-safe phase.

## THM-3818 inheritance and physical typing obstruction

The audit also maps every positive low profile to the natural control
decomposition `s=1`, `t=d`:

- `d=2`: pair `(p,q)=(1,12)`, body `2*{2,...,10} union E`;
- `d=3`: pair `(p,q)=(1,10)`, body `3*{2,...,9} union E`;
- `d=4`: pair `(p,q)=(1,4)`, body `4*{2,3,5,...,10} union E`.

All rows pass the local distinctness/primitivity/gcd typing checked by the
script.  None passes the full THM-3818 filter:

| `d` | physical passes | pair-atlas failures | crossing-height failures | largest crossing height in bank |
|---:|---:|---:|---:|---:|
| 2 | 0 | 759 | 759 | 79 |
| 3 | 0 | 0 | 312 | 30 |
| 4 | 0 | 0 | 67 | 19 |

Here the inherited crossing threshold is
`Q=91^6=567869252041`, and each checked crossing height is
`max(tp,s*u_i)/gcd(tp,s*u_i)`.  The `d=2` pair also fails the atlas because
`1+12=13`; the `d=3,4` pairs have the permitted sums `11` and `5` but fail
height overwhelmingly.

Therefore there is no finite physical profile bank in THM-3818, THM-4004,
THM-4030, THM-4032, or THM-4041 to which one can honestly promote the finite
result above.  The height condition is a lower bound, not a coefficient cutoff.

The minimum extra input needed for the intended exhaustive physical test is
one of the following equivalent data packages:

1. an exact finite producer list of tuples
   `(s,t,p,q,u_1,...,u_11,d)` after the THM-3818 atlas, gcd, distinctness, and
   crossing-height filters, together with the projection of each tuple to its
   divided profile `(d,H,E)`; or
2. a proved finite reduction bounding those projected `(H,E)` profiles and a
   generator shown to enumerate every physical survivor.

Without that input, a scan of the three theorem audit banks is a control-bank
result only.  With it, this interval engine is already parameterized by
`(H,d,E)` and can perform the physical intersection without further geometry.

## Candidate statement

**Declared-bank two-endpoint escape (`FINITE-EXACT`).**  In each of the three
finite universes in the table above, every affine-certificate-positive exception
profile `E` satisfies `R_d(H,E) != empty`; indeed one of the displayed two
rational endpoints lies in `R_d(H,E)`, and the exact bankwise measure minima
are the three fractions displayed above.  Consequently every declared row
`dH union E` has a `1/14`-safe phase.

This is suitable as a finite computation theorem or hostile-control lemma.  It
does not yet reduce or close the physical LRC14 branch.

## Replay

From the worktree root:

```powershell
& 'C:\Users\Eliott\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' `
  '.scratch\lrc14_d234_pack_safe_intersection_20260824\audit.py'
& 'C:\Users\Eliott\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' `
  -O '.scratch\lrc14_d234_pack_safe_intersection_20260824\audit.py'
```

Both modes return `RESULT=PASS` with identical semantic digests.  The final
summary digest is
`fc375ff65bc4e5c36b4e60a906291fa55088a0e149b130867da1a9341e4f78b8`.
The LF-normalized SHA-256 of `audit.py` is
`3932b4b25838dcf3b239ca0285bbb8f34362a44a31371a1ef28ce7be5ee7c09b`.
