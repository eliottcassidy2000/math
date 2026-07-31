# The half-cycle debt cylinder has every moving 13-character but no left residue

**Status: FINITE-EXACT positive partial result + exact typing obstruction.**
This note executes the paired blocker--graft test proposed in equation `(18)`
of
[`lrc14-guard-source-debt-cone-is-target-null-20260728.md`](lrc14-guard-source-debt-cone-is-target-null-20260728.md)
on the two proved THM-2698 half-cycle packets.  Every admissible unresolved
`k_a` pivot gives a nonconstant paired-shift table, and in fact all twelve
primitive `13`-characters are nonzero.  The calculation is on one positive
physical ancestry interval at a time.  It is not a completed THM-2334 target
current: THM-2698 retains no atomic left relation index `u`, so `eta.u` is
undefined rather than known to vanish.

## 1. Inherited carrier and the unresolved pivot

Use the canonical typed row

```text
(H,q1,...,q5,c1,c2,c3)
 =(1,14,27,40,53,66,13,13^3,2*13^5).                    (1)
```

For the first target dipole put

```text
a=c2=2197,                 c=c3=742586.                  (2)
```

THM-2445 assigns the unique maximal `7`-adic ordinary role `14` to the
deepest-target graft `k_b`.  Guard `H=1` and source `j=c1=13` are the two
target-neutral debt roles in THM-2461/THM-2701.  The retained THM-2698 labels
do not select the other graft pivot, leaving exactly

```text
k_a in {27,40,53,66}.                                    (3)
```

Uniformity over `(3)` will make this ambiguity harmless for numerical
nonvanishing, but it remains a typing ambiguity for a named target dipole.

The two THM-2698 event points and their selected left roots are

```text
(x0,r0)=(55232507/115843416,7),
(x1,r1)=(58313459/115843416,3).                           (4)
```

Let

```text
epsilon=1/301082946198216                                 (5)
```

be the proved symmetric packet slack, and let `I_i` be the open interval
`|x-x_i|<epsilon`.  Its indicator is a legitimate refinement of the whole
physical THM-2698 packet at event `i`; in particular every inserted factor
below lives on the same `x`-ancestry.

## 2. Exact paired tables on the open packet cylinders

For a provisional choice `k` from `(3)`, define

```text
K_(i,k)(s)
 = integral_(I_i) d_1(c x-r_i/13)
     u_1(a x-s/13) u_1(k x+s/13) dx,       s in F_13.     (6)
```

At both event centres,

```text
centered(c x_i-r_i/13)=-1/156.                            (7)
```

Thus the deep factor is strictly dangerous.  Its exact `x`-distance to the
nearest tooth wall is

```text
(1/14-1/156)/c
 =71/810903912
 =26361803 epsilon.                                      (8)
```

An exact audit of all `2*4*13` paired cells shows that `(8)` is also the
minimum distance to any of the three factor walls.  Hence every factor in
`(6)` is constant throughout `I_i`, for every shift, and

```text
K_(i,k)(s)=2 epsilon * 1_(s in S_(i,k)).                  (9)
```

The support sets are:

| event | `k_a` | `S_(i,k)` | size | `s=0` |
|---:|---:|:---|---:|:---:|
| `0` | `27` | `{0,3,4,5,8,9,10,11,12}` | `9` | yes |
| `0` | `40` | `{0,1,2,3,4,5,8,9,10,11}` | `10` | yes |
| `0` | `53` | `{0,1,2,3,4,5,8,11,12}` | `9` | yes |
| `0` | `66` | `{0,1,2,3,4,5,8,9,10,11,12}` | `11` | yes |
| `1` | `27` | `{1,2,3,4,7,8,9,10,11}` | `9` | no |
| `1` | `40` | `{1,2,3,4,5,6,7,8,9,10}` | `10` | no |
| `1` | `53` | `{1,2,3,6,7,8,9,10,11}` | `9` | no |
| `1` | `66` | `{1,2,3,4,5,6,7,8,9}` | `9` | no |

This clears the first two numerical gates of the proposed test on the
half-cycle carrier: the paired operation is on common local ancestry, and
its moving shift dependence is nonconstant for every possible pivot.

## 3. Every primitive 13-character survives

Let `zeta` be a primitive thirteenth root and use the normalized transform

```text
Khat_(i,k)(h)
 =1/13 sum_(s in F_13) K_(i,k)(s) zeta^(hs)
 =(2 epsilon/13) sum_(s in S_(i,k)) zeta^(hs).            (10)
```

For any row in the table, put

```text
P_S(z)=sum_(s in S) z^s.                                 (11)
```

Every displayed `S` is a nonempty proper subset of `F_13`.  If
`P_S(zeta^h)=0` for some `h!=0`, then `zeta^h` is again primitive, so its
minimal polynomial

```text
Phi_13(z)=1+z+...+z^12                                   (12)
```

divides the rational polynomial `P_S`, whose degree is at most twelve.
Therefore `P_S` would be a rational multiple of `(12)`, forcing all thirteen
of its zero--one coefficients to be equal.  That would make `S` empty or
all of `F_13`, contrary to the table.  Consequently

```text
Khat_(i,k)(h)!=0
 for i in {0,1}, k in {27,40,53,66}, h in F_13^*.        (13)
```

The companion independently reduces all `96` polynomials in `(13)` modulo
`Phi_13`; each remainder is nonzero.  The result is stronger than the
requested existence of one nonzero shift colour.

## 4. The first endpoint-current obstruction is `eta.u`

For `eta=e_a-e_(k_a)`, the `s`-transform in `(10)` sees the moving-endpoint
residue, with the sign convention of THM-2563,

```text
-eta.v.                                                   (14)
```

A completed THM-2334 target character must instead see

```text
eta.(u-v).                                                (15)
```

The THM-2698 packet state retains physical rail, dynamically typed present
factor, delayed word, predecessor carry, half digit, private root, and
primitive unit.  It does **not** retain an atomic left-present relation index
`u`.  The source indices used to construct its rail/root polynomials are not
THM-2334 relation indices.  Therefore `eta.u` is undefined on this state,
not zero and not recoverable from `r_i`, the owner clock, or the primitive
unit label.  Multiplying `(10)` by a guessed left phase would assume the
missing endpoint map.

This is the first exact obstruction to promoting `(13)` from a moving
boundary signal to a target current.  Uniform nonvanishing over all four
possible `k_a` values does not repair it: the missing datum is a left/right
difference coordinate, not another right factor.

## 5. Two secondary boundaries become exact

First, all four candidates in `(3)` are safe at both old event heads.  The
THM-2563 mechanism instead assumes `k_a` dangerous on the old head, which
forces the entire `s=0` plane to vanish after moving repair.  Here event zero
has `s=0` positive for every pivot.  Event one has `s=0` zero only because
the `a` factor is dangerous there.  Hence there is no inherited uniform
`k_a`-danger zero plane; the positive characters in `(13)` come from the
explicit support geometry, not the THM-2379/THM-2563 anchored-plane argument.

Second, this calculation does not identify two different recurrent objects:

```text
THM-2701:  B_0(y)={13y}, BABA debt orbit of denominator 170;
THM-2698:  B_1(y)={13y+1/2}, fixed delayed phase 11/24
           with a physical two-event lift.                (16)
```

No canonical ancestry map from the first line of `(16)` to the second is in
canon.  The phrase "debt cylinder" in this note refers to applying the
paired operation requested by the debt-cone test to the THM-2698 packet; it
does not assert that the half-cycle is a lift of THM-2701's BABA SCC.

## 6. Strongest survivor and next decisive object

The strongest justified statement is:

> Each canonical THM-2698 event packet contains an explicit positive open
> interval on which every provisionally lawful `a--k_a` paired table has all
> twelve nonzero moving-endpoint characters.

It proves neither a semantic endpoint nor a scalar-row exclusion.  The next
decisive object is not another one-sided shift table.  It is an atomic
common-ancestry refinement

```text
K_i(u,v)=integral w_i(x)
  1_(left relation index=u)
  1_(moving relation index=v)
  [paired blocker--graft factors] dx,                    (17)
```

together with the actual pivot selection.  Only then can one transform in
`eta.(u-v)` and test whether full-`X` aggregation preserves or cancels the
signal.  Alternatively, the THM-2701 BABA SCC must be given its own physical
endpoint ancestry rather than identified with `(16)` by analogy.

No scalar row is excluded and LRC(14) remains open.

## Reproduction

```bash
python3 04-computation/lrc14_half_cycle_paired_debt_cylinder_probe_20260728.py
python3 -O 04-computation/lrc14_half_cycle_paired_debt_cylinder_probe_20260728.py
```

Both runs byte-match

```text
05-knowledge/results/lrc14_half_cycle_paired_debt_cylinder_probe_20260728.out
```

The script is dependency-free, audits exact wall margins and all eight
support rows, and performs explicit cyclotomic reduction for all `96`
primitive-character evaluations.  Its logical guards remain active under
optimized mode.
