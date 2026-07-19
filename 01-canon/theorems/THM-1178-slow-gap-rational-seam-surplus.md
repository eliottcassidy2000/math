---
id: THM-1178
title: Rational seam surplus for every covered slow gap
status: PROVED quantitative all-cardinality implication.  If r active faster danger combs cover a complete c-slow gap, their harmonic excess above 7-r pays an exact spanning-tree sum of rational handoff quanta.  This strengthens strictness but is scale-sensitive and does not prove universal six-comb noncoverage or LRC(14)
source: codex-2026-07-18-S76
depends_on: [THM-1094, THM-1156, THM-1176]
related: [THM-1166, HYP-7678]
script: 04-computation/lrc14_slow_gap_rational_seam_surplus_codex_20260718.py
output: 05-knowledge/results/lrc14_slow_gap_rational_seam_surplus_codex_20260718.out
---

# THM-1178 -- rational seam surplus

At radius `1/14`, put

```text
D_s={t in R: ||st||<1/14}
```

and let

```text
G=[(14k+1)/(14c),(14k+13)/(14c)]                      (1)
```

be a complete safe gap of the speed-`c` comb.  THM-1176 proves that `r`
faster combs covering `G` force

```text
c sum_i 1/d_i>7-r.                                    (2)
```

The proof there makes equality impossible.  The result below quantifies how
far above equality an integer-endpoint cover must lie.

## 1. Connected nerve and rational handoff quantum

Let `c<d_1<...<d_r`, where `2<=r<=6`, and put

```text
A_i=G intersect D_(d_i).                              (3)
```

Assume every `A_i` is nonempty and the `A_i` cover `G`.  Their intersection
graph is connected.  Indeed, a disconnected vertex partition would express
the connected interval `G` as the union of two nonempty disjoint relatively
open sets.

Every edge `{i,j}` of this graph carries an exact positive quantum:

```text
|A_i intersect A_j|>=1/(14d_i d_j).                   (4)
```

To see this, choose one nonempty connected component of the relatively open
intersection.  Its two boundary points are selected from

```text
(14k+1)/(14c), (14k+13)/(14c),
(14m+-1)/(14d_i), (14n+-1)/(14d_j).                   (5)
```

The positive difference of any two such endpoints is an integer divided by
`14uv`, where the endpoint owners `u,v` lie in `{c,d_i,d_j}` (with the
obvious larger quantum when the owners agree).  Since `c<d_i,d_j`, this is at
least (4).  More precisely, the owner-aware quantum is

```text
gcd(u,v)/(14uv).                                      (6)
```

Thus (4) is a phase-free lower envelope of a stronger endpoint-labelled law.

## 2. The seam-surplus theorem

Choose any spanning tree `T` of the connected intersection graph and define

```text
delta=c sum_(i=1)^r 1/d_i-(7-r).                      (7)
```

> **Theorem (rational seam surplus).** Under the hypotheses above,
>
> ```text
> delta >=(7c/12) sum_({i,j} in T) 1/(d_i d_j).       (8)
> ```

**Proof.**  Let `C(t)=#{i:t in A_i}`.  Coverage gives `C(t)>=1` on `G` and

```text
sum_i |A_i|-|G|=integral_G (C(t)-1)dt.                (9)
```

At a point contained in `q` of the sets, the tree induced by those `q`
vertices has at most `q-1` edges.  Hence, pointwise and after integration,

```text
sum_({i,j} in T)|A_i intersect A_j|
 <=integral_G(C(t)-1)dt.                              (10)
```

The sharp one-comb discrepancy law from THM-1094/1176 is

```text
|A_i|<=|G|/7+6/(49d_i),       |G|=6/(7c).             (11)
```

Combining (9)--(11) gives

```text
integral_G(C(t)-1)dt
 <=sum_i(|G|/7+6/(49d_i))-|G|
 =6 delta/(49c).                                      (12)
```

Finally (4), summed over `T`, and (10)--(12) yield

```text
sum_T 1/(14d_i d_j)<=6 delta/(49c),
```

which is exactly (8). ∎

If a larger packet covers `G` but some comb never meets it, discard every
inactive comb first.  The theorem applies to the active set.  One active comb
cannot cover a slow gap, and THM-1176 already rules out active cardinality at
most three, so only `r=4,5,6` can survive.

## 3. Ordered and H-drift forms

For every tree on the ordered vertices, root it at `d_r`.  If `p(i)` is the
parent of `i<r`, then `d_(p(i))<=d_r`, whence

```text
sum_({i,j} in T)1/(d_i d_j)
 >=(1/d_r)sum_(i=1)^(r-1)1/d_i.                       (13)
```

Therefore every cover obeys the tree-free corollary

```text
c sum_i 1/d_i-(7-r)
 >=(7c/(12d_r))sum_(i=1)^(r-1)1/d_i.                 (14)
```

For the six-comb slow-gap branch of THM-1176, write

```text
H_6=c sum_(i=1)^6 1/d_i.
```

Then

```text
H_6-1 >=(7c/(12d_6))sum_(i=1)^5 1/d_i.               (15)
```

In normalized tooth coordinates `x_i=c/d_i`, this is

```text
H_6-1 >=(7/(12c))x_6 sum_(i=1)^5 x_i.                (16)
```

Equation (16) is the discrete functional form hidden behind the qualitative
strictness in THM-1176.  It is deliberately scale-sensitive: under a common
dilation the right side decays like `1/c`.  If additionally `d_r<=Kc`, then
order gives the explicit finite-carrier floor

```text
delta>=35/(12K^2 c)                    (r=6).          (17)
```

Thus a covered packet cannot approach the harmonic equality surface faster
than the rational seam lattice permits.

## 4. Relation to chi7, Fano, and Kakeya needles

The protected slow gap is the Kakeya needle.  Connectedness forces a tree of
positive handoffs along it; exact zero seams cannot cover a closed point with
two open teeth, agreeing with THM-1156's third-support law.  Formula (6) keeps
the endpoint owners and their gcds, so it is a local counterpart to
THM-1166's global pair/Fano gcd budget.

The ordinary runner-order tournament, with observable
`sign(c/d_i-c/d_j)`, is transitive: score histogram `0,...,r-1`, no directed
cycles, singleton SCCs, and one Hamiltonian path.  It retains (13) but loses
which pair actually hands off, the component length, endpoint owners, gcd,
and triple support.  The faithful object is instead

```text
connected nerve tree on active combs
 + one endpoint-labelled overlap component per tree edge
 + chronological position on G.                      (18)
```

Orienting the tree edges by the first overlap position supplies a chronology,
but completing it to a tournament destroys the multi-owner events.  This
records the challenged assumption explicitly: the proof vertices are comb
families only as a first layer; the quantitative vertices are rational
handoff components.

## 5. Scope and remaining frontier

The theorem is uniform in `c`, the phase `k`, the active cardinality, and all
faster integer speeds.  It supplies a genuine positive correction to the
zero-surplus equality surface and a new exact consumer for endpoint/Fano
data.  It does **not** give a scale-invariant correction.  Consequently it
does not by itself exclude covers with large common dilation, prove the
universal six-comb slow-gap conjecture, close the sporadic `n=12` branch, or
prove LRC(14).

The next useful composition is to combine (6)/(8) with a phase-independent
gcd or beat constraint that prevents every spanning-tree quantum from
shrinking simultaneously.  THM-1175's defining beat pair redundancy and
THM-1166's Fano line periods are natural suppliers; a bare speed-order
tournament is not.

The dependency-free exact replay checks 13,035 phase-labelled pair banks,
7,185 nonempty intersection components (the minimum component/quantum ratio
is exactly one), all 1,441 labelled trees through six vertices, and 85,787
tree/activity masks.  Normal and optimized runs are byte-identical to the
stored output.  Frozen SHA-256 hashes are

```text
source  fd0997cfa4c9ee90b1fe4c1496a276bb64d9a850803ec64fba61aabdfbcf750a
output  8ca6026ec0b1a1ba7d348400c24eb7fc590774deda8ca2cdd0313593c24fb6e5
```
