---
id: THM-3878
title: "LRC(14) eleven-plus-two seams collapse by harmonic absorption"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In THM-3818's
  exact rank-eleven two-component 11+2 branch, THM-668 harmonic absorption
  closes every seam of scale at least three and every scale-two seam with an
  even pair coordinate.  Exact two-lift geometry closes (1,3,2).  The
  necessary ledger falls from 52,692 to 7,505; LRC(14) remains open.
source: root / THM-3818 cyclic seam and THM-668 dispatch join, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS.  A separate 93,524-gate checker rebuilds
  the atlas through Euler-totient blocks, counts seams by divisor formulas,
  constructs every obstructed odd-pair phase from the cross-centre gap
  1/(2pq), and computes the control row by exact pair-sum events.  It confirms
  45,186 absorption closures, the unique universal scale-two atlas pair
  (1,3), and the final 7,505 necessary triples.  It also proves the (3,7,2)
  control has maximum 1/10 and connected decoder graph, so it is only a
  pair-selector hostile.  Normal, optimized and frozen streams byte-match.
depends_on:
  - THM-3818-scaled-inert-cubeclass-support-two-pair-packet
  - THM-668-detuned-harmonic-dispatch
related:
  - THM-3743-lonely-runner-polyhedron-khinchin-flatness-relation-reduction
script: 04-computation/lrc14_eleven_plus_two_harmonic_absorption_seam_collapse_thm3878.py
output: 05-knowledge/results/lrc14_eleven_plus_two_harmonic_absorption_seam_collapse_thm3878.out
script_sha256: 246dcb77753616aa399300daad62adaedfa838a148ea1b63edf5f75e4f4eae69
output_sha256: 2c4ad16022fddbb20bcd3a407bdd32cececc3e8e52ca6201f132a92d05767500
semantic_sha256: 58a950770c4984d4a1a3f4a4031a360042ca471ed406d529dc28aad7812c1de6
independent_audit_script: 04-computation/lrc14_eleven_plus_two_harmonic_absorption_seam_collapse_independent_audit_thm3878.py
independent_audit_output: 05-knowledge/results/lrc14_eleven_plus_two_harmonic_absorption_seam_collapse_independent_audit_thm3878.out
independent_audit_script_sha256: b9760e34de5779e5ccf328b12058feb1769965c0e7d614b5fd89f779141a9143
independent_audit_output_sha256: 6dfca0670072ca39c6cba48e895cc79abcb4548b17223989518bbd7c657450a3
independent_audit_semantic_sha256: 22bffbd1a61a3e19596408c6c73bc1db6562896c3e314e327b285e35dfe73b86
hash_basis: raw LF bytes
---

# THM-3878 -- harmonic absorption collapses the 11+2 seam ledger

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This is a
bounded theorem inside THM-3818's exact oriented rank-eleven two-component
branch.  It proves no statement about other component shapes, physical entry,
owner/arrival, or arbitrary 14-runner rows.  In particular, **LRC(14) remains
OPEN.**

## Bounded theorem

Assume THM-3818's exact `11+2` equality branch

```text
n = s u direct-sum t(p,q),
s,t>0, gcd(s,t)=1,
p<q, gcd(p,q)=1,
```

where `u` has eleven coordinates and `(p,q)` belongs to the `5,855`-ratio
cube-class atlas.  Then a hypothetical counterexample must satisfy exactly
one of the following necessary alternatives:

```text
s=1; or
s=2, p and q odd, and (p,q)!=(1,3).
```

Equivalently, all `s>=3` seams, all even-coordinate `s=2` seams, and the
odd cell `(1,3,2)` are lonely.

### Harmonic absorption

Suppose `s|p`.  Rewrite the row as

```text
n = s( u direct-sum (tp/s) ) direct-sum {tq}.          (1)
```

The parenthesized pack contains twelve nonzero speeds.  Since `gcd(p,q)=1`,
`s|p` implies `gcd(s,q)=1`; together with `gcd(s,t)=1`, this gives
`s` not dividing `tq`.  THM-668 applies to (1) and supplies a common time of
clearance at least `1/13`, strictly above `1/14`.  The case `s|q` is
symmetric.

THM-3818 already proves that every surviving tariff seam with `s>=3` divides
one of `p,q`, so harmonic absorption eliminates all `40,982` such triples.
At `s=2`, exactly one coordinate is even unless both are odd; absorption
eliminates another `4,204` triples.

### Exact odd scale-two geometry

It remains to understand `s=2` with `p,q` odd.  Choose any `y` at which the
eleven-speed component `u` is `1/14`-safe.  Since `t` is odd, its two cyclic
lifts present the pair with phases

```text
z=t y/2,             z+1/2.                            (2)
```

Let

```text
D_w={z in R/Z: ||wz||<1/14}.
```

Both lifts fail precisely when

```text
z in (D_p union D_q) intersect
     ((D_p union D_q)-1/2).                            (3)
```

For odd `w`, `D_w` and `D_w-1/2` are disjoint: their centres are at least
`1/(2w)` apart, while the two open radii total only `1/(7w)`.  A cross term,
say `D_p intersect (D_q-1/2)`, compares centres

```text
a/p,                  (2b+1)/(2q).
```

Because `p,q` are coprime and odd, the numerator

```text
2qa-p(2b+1)
```

is always odd and can equal `1`: Bezout solves
`2qa-2pb=p+1`.  Hence the minimum circular centre gap is exactly

```text
1/(2pq).
```

The two open radii sum to `(p+q)/(14pq)`.  Therefore (3) is empty exactly
when

```text
p+q<=7.                                                (4)
```

Among the THM-3818 atlas's odd pairs, only `(1,3)` satisfies (4).  This proves
that every phase has a pair-safe lift for `(1,3,2)`, while retaining all
eleven `u`-clearances from the chosen `y`.

## Exact census

Two independent atlas builders—trial factorization with a pair loop and an
SPF sieve with a sum/totient loop—give identical objects.  Two independent
seam builders—direct scale scanning and divisor-set union—also agree.

```text
primitive cube-class ratios                         5,855
original necessary seams with s>=2                46,837
  s>=3                                             40,982
  s=2                                               5,855

closed by harmonic absorption                     45,186
  all s>=3                                         40,982
  s=2 with one even coordinate                     4,204

odd-coordinate s=2 residual                        1,651
closed by universal two-branch geometry: (1,3)         1
unresolved s=2                                     1,650
unresolved s=1                                     5,855
combined exact residual                            7,505
```

The census is **FINITE-EXACT**; the absorption implication and odd-pair
criterion are proved for all tuples in their stated symbolic scopes.

## Minimal hostile and controls

The first atlas pair for which the universal two-branch selector fails is

```text
(p,q,s)=(3,7,2),             z=13/20.
```

Indeed,

```text
min(||3z||,||7z||)=1/20,
min(||3(z+1/2)||,||7(z+1/2)||)=1/20,
```

so a different pair runner kills each lift.  This is a hostile to the
**universal pair-only selector**, not a counterexample to LRC(14).

The phase-uniform selector obstruction can be realized at an actual
eleven-speed good point.  This is a method control, not an example in the
rank-eleven two-component equality branch: its full THM-3818 decoder graph
is connected, with 33 decoder edges.  Take

```text
u=(1,2,3,4,5,6,7,8,9,11,12),       y=3/10,
n=2u direct-sum (3,7).
```

Every coordinate of `u` is safe at `y`; the selected lifts `3/20,13/20`
both have row clearance `1/20`.  But the same full row is safe at `1/20`
with clearance `1/10`.  This explicitly prevents the selector hostile from
being misread as an LRC hostile or as a realization of an unresolved
two-component seam.

The off-atlas pair `(1,5)` is a positive geometry control: it also satisfies
`p+q<=7` and is universally two-branch safe, but it is absent from the cube
atlas because its sum has a factor `3`.  This keeps geometric classification
separate from atlas membership.

## Typed connection ledger

```text
source:      THM-3818's oriented 11+2 quotient and cyclic tariff
target:      THM-668's scaled twelve-pack plus one detuned runner
map:         s|p sends tp into the pack as s*(tp/s), symmetrically for q
preserved:   all eleven u-clearances and the absorbed pair-runner clearance
destroyed:   the chosen LRC(<=13) time, winning branch, owner, and arrival
sidecar:     gcd(s,t)=gcd(p,q)=1, ensuring the other runner is detuned and
             its branch orbit is complete
decisive test: exact 46,837-seam census plus rational scale-two wall cells
```

There is a real but limited connection to the factorial local-rule audit.  In
both cases, a local divisor condition is a **selector for the next proof
mechanism**, not a stable global barcode.  At factorial `d=9996`, divisor
places leave degree `3998` and an adaptive nondivisor `p=19` sidecar is
needed.  Here `s|p` does something stronger than filtering: it changes the
decomposition itself by absorbing one runner into a coherent harmonic pack,
after which THM-668 is proof-forcing.  There is no claimed arithmetic map
between factorial primes and LRC scales; the shared content is only this
typed selector/sidecar architecture.

## Remaining frontier

The bounded result leaves exactly:

- all `5,855` scale-one triples, where no cyclic branch orbit exists; and
- `1,650` odd-coordinate scale-two triples beginning with `(3,7,2)`.

The hostile shows that the latter require information about which points of
`G(u)` are available, or an owner/arrival sidecar; pair geometry alone cannot
choose a safe branch.  Component shapes `3+10` through `6+7`, the rank-twelve
terminal, physical entry, and LRC(14) are untouched.

## Exact replay

Run from the repository root:

```bash
python3 -B 04-computation/lrc14_eleven_plus_two_harmonic_absorption_seam_collapse_thm3878.py
python3 -B -O 04-computation/lrc14_eleven_plus_two_harmonic_absorption_seam_collapse_thm3878.py
python3 -B 04-computation/lrc14_eleven_plus_two_harmonic_absorption_seam_collapse_independent_audit_thm3878.py
python3 -B -O 04-computation/lrc14_eleven_plus_two_harmonic_absorption_seam_collapse_independent_audit_thm3878.py
```

The 96,102-gate primary companion independently builds the atlas and seam
union twice and checks all 1,651 odd pairs by exact rational wall cells.  The
93,524-gate audit uses different atlas, divisor-count, modular-inverse and
pair-sum-event routes.  Normal and optimized executions byte-match their
raw-LF frozen outputs.  The finite census controls the ledger; the absorption
and odd-pair arguments above prove their stated symbolic scopes.  **QED.**
