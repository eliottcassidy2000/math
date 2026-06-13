# HYP-2458 - Faulhaber odd moments are the scalar odd-atom layer of a hidden OCF-style compatibility lift

**Status:** OPEN synthesis; exact identity verified and asymptotic scout supports the proposed anchor.
**Source:** codex-2026-06-13.
**Companions:** HYP-2457, HYP-2456, HYP-2455, HYP-2454, HYP-2453, HYP-2444,
HYP-2443, HYP-2425, HYP-2426, THM-466, OCF/`Omega(T)`.
**Computation:** `04-computation/faulhaber_odd_moment_ocf_bridge_codex.py`;
stored output `05-knowledge/results/faulhaber_odd_moment_ocf_bridge_codex.out`.

This is an addendum to HYP-2457, not a replacement for it.  HYP-2457 carries
the stronger formal anchor expansion through `gamma_p/u^2`; this note asks what
OCF-style compatibility object should be attached once the odd Faulhaber atom
inventory has been isolated.

## Statement

Let `a_p(n)` be the real anchor solving

```text
sum_{j=0}^n (a+j)^p = sum_{j=1}^n (a+n+j)^p.
```

Set `c=a+n`, so the left interval is `c-n,...,c` and the right interval is
`c+1,...,c+n`.  Then the exact balance defect is

```text
D_p(c,n)
 = c^p + sum_{j=1}^n ((c-j)^p-(c+j)^p)
 = c^p - 2*sum_{r odd<=p} binom(p,r)c^(p-r)S_r(n),
```

where `S_r(n)=sum_{j=1}^n j^r`.  Thus every even Faulhaber moment cancels from
the balance equation.  Only odd moments appear.

The working hypothesis is that these odd moments are the scalar "odd atom"
layer of a larger hidden-lift object, in the same sense that the OCF starts
with odd directed cycles but only becomes `H(T)` after retaining compatibility:

```text
H(T) = I(Omega(T),2) = sum_k 2^k alpha_k(T),
```

where `alpha_k` counts compatible `k`-packets of vertex-disjoint odd cycles.
The power-balance side has the odd moment inventory `S_1,S_3,...`; the missing
proof object should be an OCF-style compatibility ledger telling which odd
moment atoms can coexist across shell, support, carry, and boundary constraints.

## Exact Low-Power Anchors

Let `u=n(n+1)`.  The formula gives:

```text
p=1: D_1(c,n)=c-u, so c=u and a=n^2.
p=2: D_2(c,n)=c^2-2*c*u, so the nonzero root is c=2u and a=2n^2+n.
```

These recover the two exact towers:

```text
p=1 anchor: a=n^2, center c=n(n+1)=2T_n.
p=2 anchor: a=2n^2+n, center c=2n(n+1)=4T_n.
```

The square-pyramid packing sits one layer away from the balance equation:

```text
6*S_2(n)=n(n+1)(2n+1).
```

So the p=2 tower has a geometric square-pyramid/cuboid packaging, while its
antisymmetric balance is driven by the odd first moment `S_1(n)`.  This is the
new warning signal for the broader project: visible volume and proof-driving
imbalance may be adjacent layers rather than the same layer.

## Asymptotic Carrier

The prompt's expansion is best written at the center as

```text
c = p*u
  + (p-1)(p-2)/(12p)
  - (p-1)(p-2)(2p^2-4p-1)/(180*p^3*u)
  + O(n^-4),
```

and therefore

```text
a_p(n)
 = p*n^2 + (p-1)*n
 + (p-1)(p-2)/(12p)
 - (p-1)(p-2)(2p^2-4p-1)/(180*p^3*n(n+1))
 + O(n^-4).
```

The stored computation verifies the odd-moment identity exactly for `792`
small checks, recovers the exact `p=1,2` anchors, and numerically compares the
true bisection root against the expansion for `p=3..8` and `n=10,30,100`.
The error times `n^4` stabilizes in each row, matching the claimed scale.

This makes `u=n(n+1)` the shared triangular carrier:

```text
p=1,2: Bernoulli correction vanishes and the address is integral.
p>=3: Bernoulli/Faulhaber correction produces a fractional address.
HYP-2456: the same carrier has a Beatty/Pell endpoint-carry layer.
```

## OCF Analogy And Its Limit

The analogy to OCF should be kept sharp:

```text
Faulhaber side: odd moments S_r are single-layer odd atom totals.
OCF side: odd cycles are single-layer odd atom totals.
OCF lift: alpha_k records compatible independent packets in Omega(T).
Missing lift: an odd-moment compatibility object should record which moment,
              shell, support, carry, and boundary atoms can coexist.
```

It would be a mistake to claim that ordinary power sums are literally the
tournament OCF.  The better transfer is procedural: first identify the odd
atoms, then refuse to collapse compatibility.

The script illustrates this with two disjoint directed triangles and all cross
edges one way.  There are `2` odd cycles and `1` disjoint pair, so

```text
I(Omega,2)=1+2*2+4*1=9=H(T).
```

The single-cycle inventory alone would only see `2`; the compatibility packet
is what turns it into the Hamiltonian-path count.  The analogous next step for
Faulhaber/LRC/code72 is to construct the compatibility packet ledger.

## Tournament Analysis

The session's carrier tournament deliberately uses proof carriers, not runners
or arcs.  Alternate vertex sets considered: runners, gaps, fixed circle
sections, boundary events, residues, cover arcs, Fourier modes, moment atoms,
odd cycles, hidden lifts, and proof obligations.  The chosen vertices were:

```text
faulhaber_odd_moments
ocf_odd_cycles
omega_alpha_packets
beatty_pell_addresses
convolution_hidden_lifts
lrc14_q27_support
code72_support_design
unit_distance_nonproduct_fibers
```

Pairwise observable: majority comparison over seven retained-channel criteria:
odd atom visibility, compatibility layer, address/carry layer, support transfer,
exactness, computation readiness, and hidden-lift retention.

Switch/gauge: higher retained-channel score wins; ties use the displayed
Hamiltonian path order.

Fingerprint:

```text
score_hist={0:1, 1:1, 4:3, 5:3}
directed_3cycles=8
scc_sizes=[6,1,1]
hamiltonian_paths=45
```

For LRC, this quotient preserves blocker/support address, carry class, and
witness target.  It destroys raw time geometry and exact interval endpoints.
The challenged assumption is that the useful tournament must live on runners or
runner pairs; here the live tournament is on odd atoms and proof obligations.

## Transfer Targets

1. **LRC14 Q27.**  Upgrade the HYP-2456 address ledger by adding odd-moment or
   antisymmetric imbalance atoms.  A scalar "q blocked" row should retain
   owner support, divisor fiber, carry residue, endpoint atom, and the odd
   imbalance atoms that force or prevent compatibility.

2. **The `[72,36,16]` code.**  Treat the weight enumerator as a boundary total
   and the `5-(72,16,78)` minimum-word design as the first support lift.  The
   recurring `78/90` triangular shadow should become a compatibility constraint
   on support packets, not a scalar numerology claim.

3. **Polynomial irreducibility.**  HYP-2452 gives the convolution lift.  This
   hypothesis suggests adding an "odd atom" side channel: factor-capture
   values, residue classes, Newton slopes, and coefficient moments should be
   retained only insofar as they constrain compatible factor grids.

4. **Unit distance.**  Product-reducible unit-distance constructions are the
   baseline hidden lift; non-product Moser fibers are the compatibility layer.
   The OCF lesson is to measure packet coexistence, not only individual norms.

## Open Problem

Build an explicit odd-moment compatibility lift `M_p(n)` whose one-particle
shadow is the odd Faulhaber list `S_1,S_3,...`, and whose packet terms predict
when a scalar shell balance can be lifted to an LRC support ledger, polynomial
factor obstruction, or code72 support-design incidence constraint.
