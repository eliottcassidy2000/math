# Session 109: Pisano clocks, ramified symmetry, and exact attachment

2026-10-30. **Completed mathematical packet: elementary proofs, independent audits,
and reproducible exact checks.** Classical inputs and recovered repository results
are credited; no global originality claim is made.

## Main results

**1. The ramified gap becomes an exact recursive tournament bank.**
Let mu²−a mu−1=0 and let p be an odd prime with v_p(a²+4)=1. On the additive
group of Z[mu]/p^k, the Cayley tournaments invariant under multiplication by
mu⁴ have exactly B_k=p^k−1 independent orientation bits:

    B_0=0,   B_k=p B_(k−1)+(p−1),
    T_k=2^B_k=2^(p−1) T_(k−1)^p.

For Fibonacci (a=1,p=5), B_k begins4,24,124. Each parent orbit pair has five
children; the kernel adds four pairs. A fixed parent tournament has16 coherent
lifts, each an ordinary lexicographic product with a25-vertex kernel tournament.
Unrestricted choices need not descend. Fully coherent towers have2^(4k) choices,
against2^(5^k−1) unrestricted choices at the top.

These are **labelled, specified-symmetry counts on p^(2k) vertices**, not A000568
or tournament-isomorphism counts. This extends the explicit ramified exclusion
in the earlier all-depth unramified theorem U1.

**2. The same ramification identifies exactly what a trace forgets.**
For Q the Fibonacci matrix and A=Q², adjacent traces have pairing determinant5.
At5^e the matrix/trace clocks are4·5^e versus4·5^(e−1) for Q, and2·5^e versus
2·5^(e−1) for A. Another CRT component can mask this loss. In the general
ramified algebra, an element of remaining uniformizer depth r has full clock
p^floor(r/2), trace clock p^floor((r−1)/2). Only even layers lose a factor p.

One extra p-valued **snapshot** restores a lost p-state phase. For p=5, a
binary sensor with cyclic word00011 can instead recover it in three serial
readings. This trades observation time for alphabet size; it does not create
information or supply a free implementation of the sensor.

**3. Affine clocks retain a different missing coordinate.**
We give the Smith factors of I−Q^r for every positive stride. At strides5,6,7,9,10,
the obstruction groups are Z/11, (Z/4)², Z/29, Z/2×Z/38, (Z/11)².
For x↦Ax+b, if t=ord(A), then

    ord(x↦Ax+b)=t·ord_add((I+A+...+A^(t−1))b).

Stride5 modulo11 has the same linear order2 but affine orders2 or22. Conversely,
stride4 modulo5 has nonremovable carries that do not enlarge the clock. Neither
an obstruction's existence nor its cardinality alone determines the period.

**4. Preserving, reversing and quotienting symmetry are separated.**
An automorphism of a finite odd abelian group preserves some Cayley tournament
iff it has odd order. It reverses one iff its order is2m with m odd and its
midpoint is−I. Applying this to Pisano powers gives an exact composite-modulus
criterion. Quadratic conjugation and reversal act freely as a specified Klein
observer quotient. That quotient is not full graph-isomorphism classification.

## Read the packet

- [Synthesis, composition, 3/5/19, and temporal bits](SYNTHESIS.md)
- [Prior Pisano, golden-ring, cat-map and Zeckendorf work](history/RECOVERED_PISANO_ADJACENT_WORK.md)
- [Observer clock proofs](clocks/PROOFS.md)
- [Affine Smith and carry proofs](affine/README.md)
- [Tournament preserving/reversing and ramified bank proofs](tournaments/results.md)
- [Independent ramified attachment and lexicographic-product audit](audit/RAMIFIED_ATTACHMENT_AUDIT.md)
- [Independent clock/carry audit](audit_clocks/AUDIT.md)
- [Historical missing-hypothesis correction](history/PISANO_LIFT_CORRECTION.md)
- [Verification receipt](verification.json)

## Reproduce

From this directory run:

```sh
python3 verify_packet.py
```

Only the Python standard library is required. The runner verifies the SHA256
manifest, copies the packet into isolated temporary directories, and runs
seven checkers both normally and with Python optimization. All14 runs must
match their recorded JSON outputs and their paired runs. Validation checks
remain active under `-O`.

Independent evidence includes151,516 clock/carry checks,14,868 exhaustive carry
cases, all ramified points through15,625 vertices,390,000 ordered adjacency
comparisons for the lexicographic bridge,728 preserving/reversing constraint
comparisons, and general ramified controls at13,17,29,41. These scopes overlap
and are not added into one misleading total. The general claims rest on the
written proofs, not on extrapolation from the finite checks.

## Publication and limits

First verified main checkpoint: `94f5716c76a1056590bd70acd8cefc15f41654f2`.
The final checkpoint completes these branches and adds a narrow warning to the
historical THM-491: the generic prime-power Pisano lift formula needs the
return-depth parameter h. No exceptional Fibonacci prime is asserted.

The exact common language with Collatz is finite affine composition. Integer
sources, exact valuation guards, ordered words and original-source descent
remain necessary. No Collatz convergence, Wall–Sun–Sun classification, modular
category saturation, LRC closure, or unrestricted A000568 enumeration is claimed.
