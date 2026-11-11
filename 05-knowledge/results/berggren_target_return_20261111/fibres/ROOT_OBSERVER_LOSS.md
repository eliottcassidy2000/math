# The full seven-state root observer still loses arbitrary-star content

**PROVED for an infinite family, with independent literal controls.** This
is an operation-scope separation, not a defect in the inherited exact
seven-state formulas for their declared consumers.

Let P and P' be two triangle-chain sources with the same r>=2 triangles
and the same g>=r-1 singleton vertices. In P separate every consecutive
triangle pair by at least one singleton. In P' move all singletons to the
left. Append a common new global sink root rho to each graph. Then their
entire rooted response vectors `(H,E,S,U,V,d,K)` agree:

```
(3^r, 3^r, 0, 3^r, 12^r*2^g, 3^r, 12^r*2^g).          (1)
```

Here E,S count full paths ending/starting at the root; U,V count ordered
two-bin covers with the root bin starting/ending there; d is the root-
deletion H; K is the deletion graph's full ordered two-bin count.

Proof: every full path ends at rho, so E=H=3^r and S=0. A path starting
at rho must be the singleton rho, giving U=H(P)=3^r. A path ending there
can use any old subset, giving V=Q_P(2). Also d=H(P), K=Q_P(2). The
ordered-cover polynomial multiplies under ordinal joins, with Q_C3(2)=12
and Q_K1(2)=2, so Q_P(2)=12^r*2^g independent of the gaps. The same proof
applies to P'.

Nevertheless their first extension contents are different:

```
C_1(P->rho)=3^(r-1),
C_1(P'->rho)=3^(r-2).                                  (2)
```

This follows from the complete first-step theorem in
[TRIANGLE_FILTRATION.md](TRIANGLE_FILTRATION.md), since appending rho only
changes exterior padding. Thus no function of the full seven coordinates
alone computes the gcd of ALL possible one-vertex extension counts.

The smallest displayed pair is

```
C3 -> K1 -> C3 -> rho       versus       K1 -> C3 -> C3 -> rho,
state=(9,9,0,9,288,9,288),   C_1=3 versus1.
```

`gap_transport.py` builds both actual eight-vertex graphs, independently
counts every root coordinate by subset paths, and computes each complete
star polynomial. The two graphs are also linked by the named neutral
macro, with every original vertex retained.

By the inherited exact formulas, these equal states remain indistinguishable
to any common O/Z/J continuation and to each fixed two-sided rooted
amalgamation context. Those consumers are developed in
[dynamics/PROOF.md](../../h13_paired_descent_20261111/dynamics/PROOF.md) and
[ADDITIVITY.md](../../h13_paired_descent_20261111/modular_cycles/ADDITIVITY.md).
Arbitrary vertex attachments are a larger operation language, and they
distinguish the pair by (2). The missing sidecar is actual ordered SCC
placement or the full attachment polynomial, not another scalar label
chosen merely because four or five coordinates have appeared elsewhere.

This connects the Berggren H9 collision to the return problem precisely:
source information must be adequate for the NEXT allowed operation. The
four-dimensional append compiler is minimal in its own linear language;
seven rooted coordinates are exact for their own contexts; neither claim
implies completeness for all metagraph edges or proves Chan–Pak.
