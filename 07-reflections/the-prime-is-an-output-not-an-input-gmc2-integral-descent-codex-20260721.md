# The prime is an output, not an input: direct integral descent for GMC(2)

The original proof of THM-2022 chooses a large rational prime, chooses a
prime of a number field above it, clears denominators, and reduces.  That is
mathematically sound, but it is not the shortest formal route.

The decisive change is that the formal multinomial-isolation lemma is
stronger than the elementary two-digit Kummer narration: for a prime `p`, a
multinomial channel of total mass `p*m` which is not componentwise divisible
by `p` already vanishes modulo `p`.  It does **not** require `p > m`.

Consequently the proof need not choose its prime.  Starting from a complex
torus point `c` and a nonzero lowest-face seed `Q(c)`, form

```text
A = Z[c_i, c_i^(-1), Q(c)^(-1)] subset C.
```

This is a nontrivial finite-type `Z`-algebra.  A maximal quotient of `A` is a
finite field of some prime characteristic `p`; all coefficient coordinates
and `Q` survive because they were made units.  The prime is now an output of
the specialization.

There is a quantifier subtlety here.  The normalized moment polynomial used
in the contradiction depends on the resulting `p`, so it cannot be placed in
a finite relation list chosen before specialization.  The correct interface
is universal:

```text
for every f in Z[X_i],  f(c)=0 implies f(c_bar)=0.
```

After the quotient reveals `p`, choose the integral normalized relation at
mass `p*m0` and minimum height `p*A0`.  Complex moment-nullity makes that
relation zero before specialization; the universal property makes it zero
afterward.  The residue decomposition should make the same value `Q_bar^p`,
which is nonzero.

This rerouting removes from the main formal proof:

- number-field construction as a required step;
- denominator and integral-model bookkeeping in a number field;
- selection of a good prime ideal;
- the artificial largeness condition `p > m0`.

The number-field/good-reduction modules remain valuable as an independently
checked route, but they are no longer on the critical path.

The remaining formal bottleneck is now sharply combinatorial: reindex the
surviving normalized channels by componentwise dilation, identify the
equality-face image with the undilated constant-term channels, and feed the
three cases into the generic residue-assembly theorem.  The face-height
floor, normalized-moment cancellation, Lucas congruence, off-face factorial
gap, and Frobenius power law are already separate checked interfaces.

## Assumption challenge

The discarded assumption was that arithmetic control requires choosing a
prime first.  It does not: making the finite set of quantities that must
survive into units lets an arbitrary maximal quotient choose the prime.

A second tempting compression is to orient channels pairwise and treat them
as tournament vertices.  That quotient preserves comparisons such as
"lower face versus off face", but destroys coefficient magnitudes and tied
whole-face cancellation.  The Frobenius step acts on the complete face sum,
not on a winner selected by pairwise orientation.  Channels are useful as
vertices only while the tied face remains a retained hyperedge/finite sum.
