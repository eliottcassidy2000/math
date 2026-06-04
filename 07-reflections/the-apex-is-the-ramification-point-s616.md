# The apex is the ramification point (S616)

HYP-2101 gave the apex-lift certificate sheaf a name and a stalk, and stopped at a bridge: "if the section can't
glue, the failure forces a positive-measure interval." I spent this session asking what the *sheaf* actually is —
not the germ, the sheaf — and the answer reorganizes everything I thought was separate.

A certificate is a point that avoids every runner's forbidden line. "Avoid every line indexed by a set of lanes" is
visibly a sheaf condition: the locus over a union of lanes is the **intersection** of the loci. Nothing to prove
there — gluing is automatic, because avoidance is a conjunction. The whole drama is in one word: *nonempty*. The
sheaf always glues; the question is whether what it glues to is inhabited. So H⁰ is the global certificate locus,
and the obstruction — HYP-2101's positive-measure interval — is exactly H⁰ being empty. The resonance, the thing the
user keeps telling me to sidestep, is the H¹ of this sheaf.

Then the apex stops being a special case and becomes a *position*. Over a transverse lane the forbidden set is an
honest line — codimension 1, leaves room. The apex lane forbids the **whole plane**: its stalk has no sections at
all. A single empty stalk empties the intersection. That's the entire obstruction, localized to one lane. And the
r/p lift is not a trick — it's the statement that adding the apex-speed-mod-q coordinate (a unit) turns the
whole-plane forbidder back into a proper hyperplane, refilling the stalk. The lift restores codimension 1 exactly
where the projection had collapsed it. I formalized both halves: the apex empties H⁰, the lift makes its stalk
nonempty again.

The part that made me sit up is the symmetry. Antipodal speeds `v` and `−v` give the *same* forbidden hyperplane —
the equation is homogeneous, so `a x+b y=c` and its negation are one line. So the involution `σ : v ↦ −v` doesn't
permute distinct lanes into new ones; it *identifies* them in pairs. The six nonzero lanes are three σ-orbits. And
the apex lane — residue 0, the multiple-of-n runner — is the one lane with `0 = −0`: the **fixed point** of σ. The
certificate sheaf is a ℤ/2-equivariant sheaf on a double cover, and the obstruction lives precisely at the branch
point. That is the perspective key the user has been circling: it was never chirality, it's the fixed locus of the
antipodal automorphism. Even `n = 2q` forces a `⟨−1⟩` fixed point; the fixed point is an apex; the apex is where the
sheaf can't see a section. The dynamical face fragmenting at 2q (HYP-2150), rigidity leaking through the 2-block
(HYP-2145), the multiple-of-n residual (HYP-2170) — they are one statement read at the ramification point.

So the n=14 proof, in this language, is two clean tasks instead of sixty-four. Show the transverse lanes always glue
to a nonempty section (a line arrangement can't cover the plane unless its slopes cover ℙ¹). Then show the single
apex stalk, once lifted, doesn't re-empty it. The sheaf turned a class table into a fixed-point computation.
