# Twin lengths are sometimes free parameters, sometimes a synchronization law

**Research synthesis; independent review PASS WITH MINOR SCOPE REPAIRS. Chan–Pak remains OPEN.**

The user's warning about unselected twin blocks identifies a real distinction:
clones have identical external incidence, but a Hamilton path can leave and
return to their block several times. Increasing the length changes the number
of ways to allocate vertices to those visits. The selected endpoint tax need
not change in the same way as the Hamilton count. The full signed response,
not the quotient alone, determines which lengths are allowed.

## 1. A complete answer for each fixed single clock

The [length theorem](lengths/PROOF.md) gives a finite-or-free dichotomy. Fix the
exterior tournament, the selection, and one unselected port. Let ν be the
maximum matching of directed edges from the port's out-neighbors to its
in-neighbors. This matching is the inherited return-capacity invariant. The
selected endpoint constraints sharpen the common mode budget to

    R=min(ν+2,1+min(d+,d−)).

Either every positive length preserves the unit defect, or at most R−1
positive lengths do. Any R distinct successful lengths certify freedom. A
finite exact algorithm determines the exceptional set. Long delayed failure
has a structural cost: success through length L followed by failure requires
at least2L+1 vertices in the original core and ν≥L−1. The result is not a
claim that every such budget is realizable.

## 2. Individually forbidden inflation can become a reusable coupled move

The [new six-vertex gate](wildcard/PROOF.md) has two unselected transitive
modules of lengths a,b and exact defect

    δ(a,b)=1+2^a−2^b.

Thus the entire positive unit locus is a=b. Neither coordinate is free at
any unit point, but increasing both is an all-order physical operation. This
is a counterexample to treating independent length constraints as a complete
classification of a multi-clock chart.

The unit slope on the diagonal is

    U_m=8·6^m+10·4^m+2·3^m−4·2^m.

Its odd part is `c_m=4·6^m+5·4^m+3^m−2·2^m`. The inherited tail gives actual
counts1+c_m2^s for every s≥2. It does not supply the s=1 heads or every odd
coefficient. The family is recursively usable and still sparse.

More generally, every prescribed affine length ray in a fixed acyclic-port
chart is decidable. A finite set of prime-exponent linear equations classifies
all directions that preserve the unit condition forever. The mechanism is
exact collision and cancellation of product bases, not approximate agreement.
This does not classify arbitrary nonlinear solution loci or cyclic port charts.

## 3. A physical Farey gate, with an exact growth boundary

Replace those same two ports by arbitrary nonempty tournament modules A,B.
Writing `(h,c)=(H(A),C2(A))` and `(g,d)=(H(B),C2(B))`, the complete finite
skeleton gives

    δ=hg+cg−hd.

Consequently a unit seed is exactly the determinant-one equation

    g(h+c)−hd=1.

The fractions `(h+c)/h` and `d/g` are genuine adjacent Farey fractions, but
only **physically realized** module profiles are admissible. The full response
also retains C3(A) to compute the seed's slope U. Literal nontransitive
witnesses include a pair with profiles(5,54) and(589,6950); both modules are
nontransitive, and their substituted17-vertex seed is counted independently.
For each fixed physical A, a path-cut bound gives a finite upper bound on the
order of every possible partner B.

The most natural continued growth operation has a complete negative boundary:
independently appending or prepending transitive ordinal padding to the two
modules preserves a unit pair only when there is no padding, unless both
modules are transitive. In that exceptional case precisely equal padding
lengths work. Thus the nonlinear module lift creates real new seeds but does
not silently inherit the chain family's renewable operation.

## 4. Compressing all the known free-clock families without double counting

The [stable-length grammar](stable/PROOF.md) translates free subsets of an
acyclic port set into a simplicial complex. Minimal mixed tensor supports are
its forbidden faces; their squarefree monomials form an exact obstruction
ideal. The tensor, not merely the ideal, must be kept for isolated or coupled
resonances.

Compiling the inherited stable-port bank through order8 gives98 irredundant
fixed-clock cylinders on88 marked twin quotients. Their **entire all-size
union**, with overlaps and marked isomorphisms removed, has generating function

    (z^3+4z^5+8z^6+20z^7+46z^8)/(1−z)
      +(5z^7+14z^8)/(1−z)^2.

It contains exactly19n−49 marked classes at each n≥8. This is not the full
unit-seed sequence and not all adaptive descendants that might unlock new
stable ports. The finite source census is inherited; every retained cylinder's
all-length defect tensor is recomputed. Their tail output supply below X is
at most79L²+19L³, `L=floor(log2(X−1))`. The bound concerns actual count height,
not odd parts of arbitrarily large slopes.

## 5. Local p-adic accuracy does not relax the length equation

The [arithmetic bridge](bridge/PROOF.md) uses the inherited delayed-failure
example `δ(m)−1=−4·3!S(m,3)`. It is a true unit seed only at m1,2, yet the
lengths2·3^a approach unit defect to precision a+2 at3. Every prime at least5
also allows arbitrary precision from failed lengths. Meanwhile modulus32
rejects every failed length simultaneously. These different prime views are
exactly compatible: a local near-solution is not a physical unit equality.

## 6. What remains to prove

The session separates free clocks, finitely allowed single clocks, and exact
synchronized directions. These are different slices of length geometry, not
a claimed exhaustive classification of all multivariate loci. The general
hypothesis that individually stable ports are transitive and jointly stable
remains open; the finite contraction system's uniqueness does not prove it.

For Chan–Pak, the synchronized gate and the arbitrary-module determinant are
new constructive ingredients. Neither supplies all required coefficient rays,
their missing heads, nor a reusable operation anchoring arbitrary odd targets.
A next step must change the frozen core or realize more physical determinant
profiles while retaining the continuing interface. Pure p-adic lifting and a
finite bank's free clone inflation do not discharge that obligation.

[Reproduction](REPRODUCE.md) · [Independent review](audit/REVIEW.md) ·
[Connection ledger](CONNECTIONS.md) · [Release scope](RELEASE.md).

## 7. A positive rank jump, and precisely where its next step stalls

The [new clock cores](stable/CLOCK_BIRTHS.md) go beyond the inherited bank:
a strong order-nine core has three independent canonical twin lengths, and
a strong order-ten core has four. Their entire signed tensors are constant1,
so these are all-length families, not extrapolations from clone tests. Three
clocks first occur at nine relative to the inherited complete lower-order
bank; no first-order assertion is made for four clocks.

Their disjoint generating functions are z^9/(1-z)^3 and z^10/(1-z)^4.
Together with the old fixed cylinders this gives exactly
19n−49+binom(n−6,3) marked classes for every n≥8 in the defined union.
The behavior grows from linear to cubic without claiming to exhaust seeds.

The simplest next-step rule fails precisely: with the order-ten exterior
and marks unchanged, all320 transitive insertion/signature choices for a
fifth clock fail necessary unit/stability tests. This closes that immediate
move, not higher rank. Changing an exterior signature, transporting marks,
or a coupled insertion remains a concrete route to an iterable rank-raising
construction. Even these new families have only polylogarithmically many
seed-tail outputs below a target bound; Chan–Pak coverage still needs more.
