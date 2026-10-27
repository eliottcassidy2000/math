# Integral slack, exceptional norms, and a finite catalyst threshold

Let `Q(a,b)=a²−ab+b²` on `Z[ω]`. A binary merge of charges `x,y` has
polarization `β=Q(x+y)−Q(x)−Q(y)` and paid cost `c=max(β,0)`.
Set `h(z)=floor(Q(z)/3)`. Since Eisenstein norms are zero or one modulo
three, write `Q(x)=3r+e`, `Q(y)=3s+f`, `e,f∈{0,1}`. Then

    h(x+y)−h(x)−h(y)=floor((e+f+β)/3)≤max(β,0).

For `β≤0` the numerator is at most two; for `β≥1`, it is at most `β+2≤3β`.
Thus the slack `σ=c+h(x)+h(y)−h(x+y)` is a nonnegative integer.

## A general exceptional-set lifting lemma

For any integral-slack binary operation, `φ=h+1_E` remains a subsolution
if and only if every zero-slack merge with parent in `E` has a child in
`E`. Indeed its new slack is

    σ+1_E(x)+1_E(y)−1_E(x+y).

Only a zero old slack and an exceptional parent can make this negative.
For any binary source-preserving tree, telescoping therefore yields

    total paid cost ≥ φ(total charge)−Σ_original leaves φ(leaf).

The original indivisible leaves matter; one may not replace a leaf by
convenient unit pieces for free.

## A parity-shell exception containing eleven

Call `z` exceptional when `Q(z)>1`, `Q(z)≡1 (mod3)`, `Q(z)+2` is not an
Eisenstein norm, and the norm-`Q(z)` shell in the parity class of `z`
contains only `z,−z`.

For a zero-slack parent with norm residue one, the only residue and
polarization possibilities are

    (e,f,β)=(0,1,0),(1,0,0),(1,1,−1).

For `v=x−y`, we have `Q(v)=Q(x+y)−2β` and `v≡x+y (mod2)`.
The `β=−1` case is excluded by the nonnorm condition. In the `β=0`
cases, parity-shell rigidity forces `v=±(x+y)`, so one child is zero
and the other is the exceptional parent. The lifting lemma applies.

The integer charge eleven is exceptional. Its norm is121; that shell
consists of its six unit multiples, and only `±11` have its parity.
The number123 is not a norm: `123=3·41`, and the inert prime41 occurs
to odd valuation. For completeness, if `q≡2 (mod3)` divides
`a²−ab+b²`, reduction modulo `q` forces both coordinates divisible by
`q` (a nontrivial root would give an element of order three). Consequently
its valuation in a nonzero norm is even.

The same argument gives an infinite family `11^(1+820k)`, `k≥0`.
The congruence `11^820≡−1 (mod41²)` implies its square plus two is
`123 (mod41²)`, with exact41-valuation one. Inertness of11 makes the
norm shell consist of six unit multiples; odd parity picks only the
opposite pair. The checker verifies the elementary congruence and the
finite shell at eleven, not an unbounded search.

For unit leaves summing to11 the theorem gives total paid cost at least41.
A two-unit exceptional bump is false: `x=(−1,−4)`, `y=(12,4)` sum to11
but the base slack is only one. A single indivisible leaf11 has no merge
cost, illustrating why its own potential must be subtracted.

## The catalyst threshold: finite verification, not an infinite search claim

The separate exact finite audit established, for leaf multiplicities
`(11+k,k,k)` of `1,ω,ω²`,

    minimum paid cost =41 for k≥36;
    minimum paid cost ≥42 for k<36.

The latter is a lower bound, not a claim that every smaller case has cost42.
At `k=36` the witness has119 original unit leaves, positive invoice41
and negative invoice39, consistent with `121−119=2`. Adding a
zero-total triple is a source-preserving extension with three extra units
of source norm, not a free rewrite of the original instance. Such a triple
can itself be merged at zero paid cost: first `1+ω`, then `+ω²`; its zero
root can be appended to the existing witness at zero cost. This proves
the extension from the finite witness at36 to every larger `k`.

Here is the scope of that finite proof. Every unit-leaf subtree satisfies
`Q≤3F+1`, where `F` is its paid cost. Any competitor of cost at most41
therefore lies in the finite disk `Q≤124`, containing451 lattice states.
The audit checked the resulting Bellman inequalities and an acyclic
witness, including multiplicity constraints, to exclude all smaller
catalysts. Those certificate tables are intentionally not bundled in
this explanatory package. The included small script verifies the general
potential, shell argument and hostile examples; it does **not** reconstruct
the full catalyst certificate or provide an independent replay of this
finite threshold. This result is reported with that finite-verification
status, rather than promoted to a new proof from the small script alone.

Finally, a sum of positive node costs `F` is not the ROOT invoice `D`
formed by first summing a whole level and only then taking its positive
part. The earlier exact cutoff `N=4096` gives `F=11419`, `D=460`, source
norm9949. Bounds for the first quantity cannot silently be used for the
second, and no all-scale ROOT debt theorem follows here.
