# Exact first gap for a factored tournament

Follow-up to the [triangle budget](CORE.md). The obstruction gets
much stronger when the old witness factors. It is possible to compute
the first positive induced gap from the individual strong components,
with no joint search across components.

Write T = C₁ → ... → Cₖ for the ordinal sum of its strongly connected
components, hᵢ = H(Cᵢ), H = ∏hᵢ. For each nonsingleton Cᵢ let gᵢ be the
least positive difference H(Cᵢ+x)−hᵢ over its one-vertex extensions.
This exists: let x beat only one old vertex u and pick an outneighbor
of u inside the strong component to make a triangle through x.
Ignore singleton components in the inner minimum.

**Exact induced-gap formula.** The least positive value of H(U)−H(T)
among all U containing T as an induced subtournament, with any number
of extra vertices, is

$$
\min\left(2H,\ \min_{i:\,|C_i|>1}\frac{H}{h_i}g_i\right).
$$

**Lower bound.** Consider any single insertion x which strictly
increases H. For each component Cᵢ, take an odd-cycle packing using x
whose cycle through x is contained in Cᵢ+x; take arbitrary odd-cycle
packs in every other old component. The families for distinct i are
disjoint because the cycle through x is unique. Using the positive
weighted cycle-pack expansion gives

$$
H(T+x)-H(T)\ \ge\ H\sum_i \frac{H(C_i+x)-h_i}{h_i}.
$$

If x has nonuniform contacts with a nonsingleton Cᵢ, then some arc from
its outneighbors to its inneighbors lies in Cᵢ, by strong connectivity;
it makes a triangle with x. The indicated term is at least gᵢ/hᵢ. If
all contacts are uniform, collapse each Cᵢ to one vertex. The quotient
is a tournament. Each of its Hamiltonian paths lifts in H ways by
traversing each old component as a single block. If the quotient is
nontransitive, its count is at least three by the cycle-pack formula.
Otherwise the whole extension is an ordinal sum and H is unchanged.

For an arbitrary induced extension, inspect its first strict
increase. All previous neutral insertions added singleton components
without changing any old strong component: by the same expansion a
neutral vertex is in no directed triangle, and in a tournament every
vertex of a nonsingleton strong component lies in a directed triangle.
The bound therefore still applies.

**Attainment.** Add a new directed triangle as a last ordinal block
to attain 2H. For any i choose a minimizing one-vertex extension of Cᵢ
and make all contacts from the new vertex to other components agree
with the ordinal order of Cᵢ. The new count is exactly H+(H/hᵢ)gᵢ.
In particular minimum positive gap is an operation on one component
or the addition of a triangle, not a cooperative effect between distant
components. This does *not* give the later points of the induced spectrum.

## Consequence for Chan–Pak

If M=max hᵢ≥3, then with the triangle budget,

$$
H(U)-H\ \ge\ \frac{2H\log_3 M}{M}\qquad(H(U)>H).
$$

Indeed gᵢ≥2⌈log₃ hᵢ⌉, and log₃(t)/t decreases for t≥3.
A family with H→∞ and an induced repair of cost O(log H) must
have M=Ω(H): when M<√H the lower bound already exceeds
√H·log₃ H, and otherwise log₃ M≥(log₃ H)/2. In particular even
logarithmic induced repairs are impossible when every component
carries a vanishing fraction of H. If T is an ordinal sum of
k≥1 directed triangles and any singletons, the exact least positive
gain is 2·3^(k−1), since gᵢ=2 (triangle extended to the strong
four-vertex tournament). These are precisely the witnesses of
H=3^τ from the triangle code, yet their attainable first increments
are exponential in τ, not linear.

The proof suggests measuring how a candidate repair changes the
dominant strong component, not just improving the precision of an
arithmetic offset across a bank of small frozen factors. It does not
constrain the one-arc rewiring pairs in [BALANCE.md](BALANCE.md);
those deliberately do not keep the old graph induced.
