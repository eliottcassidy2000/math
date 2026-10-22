# Blue tilings as an endpoint gate: completion boundary flux

Status: PROVED transport identity; FINITE-EXACT complete order8 control.

Assume n=2m≥2; the empty object has no half-path endpoint and is excluded.

Let T be a self-converse almost-regular tournament of even order2m, with
score defect d=S1 in{−1,+1}. Let R be its unique regular completion, and
choose an anti-involutionσ ofT. It has m swapped pairs and no fixed vertices.
Extendσ toR by fixing its new vertex*. Write c for its centralizer order
inAut(T).

A signed half-path P chooses one vertex from each swapped pair, in an order
whose consecutive arcs are directed forward inT. Let h_v count such half-paths
ending atv. The same half-path state space has two distinct endpoint gates:

- T gate: append the reversedσ-image to obtain a symmetric Hamilton path inT.
  This is possible exactly when S_(v,σv)=+1.
- R gate: insert* before the reversedσ-image. This is possible exactly when
  S_(v,*)=+1, or equivalently d_v=−1.

Therefore the blue-orbit counts satisfy

    b(T) = (1/c) Σ_v h_v (1+S_(v,σv))/2,
    b(R) = (1/c) Σ_v h_v (1−d_v)/2.

In particular define the completion boundary flux

    F(T)=b(T)−b(R)=(1/(2c)) Σ_v h_v [S_(v,σv)+d_v].

This is an exact endpoint-weighted transport of blue count, not a spectral
approximation or an equality of total Hamiltonian counts.

Why the denominator is shared: Aut(T) extends to the stabilizer of* inAut(R).
Any automorphism commuting with the extendedσ fixes its unique fixed vertex*.
Thus the two centralizers coincide. In a self-converse tournament all
anti-involutions are conjugate under its odd-order automorphism group. Counting
symmetric Hamilton words by their induced anti-involution and dividing the free
Aut action yields b=f_σ/c on both sides. The argument includes the orbit
multiplicity that would be lost by simply comparing labelled fixed paths.

## Tests and a rejected shortcut

For the two order6 isospectral witnesses, T hasblue9 and its completion also
hasblue9. Their completions nevertheless have H175 and171. Those examples
might tempt a claim that regular completion always preservesblue count.

The complete order8 test rejects that claim: among the15 self-converse classes,
seven havezero flux andeight have nonzero flux. Literal inputs, anti-involutions,
centralizers, endpoint vectors and gate counts are inblue-flux-census.json.
The formula, rather than an invariant-blue shortcut, is the surviving statement.

This explicitly connects the signed half-path methods fromSession29 with the
score/endpoint corrections ofSession30. The preserved object is the full
half-path transfer state; the changed object is its terminal acceptance vector.
Neither score variance nor spectrum records this endpoint vector by itself.
No statement about all black-edge wiring follows fromblue flux alone.

Reproduce: python3 check_completion_blue_gates.py and
python3 check_blue_flux_census.py, with the sibling order8 census bank present.
