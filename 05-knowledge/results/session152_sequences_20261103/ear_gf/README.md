# Transfer for ear-response enumerators

Status: an exact finite generating-function identity, not a proof of a global spectrum theorem. This repackages the integral directed-cut carrier in [`cut_induction`](../cut_induction/README.md) and the self-converse reciprocity in [`recursions`](../recursions/RESULT.md). The quadratic response, explicit order-five polynomial, and the insufficiency of H alone were already recorded in THM-4118 and the session's `ear/README.md`; no literature priority is claimed.

Let c be any integral carrier for a parent tournament T of order n: with x -> i iff i in S, F(S)=H+2 sum_{i in S,j not in S} c_ij. For arbitrary weights c_ij and vertex fugacities t_i, define

  Z_C(t_1,...,t_n;q) = sum_{S subset [n]} (product_{i in S} t_i) q^{sum_{i in S,j not in S} c_ij}.

For m=0 define Z_empty=1. Eliminating vertex m from an m-vertex matrix gives the exact transfer

  Z_C(t;q) = Z_{C-m}((t_i q^{c_im})_{i<m};q)
             + t_m q^{sum_{j<m} c_mj} Z_{C-m}((t_i q^{-c_mi})_{i<m};q).    (1)

Proof: split the subsets according to whether m is absent or present. If absent, every selected i contributes c_im across the cut; if present, m contributes sum_j c_mj minus the terms c_mi for selected i. This is an identity in Laurent polynomials, even though negative powers may arise during the computation; for nonnegative c the final cut enumerator is an ordinary polynomial. In any carrier, the response enumerator is

  E_T(u,z) = z^H Z_C(u,...,u; z²).

The proper, nonempty ear-image is the support of E_T(1,z)-2z^H, with coefficients counting labeled subsets. For strong parents nontrivial cuts are strictly positive, so the subtraction removes exactly the two constant cuts. If C differs among feasible max-flow allocations, Z_C is unchanged: all its cut exponents are fixed by F(S). If T is self-converse, E_T(u,z)=u^n E_T(u^{-1},z), independently of whether C itself has a manifest matching symmetry. At odd order, each response coefficient is even by the anti-automorphism argument in `recursions/RESULT.md`.

The identity is a valid induction in the size of a fixed weighted carrier, not a closed induction on tournament order. It can branch exponentially and requires retaining the matrix and shifted fugacities. Removing a carrier vertex does not generally give a carrier for the subtournament, and H alone is insufficient: the two H=9 parents code 8 and code 759 already have different images. Thus (1) gives neither a solid interval nor the global coverage conjecture. The observed holes in code 8 provide a direct counterexample to inferring interval coverage from a single parent.

## Exact example and comparison

For the code-8 carrier in `cut_induction/README.md`, (1) gives

 E_T(u,z)=(1+u^5)z^9+(u+u^4)(z^15+z^17+2z^23+z^33)
  +(u²+u³)(z^17+z^19+z^23+z^25+z^27+z^29+z^33+2z^37+z^41).

The proper image is A={15,17,19,23,25,27,29,33,37,41}. For the separately supplied B={2,3,5,7,11,13,17,19,23,27}, termwise support intersection of its indicator polynomial with the ear polynomial is

  4z^17+2z^19+6z^23+2z^27.

This counts 14 labeled cuts, of which 12 have prime counts 17,19,23 and two have count 27=3³, which is composite. This is a set/filter comparison only; neither the transfer nor this shared support provides a rule generating B or a causal or number-theoretic relation between the two sets. The distinct additive-semigroup consequences are in `number_theory/NOTE.md`.

## Verification

Run `python3 check.py`. It recursively computes the full transfer with integer field shifts, compares every bivariate coefficient with independent enumeration of all 32 cuts of the published carrier, checks the claimed support and multiplicities, and prints the comparison with B. The carrier's separate script verifies its cuts against literal child Hamiltonian-path dynamic programming; the transfer check does not independently reverify that tournament calculation. See `out.txt` for the transcript.
