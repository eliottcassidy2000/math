# What an unmarked tournament module can tell any surrounding tournament

Let c_r(U) be the number of ordered covers of a nonempty tournament U by r **nonempty** directed paths, with c_r=0 above its order. Write H=c_1 and P_U(k)=Σ_r binom(k,r)c_r(U) for the number of k ordered slots, now allowing empties. Put
```
d_r(U) = (r−1)! c_r(U).
```

**Theorem.** The sequence d(U) is an exact arithmetic interface for unmarked substitution. For any tournaments U,V, the gcd of all differences
`H(S[z←U])−H(S[z←V])`, over every marked tournament host (S,z), is exactly the gcd of the coordinates d_r(U)−d_r(V). In particular modular equivalence in every host is decided by these coordinates, and finitely many physical probes suffice to recover them integrally. The same statement without differencing gives the common divisor of all counts of hosts containing U as an unmarked module.

This is an arithmetic analogue of an interface in a modular-decomposition tree. It is stronger than treating individual counts or orbit sizes as the interface. The integral span in the theorem is **not** an assertion that integer linear combinations of counts can be physically realized.

## 1. The first factorial comes from excursions; the other comes from order of slots

For a fixed marked host (S,z), delete U from a Hamiltonian path of S[z←U] and replace each maximal U-run by z. The resulting word uses outside vertices once, z a number r of times, and satisfies the ordinary arc rule between consecutive symbols. Conversely any legal outside word with r copies of z accepts exactly c_r(U) fillings by its ordered U-runs. Thus
```
H(S[z←U]) = sum_r a_r(S,z) c_r(U).
```
The r−1 interior excursions from z to z are nonempty directed lists of outside vertices. They can be permuted arbitrarily, while preserving the initial and final outside pieces. This is free because outside vertices are distinct. Hence (r−1)! divides a_r. Independently, r! divides c_r because the orders of the r nonempty slots act freely. These two factorials come from different operations and are not inferred from any vertex automorphism of U.

For an exact physical probe at level r, use z and r−1 vertices O_i that it beats and r−1 vertices I_i that beat it. Choose O_i→I_i and no other O→I arcs; orient the remaining pairs in either side transitively. Every interior excursion visits at least one O and one I, so at most r runs are possible. A word using exactly r runs has no exterior pieces and its excursions are exactly the r−1 matched pairs (O_i,I_i) in arbitrary order. Its coefficient a_r is therefore (r−1)!, exactly. The level-one probe is a singleton.

In the d basis the physical probe matrix is integral lower triangular with diagonal ones. Solve successively to recover all d_r, including modulo any integer, without division. For U,V up to order n only probes r≤n are required, with hosts of orders 2r−1 before substitution. The gcd assertion follows both ways. The normalized probe matrix explains why treating only maximum-return configurations as a supply of small shifts incurs a large factorial cost.

**Probes can require the final tournament to be strong.** All r≥2 selector hosts are strong, and inflating a vertex of a strong host of order at least three yields a strong tournament even if U is not strong (paths can leave and re-enter the module). To avoid the singleton probe, add the nearly transitive order-four host in natural order 0,1,2,3 with 3→0 reversed, marking vertex 1. It and the triangle probe each allow at most two visits. Their first two coefficients are (5,2) and (3,1) respectively. In the order-four host a two-visit excursion must end at the unique in-neighbor 0, and before 0 pass through 3; it is either (3,0) with 2 afterward, or (2,3,0). Hence the coefficients have determinant one (triangle row first), and both ordinary H and P(2) can be recovered integrally from strong outputs alone:
```
H(triangle probe[U]) = H + P(2),
H(order-four probe[U]) = H + 2 P(2).
```
Continue with the selectors at r≥3. The same gcd theorem therefore holds even when only strong output probes are allowed.

## 2. Recursive modular arithmetic, and rigidity without graph symmetry

For a modulus M only finitely many r can matter: d_r(U) is divisible by r!(r−1)!, so all sufficiently large coordinates vanish modulo M, uniformly in U's order. Contextual equivalence is a congruence of tournament substitution: if U,V are equivalent in every host, wrap them in any further marked host and use associativity (flatten the surrounding substitutions while keeping its final hole marked). Thus these finite states compose to **all** depths. For a fixed finite palette of substitution skeletons, closure of the reachable states is a finite computation, certifying necessary residue conditions at unrestricted depth. It does not provide a count of a specified size.

For a prime p, the universal p-divisibility condition is exactly
`P_U(k)=0 mod p` for every k in F_p. Indeed for r<p the factor (r−1)! is invertible, and r≥p gives automatic p-divisibility of c_r; the P_U values at k=1,...,p−1 and these coefficients are related by an integral unit triangular matrix. Values at other integer k agree modulo p: write binom(k,r)c_r as the falling factorial times c_r/r!, an integer-coefficient polynomial. This also avoids using a binomial with a zero denominator in F_p.

For a tournament, the [Grinberg–Stanley expansion](https://arxiv.org/abs/2307.05569), evaluated at k variables equal to 1, gives
```
P_U(k) = sum over disjoint nontrivial odd cycle packings C
               2^|C| k^(|U|−sum_{gamma in C} (len(gamma)−1)).
```
The principal specialization counts exactly ordered path slots, by reading the permutation word collapsed into its monochromatic directed runs (with empty colors allowed). In particular P_U(−k)=(−1)^|U|P_U(k). For odd p only (p−1)/2 nonzero residue tests are needed. For **p=3** the test reduces to H(U)=0 mod 3. Hence every module with 3 dividing H forces 3 to divide the Hamiltonian count of **every** surrounding tournament, even for a rigid module having no nonidentity vertex symmetries. Example: the literal strong order-five tournament with outgoing rows (16,1,3,7,14) is rigid, has c=(9,66,186,240,120), and has universal common divisor exactly 3. This is arithmetic divisibility inherited through composition, not a count of automorphisms.

At p=5 the two tests really differ. The nearly transitive order-four module has H=5 and P(2)=32, since its two directed triangles contribute 4k² to k⁴. The r=2 probe is the directed triangle; its output is 3H+c_2=H+P(2)=37. It does not preserve divisibility by 5. Actual vertex automorphisms always divide every c_r (they act freely on ordered full lists), but do not explain the rigid mod-3 example.

## 3. An all-size obstruction to unit tuning with unmarked modules

For every U of order at least two, P_U(2) is a multiple of 4: the empty packing contributes 2^|U|, and any nonempty packing has at least one cycle weight 2 and at least one power of 2 from the specialization. Thus c_2=P_U(2)−2H is 2 modulo 4. In particular if U and V have the same H>1, their H counts after replacement in any unmarked host differ by a multiple of 4. If they also have the same P(2), the difference is a multiple of 12 by the two factorials starting at r=3. An exact adjacent pair produced by a **marked** port change is not subject to this obstruction.

There is also a symmetry law for simultaneous substitutions. Let S have vertices i and substitute U_i of orders at least two. Modulo 4, schedules with at least two repeated module labels contribute zero, as each repeated-label factor c_r (r≥2) is even. For only one repeated label, r=3 contributes zero because (r−1)! divides the outside coefficient; r≥4 is zero already from r!. The remaining terms are
```
H(S[U_i]) = H(S) product_i H(U_i) + 2 sum_i a_2(S,i)  (mod 4).
```
In the corrective term the other H factors are all odd. Thus permuting modules of orders at least two in a fixed quotient preserves the count modulo 4, even if the modules are nonisomorphic and the quotient has no symmetry exchanging their positions. Singletons are a genuine exception: their c_2 is zero instead of 2. Consequently swapping whole nontrivial modules or preserving ordinary H between unmarked modules cannot execute a +2 repair. Root/port information can break this modular constraint.

There is a related exact cost for modules that differ only at one run. If c_r(U)=c_r(V) for all r≥2 and H(U)−H(V)=2, then in every surrounding host the difference is exactly 2H(S): the one-run exterior coefficient is the Hamiltonian count with the module replaced by a singleton. Thus even the physical one-visit-only switches in [COMPOSITION.md](COMPOSITION.md) cannot give a unit change of the odd count inside a strong host, since H(S)≥3 there. Marked switches, or simultaneous changes that compensate across genuinely different host coefficients, are still possible.

## 4. Scope and an exact next test

A round of substitution can be analyzed as a finite congruence state by recording d_r rather than the ordinary count, and integer-linear probe differences can be recovered with no loss. Neither fact establishes exact **positive** realizations of integer-linear combinations. To exploit the probe family toward Chan–Pak one needs a marked or multi-interface compensation operation that implements a controlled difference physically, with a state that survives the next round. The all-size congruence above is a fast way to rule out candidate unmarked swaps before searching for such an operation.

[Replay](verify_observables.py) checks the probe coefficients, the explicit rigid example, and the modular laws through literal independent path/cover enumeration on small hosts.
