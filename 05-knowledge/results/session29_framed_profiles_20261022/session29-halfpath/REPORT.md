# Mirror half-paths, ternary state compression, and the centralizer quotient

Status: exact proofs and independent finite controls, Session29,2026-10-22. This packet changes no repository files. Its algorithms compute symmetry-refined counts, not the complete Hamiltonian spectrum or A000568.

## 1. Setup and prior work

Let T be a self-converse tournament on n vertices and let sigma be an involutory anti-automorphism: u->v iff sigma(v)->sigma(u). Put G=Aut(T) and C=C_G(sigma). Every automorphism group of a tournament has odd order. Otherwise Cauchy's theorem gives an involution, whose transposed vertex pair would have its edge reversed by an automorphism. For n>=2, the extended automorphism/anti-automorphism group has order2|G|; alternatively, taking the appropriate odd power of any anti-automorphism yields an anti-involution. The n=1 case uses the identity directly. An anti-involution has at most one fixed vertex, since an arc between two fixed vertices would reverse. Thus n=2m has m exchanged pairs and n=2m+1 has m pairs plus one fixed vertex v.

All anti-involutions are conjugate under G. An elementary proof avoids an assumption about the automorphism group being abelian: if a,b are two anti-involutions, g=ab has odd order d and g^((d-1)/2) a g^(-(d-1)/2)=b. There are therefore |G|/|C| anti-involutions.

The repository already contains the half-path principle, explicitly for Paley examples, in `07-reflections/the-reversal-involution-and-the-half-principle-opus-20260629.md`. It reports the fixed Hamiltonian-path counts1,9,185 at p=3,7,11. The present packet does not claim the underlying half-path bijection is new. It provides a uniform signed-pair dynamic program, handles both parities and fixed odd cycles, and proves the automorphism-centralizer conversion needed for actual blue tile counts. Only the relevant tournament statements from that historical reflection are used; its broader cross-topic analogies are not adopted as theorems.

## 2. Exact half-path bijections

Write each exchanged pair as {a_i,b_i}, with sigma(a_i)=b_i.

A Hamiltonian path fixed by P -> sigma(reverse P) is uniquely determined by a first half x_1,...,x_m, with exactly one vertex chosen from each pair.

- Odd order: P=(x_1,...,x_m,v,sigma(x_m),...,sigma(x_1)). It is directed iff x_1->...->x_m->v.
- Even order: P=(x_1,...,x_m,sigma(x_m),...,sigma(x_1)). It is directed iff x_1->...->x_m and x_m->sigma(x_m).

The remaining arcs follow from anti-symmetry. These are bijections of directed paths, not merely equal cardinalities. For n=1 the empty half gives the unique path.

A directed odd cycle fixed by reversal followed by sigma must have a sigma-invariant odd vertex set, hence contain the unique fixed vertex v. Thus there are none at even order. At odd order it has the unique form, up to its conventional cyclic rotation,

(v,x_1,...,x_j,sigma(x_j),...,sigma(x_1)), 1<=j<=m,

where the x_i choose at most one vertex from each pair. The necessary and sufficient half conditions are v->x_1->...->x_j and x_j->sigma(x_j). Starting the cycle at v makes this encoding unique. Its length is2j+1. Fixed cycles and fixed full paths are different statistics.

## 3. Ternary dynamic program

A state q in {0,+,-}^m says whether pair i is unused or which of its vertices has been selected. Let S(q) be the selected vertices. For x in S(q), define

D(q,x) = number of directed paths using exactly S(q), ending at x.

For singleton q selecting x, D(q,x)=1. For larger states,

D(q,x)=sum_(y in S(q)-{x}, y->x) D(q with x's pair unused,y).

Equivalently extend every state by either sign of an unused pair when the next arc exists. The fixed full-path count f is the sum over full states (no zero digits), using the odd/even terminal tests from Section2.

For fixed odd cycles maintain a second array D_v with singleton initialization D_v(q,x)=1 if v->x, and0 otherwise, using the same transitions. Sum D_v(q,x) over all nonempty q and endpoints satisfying x->sigma(x). Grouping by |S(q)| records the full fixed-cycle length distribution.

There are exactly2m3^(m-1) possible nonempty (q,x) pairs. Each has at most m-1 predecessor endpoints, so both counts use O(m^2 3^m) arithmetic operations and O(m3^m) integers of storage. The straightforward supplied implementation allocates two rectangular3^m by2m arrays and skips invalid endpoints. It is not optimized to the exact valid-state count. Counts have O(m log m) bits, so the arithmetic-operation bound is not a unit-cost bit-complexity claim. Layered storage can discard states after their cardinality layer has been used, but that optimization is not implemented here.

This improves the exponential base from the full vertex-subset path algorithm's2^(2m) to3^m for these symmetry-fixed observables. It does not compute H(T), Aut(T), or a canonical labeling for free; sigma must be supplied or found separately.

## 4. Blue spine tiles require the centralizer quotient

The automorphism group G acts freely on full Hamiltonian paths: an automorphism fixing an ordered path fixes every vertex. Its orbits are the distinct fixed-spine tournament tilings belonging to T, so there are H(T)/|G| of them.

A tile orbit is mirror-fixed (blue) precisely when some anti-involution tau fixes its representative path under reversal. For each ordered Hamiltonian path, at most one such tau exists, since tau must send its i-th vertex to its (n+1-i)-th vertex. Consequently each blue tile orbit contributes exactly |G| pairs (P,tau).

On the other hand, the |G|/|C| anti-involutions are conjugate, and each fixes f paths. Double counting gives the exact formula

b(T)=f/|C_G(sigma)|.

In particular the centralizer divides f. The path involution has f congruent to H(T) mod2; by Redei parity, f is odd. Since the centralizer is odd, b(T) is a positive odd integer. The number of non-blue tile incidences is H(T)/|G|-b(T), an even integer.

Pure-blue status is equivalent to H(T)|C|=f|G|. No characterization of primitive exceptional tournaments is being imported into this statement.

## 5. A counting bound and an exact ensemble mean

For m>=1 there are2^m m! signed half-orderings. Once the last pair is specified, exactly one of its signs satisfies the terminal fold condition: x->v in odd order or x->sigma(x) in even order. Hence

f <= 2^(m-1)m!.

For m>=2 the right side is even while f is odd, improving it by1. These bounds ignore the internal consecutive-arc tests and need not be sharp. They imply a simple sufficient obstruction to pure-blue status:

H(T) > (|G|/|C|) 2^(m-1)m!  =>  T is not pure blue.

There is also an exact mean, under a specified measure. Fix sigma and choose uniformly among all labelled tournaments admitting it. There are2^(m^2) choices at n=2m and2^(m(m+1)) at n=2m+1: these are the independent unordered-edge orbits under sigma. Every signed half-ordering imposes m independent orientation bits. The m-1 internal half-path edges lie in distinct edge orbits because only one vertex from each pair has been chosen; the terminal fold edge belongs to a further distinct orbit. Therefore

E_sigma[f]=2^m m! * 2^(-m)=m!.

This is not the mean over unlabelled self-converse classes, nor the mean blue tile count. The centralizer correction cannot be dropped.

## 5b. Signed-pair isomorphism and two global checksums

The centralizer of sigma in S_n is the signed-pair permutation group of order2^m m!. Two sigma-invariant tournaments are ordinarily isomorphic iff they are isomorphic by an element of this group. Indeed, an arbitrary isomorphism transports sigma to an anti-involution of the second tournament; an automorphism of that tournament conjugates this anti-involution back to sigma. Their composite then commutes with sigma. Thus fixed-sigma tournaments modulo signed-pair permutations enumerate self-converse isomorphism classes, with no duplication from different anti-involutions. This is a symmetry reduction, not a claim of a new general tournament-enumeration algorithm.

Let o=m^2 at even order and o=m(m+1) at odd order. A self-converse class T contributes exactly2^m m!/|C_G(sigma)| distinct labelled tournaments in the fixed-sigma ensemble. Consequently, summing once over self-converse classes,

sum_T 1/|C_G(sigma)| = 2^(o-m)/m!,

sum_T f(T)/|C_G(sigma)| = 2^(o-m).

The second is also the global blue-tiling count, because f/|C|=b(T). It equals2^(m(m-1)) for n=2m and2^(m^2) for n=2m+1. These weighted sums are exact consistency checks for a class census; neither is the unweighted number of classes. They agree with counting free half-tiling bits directly.

## 6. Exact examples and an important distinction

The independently verified tuples (H,|G|,|C|,f,b,total tiles) are:

- C3: (3,3,1,1,1,1).
- C5{1,2}: (15,5,1,3,3,3).
- Paley P7: (189,21,3,9,3,9).
- Session29 seed A: (15,3,1,1,1,5).
- Session29 seed B: (13,1,1,3,3,13).

Thus P7 is mixed, with3 blue and6 black tile incidences. Its fixed full-path count is9, not27, and it must be divided by the centralizer3. Even the coarse bound proves it cannot be pure blue:189>7*24=168. A different primitive exceptional status of P7 cannot be identified with this metagraph property without a separate theorem.

For P7 the fixed odd-cycle count is12: six fixed5-cycles and six fixed7-cycles, with no fixed triangles. This is not its fixed full-path count9. C3 has1 fixed odd cycle; C5{1,2} has3; seed A has3; seed B has2.

## 7. Verification and limits

`check.py` uses only Python's standard library and explicit exceptions, not removable assertions. It independently enumerates all vertex permutations for Hamiltonian paths and automorphisms; constructs distinct normalized fixed-spine tiles; tests their mirror symmetry directly; enumerates directed odd cycles up to rotation; and compares these with the ternary DP and centralizer formula.

All599 tournaments invariant under the standard reversal involution at orders1 through6 are checked. The ensemble mean identity is checked over each full labelled ensemble. Canonical normalized-spine representatives recover1,1,2,2,8,12 self-converse classes at orders1 through6, and both global centralizer-weighted checksums are independently verified on those class sets. The named examples above and two additional examples each at orders7 and8 are also checked. Every bound, parity, quotient and count passes. Normal and optimized-mode outputs are identical.

These finite controls support the implementation, while the bijections, group action, and recurrence prove the general statements. The construction does not settle the general Hamiltonian-count spectrum, the {7,21} exclusion proof, complete metagraph wiring, or the identification of every primitive exceptional class. Its direct use is faster computation of blue incidences and fixed-cycle congruence data once a converse symmetry is known.
