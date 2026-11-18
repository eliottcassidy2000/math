# Four centers cover every odd residue modulo 8·7^k in the strong blue sector

For each k>=1 and odd integer N, the literal cyclic mirror construction
below gives a strong self-converse tournament with H congruent to N
modulo 8·7^k. Unlike the accompanying two-vertex-center theorem, this
allows choosing the mod-8 sheet and *any* mod-7^k residue simultaneously.
Unlike that theorem, no arbitrary induced source is required. This is
a modular coverage theorem, not an exact count or prime theorem.

Take X to be j ordinal copies of D=(2,4,9,3) followed by a transitive
chain of length m+1, and use C3[X,X^op,Z]. The center Z is selected by
the desired mod-8 residue using the center law in PROOF.md:

| H mod 8 | Z |
|---:|:---|
| 1 | transitive triple |
| 3 | out-masks (14,20,24,18,1) |
| 5 | transitive pair |
| 7 | directed triangle |

The mask center is reversed by index reflection; the other centers
are self-converse as well. All resulting hosts have nonempty cyclic
modules, so are strong and self-converse. Their order is
2(4j+m+1)+|Z|.

For fixed center, classify directed Hamilton paths by their cyclic
sequence of nonempty runs in the three modules. For each possible
run length L from 3 through 3|Z|+2 and each possible starting module,
let (a,b,c) be the individual multiplicities. Then the Hamilton total is

    sum Q_a(X) Q_b(X) Q_c(Z),

where impossible runs vanish. This is exact since the run sequence and
ordered covers uniquely determine the path. By inclusion-exclusion on
empty slots, Q_r(X)=sum_i (-1)^(r-i) binom(r,i) P_i(X); only r<=6 enter.
Here P_i(X)=P_i(D)^j i^(m+1), so for fixed j the expression is a polynomial
with integer coefficients in 1^(m+1),...,6^(m+1).

The pairs (j,m0) in crt.py are an explicit table of just 28 seeds:
one for each center and each H mod 7. The verifier checks their
exact residues, and that H(m0+6)−H(m0) is nonzero mod 49. That
condition certifies every higher digit, not just mod 49. Indeed,
for i=1,...,6 write i^6=1+7u_i mod 49. Binomial expansion gives

    i^(6·7^e·d)=1+7^(e+1) d u_i mod 7^(e+2).

Taylor expansion of the integer count polynomial modulo 7^(e+2)
then has exactly the same first-order direction and nonzero slope
mod 7 as at the seed. Every step m -> m+6·7^e·d for d=0,...,6 thus
selects the next base-seven digit uniquely. All these steps preserve
the mod-8 residue by the center law. The final width satisfies
m<6·7^(k−1) and j<=4; this is still a large absolute-height cost.
Run python crt.py for the seed check, independent literal DAG-subset
counts, and sample lifts through 7^4. Further improvement needs an
independent way to set the high part of the count, not more residue
tests.
