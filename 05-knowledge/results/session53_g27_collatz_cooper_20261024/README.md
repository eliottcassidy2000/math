# Session53: three meanings of 27, and the data each needs

24 October 2026. First independently audited checkpoint. The anchor is the
Cooper s7 square-root congruence; the graph-growth response is the niche;
the Collatz/Paley character comparison is the wildcard. Coinciding indices
are not treated as a map between problems.

## 1. The prime dichotomy is recovered, not newly discovered

Normalize F(t)=sum s7(n)t^n=1+4t+48t^2+760t^3+..., and G(t)=sqrt(F(t)) with
constant term1. Its coefficients start1,2,22,336,6006. The exact statement
requested by the user is already Theorem D of the October18 repository
note: for every n>=0, c(pn)=c(n) modulo p^2 exactly for p=7 or p=1,2,4 modulo7.
For each other prime there exists a failure already modulo p. The proof
does not assert that n=1 is the failing index at every inert prime.

The [recovery packet](cooper-prime-dichotomy/REPORT.md) pins the normalization,
literature and exceptional primes. The [independent audit](audit-inherited-s7/REPORT.md)
accepts the inherited proof under its explicitly stated classical modular
inputs, reconstructs the exceptional-prime certificates, and independently
tests all primes up to500 with pn<=20000. It does not certify unrelated CM
claims, all the old GP transcript, higher-power conjectures or novelty.
The older script is tracked upstream; its initial absence was only sparse
materialization, not missing repository evidence.

Cooper2024 explicitly lists the conjectural question. We found no later
primary proof in this search, but do not claim exhaustive literature coverage
or that this repository proof has received external peer review.

## 2. A literal role for 27 in the coefficient operator

In characteristic3 the known Lucas identity for F gives F=(1+t)F(t^3).
Thus G^4=(1+t)^(-1), and the unique normalized fourth root is

    G(t) = product_(j>=0) (1+t^(9^j))^2  (mod 3).

The proof is a telescoping product, using (1+x)^9=1+x^9. Consequently c(n)
is nonzero modulo3 exactly when every base9 digit of n is0,1 or2. On this
support, c(n)=(-1)^(number of digits equal to1) modulo3.

Put s^2=(1+t)(1-27t), with constant term1, and H=G/s. The Cartier operator
C3 extracts every third coefficient. Modulo3 it swaps G and H:

    C3(G)=H, C3(H)=G; C3^2(G)=G, C3^3(G)=H.

Therefore c(9n)=c(n), while c(27n)=c(3n)=[t^n]H, modulo3 for all n.
G is not a negative eigenvector: already c(3)=0 modulo3 differs from -c(1).
This is an exact two-state arithmetic mechanism, not a statement about
Collatz convergence. Both the recovery and independent audit include it.

## 3. G27 creates compatible joins, not the first Hamiltonian path

The [graph theorem](g27-growth/THEOREM.md) retains the actual new neighbors:
9 and22. Among the12 Hamiltonian paths of G26, counted up to reversal,
exactly3 extend to an endpoint at27. The other32 Hamiltonian paths of G27
join two different old paths through9--27--22. G26 has1116 spanning two-path
covers, and exactly32 expose those two neighbors on distinct components.
Each compatible deletion fibre has size1. Since9+22=31 is not a square,
these joins are not subdivisions of an existing edge.

The adjacent births25,26,27,28 split respectively as0+10,2+10,3+32,5+47.
The old deletion correspondence and known35 count are explicitly inherited;
the response split and full endpoint-state law retain the missing detail.

For any graph H and new vertex with neighbor set A, define L_H(t) as the
generating polynomial of spanning linear forests by component count. For
each old forest F with k components, let e_A(F) count exposed vertices in A
(an isolated vertex counts once), and j_A(F) count unordered pairs of such
vertices in different components. Then

    L_(H+v)(t)=sum_F [t^(k+1)+e_A(F)t^k+j_A(F)t^(k-1)].

Exact endpoint-pair states suffice to transfer this response. Merely knowing
the minimum-cover number does not: a near-minimum old cover can join into a
new optimum even when the old optimum cannot extend.

## 4. Collatz27: an infinite first-descent family, not global coverage

Use the shortcut map n/2 when n is even, (3n+1)/2 when n is odd. The
[Collatz packet](collatz27-return/REPORT.md) proves that every integer

    n=27+2^59 t, t>=0,

first falls below its own starting value at shortcut step59, ordinary
Collatz step96. Its endpoint is23+3^37 t. The proof checks every earlier
prefix multiplier and keeps the actual clock, not just a formal affine map.
Sharing27's full70-step shortcut prefix does not mean that every lift reaches
1 at step70. No universal convergence result follows.

Every finite shortcut parity word occurs in every residue class modulo7,
by its unique dyadic cylinder and CRT. Thus quadratic character modulo7
cannot exclude or determine any finite parity prefix. On the other hand,
a *fixed* word with s odd steps and k total shortcut steps acts modulo7 by
an affine map with slope3^s/2^k. Its quadratic character is(-1)^s, so it
preserves or reverses Paley7 orientation. The59-step and70-step words of27
induce1-x and-x, respectively. This requires retaining the chosen word.
G27 is a square-sum graph, not Paley7; no direct identification is claimed.

The [joint independent audit](audit-birth-return/REPORT.md) checks both graph
and Collatz proofs, relocatable normal/optimized executions, and hostile
controls. It retains the full59-row clock/affine certificate. The graph
provenance documents are hash-referenced; the proofs and executable runtime
are self-contained.

## Boundaries

The matching number27 points to different objects: a vertex birth, a long
actual trajectory, and three iterations of a coefficient-extraction
operator. The productive connections preserve those operations and their
state. This checkpoint proves no equivalence of the three problems, no
all-size tournament wiring theorem, and no stronger uniform prime-power
supercongruence beyond its stated moduli.
