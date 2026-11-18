# Folded exponent boxes, triangle refinement, and the Goldbach boundary

## Provenance

This note synthesizes, rather than reclaims, the exact divisor result in
`math/05-knowledge/results/arithmetic_braids_20260917_divisors.md` (DB1),
the folded-cube comparison THM-4733, `perpendicular_refinement/DIVISOR.md`,
and the warning about order/count/slot readings in `pair_birth/root/PROOF.md`.
The new contribution is an explicit conditional compatibility criterion,
not a construction proving Chan–Pak or Goldbach.

## 1. What F=U+S really measures

For N>=2, F counts proper divisors other than 1, U counts prime proper
divisors, and S counts squarefree proper divisors other than 1.
For N=product p_i^a_i with r distinct primes,
F=product(a_i+1)-2, S=2^r-1-[N squarefree], U=r-[N prime].
Consequently F=U+S iff N=p, p^3, or p^2qr (three distinct primes).
Proof: primes give 0=0+0; squarefree composites fail. Otherwise
product(a_i+1)=2^r+r+1. For r=1 this gives exponent 3; for r=2
it asks two integers >=2 to multiply to 7; for r=3 only 3*2*2
works. For r>=4, product>=3*2^(r-1)>2^r+r+1.
Thus 10=3+7 at p²qr is a genuine dimension cutoff, not a
Pythagorean identity or an additive prime theorem.

## 2. The precise 3-to-5 cross and what it forgets

On exponent triples (a,b,c) in {0,1,2}x{0,1}² with the two
extremes deleted, let Q(a,b,c)=(a,c,b) and
C(a,b,c)=(2-a,1-b,1-c). The coarse coordinate x=b-c
has three values with fibre sizes 3,4,3. Refine it by keeping
(x,2|a-1|-1) when x!=0 and one central symbol when x=0.
This has five values: a fixed center and two exchanged pairs.
Both C and Q act as the SAME horizontal reflection. Their
product is therefore invisible. Keeping a-1 rather than its
absolute value gives seven values and distinguishes C from Q.

Proof is substitution. The five-state quotient merges a=0 with a=2
on each wing. C reflects a about 1 while Q leaves a unchanged,
so their difference is precisely the information discarded.
This is an exact analogy to a continuation profile that preserves an
answer but forgets which future move is possible. It does not preserve H.

The independent tournament example: append a fourth vertex to a directed
triangle. Of its eight contact patterns, the two uniform patterns give
H=3; all six mixed patterns give H=5, and form the unique strong order-four
isomorphism class. Thus changing one contact from uniform to mixed adds
exactly two paths. Proof: the uniform vertex must be first or last,
leaving the three cyclic paths. The six other contact patterns are equivalent
under triangle rotations and converse; directly enumerate one. This 3-to-5
counts paths, whereas the exponent cross counts observer states. They
share a two-sided refinement pattern but are not the same construction.

The 3-4-5 triangle obeys 3²+4²=5². Here (order,H)=(3,3) then (4,5).
The square identity supplies neither endpoint counts nor a physical gluing.
THM-4733 does give a genuine equivariant map from four mirror-moved bits:
(u,v,b,c) -> (u+v,b,c). Its missing mirror-fixed tiles obstruct extending
that comparison to the full marked tournament cube. It transports the
Klein action, not Hamilton counts or isomorphism classes.

## 3. Exact conditional Goldbach compatibility criterion

Let B be any set of odd primes for which a specified tournament witness
and its interface data have been supplied. For each unordered pair p<=q
in B, let K(p,q)=1 exactly when a CERTIFIED physical operation on the
chosen witnesses produces H=p+q-1. Interface data, not just p,q, belong
in this definition. Define

R_B(E)=#{p<=q in B:p+q=E},
D_B(E)=#{p<=q in B:p+q=E and K(p,q)=0}.

**Theorem (compatible Goldbach cover).** Every odd N with
R_B(N+1)>D_B(N+1) is a tournament count. Thus eventual strict
inequality plus separately supplied small witnesses proves cofiniteness;
strict inequality for every admissible N supplies full coverage.

Proof: subtracting the bad-pair count leaves exactly the number of
certified constructions of N. Positivity selects at least one.
This is a counting certificate for physical compatibility, not a claim
that numerical Goldbach representations are sufficient.

If B omits at most b odd primes, and every even E under consideration
has more than b+D_B(E) unordered odd-prime representations, the same
conclusion follows: each excluded prime removes at most one pair.
This yields a quantitative route: bound the incompatible interface pairs,
rather than insist that every pair glue. A target-dependent choice can
switch BOTH primes and BOTH interfaces.

**Necessary obstruction.** A universal p+q-1 rule on all realizable prime
counts cannot exist: p=q=11 would give the forbidden count 21. Excluding
only 7 as a source does not fix this, since 11 is allowed. A local rule
which accepts the same numerical inputs cannot cure the contradiction;
it must reject every 11/11 interface, change the output formula, or operate
in a different class of graphs. Even assuming ordinary binary Goldbach
does not prove the displayed strict inequality: its one promised pair
might be excluded or incompatible.

No new unconditional prime-placement result is claimed. The useful
research target is a density bound on failed *physical pair interfaces*
along an absolute sum diagonal, with a source supply theorem.
