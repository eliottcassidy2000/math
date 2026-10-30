# The unifying transport theorem, and the oddness firewall

**PROVED elementary synthesis; classical mechanisms, no priority claim.**
The word "odd" alone transfers no theorem. The meaningful shared structure is a
specified action, its accumulated transport around a closed orbit, and the value
group in which the resulting obstruction must be solved.

## 1. Affine transport around a cycle

Let A be an abelian group, let alpha_i be automorphisms of A, and let b_i in A.
Consider an oriented cycle of length k with equations

    x_(i+1) = alpha_i(x_i) + b_i,   i modulo k.

Compose its affine edge maps in their actual cyclic order. Write the composition
as t -> M(t)+B. Then solutions correspond bijectively to

    (id-M)(t)=B.

If a solution exists, its set is a translate of ker(id-M); every other x_i is
reconstructed uniquely by forward transport. This proves the claim, including
inhomogeneous constants. In particular order cannot generally be replaced by a
multiset of edge labels. If all alpha_i are scalar multiplications by commuting
units a_i, M is multiplication by their product, but B still depends on order.

For a connected graph choose a spanning tree and transport one root value along
it. Each remaining edge imposes one affine fixed-point condition on that root.
All conditions must hold simultaneously. This is a constructive equivalence,
not a claim that one cycle or a determinant alone checks every right-hand side.

Three exact specializations explain this session's bridges:

1. Signs alpha_i=+1 or -1 on integer values: negative total transport imposes
   2t=B; positive transport imposes B=0 and leaves t free. The signed-cycle and
   signed-graph companions classify the full integer and modular cokernels.
2. Multiplication on an odd residue domain, with reflection-odd values: quotient
   by x~-x retains a sign on each orbit. Dropping that sign loses the integral
   obstruction. The odd-function companion gives the exact arithmetic cycle count.
3. Weighted arithmetic transport x_(i+1)=(u*x_i+b_i)/p^(a_i), first over Q:
   clearing denominators produces the obstruction p^A-u^k. The arithmetic
   companion proves the integral Smith form and ordered-word condition. The
   signless odd-cycle equation is the p=1,u=-1 instance; accelerated Collatz is
   p=2,u=3,b_i=1. This shared matrix mechanism does not by itself supply a nontrivial positive
   cycle, exclude all such cycles, or prove convergence of arbitrary trajectories.

For character-covariant signs, transport around a stabilizer loop multiplies the
value by chi(g). A negative stabilizer loop forces a value equal to its negative.
It permits only zero in a torsion-free target and forbids a nonzero +/-1 sign.
This is the obstruction behind the complete composite-modulus classification.

## 2. Odd Boolean weight is a different condition

For n>=2 let f:F_2^n -> F_2 have odd truth-table weight. Then for every h!=0,
its Boolean derivative D_h f(x)=f(x+h)+f(x) is **not constant**.

Proof: translation by h partitions the cube into 2^(n-1) pairs. If D_h f=0,
each pair contributes zero or two ones, making the total weight even. If
D_h f=1, each pair contributes one, giving weight 2^(n-1), again even. Both
contradict odd weight. The n=1 exception is real: f(x)=x has odd weight and
constant nonzero derivative.

Moreover no permutation pi of the input set can satisfy f(pi x)=1+f(x): weight
preservation would force wt(f)=2^(n-1). Thus odd-weight Boolean phases have no
input-permutation symmetry that negates their sign phase. This is not the same
as a reflection-odd real function. On F_2^n the literal map x->-x is the identity;
its interpretation must not be silently replaced by coordinate complementation.

This strong absence of constant derivatives is not circuit hardness. The inherited
Session 117 point-defect family has odd weight for n>=3 and exact full-B2 size
n-1 or n+k-2. In particular it can have this property with only n-1 gates.
The exact two input-affine orbits in that earlier work must not be merged by
pretending output complementation is an input symmetry.

## 3. Recovering the four odd eighths through actual characters

On U(8)={1,3,5,7}, define

    chi4(a)=(-1)^((a-1)/2),
    chi8(a)=(-1)^((a*a-1)/8).

Together with 1 and chi4*chi8 these are the four real multiplicative characters.
Their four-by-four value matrix is orthogonal: summing chi(a)psi(a) over a in
U(8) gives four when chi=psi and zero otherwise. Hence every function on the four
odd residues has a unique character expansion

    f(a) = sum_chi [ (1/4) sum_r f(r)chi(r) ] chi(a)

over any coefficient field where four is invertible. This coefficient-ring
condition matters; the formula does not transfer unchanged to characteristic two.
The exact residue projector is (1/4)sum_chi chi(a)chi(r).

Reflection a->-a gives two even character directions, 1 and chi8, and two odd
ones, chi4 and chi4*chi8. The four complex eighth phases satisfy

    exp(2*pi*i*a/8) = [chi8(a)+i*chi4(a)*chi8(a)]/sqrt(2).

Thus the familiar fractions 1/8,3/8,5/8,7/8 have an explicit arithmetic-character
connection to odd functions, with both real/even and imaginary/odd parts retained.
This is a classical Fourier identity, not evidence that all four uses of "odd"
are interchangeable or that a new analytic number-theory bound follows.

## 4. A decisive eighths-to-odd-cycles hostile

An odd sign on the unit differences of Z/8 defines an orientation of K_(4,4),
not a tournament on eight vertices: every retained difference is odd, so every
edge joins opposite parity classes. There are no undirected or directed odd
cycles. More generally unit differences on Z/(2^a), a>=2, give the complete
bipartite graph K_(2^(a-1),2^(a-1)). At a=1 the sole unit is fixed by negation,
so no nonzero odd sign exists.

Consequently odd numbers used as differences, an odd sign function on those
numbers, and odd-cycle existence are independent requirements. The full
circulant-tournament dictionary requires every nonzero difference of an odd-order
group; restricting to units of an even-order group destroys that target.

## 5. Exact checks

`check_transport_eighths.py` tests the character table, exact projectors and phase
coordinates; all odd-weight Boolean truth tables through n=4 against every
nonzero translation; and the bipartite unit-difference graphs through2^6.
It records explicit boundary counterexamples. These finite checks corroborate the
short uniform proofs; they do not replace the separately stated maps and domains.
