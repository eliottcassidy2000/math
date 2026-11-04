# A finite decision criterion for one homogeneous exponential clock

2026-11-03. **PROVED corollary** of the inherited positive clone-run formula
and the cited strong-order bound, independently reviewed in this session.
This is a decision criterion for a specified operation class, not for the
global Hamiltonian spectrum or an unrestricted multi-port grammar.

Fix integers d>=2, a>0 and c, with h=ad+c>0. For a fixed marked nontrivial
strong tournament (R,x), write A_r for its clone-run coefficients as in
[native/clone_escape/PROOF.md](native/clone_escape/PROOF.md). Then

    H(R[x := TT_t]) = a*d^t+c for every sufficiently large integer t

holds if and only if

    A1=h,
    A_r=a*binom(d,r) for 2<=r<=d,
    A_r=0 for r>d.                                      (1)

If it holds eventually, it holds for every integer t>=1. In particular
there is a terminating finite decision procedure for whether any such
marked strong base exists for the given (a,d,c).

## Proof and exact scope

For t>=1, classify functions from a t-element set to a d-element set by
the number r of used values:

    d^t = sum_(r=1)^d binom(d,r) r! S(t,r).

Also 1! S(t,1)=1, so the constant c changes only the r=1 coefficient.
This proves sufficiency of (1). For necessity, the largest nonzero clone
coefficient A_m gives a positive coefficient of m^t in the ordinary
exponential expansion. An eventual identity with leading base d forces
m=d. Linear independence of the finite exponential sequences, followed
by the triangular change to the surjection basis r! S(t,r), gives (1).
This argument does not substitute t=1 into an assumed eventual identity.

Equation (1) forces H(R)=h. The Moon--Busch lower bound for strong
tournaments tends to infinity with their order. Therefore only finitely
many base orders can have H=h. Enumerate all tournaments through that
bound, retain the strong ones with H=h, and test every physical mark and
its complete run vector against (1). Isomorphism reduction is optional;
it does not affect termination or completeness. The singleton cannot
produce a growing target, and is excluded separately. If h is even,
Hamilton-count parity already gives a negative answer.

This uses no assumption that finite tests prove a uniform rule. The
order reduction and coefficient identity are the uniform proof; the
finite test checks their exact remaining predicate. It offers no useful
complexity guarantee at large h and no smaller unrestricted recurrence
for A000568.

## Cheap arithmetic and size obstructions

The run formula has two additional inherited restrictions:

    (r-1)! divides A_r,
    A_r=0 for r>floor((|R|-1)/2)+1.

At r=d these imply

    (d-1)! divides a,        |R|>=2d-1.                 (2)

The factorial condition also implies every lower necessary divisibility
in (1), since (r-1)! divides (d-1)!. Combining the size condition with the
strong lower bound can reject a proposed clock before enumerating a base.
For example, an odd coefficient of 3^t is impossible in a pure family
a*3^t+c arising from this single homogeneous transitive inflation. This
does not obstruct iterated tripling by changing graphs or adding different
modules: those are different operations.

Conditions (2) are necessary, not sufficient. The exact positive control
in the clone packet has A=(89,63,0,0), and realizes 63*2^t-37. The higher-run
hostile has A=(9,10,2), so although its leading coefficient satisfies
2!|2, it gives 2*3^t+4*2^t-5, not a pure 2*3^t+c clock. Its second-run
coordinate would have to be 2*binom(3,2)=6, rather than 10. This is the
smallest missing coordinate in that attempted simplification.

For the present 24/44 families, d=2,a=63,c=-35 or -29 gives h=91 or97.
The complete marked-base tests and their independent audit establish
failure of (1). A changing reader escapes the fixed-base premise, while
retaining a new acceptance obligation at every stage. The next useful
analogue for several growing blocks is a joint run tensor, together with
its realizability constraints; a collection of marginal run vectors is
not asserted sufficient.
