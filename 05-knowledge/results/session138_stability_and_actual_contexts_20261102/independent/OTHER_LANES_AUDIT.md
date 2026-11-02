# Independent audit of Session138 structural and exact-resource lanes

Status: ACCEPT, in the stated scopes. All-order scalar zero-rigidity and the
unrestricted odd Hamiltonian-count spectrum remain OPEN.

## Homogeneous pair and actual-kernel compiler

The proof in `../zero_rigidity/README.md`, sections1--3, is valid.
Swapping the homogeneous pair in the cyclic completion changes exactly its
internal edge relative to the middle completion. Reversing an adjacent pair
preserves both exterior attachments because its vertices are twins. The only
unmatched words contain exactly one of `x,w,u` or `w,u,x` consecutively.
Contracting them gives the stated endpoint and free-gap coefficients.
Polarization therefore gives

    H(C)-H(M)=2F_D(r)+(s-r)·grad F_D(r)
             =2H(D)+sum_(i->j)c_ij[r_i(1-s_j)+s_i(1-r_j)]
             >=2H(D)>0,

including the empty case, where exactly two unmatched words remain. The
quadratic gradient sign and resulting capacity recursion both check.
For strong D of order>=3, the stated partition argument proves lower-bound
equality only at R=S=empty or R=S=V(D).

The homogeneous-C3 identity is also correct: when the third triangle vertex
is separate, the marked pair has two orders; when adjacent, one marked-pair
order admits both placements of the third vertex and the other admits none.
Either way the count is twice H after one triangle-vertex deletion.

The quotient run factors C3:(3,6,6), TT2:(1,2), K1:(1) are correct ordered
path-cover counts. Differentiating the unrestricted transition polynomial
at a binary quotient counts a single free marked transition, even if the
marked adjacency variable is zero. Independent intact-block rotations give
the factors9 and3. TT2 endpoints are only summed, not falsely identified;
the pair-gradient difference reconstructs them individually. These statements
relate actual induced kernels and do not alone prove scalar zero uniqueness.

The strengthened formal firewall was replayed with its active checks. An
additional exact check found actual zero cuts `[1793]`, while the formal zero
cuts are `[480,736,1248,1793]`. All three escaped cuts have actual value239643,
not target241137. The modified minimum capacity is2538>0. The advertised
capacity, symmetry, layer-mean, gcd, pair-deletion and twin constraints are
preserved, but the actual cross-block quotient-derivative identities are not.
This is a valid formal-network countercontrol, not an actual tournament
counterexample. Receipt: `formal_firewall_replay.log`.

## Literal one-edge reflector

The graph maps and all polynomial identities in `../exact_resources/PROOF.md`
check. Both quotients contain a directed Hamilton cycle and are strong;
substitution of any nonempty U preserves strongness. The two quotient maps
change exactly exterior edge{1,3}. Their run responses are11h+12p2 and
13h+12p2, giving6E-h and6E+h. The sole exterior outneighbor bounds the module
run count by2, so the finite response-row calculation proves a universal
identity rather than extrapolating from sampled modules.

The transitive-module cycle proof has at most one odd cycle: every cycle
meeting U contains exterior vertex4 and at least one of1,2,3, while the only
possible exterior triangle avoiding4 is123. The external path-count vectors
(2,3,1) and(2,2,2) are correct. Algebra simplifies to the stated shift operators
and difference operator. Linearity in the full path-cover polynomial extends
these identities from the basis z^k (k>=1) to every actual nonempty module.
The negative p2 change at k4 is correctly -22 despite scalar increment+2.

Ordered joins multiply both H and E, including the empty-path convention.
This proves the unbounded strong fixed-gap families without making them
induced extensions. The source-certified inverse retains an actual core and
its E, never admitting a formal pair. Given the inherited exact core-order
bound K(h), the height-limited enumeration is finite and complete for this
constructor image. The inequality E>=(|U|+1)H(U) is an injective Hamilton-word/
cut count; at h>1 it forces minus output>=23h and plus>=25h. Consequently
all outputs<=200 have h=1 or h<=8, where core order<=4, yielding exactly the
listed16 outputs. The minimal-context claim has the correct finite universe:
all port0 labelled contexts through order5 and all six exterior edge flips.
It does not assert uniqueness at larger context orders.

The output-counting bound over labelled cores and padding exponents is valid.
Using K(B)<=3log_3 B gives exponent9c^2/[2(log_2 3)^2], hence the stated
zero-density resource threshold. The self-contained verifier was rerun both
normally and with Python -O; receipts `exact_resources.normal.json` and
`exact_resources.optimized.json` agree byte-for-byte.

## Bounded-source profinite obstruction

`../exact_resources/RESOURCE_BUDGET.md` is valid. The inherited source normal
form retains every nonsingleton SCC, so replacing its coarse retained-core
bound h by K(H) gives precisely N(H,0)=K(H)+1 and
N(H,I)=K(H)+3+floor(log_2 I) for I>=1. Labelled n-vertex parents times all cuts
number2^(n(n+1)/2), hence C_N is a deliberate complete overcount of bases.
Every base supplies at most one exceptional native value and one dyadic ray.
Modulo2^r-1 each ray has at most r residues, regardless of signed coefficients,
offsets or repeated residues. Thus C_N(r+1)<2^r-1 certifies an omitted residue
for the entire bounded-budget class, not for the unrestricted spectrum.

The quantitative necessary inequality and asymptotic lower bound
N>=(sqrt2-o(1))*sqrt r follow from C_N<=N*2^(N(N+1)/2). Substituting the
upper bound on N gives the stated necessary height/interior resource growth.
The H1,I0 residue language modulo15 is exactly{1,2,3,5,9}. No substitution-width
bound is silently identified with a bound on these actual-source resources.
`budget.py` was rerun normally and with -O into this audit directory; matching
receipts retain its10000 exact core-budget comparisons and1331 signed ray
period checks.
