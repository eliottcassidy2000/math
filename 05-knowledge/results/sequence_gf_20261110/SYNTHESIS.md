# Sequence discovery by choosing the right carrier

**Session249, 10 November 2026. COMPLETE.** New proved-in-note identities,
structural criteria and runnable exact examples. Inherited results are
credited rather than re-audited. No independent canon audit or external
mathematical-priority claim. **Chan–Pak remains OPEN.**

The main improvement is a usable way to choose a generating-function class
before guessing coefficients. Four seemingly different objects now share
one workflow: retain the operation, expose its connected pieces or boundary
state, then choose the variable in which the operation is simple.

## 1. What was actually added

### A. A connected-interaction compiler for A000568-related sequences

A000568(n) counts ordinary tournament isomorphism classes on n vertices.
Its inherited Burnside formula sums over odd cycle partitions. Retaining
the edge-orbit variable q, expand at q=1+s rather than immediately q=2.
At the soluble point q=1 the series is

    A(x,1)=sqrt((1+x)/(1-x)).

The coefficient of s^r in log A requires only connected graphs on at most
r+1 permutation cycles. An edge carrying b marks has arithmetic weight

    Psi_b(d)=sum_(c|d) mu(d/c) binom(c,b),

the number of primitive b-subsets modulo d. Divisor constraints merge to
an lcm at every graph vertex. This produces exact Lambert-type generating
functions and convolution recurrences, with polynomial complexity in the
requested size horizon when r is fixed.

For instance, with D=x d/dx and
V[a,d]=binom((D-1)/2,a) artanh(x^d)/d,

    [s] log A(x,1+s)=V[1,1]+(1/2)sum_(d odd) phi(d)V[0,d]^2.

The first factorial-normalized marked sequence starts

    0,0,1,5,22,138,927,7431,57564,...

These count odd-order permutations with a marked pair-orbit. They are
useful A000568-related sequences, not new Hamiltonian-path values. All
Taylor layers are nonnegative and their full sum gives A000568, but every
fixed finite jet captures an asymptotically vanishing fraction of that
superexponentially growing census. This precise boundary prevents a false
fast-count claim while retaining a reusable fixed-interaction engine.

[Full connected-diagram theorem, formulas and compiler](burnside/README.md).

### B. An exact criterion for rational growth-family generating functions

Fix a quotient tournament Q. Some modules grow as repeated ordered joins
A_i^{▷r}; all other modules remain fixed. Then

    H_r is C-finite  iff  Q restricted to the growing modules is transitive.

If the fixed context has b actual vertices, each growing module appears
in at most b+1 maximal runs of a Hamiltonian path. The path-cover/Newton
transform therefore gives an explicit finite exponential sum with bases

    product_i Q_(A_i)(j_i),   1<=j_i<=b+1.

Equal bases combine; the surviving distinct bases give the exact minimal
characteristic polynomial. The result includes multivariate repetition
parameters, not just one diagonal sequence.

Examples:

    H(C3[C3^{▷r},1,1]) = 3^r+12^r,
    H(C3[P_r,P_r,1]) = 4^r-2*2^r+3.

If three growing modules contain a directed triangle, they already support
at least3(r!)^3 paths. This rules out any eventual constant-coefficient
recurrence. It does NOT rule out polynomial-coefficient recurrences.
The boundary depends on which modules grow, not whether all of Q is cyclic.

A second compiler converts a nonnegative finite matrix grammar into
positive multivariate functional equations that retain the actual output
value. For a fixed target h, clipping every state coordinate at h+1 is
exact: it produces a finite automaton, rational word-length GF and a
terminating witness search. Signed or arbitrarily minimized matrix
representations do not support this clipping argument.

The one-apex example has explicit catalytic equations, all-length fixed-H
word GFs C_h(z)/(1-z)^2, and an all-order certificate that this family misses
13. That failure is useful:13 occurs in the full spectrum, so richer
branching is necessary. A finite bank of the one-parameter rational
families above yields only O(log X) distinct values below X.

[Growing-quotient theorem, support compiler and examples](grammar/PROOF.md).

### C. Two across-order recurrences from heat flow and coherent birth

S248 supplies a vertex-support spectral coordinate; the even-path forest
expansion supplies arc degree. Retain both. For normalized path count
kappa_n=2^(n-1)H/n!, let arc noise have multiplier sqrt(y) and vertex
retention probability z. Then

    F_n(y,z)=n! E[kappa_n · N_sqrt(y) K_z kappa_n],
    sum_n F_n(y,z)x^n=sum_k k! f(x;y,z)^k,
    f=x+2yz^3 x^3/(1-yz^2x^2).

The rational pullback of the factorial-series ODE yields an explicit
seven-step polynomial recurrence in n. The ODE technique and unweighted
specialization were already known; the jointly support-refined kernel is
new here. For example,

    F_6=720+192yz^3+8y^2z^5+8y^2z^6.

The two equal-edge-degree terms with support five and six are now distinct.
These are uniform labelled-tournament covariances, not unlabelled counts.

There is also an exact coherent-birth martingale in the unnormalized
signed forest coordinates. Its two-time energy B_n(a) has EGF

    sum_n B_n(a)x^n/n! = exp(x+2ax^3/(1-ax^2)),

and a five-step polynomial recurrence. At a=1 the sequence begins
1,1,1,13,49,361,3121,... . Fixed Fourier-degree energies are polynomials
in n and have finite-difference recurrences. Thus ordinary factorial
composition and exponential assembly yield different, explicitly related
sequence engines from the same small connected path kernel.

For higher moments, the diagonal forest pair is replaced by connected
colored overlap graphs, followed by factorial extraction. At moment order
three, every edge has one of three pair-colors and the union graph is
subcubic. The first genuine three-layer terms are a positive triangle and
a negative three-star. The sign is essential: positivity of the Burnside
connected diagrams cannot be exported to this different carrier.

Finally, write Y=(H-1)/2 and beta_r=E binom(Y,r). For odd L,

    sum_(r=j)^(j+L) (-1)^(r-j) binom(r,j) beta_r <= Pr(H=2j+1).

A positive finite sum certifies an actual tournament with that exact path
count. The full hierarchy is necessary and sufficient for support, but
there is no uniform positivity bound yet. This is the precise coefficient
obligation that a higher-moment recurrence would need to discharge.

[Covariance, birth energy, replica kernels and support certificates](heat/README.md).

### D. A Mahler compiler for certified values, and why11 is exceptional

THM-5126 already proves65 dyadic lines with odd part c<=129, omitting only
7 and21. The new step is to convert this set into a positive functional
equation. With one seed onset s_c per odd part,

    B(z)=R(z)+B(z^2),   R(z)=sum_c z^(c2^s_c),

and zB marks the values1+c2^s. This yields a regular binary language,
exact integer counting formulas and optional witness-size weights.
The finite-ray indicator is2-automatic but is NOT P-recursive: it has
arbitrarily long zero gaps and infinitely many nonzero terms. The same
selected witnesses counted by graph order have a rational series.
The sequence index changes the lawful recurrence class.

The full conjectural target also has a positive Mahler source, but only
after separating11, since2*11-1=21 is forbidden:

    P_*(z)=z+z^11+z sum_(k>=0) R_*(z^(2^k)),
    R_*(z)=z^2/(1-z^4)-z^6-z^10+z^12+z^40.

R_* has only nonnegative coefficients. The formula is an exact description
of the conjectured target SET, not a realization theorem. Its unbounded
source still needs graph certificates, and a root giving a late ray still
leaves the points below that ray to settle. This is a concrete loss ledger
for attempts to prove completeness using dyadic patterns alone.

[Positive source theorem,65-line equation and automaticity](rays/PROOF.md).

## 2. A practical recurrence-discovery workflow

1. **Specify the coefficient.** Class count, marked permutation count,
   family path count, labelled moment, fixed-value multiplicity and value
   indicator are different sequences. Specify its index as well.
2. **Identify a closed operation.** Odd-cycle interactions suggest log and
   connected diagrams. Fixed-context substitution suggests Newton inversion.
   Nonnegative boundary matrices suggest catalytic monomial equations.
   Affine value doubling suggests a Mahler source.
3. **Test growth before fitting.** A cyclic growing quotient forbids
   C-finiteness. Long zero gaps forbid P-recursiveness of finite-ray support.
   A000568's superfactorial census cannot be inferred from a fixed marking
   jet. None of these obstruct the correctly typed equations above.
4. **Preserve the target-bearing coordinate.** Vertex support separates
   heat modes; nonnegative states permit exact clipping; H must remain an
   exponent or a complete moment hierarchy for support questions.
5. **Derive and then check a small example.** The four compilers implement
   proofs rather than using finite-prefix fits as evidence for a theorem.

The connected-graph commonality is a real mechanism, but has different
signatures. Burnside log diagrams are positive; replica moment diagrams
can be signed; positive support languages retain actual realizability.
The graph decomposition transfers, not automatic coefficient positivity.

## 3. Strongest next research moves toward Chan–Pak

- Apply the positive support compiler to an inherited two-apex or
  cutwidth-two grammar. Produce its catalytic equations and seek one
  target-uniform reduction for accepting states. The fixed-target
  clipped automata are exact but are not that induction.
- Develop a finite rooted grammar for the TREE sector of the signed
  three-replica overlap kernel, retaining unfinished-path parity and
  orientation. Then add controlled cycle excess. A tree-only recurrence
  would still not enumerate the complete third moment.
- Lift the positive all-line source R_* to actual rooted tournament
  states, with explicit thresholds and finite below-ray certificates.
  Do not identify its arithmetic labels with realizable seeds by fiat.
- Retain marked path information in the Burnside connected compiler.
  Scalar census marks alone cannot see whether a particular H occurs.

These are four concrete next objects rather than four requests for a
larger verification sweep. No all-order coverage or conjecture closure is
claimed in this session.

## 4. Reproduction and publication scope

From this directory run `python3 run_all.py`; all scripts use the Python
standard library. It runs the four new-formula examples. The heat script
refreshes its exact `output.json`; other outputs are printed.

Recorded controls cover marking layers through order4 and size24, exact
context/substitution cases, the new covariance/birth recurrences through
n=40, signed replica components through five vertices, finite Bonferroni
bounds, and the ray functional equations through100,000. These checks are
not independent re-audits of inherited tournament theorems. Individual
notes list exact universes, weights, commands and stopping boundaries.

No new canon ID is assigned. The notes distinguish proved formulas,
inherited membership, finite controls and open coverage. The previous
S248 symmetry and S239 sequence engines remain the starting points,
not results rebranded as discoveries of this session.
