# Local lattice universality and coupled valuation certificates

Session 188 extension of `../primeflow/PROOF.md`.

**Proved:** every finite lattice occurs as an actual interval of the local
prime-valuation closure lattice; finite coupled row/column propagation is sound
and terminates; a two-source arithmetic branch defeats all separate row and
column feasibility checks. A four-source arithmetic control proves that even
exact interval propagation remains incomplete, and is rejected by a parity
certificate. These are local and partial-branch results, not global OPN bounds.

## 1. Recall the exact local compiler

For an odd prime source p, a finite set of other odd prime targets R, and
ordinary or Euler exponent parity, put U_p(n)=(p^n-1)/(p-1). Let
Phi_p(n)_r=v_r(U_p(n)). For allowed targets let d_r=ord_r(p) and

    rho_r(k)=r^k                         if d_r=1;
    rho_r(k)=d_r r^max(0,k-c_r)           if d_r>1,

where c_r=v_r(p^d_r-1), k>=1, and rho_r(0)=1. With b=1 for ordinary and b=2
for Euler sources, define L_p(u)=lcm(b,{rho_r(u_r)}) and C_p(u)=Phi_p(L_p(u)).
Positive demands at parity-disallowed targets are impossible.

The preceding packet proves that u<=Phi_p(n) exactly when L_p(u)|n, so C_p(u)
is the least realizable profile above u. A valuation box [u,v] is feasible
exactly when C_p(u)<=v. The ordinary all-zero corner has a positive physical
witness at a suitable odd prime length. Exact profiles are the fixed points
of C_p; their meet is coordinatewise minimum and their join is closure of the
coordinatewise maximum.

## 2. Every finite lattice is an arithmetic interval

**Theorem 1 (interval universality).** For every finite nonempty lattice K,
there is an odd prime p, a finite set R of other odd primes, and a realizable
finite profile c such that the interval [0,c] in the ordinary-source profile
lattice is isomorphic to K. The meet and join are those of the ambient
valuation-profile lattice, not merely newly defined operations on the image.

**Proof.** The one-element lattice is realized with R empty. Otherwise let
0_K and 1_K be the bottom and top. Assign distinct odd primes ell_t to each
t in K other than 1_K, and set

    M = product_{t != 1_K} ell_t,
    d_x = product_{t != 1_K, x not<=t} ell_t.

Then d_0=1, and d_x>1 for x!=0 because ell_0 occurs. Directly,

    d_x divides d_y  iff  x<=y,
    d_(x join y) = lcm(d_x,d_y).                              (5)

For the reverse implication in the first statement, when x not<=y take t=y;
then ell_y occurs in d_x but not d_y. The exceptional t=1_K cannot occur
in this argument because x<=1_K always. The second statement follows from
x join y<=t iff both x<=t and y<=t.

Choose a primitive prime divisor p of 2^M-1. Bang–Zsigmondy supplies it because
M>1 is odd, so neither the exponent-two exception nor (base,exponent)=(2,6)
is possible. Thus ord_p(2)=M, p is odd, and p=1 mod M.

For each x!=0 choose a primitive prime divisor r_x of p^d_x-1. Again d_x>1
is odd, so the exceptions are absent. Then ord_(r_x)(p)=d_x, so r_x is odd,
different from p, and different for distinct x. Moreover r_x divides no label
d_y: all label primes ell_t divide p-1, whereas a primitive target has order
greater than one. This is the crucial control on unintended LTE lifting.

Set R={r_x:x!=0} and c_x=v_(r_x)(p^d_x-1). For every y in K, the exact
valuation law now gives

    f(y)_x := Phi_p(d_y)_(r_x) = c_x * 1[x<=y],              (6)

because r_x does not divide d_y. Thus f(0)=0 and f(1)=c. The mathematical
length 1 used for f(0) is harmless: a fresh odd prime length realizes the same
zero profile at a positive even exponent.

For meets, min(f(y)_x,f(z)_x)=c_x exactly when x<=y and x<=z, equivalently
x<=y meet z. Hence the **ambient coordinatewise meet** is f(y meet z).
For joins, the lower rank of a positive coordinate at most c_x is d_x.
Therefore the LCM of the demands in max(f(y),f(z)) is
lcm(d_y,d_z)=d_(y join z), so the **ambient closure join** is f(y join z).

Finally take any realizable profile w with 0<=w<=c. If its x-coordinate is
positive it is at least c_x, so it equals c_x. Its lower length is the LCM
of the corresponding d_x, hence d_y for their join y (empty join means 0).
Since w is fixed by closure, (6) gives w=f(y). Thus the interval contains
exactly the embedded image, not extra profiles. The map is injective by (5),
or by coordinate y when comparing y not<=z. QED.

This construction gives no bound useful for computation: its primes and
labels can be enormous. Its force is structural. No lattice identity that
fails for some finite lattice can hold for every local prime-flow closure
lattice. In particular neither distributivity nor modularity is available as
a universal local obstruction. This does not assert that these configurations
occur inside any perfect number.

### A small explicit pentagon, in addition to the earlier diamond

The nonmodular pentagon N5 has 0<a<b<1 and 0<c<1, with c incomparable to a,b.
Choose label primes (3,5,7,11) for (0,a,b,c). The labels for (a,b,c,1) are
(33,165,105,1155). An explicit source and targets are

    p=18481,
    (r_a,r_b,r_c,r_1)=(67,331,50821,4621).

All five displayed numbers are prime. Their orders are (33,165,105,1155),
their base valuations are all one, and p=1 mod 1155. Thus the interval below
1111 consists exactly of

    0000, 1000, 1100, 0010, 1111.

The first profile has physical witness length 3; the others have lengths
33,165,105,1155. In particular a<=b but

    a join (c meet b)=a,     (a join c) meet b=b,

so the ambient local profile lattice is not modular. Unlike the general
existence proof, every prime, order and valuation in this example is checked
by deterministic integer arithmetic in the reproduction script.

## 3. Coupling rows by exact column totals

Let P be finitely many prime-power sources with specified roles, and let R be
finitely many tested target primes. If a target equals a source, set that
entry permanently to zero. Suppose a branch supplies finite integer boxes

    ell_pr <= w_pr <= u_pr,

and exact residual column totals A_r. Every row must be a locally realizable
profile, and the desired allocation satisfies sum_p w_pr=A_r. Residual totals
can account for already known contributions from sources outside P.

These constraints are only a relaxation of an OPN branch: source exponents
are not required to equal the column sums at their own primes, and outside
factors are not excluded. Therefore infeasibility rejects the original branch,
while feasibility is inconclusive.

**Theorem 2 (finite propagation certificate).** Repeatedly perform full sweeps
of the following rules, stopping on contradiction or when a sweep makes no
change:

1. Replace each row lower bound ell_p by C_p(ell_p). Reject if it exceeds u_p
   or asks for a disallowed positive valuation.
2. Reject if A_r lies outside [sum_p ell_pr, sum_p u_pr].
3. Replace ell_pr by max(ell_pr, A_r-sum_(s!=p) u_sr).
4. Replace u_pr by min(u_pr, A_r-sum_(s!=p) ell_sr).
5. Optionally replace u_pr by the greatest k<=u_pr such that
   C_p(max(ell_p,k e_r))<=u_p. If the row was feasible, this set is nonempty.

Every step preserves every feasible allocation. With finite starting boxes,
this full-sweep procedure terminates in a contradiction or a fixed point.
Under any other scheduling, there can be only finitely many strict updates.
Every rejection
has a finite replayable certificate using only integer sums, prime-power
ranks, LCMs, and valuations at the resulting lengths.

Proof. Step 1 is exact local least completion. Steps 2–4 are consequences of
the column equation and the other entries' bounds. For Step 5, every feasible
row with r-coordinate k dominates max(ell_p,k e_r), and its least completion
must stay inside u_p. Conversely the closure supplies a row attaining at
least each feasible threshold k. Threshold feasibility is downward closed,
so the maximum may be found by binary search. Each strict nonrejecting update
decreases the nonnegative integer sum_(p,r)(u_pr-ell_pr), proving termination.
Recording the previous boxes and each displayed inequality gives direct
certificate replay. QED.

### A contradiction invisible to separate rows or separate columns

Use ordinary source 3 and Euler source 5, targets (11,23), lower bounds zero,
and row upper bounds

    u_3=(2,1),    u_5=(1,1).

Require exact column totals (A_11,A_23)=(2,2). Both targets have ordinary
even exponents in the corresponding OPN branch. The two rows separately are
feasible, and each column equation separately is feasible even while retaining
the full row restrictions:

* Total two at 11: source lengths (n_3,n_5)=(5,2) give profiles (2,0),(0,0).
* Total two at 23: lengths (11,22) give profiles (0,1),(0,1).

But jointly the first column forces w_3,11>=2-u_5,11=1; the second forces
w_3,23>=2-u_5,23=1. Exact local closure at source 3 gives

    L_3(1,1)=55,    C_3(1,1)=(3,1),

contradicting u_3,11=2. This is a two-column, two-source arithmetic certificate,
not a single-row contradiction already present in the initial bounds. It
does not assert a complete support or hypothetical perfect-number witness.

## 4. A modular obstruction beyond the propagation fixed point

Interval propagation is not complete, even if Step 5 enforces exact attainable
row endpoints. Consider four ordinary sources

    P=(191,653,991,1543)

and the single target 13. Every source has order three modulo 13 and
v_13(p^3-1)=2. Give every row the interval [0,2] and require column total five.
Each row's exact feasible values are {0,2}: length 5 witnesses zero and length
3 witnesses two. Every interval endpoint is attainable. The column interval
is [0,8], and each residual lower bound is max(0,5-6)=0; its upper bound remains
two. Thus all five propagation rules reach a nonrejecting fixed point.

Nevertheless each contribution is even, so their sum cannot be five. The
one-bit modular certificate rejects the allocation. Target 13 could have the
Euler exponent five; the parity issue is not an artificially invalid target
role. In OPN language the conditional branch assumes that these four sources
are the only positive contributors at 13 and that each contributes at most
two. No such branch can occur.

This connects weighted local closure to S184's modular selection layer using
an actual arithmetic partial-row example. It is not a claim that an entire
support survives S179 exact closure or S184's complete algorithm. Global
exponent compatibility, exterior factors, and unrestricted support coverage
remain separate requirements.

## 5. Sources and verification

The classical existence input for Theorem 1 is Bang–Zsigmondy, with exceptions
explicitly excluded in both uses. See [Gyory and Smyth, *The divisibility of
a^n-b^n by powers of n*, Integers 10 (2010), Theorem 13](https://math.colgate.edu/~integers/k27/k27.pdf).
Its role is primitive prime existence, not novelty of the present lattice
formulation. No historical-priority claim is made.

Run `python verify.py` and `python -O verify.py`; the receipts are deterministic.
The script checks the order-label construction over finite lattices extracted
from all naturally ordered relations through six elements. Those are symbolic
encoding checks, not claims to construct all their enormous arithmetic
witnesses. The pentagon and both multi-source examples are checked by direct
integer arithmetic, including direct geometric-sum profiles. The proofs,
not finite enumeration, establish the universal statements.

Both execution modes pass identically. The symbolic census contains 51
naturally ordered lattice relations, with 1,625 pair comparisons and 1,383
binary interval checks. The propagation control covers 400 cases, including
176 feasible ones and 206 preserved witness allocations. The explicit
four-source parity example remains nonrejecting under interval propagation
and is rejected by its separate modulo-two certificate, as claimed.
