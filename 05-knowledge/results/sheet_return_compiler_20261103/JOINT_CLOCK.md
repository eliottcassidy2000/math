# A scalar cubic, a joint cubic, and the exact retained phase

**PROVED with exact transfer certificates.** Operation scope: grow the
specified block of A_n or B_n; observe the fixed old root r=(0,2), including
any of the five inherited Q4 presentations as readouts. This is not closure
under arbitrary re-marking, complements, or arbitrary vertex attachments.

## 1. Two physical probes recover the whole integer observer

Use the inherited definitions

    v_n=(H_n,X_n,Y_n),
    E_n=Y_n-X_n = Hamilton paths ending at r,
    S_n=X_n+Y_n = ordered two-bin covers whose first path starts at r.

Empty bins follow the inherited convention; the first path is nonempty.
These are joint path-cover counts, not products of start/end marginals.

Attach the actual q2 and q5 four-vertex gadgets at this *same physical root*.
Let their Hamilton counts be a_n,b_n. The inherited Q4 matrices act on row
vectors, giving

    a=5H+8X+6Y,    b=5H+10X+8Y,
    S=(b-a)/2,    E=(7b-9a+10H)/2,
    X=(4a-3b-5H)/2, Y=(4b-5a+5H)/2.               (1)

The integer parity guards are retained; division by2 is not licensed after
reduction modulo an even modulus. The full relation, root address, and
ordered gadget attachment are in `bulk.graph`. No class-valued reset occurs.

The two probes share the same17-state paired-boundary bulk matrix and initial
row. They have different terminal columns. `phase_bulk.json` saves them
separately. The exact characteristic polynomial of that bulk matrix is

    X^5 (X-1)^2 (X+1) f1(X)^2 f2(X),
    f1=X^3-X^2-X-1,  f2=X^3-2X^2-X-2.

Let beta2,beta5 be its two terminal columns and alpha its actual initial row.
The exact reachable-space certificate for S is

    alpha B^j f2(B)(beta5-beta2)/2 = 0,
                    0<=j<17.

Cayley--Hamilton proves this for all j. A block-diagonal certificate with the
four-state H transfer proves `(shift-1) E=0`. These checks are over Q with
exact rationals; no floating point or inference from a fitted prefix is used.

Thus, for both actual families at every n>=2,

    H_(n+3)=H_(n+2)+H_(n+1)+H_n,
    S_(n+3)=2S_(n+2)+S_(n+1)+2S_n,
    E_n=3 on A; E_n=2 on B.                            (2)

The initial triples at n=4,5,6 are

| Family | H triple | S triple | E |
|---|---|---|---:|
|A|(873,1649,3007)|(24555,65781,174783)|3|
|B|(889,1631,3007)|(10486,27854,74078)|2|

Literal n=2,3 controls extend the same recurrences to the stated starting
range. The S Hankel determinants are -36860651424 and -53890048. Both are
nonzero. The two irreducible-over-Q cubics are coprime, and E is nonzero.
Consequently the **minimal common rational linear recurrence** of the
three-coordinate observer has order7, with polynomial `(X-1) f1 f2`.
Equivalently, any constant linear state system generating this whole vector
sequence has dimension at least7. If the family's known E is supplied as an
external affine constant, six dynamic coordinates suffice. This is a linear
realization claim, not an information lower bound against storing n itself.

At n=6 the scalar counts coincide, but the rooted observers are

    A: (3007,87390,87393),   B: (3007,37038,37040).

They have phases `(5,5,6)` and `(5,3,1)` respectively. Their next scalar
counts5529 and5527 also differ. Joint incidence has survived scalar collision.

## 2. A seven-coordinate integer carrier and its inverse guard

Use `(H_n,H_(n+1),H_(n+2),X_n,X_(n+1),X_(n+2),E)`.
The first triple advances by f1, the last coordinate is fixed, and

    (x0,x1,x2) -> (x1,x2,2x0+x1+2x2+2E).             (3)

This integer carrier reconstructs `Y=X+E` without fractions. In S-coordinates
the corresponding realization lattice is `S0=S1=S2=E (mod2)`, of index8.
Retaining X rather than arbitrary unconstrained S removes this artificial
lattice expansion. Physical tournament realizability still requires the
actual source family and its parameter; lattice membership alone is insufficient.

The transfer D has determinant2 and, for every j>=0, Smith factors

    (1,1,1,1,1,1,2^j) for D^j.                       (4)

Proof: modulo2 the H block is invertible, E is fixed, and the X block has a
two-dimensional image on which it exchanges two coordinates. Every positive
power therefore has rank6 modulo2. Since the determinant of D^j is2^j,
only one Smith factor can be even and it must be2^j. Inverting one step
requires the actual guard

    x2' = x0' (mod2),
    x0=(x2'-x0'-2x1'-2E)/2.

This is a concrete instance of the inherited integer j=2 recurrence clock,
not a new general Smith theorem and not field Frobenius. A guard cannot be
discarded because the same matrix is invertible modulo7.

## 3. F49 survives, and F343 is also needed

Modulo7,

    f1=(X-3)(X^2+2X-2),
    f2=X^3-2X^2-X-2 is irreducible.

In the already-used field `F49=F7[theta]/(theta^2-4theta-2)`, put
`beta=2+2theta`. It satisfies the quadratic factor, has order48, and
`beta^7=3+5theta`. For the H triple define

    r0=H_(n+1)-3H_n, r1=H_(n+2)-3H_(n+1),
    z=(r1+4r0)+2r0 theta.

Then one *actual bulk insertion* sends z to beta*z. This supplies an explicit
F49 transport map from the recurrence, not just a field-cardinality analogy.

It is distinct from the current golden observation

    (a,b,g)=(4(H+2X-Y),Y,H-Y) mod7.                   (5)

Equation(5) remains an invertible three-coordinate observation of `(H,X,Y)`;
the five Q4 matrices act on the full integer observer before reduction.
Frobenius in the theta chart is `(a,b)->(a+4b,-b)`, which swaps the inherited
normal-basis coordinates `(b+2a,2a)`. It has order2. Bulk transport has order48
on its quadratic factor. These operations are not conflated.

The S triple is a one-dimensional space over F343, with growth equal to
multiplication by a root eta of f2. Its order is171, not342:

    eta^57=2, eta^171=1.

Its field Frobenius has order3, again different from the growth clock.
Both actual families have nonzero quadratic and cubic components. Therefore
their full joint state has exact modulo-7 period `lcm(48,171)=2736`.
The checker also enumerates the orbit and verifies minimality independently.

The current phase and its next two growth observations recover H's three
coordinates, S's three coordinates, and E. Distinct states in either orbit
therefore have distinct future phase observations. Each actual orbit needs
exactly2736 deterministic phase states for this language. Their disjoint
union needs5472 because E=3 and E=2 are themselves observable. This does not
make 2736 the smallest carrier for arbitrary strong atoms.

For every precision a>=1 the exact period becomes

    2736*7^(a-1).                                    (6)

Here D^2736=I mod7, while `(D^2736 v-v)/7` is nonzero mod7 for each initial
integer state. In the displayed H/X/E coordinates those residue vectors are
`(0,3,0,1,2,2,0)` and `(0,1,3,3,6,4,0)`.
Writing D^2736=I+7A and expanding successive seventh powers proves exact
valuation growth by one at every lift. The period modulo7 already forces
every higher period to be a multiple of2736, proving (6).

The inherited twelve-vertex 81/27 separator is unaffected: this carrier is
defined for two explicit marked families and their bulk operation. It is not
asserted to reconstruct arbitrary tournaments from these seven numbers.
