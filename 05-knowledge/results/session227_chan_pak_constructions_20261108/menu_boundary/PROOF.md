# The bounded-intercept obstruction survives every finite graft menu

**Status: PROVED; independent audit PASS.** The mathematical argument is
all-order. The accompanying controls are FINITE-EXACT. This extends the
construction language of Session226's obstruction, rather than claiming a
new original obstruction or a literature-priority result. The Chan–Pak
all-odd Hamilton-count problem remains OPEN.

## 1. Inheritance and exact scope

The inherited result is
`session226_structural_continuation_20261108/root_descent/PROOF.md`, Sections
2–4: bounded endpoint intercepts give only logarithmically many distinct
counts in the O/Z/J language, and a common even-g residue gap then obstructs
a finite line bank. The present extension replaces O/Z/J by O, Z, and **any
fixed finite menu of nonempty rooted tournament grafts, of arbitrary fixed
orders**. A singleton graft is J and may be included. The closing OJJ is
available irrespective of whether J belongs to the interior menu.

We inherit the graph construction and its four-state law from
`session225_graft_recovery_20261108/PROOF.md`, Sections 1–2 and 5. The
physical inequalities and their path-cut/insertion proofs are in Session226
`root_descent/PROOF.md`, Section 1, and the exact closed cone is in this
session's `graft_images/PROOF.md`, Section 1. The inequalities themselves are
older endpoint-defect facts, as those sources explain. Only cone containment
is needed for the obstruction; exactness additionally proves optimality of
the ratio constant below.

Closest mechanism: the inherited last-non-O decomposition. Hostile: even
an escaping graft state can have an already available scalar count.
Corrected near miss: a more expressive state grammar need not defeat a
bounded-intercept counting obstruction. Least-used sidecar: positivity of
the graft's D row on all four extreme rays. Concept board: physical cone;
graft ratio; finite cores; O tails; exact residue gaps.

The map takes physical rooted graph histories to endpoint coordinates and
then scalar Hamilton counts. It preserves the stated formulas and necessary
congruences; it loses graph structure and reachability. Actual graphs or
full states remain the needed sidecars. None of the following identifies
the scalar image, excludes every unrestricted finite nonunit line bank,
or treats arbitrary nonphysical seeds or alternative final constructions.

## 2. States, rays, and graft ratios

For a nonempty rooted tournament (A,a), use

    x=(D,E,U,d),

where D=H(A), E counts Hamilton paths ending at a, d=H(A−a), and U counts
ordered spanning path pairs (P,R), with P possibly empty and R nonempty and
starting at a. An allowed empty path counts once, including H(empty)=1.

The inherited physical inequalities give

    E≥0, d−E≥0, D−d≥0, U−D−d+E≥0, D≥1.

Indeed delete a terminal root for E≤d, insert the root deterministically
into every deletion path for d≤D, and cut a Hamilton path before its root
for U−d≥D−E. Insertion also gives Hamilton-path existence. Consequently

    x=E r1+(d−E)r2+(D−d)r3+(U−D−d+E)r4,

    r1=(1,1,1,1), r2=(1,0,2,1),
    r3=(1,0,1,0), r4=(0,0,1,0),                   (1)

with nonnegative coefficients. In particular U≥D≥1. Singletons and sink
roots are included: a sink root has x=(D,D,D,D), a multiple of r1.

For a fixed nonempty rooted gadget (B,r), G_B(A,a) retains both tournaments,
directs B→a and A−a→B, and marks r. Write H=H(B); Q for ordered two-bin
path covers of B with empty bins allowed; s for Hamilton paths starting
at r; alpha for two-bin covers whose first bin is nonempty and starts at
r; and T for three-bin covers whose first bin is nonempty and starts at
r and whose other bins may be empty. The inherited first and third rows
of the graft law, valid also for singleton B, are

    D'=H D−H E+(Q−H)U+(2H−Q)d,
    U'=alpha D+(s−alpha)E+(2T−alpha)U
          +(3alpha−s−2T)d.                         (2)

The fourth row needs the inherited singleton exception, but neither the
ratio argument nor (2) uses it. Since H≥1 and the two placements of a
Hamilton path in one bin give Q≥2H, the following denominators are strictly
positive:

    D'(r1),D'(r2),D'(r3),D'(r4)=(H,H+Q,Q,Q−H).
    U'(r1),U'(r2),U'(r3),U'(r4)
       =(2alpha, 2alpha+2T−s, 2T, 2T−alpha).

**Ratio lemma.** For every physical input and every nonempty rooted gadget,

    U'≤kappa_B D',

    kappa_B=max{2alpha/H, (2alpha+2T−s)/(H+Q),
                2T/Q, (2T−alpha)/(Q−H)}.           (3)

**Proof.** Substitute (1) into the two linear functionals. The quotient is
a weighted average of the four displayed ray ratios, using weights equal
to the nonnegative cone coefficients times the positive D' values. The
input is nonzero, so the total weight is positive. This proves (3).

Moreover kappa_B is the **optimal** uniform ratio over all physical inputs.
The exact closed-cone result shows that r1 and r2 are physical, while r3
and r4 are limits of positive multiples of physical states. Explicitly,
the singleton and source-rooted edge realize r1 and r2; the tournament
X→u→v→X with X transitive of size m, rooted at u, has state
(2^m+1,1,2^m+1,1), which tends projectively to r3; and source-rooted
transitive n-tournaments have state (1,0,2^{n−1},1), tending projectively to
r4. Homogeneity and positivity of each limiting denominator give every
ray ratio as an attained value or limit. Hence their maximum is the
supremum. It need not be attained by a finite physical input. ∎

The inherited operations are

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D−E+U,d,2U−d,d),
    J(D,E,U,d)=(D−E+U,D−E,D+U,D).

For Z the D and U evaluations on the four rays are (1,3,2,1) and
(1,3,2,2), so its optimal constant is 2. For singleton J, formula (3)
also gives 2. For the source-rooted edge L it gives 4; for a rooted directed
triangle C3 and the inherited B4=[8,1,3,6], root 0, it gives 2. Thus the
menu O/Z/J/L/C3/B4 has kappa=4 below.

There is no finite ratio constant for O: it has D'(r4)=0 and U'(r4)=2,
and physical transitive sources make U'/D' unbounded. This is why O tails
must be separated rather than charged to the finite-core estimate. There
is likewise no menu-independent graft constant: a transitive n-vertex
gadget rooted at its source, applied to the singleton input, has D'=1
and U'=2^n. The finite-menu assumption is substantive for this proof.

Finally every operation here retains its old tournament as an induced
subtournament, so D never decreases. Deterministically inserting added
vertices into Hamilton paths gives an injection; deleting the added
vertices recovers the original path. O in fact preserves D exactly.

## 3. Explicit bounded-intercept theorem

Fix a finite physical seed bank S, a finite rooted-graft menu M, and an
integer C≥1. Let

    kappa=max({2} union {kappa_B:B in M}),
    S_C={s in S:D_s≤C},
    B_C=C floor(kappa C)+|S_C|.                     (4)

Words may be arbitrary finite words in O,Z and G_B for B in M, chosen
adaptively if desired. Given such a word W, form its endpoint by applying
O then J, and apply one additional J to obtain the closing tournament.
For the state (D,E,U,d) immediately after W,

    h=D+2U, g=2U, e=h−g=D,
    H(closing tournament)=N=h+g=D+4U.              (5)

Indeed the state after OJ is (D+2U,D,D+2U,D), and its next J has Hamilton
count 2(D+2U)−D. In particular g is even, h>g>0 and g/h≥2/3.

**Theorem.** All such closing counts with e≤C belong to the union of the
dyadic sequences

    D+4u·2^k, k≥0,                               (6)

over the following finite list of core pairs:

    1≤D≤C, 1≤u≤floor(kappa C),
    together with (D_s,U_s), s in S_C.             (7)

There are at most B_C pairs in this envelope. For every real X≥1 the number
of distinct closing counts ≤X with e≤C is at most

    B_C (1+floor(log_2 X)).                        (8)

More precisely a core pair contributes zero if X<D+4u, and otherwise at
most 1+floor(log_2((X−D)/(4u))). For C<1 there are no such endpoints.

**Proof.** If W is O^k, its pair is (D_s,2^k U_s), giving the seed case.
Otherwise split just after its last non-O operation:

    W=P F O^k, F in {Z} union {G_B:B in M}.

Let (D_1,E_1,U_1,d_1) be the state immediately after F. Trailing O preserves
D, so D_1=e≤C. The ratio lemma gives 1≤U_1≤kappa D_1≤kappa C. These are
integers, yielding a core pair from (7). Trailing O multiplies U_1 by 2^k,
and (5) proves (6). Each pair contributes at most (8)'s per-pair factor,
since D,u≥1. Counting pairs proves the bound. ∎

Monotonicity additionally makes D≤C a safe pruning test on every prefix,
but the last-operation proof does not require a bound on word length or
the number of non-O operations. A sink-root seed has Zx=x, so arbitrarily
long stationary Z runs really occur. Singleton seeds and these runs are
not discarded. O-only seed tails must be retained separately: a finite
seed bank may have U_s>kappa C even when D_s≤C.

Thus every fixed finite set of boundary lines h−g=c contributes
O_{S,M,C}(log(X+2)) distinct counts, taking C to be the largest positive
intercept. Nonpositive boundary intercepts contribute nothing.

## 4. Inherited exact residue-gap obstruction, now for the larger language

Fix S,M as above. Restrict endpoints to a finite bank of strict-slope lines

    a h−b g=c, gcd(a,b)=1, 0<a<b, c an integer,

together with finitely many boundary lines h−g=c with integer c and finitely many
exceptional scalar counts. Remove the strict-slope lines with a/b<2/3.
Let L be the remaining lines and let

    m=lcm{2(a+b):(a,b,c) in L},

using m=2 when L is empty. **If an odd residue r modulo m satisfies**

    a r not congruent to c modulo 2(a+b)

for every line in L, then these selected closing counts omit infinitely
many positive odd integers. More quantitatively, at least

    X/m−O_{S,M,line bank}(log(X+2))                (9)

members of the progression N≡r mod m through X are missing.

**Proof, with explicit finite thresholds.** All endpoints obey g≥2h/3.
For a discarded line put t=2b−3a>0. Its equation gives

    c≤−t h/3, so h≤H_c=floor(−3c/t).

If H_c<1 it contributes no endpoint; otherwise its counts satisfy
1≤N=h+g≤2H_c−1. Thus every discarded line contributes only finitely many
counts, with the explicit ceiling 2H_c−1. On a retained line,

    g=(aN−c)/(a+b), h=(bN+c)/(a+b).

The even-g requirement is exactly aN≡c mod 2(a+b) for integer candidates.
It is necessary for physical reachability, not sufficient. Therefore the
chosen progression misses all retained strict-slope lines. Section 3
bounds boundary-line contributions by B_C(1+floor(log_2 X)), or zero if
there is no positive boundary intercept. Add the finite discarded-line
ceiling and the number of exceptions, and count progression members to
obtain (9). ∎

For an explicit lower bound, use the representative r_0 in {1,...,m} and
P(X)=max(0,1+floor((X−r_0)/m)). Subtract from P(X) the bound in (8), the
number of exceptional counts, and

    sum_discarded max(0,2 floor(−3c/(2b−3a))−1).

This may overcount exclusions, which is safe. It gives a concrete
threshold test for every specified seed/menu/line bank.

As in Session226, r=1 works if c is not congruent to a modulo 2(a+b) on
every retained line, in particular whenever |c|<a. This includes zero and
nonunit intercepts. The new ingredient is exclusively the larger allowed
interior language and its ratio-controlled finite core.

## 5. Hostiles and boundaries of the conclusion

The inherited five lines 2h−3g=c for c in {0,−2,−4,−6,−8} cover all odd
scalar residues: for every odd N exactly one such c satisfies 2N≡c mod
10. Their candidates have g/h≥2/3 because c≤0. For positive odd N≥17,
h=(3N+c)/5, g=(2N−c)/5 and h−g=(N+2c)/5 are positive for the matching
line. Thus for all odd N≥17
there is a positive integer chamber candidate with even g on this bank.
This is an arithmetic hostile, not a reachability claim. A common residue
gap can genuinely fail; the present proof does not exclude such banks.

The result does not bound all counts from the enlarged grammar. It bounds
only bounded e and the indicated residue-gap line selections. Allowing
gadget menus or seed banks to vary with X is outside the fixed-constant
statement. Strong connectivity is not asserted for arbitrary seeds, and
state escape, full rank, or hull membership is not scalar completeness.

## 6. Portable verification

The standalone standard-library `verify.py` uses explicit RuntimeError
guards rather than assert, so checks remain enabled under optimization.
It independently obtains gadget statistics by permutation/bin enumeration,
compares subset-DP states, evaluates all four ratios, checks literal graph
grafts and O/Z operations, and exercises every word through a declared
finite depth over a mixed fixed menu. Singleton, sink-root, large-U seed,
large-O-tail, stationary-Z, and no-common-residue-gap controls are included.

The precise finite universe and totals are recorded in CHECK.json. The
same checker runs without repository imports or third-party dependencies:

    python3 -I -B verify.py > replay.json
    cmp replay.json CHECK.json
    python3 -I -B -O verify.py > replay.optimized.json
    cmp replay.optimized.json CHECK.optimized.json

Finite checks support the formulas and controls; they do not replace the
proof of arbitrary gadget order or unbounded word length.
