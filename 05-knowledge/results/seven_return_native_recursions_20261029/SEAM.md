# A universal seam transfer and its native-complement obstruction

**PROVED, scoped elementary statements; same-agent proof review and independent
finite algorithms.** No independent agent proof audit or literature-wide novelty
claim. The earlier two-block formula and the earlier 18-letter transfer are
inherited; the extension here is to arbitrary nonempty directed blocks, with
arbitrary entry and exit ports, including coincident ports.

## Inheritance and the operation being represented

The two-block bridge count is already in
[Session29 ADDENDUM, section 2](../session29_framed_profiles_20261022/session29-profiles/ADDENDUM.md).
Its [REPORT](../session29_framed_profiles_20261022/session29-profiles/REPORT.md)
has an exact three-state observer for the specified 18 five-vertex letters.
THM-4181 / [ordinal-sum capacity transfer](../../../01-canon/theorems/THM-4181-ordinal-sum-capacity-transfer-and-parity-component-exchange.md)
already gives the broader path-cover capacity principle. The native
complement is kept distinct, as in Session102
[SIGN_GRAMMAR](../session102_tournament_classification_20261029/metagraph/SIGN_GRAMMAR.md).
No graph-isomorphism or complement-closure conclusion is imported from a
scalar path count.

For disjoint blocks X,Y with exit a in X and entry b in Y, form X⊙Y by
retaining all internal arcs, putting every cross arc from Y to X, and
reversing just b→a to a→b. Retain the entry of X and exit of Y.
This associative operation is a tournament when the blocks are tournaments.
Internal blocks can more generally be any finite loopless directed graphs.
The graph, the chosen ports, and ordinary vertex count are not quotiented away.

For a nonempty graph X with exit a define the column vector

    v(X,a) = (H,S,R)^T,
    H = number of Hamiltonian paths of X,
    S = number of Hamiltonian paths starting at a,
    R = sum_{nonempty U containing a} E(X[U],a) H(X\U),
    H(empty)=1.

E counts paths ending at a. Thus R=E(X,a)+Q(X,a), where Q counts spanning
two-path covers whose distinguished path ends at a, the other path nonempty.
Paths are oriented vertex sequences; singleton paths are allowed. A two-path
cover is counted once with its distinguished path first.

## Universal reduction to three coordinates

Every block Y has an integer 3×3 matrix M_Y with

    v(X⊙Y, exit(Y)) = M_Y v(X,a).                         (1)

Here is a finite construction and a proof valid for every Y, rather than
verification of a preselected letter bank. Contract paths of the old prefix
to atoms, retain Y literally, and compute the three output counts. For one
old path use types I,S,E according as a is internal, its start, or its end.
Call the output columns f_I,f_S,f_E. For two old paths, the active one ends
at a; distinguish whether it is nonsingleton (DE) or singleton (DB). The
other old path is a second atom. There are no arcs between these atoms.
All Y-to-atom arcs are present except entry-to-active when a is its start;
active-to-entry is present exactly when a is its end. Let the resulting
columns be f_DE,f_DB. Type B denotes a single singleton old path.

These columns satisfy, for every directed Y and every pair of ports,

    f_DB=f_DE,       f_E-f_I=f_DE,
    f_B=f_E+f_S-f_I.                                    (2)

For the first equality, a counted output involving two old components must
use the unique active-to-entry seam. Otherwise H and S cannot connect the
two components; for R both would need to be terminal in separate paths,
leaving no path ending at the new exit. Once the seam is used, the only arc
that distinguishes DE from DB, entry-to-active, would make a directed
two-cycle and cannot be used in a path cover.

For the second equality, f_E−f_I counts precisely the covers using the new
active-to-entry arc. In an H or S cover, append the second old atom to its
last vertex, which is in Y. In an R cover with two paths, append that atom
to the path that does not end at the new exit. In an R cover with just one
path, add the atom as a singleton second path. Each construction is reversible
by deleting the second atom. All appended arcs exist. For the third equality,
adding active-to-entry and deleting entry-to-active have zero mixed term:
no counted cover can use both opposite arcs. This also covers equal entry
and exit ports in Y.

To justify completeness, an output Hamiltonian path restricts to at most
two old paths, since only one old-to-new seam exists. An output R-cover also
restricts to at most two old paths: one output path must end in Y, so at most
one old component can terminate without using that seam. Thus there is no
unrecorded three-component input term. Partitioning the old Hamiltonian
paths into I,S,E and the two-path covers into DE,DB, then using (2), gives

    M_Y = [ f_I | f_S-f_I | f_DE ].                       (3)

For singleton X, v=(1,1,1), and the last identity in (2) gives the same rule.
The dummy sink used for f_I also shows

    M_Y e_0=v(Y,exit(Y)),       e_0=(1,0,0)^T.             (4)

The construction is implemented in `seam_transfer.py`; its finite counts
use subset dynamic programming, while independent checks enumerate literal
permutations and path cuts.

## Composition, converse and exact linear dimension

Associativity of the actual seamed graph and (1) imply

    M_(X⊙Y)=M_Y M_X.                                    (5)

Equality holds as a matrix, not merely on one initial vector: the actual
prefixes consisting of the ordered edge, ordered triangle, and directed
triangle have exit responses (1,0,2), (1,0,4), (3,1,3). Their determinant is
−2, so they span Q^3.

Let X* reverse all arcs and swap the marked entry and exit. Put
J=diag(1,−1,1). The inherited two-block formula reads

    H(X⊙Y)=v(Y*,entry(Y))^T J v(X,exit(X)),

where `v(Y*,entry(Y))` uses the exit of the reversed marked block. Hence
e_0^T M_Y=v(Y*)^T J. Applying this to X⊙Y, using
(X⊙Y)*=Y*⊙X*, and using the spanning vectors, proves

    M_(X*)=J M_X^T J.                                   (6)

The reachable vectors span three dimensions and the reversed right contexts
give a nondegenerate bilinear pairing. The associated context-count matrix
therefore has rank three. Three is minimal **for linear representations of
Hamiltonian-path counts under arbitrary seam continuations over Q**.
This is not a lower bound on arbitrary nonlinear encodings, and not a
three-state reconstruction of tournament isomorphism classes.

For any finite alphabet of marked blocks, possibly of different orders,

    A(z)=sum_B z^|B| M_B,
    sum_words H(word) z^(total vertices) = e_0^T (I−A(z))^(-1) e_0. (7)

All block sizes are positive, so the inverse is valid formally. This gives
a rational, vertex-clock-preserving recurrence for the chosen composition
language. Words are counted as presentations, with multiplicity; no
unjustified quotient by isomorphism is taken.

The singleton has matrix

    [1 -1 1; 1 -1 0; 1 0 1],

with characteristic polynomial t^3−t^2−t−1. This recovers the old tribonacci
partner, not a new discovery of its recurrence: see THM-4731 /
[tiling-line metagraph and tribonacci partner](../../../01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md),
which explicitly credits the March5 proof and records MISTAKE-645 about
rediscovery.

## A five-vertex obstruction to native complement closure

Keep the Hamiltonian spine 0→1→2→3→4. Encode the other arcs by bits in
lexicographic order (02),(03),(04),(13),(14),(24), with bit one meaning
the smaller endpoint points to the larger. The marked tournaments with
bitmasks 19 and 35 have the identical full operator

    M = [13 -2 26; 3 0 6; 15 -2 40].

They are indistinguishable by H in **every** surrounding seam context,
including contexts on both sides. But the native complement, which reverses
every nonspine edge and retains every spine edge, gives respectively

    H(K_19)=11,       H(K_35)=13.                         (8)

Both counts are checked by independent literal permutation enumeration.
Their complete adjacency matrices and responses are in `seam_checks.json`.
Even the full nine entries of the transfer matrix fail to retain an observable
that native complement tests. Merely extending the alphabet of seam blocks
cannot repair this loss.

The strongest survivor is (1)–(7). The missing information is internal
marked-path response under complement. The cheapest next extension is to
retain both M_X and M_(K_P X) with the actual path P; the next decisive test
is closure after re-marking or attaching a new block and complementing again.
No finite-dimensional closure of that enlarged operation family is proved.
