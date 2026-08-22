# The ternary tree is a circuit, a Walsh transform, and several different clocks

**Status:** research synthesis and next-operation map.  The proof source is
[THM-3357](../01-canon/theorems/THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit.md),
with the Gaussian/Lorentz carrier in
[THM-3333](../01-canon/theorems/THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md),
the square-triangular compiler in
[THM-3335](../01-canon/theorems/THM-3335-square-triangular-pell-markov-pythagorean-selector.md),
and the ancestry warning in
[THM-3345](../01-canon/theorems/THM-3345-prime-xor-ancestry-path-groupoid-and-source-dependent-berggren-cost.md).
This reflection is not an additional proof source.  `LRC(14)`, the planar
Jacobian conjecture, `FC(3)`, and skew-EW design existence remain open.

## Outcome first

The useful object was not a fourth distinguished spine.  It was the whole
three-colour branching operator viewed at four resolutions:

```text
one parent u
  |
  +-- labelled siblings Lu,Mu,Ru
  |      determinant circuit -> parent Pythagorean triple
  |      norm order          -> one transitive T3
  |      support-norm circuit-> LRC determinant Horn rule
  |      common gcd          -> one prime-five edge clock
  |
  +-- all 3^d descendants
         weighted first moment -> one 2x2 transfer
         four parity characters -> M^(2d),L^(2d),R^(2d),I
         quadratic moment       -> 3x3 transfer and shell energy
         transfer determinant   -> exact FC(3) factorial moments
```

This changed the interpretation of the old Pell coincidences.  The
square-triangular/Markov state is the **first moment of an entire ternary
level**, not a property shared by its nodes.  At depth two the nine node
parameters sum to `(29,70)`, so the cannonball `70` is literally a global
level mass.  The nonlinear Gaussian image, the sum of nodewise triples, and
the consecutive square-leg selector triple are three different objects.

The strongest new gains were:

1. a positive sibling circuit whose three determinants recover the parent;
2. a four-character Walsh collapse of every ternary level to all three unary
   rays simultaneously;
3. a strict LRC determinant-gate implication from the two outer children to
   the middle child;
4. a local prime-five shared-hypotenuse clock with limiting parent density
   `1/3`;
5. exact full-level first- and second-moment recurrences, including separate
   even- and odd-parameter square totals; and
6. an FC(3)-facing cubic whose odd factorial moments vanish and whose even
   moments have one closed formula.

## Anchor / Niche / Wildcard portfolio

| lane | inherited object | what moved | honest boundary |
|---|---|---|---|
| **Anchor — LRC(14)** | THM-2056 determinant support norm | the sibling circuit makes outer gate-pass imply strict middle gate-pass | no seeds, clocks, phase height, or open row closure |
| **Niche — all three Berggren branches** | THM-3334 unary spine and THM-3339 three selected rays | a weighted whole-tree transfer and four parity-sector formulas | Parikh data forget word order and source path |
| **Wildcard — primes / FC(3) / square sums** | Gaussian norms, factorial functional, cannonball `70` | prime-five collision DFA, factorial cubic moments, two typed even/odd square splits | none proves a prime torsor, FC(3), or design theorem |

The anchor did not overtake the session.  The niche produced both a genuine
gate-pruning theorem and a much sharper explanation of why aggregate
statistics cannot supply LRC phase information.

## Inheritance pass

| required item | controlled source |
|---|---|
| closest proved mechanism | THM-3333: Gaussian square and Lorentz/minor polarization |
| native ternary operation | THM-3339/3353: the three exact parameter branch matrices |
| canonical hostile | THM-2596: three infinite-order branches are not a `C3` action |
| corrected near miss | three Fibonacci rays select a regular language; they do not classify arbitrary tree words |
| ancestry hostile | THM-3345: a reduced tree path is source-dependent even when an external prime colour is fixed |
| least-used sidecar | branch-count parity and the signed sibling circuit, rather than one scalar norm |

The outer Hadamard reflection was retained only in its corrected form.  It has
determinant `-2`, takes opposite-parity spinors to odd/odd spinors, and
requires raw-content-two normalization.  It is a useful projective reflection,
not an integral tree automorphism.

## Live concept board

| object | representation | preserved invariant | operation | missing coordinate / hostile |
|---|---|---|---|---|
| sibling packet | `(Lu,Mu,Ru)` | three signed minors `(b,c,a)` | take all three children | absolute shells lose circuit signs |
| norm tournament | `(N_L,N_M,N_R)` | middle maximum and outer sign | reset/flip/reset automaton | order loses circuit magnitude and LRC phase |
| weighted level | `(xL+yM+zR)^d u` | all Parikh-weighted first moments | matrix powering | `LM` and `ML` differ |
| parity level | four affine `F_2^3` sectors | count-parity sums | Walsh inversion | no ancestry order within a sector |
| Pell mass | `(P_(2d+1),P_(2d+2))` | Markov/Cassini/square-triangular clock | middle matrix squared | aggregate is not a selector node |
| quadratic level | `sum vv^T` or `sum Psi(v)` | all second moments and total shell energy | 3x3 transfer | even full first+second moments on same-level subpackets miss max shell |
| local gcd graph | child hypotenuse gcds | one possible shared factor `5` | projective mod-five clock | loses norm, height, word, cofactor |
| factorial cubic | determinant of weighted triple transfer | every factorial moment | simplex/radial integral | even moment at `r=2` kills FC candidacy |

## Pull 1: three siblings form the faithful local object

For `u=(m,n)`, `0<m<n`, put

```text
Psi(u)=(a,b,c)=(n^2-m^2,2mn,n^2+m^2).
```

The three children satisfy

```text
det(Lu,Mu)=b,       det(Lu,Ru)=c,       det(Mu,Ru)=a,
a Lu+b Ru=c Mu.                                          (1)
```

This is better than starting with a tournament.  It is a positive oriented
circuit, and it decodes the parent triple exactly.  Gaussian squaring sends
its minors to Lorentz shells `2b^2,2c^2,2a^2`; the full sibling triple has
volume `-4abc`.  The three norm values then give a lawful but deliberately
coarser tournament:

```text
N_M-N_L=4b,       N_M-N_R=4a,       N_R-N_L=4(b-a).     (2)
```

The middle child is always champion.  Only the outer arc moves, and its state
is reset by `L`, flipped by `M`, and reset oppositely by `R`.  The minimal
outer gap occurs exactly on the consecutive-leg middle ray.  This is a rare
case where tournament structure is intrinsic and useful, but the circuit is
still the more faithful carrier.

## Pull 2: four signed sums see all three spines at once

The matrix quartet

```text
L+M+R=M^2,       L+M-R=L^2,
-L+M+R=R^2,      L-M+R=I                              (3)
```

turns the exponential tree into four unary states.  At depth `d`, unsigned
sum lands on `M^(2d)u`; signs by `R` count, `L` count, and `M` count land on
`L^(2d)u`, `R^(2d)u`, and `u`.  Fourier inversion then recovers every
branch-count parity sector.

At the root,

```text
M^(2d)u_0=(P_(2d+1),P_(2d+2)),
L^(2d)u_0=(2d+1,2d+2),       R^(2d)u_0=(1,4d+2).        (4)
```

So the three spines contribute hyperbolic/Pell, left-parabolic, and
right-parabolic terms to one exact parity formula.  This is much more
informative than studying the spines independently.

The quotient remains commutative.  Words with the same branch counts can be
different nodes, and their unique ancestry paths can require unbounded source
memory.  The next refinement should therefore be a noncommutative rational
series or an ordered automaton, not another scalar recurrence.

## Pull 3: the square-triangular and square-sum clocks are level masses

Write

```text
alpha=P_(2d+1),       s=P_(2d+2),       beta=P_(2d+3).
```

Then

```text
alpha beta-s^2=1,
alpha^2+beta^2+4=6alpha beta,
n=(alpha+beta-2)/4,       q=s/2,       T_n=q^2.          (5)
```

This is precisely the THM-3335 clock at selector depth `d+1`, but its new
meaning is global: `(alpha,s)` is the sum of all `3^d` parameter nodes.
At `d=2`,

```text
sum_(nine nodes) u=(29,70),
70^2=sum_(r=1)^24 r^2=2600+2300.                        (6)
```

The last split is the even/odd split of the classical range `1,...,24`.
The actual full-tree parameter-square split is a different exact sequence.
If `C_d=sum(m^2+n^2)` over the level, then

```text
Eeven_d={C_d+(-1)^d(4*3^d-1)}/2,
Eodd_d ={C_d-(-1)^d(4*3^d-1)}/2.                        (7)
```

At depth two this is `364+329=693`, not `2600+2300=4900`.
The shared phrase “even and odd squares” had hidden two different index sets;
typing the range resolves the apparent coincidence.

## Pull 4: the circuit creates a Horn rule, not an LRC proof

For a fixed THM-2056 column deck let

```text
Delta(v)=max_i |det(v,c_i)|,       G(v)=||v||^2-kappa Delta(v).
```

Subadditivity applied to (1), plus the norm gaps, gives

```text
cG(Mu)>=aG(Lu)+bG(Ru)
        +2m(n-m)(3n^2+4mn-m^2).                         (8)
```

The correction is strictly positive.  Hence two outer gate-pass seeds imply
a strict middle gate pass, for every `kappa>=0`, including `kappa=91`.
Contrapositively, a bad middle child forces at least one bad outer child.

This is a directed Horn hyperedge.  It is not a pairwise tournament relation,
and it supplies no clocks or phase witness.  The same spinor/Kelvin scalar
state can occur in THM-3335's labelled tail planes with phase maxima
`1/7,...,1/12`; the missing coordinate is the labelled clock placement.

## Pull 5: a local prime clock singles out five

For primitive parents the child-hypotenuse gcd graph is empty or one edge:

```text
gcd(N_L,N_M)=5 iff 5|m,
gcd(N_M,N_R)=5 iff n=m mod5,
gcd(N_L,N_R)=1.                                         (9)
```

No higher power of five is common.  Projective slopes modulo five collapse,
for level counts, to three equitable classes.  The number `c_d` of parents
with a shared-five edge satisfies

```text
c_(d+3)=c_(d+2)+3c_(d+1)+9c_d,
c_d: 0,0,6,6,24,96,222,726,...,
c_d/3^d -> 1/3.                                        (10)
```

This is not the Gaussian prime-XOR collision torsor.  It is a local sibling
support event, has only one possible edge, and forgets exact hypotenuse and
ancestry.  The six raw projective residues are still needed for labelled
wordwise prediction; the three classes suffice only after summing outgoing
branch labels.

The natural next prime experiment is to keep the full projective automaton
for each `p=1 mod4`, then ask which outputs survive simultaneously with a
fixed determinant owner.  Incoming
[THM-3356](../01-canon/theorems/THM-3356-primitive-affine-determinant-shells-parabolic-orbits-and-prime-clock-resultants.md)
now proves a complementary clock: primitive affine determinant shells split
into parabolic residue/content orbits and their two U-spine channels have a
sharp rank-five resultant fingerprint on the current incoming rays.  Its
fixed shell/resultant state and the present local sibling-gcd state are not
the same quotient.  A useful next pull is their fibre product, retaining both
the shell residue and the child label.  Incoming
[THM-3358](../01-canon/theorems/THM-3358-admissible-composite-parabolic-compiler-and-hensel-normal-offset-atlas.md)
adds a second complementary operation: every exact-grade normal lift for an
admissible composite norm receives a fixed-prefix/unary/fixed-suffix Berggren
compiler.  This makes a composite version of the fibre-product question
concrete—can a sibling-output automaton retain the centered Hensel normal unit
without restoring the whole ancestry word?

## Pull 6: the all-branch determinant touches FC(3) exactly

The determinant of the weighted triple transfer is

```text
P=(x-y-z)(x+y-z)(x+y+z).                                (11)
```

For the factorial functional `calL(x^i y^j z^k)=i!j!k!`, simplex coordinates
turn the normalized integrand into `(uv)^r` on `-1<=u<=v<=1`.  Consequently

```text
calL(P^r)=0                         for odd r,
calL(P^r)=(3r+2)!/[2(r+1)^2]       for even r.           (12)
```

This is the cleanest FC(3) connection in the session.  Half of an infinite
moment sequence vanishes for a structural symmetry reason, while the other
half is explicitly positive.  The first detector is

```text
calL(P^2)=2240.
```

So `P` is a benchmark for parity-blind or bounded-depth factorial arguments,
not an FC counterexample.  The mechanism also suggests a procedural search:
factor a transfer determinant, move to simplex coordinates, and inspect
whether an involution kills one parity of moments before doing large symbolic
expansions.

## What the quadratic level state does and does not know

The sum of nodewise triples is one 3x3 recurrence, and Cauchy--Binet turns it
into the total pairwise determinant energy.  But even full first and second
moments can miss the largest shell.  Inside Berggren depth three, the
subpackets

```text
{LLM=(4,11),MLR=(5,18),MRL=(9,16)},
{LMM=(8,19),LMR=(3,14),LRL=(7,12)}
```

have the same vector sum and the same matrix `sum vv^T`, but determinant
patterns `(17,35,82)` and `(37,55,62)`.  They do not give two colliding
complete levels; they show that the moment tuple alone is insufficient even
on same-level subpackets.  Moment methods are suited to averages and energies;
LRC owner selection is an extremal, labelled problem.

## Procedure for another ternary tree

The reusable test is now concrete.

1. Put all three branch maps on one lattice and fix action order.
2. Compute the four signed branch sums.  If they close to simple powers,
   build the Walsh compiler; if not, stop that transfer.
3. Compute the three sibling minors.  A signed circuit is more faithful than
   an ad hoc tournament.
4. Separate first moment, quadratic moment, maximum shell, and labelled owner.
5. Transport the metric explicitly under conjugation.
6. Run a same-Parikh/different-word hostile and a same-moment/different-max
   hostile.
7. Only then connect a target frontier through the predicate it actually
   consumes.

Simultaneous conjugates of the Berggren quartet pass step 2.  An arbitrary
unimodular positive ternary triple need not.  The open classification problem
is to describe all integral positive-chamber triples satisfying the quartet,
and determine which are genuinely new trees rather than coordinate changes.

## Frontier ledger

| frontier | exact gain | missing coordinate / stopping certificate |
|---|---|---|
| LRC(14) | outer-to-middle determinant Horn rule | clocks, labelled plane, phase height, seed/exit |
| tournament structure | intrinsic transitive sibling `T3`, two-state outer arc | circuit weights; no global level tournament |
| sequences | order-2/3/4 recurrences and `O(log d)` fixed-weight evaluation | arbitrary coefficient extraction has quadratic support; word list exponential |
| primitive triples | parent decoded from sibling minors; all branches used | absolute shells lose orientation |
| square-triangular / EW | full-level mass compiles both typed order sequences | EW arithmetic necessity is not a design |
| square sums | cannonball range split and internal parameter split separated | different index sets, no general square-pyramid law |
| FC(3) | exact factorial moments of transfer determinant | even moments are positive |
| planar JC | concrete projections have nonconstant Jacobian | no Keller map or descent circuit |
| other ternary trees | quartet test and conjugacy-covariant Walsh port | arbitrary arity/unimodularity is insufficient |

## Cheapest next pulls

1. Replay the Horn rule on actual THM-2056 residual decks and count how often
   middle directions become redundant after outer certification.
2. Classify integral matrix quartets satisfying the four signed identities,
   first up to simultaneous `GL_2(Z)` conjugacy and branch reflection.
3. Replace the commutative Parikh transfer by a small noncommutative automaton
   that retains enough suffix order to talk to THM-3345 without storing the
   full word.
4. Build full `P^1(F_p)` sibling-output automata for split primes and compare
   their languages with fixed-hypotenuse Gaussian toggles.
5. Search determinants of other ternary transfer representations for simplex
   involutions giving exact factorial-moment formulas.

The session's main procedural lesson is that “use all three branches” should
not mean merging their scalar sequences.  It should mean finding the smallest
simultaneous representation in which addition, orientation, parity, and
nonlinearity can be audited separately.
