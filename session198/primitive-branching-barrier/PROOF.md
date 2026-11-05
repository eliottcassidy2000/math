# Primitive sources still require target-preserving rebasing

Session 198. **PROVED:** finite homogeneous source-extension menus cannot
give nonzero-offset affine children for every supplied pair of primitive
sources. Moreover a specific primitive representation of H=41 cannot give
either binary child 81 or 83 by source extensions of any length, even after
switching freely among all ten C2/C3 bilinear families.

**Scope:** these are obstructions to extending a supplied representation.
They are not scalar nonrepresentability results. Both 81 and 83 are
physically realized below after changing the H=41 source presentation.
Universal odd-target coverage remains open.

## 1. Inherited physical coordinates

Use the audited complete positive source grammar in
`dependencies/POSITIVE_SOURCE_GRAMMAR.md`. Its nontrivial source seed is
t=(D,E,U,d)=(1,0,2,1). Every word in O,Z,J is physical, and every nontrivial
source has such a word. The positive coordinates and their inverse are

    x=(a,f,e,p)=(D-E-d,U-D+E-d,E,d-E),
    (D,E,U,d)=(a+2e+p,e,a+f+2e+2p,e+p).

The physical transformations are

    O(a,f,e,p)=(a+e,a+2f+e+2p,0,e+p),
    Z(a,f,e,p)=(2a+f+e+p,f,e+p,0),
    J(a,f,e,p)=(f+p,0,a+e+p,e).

All source coordinates are nonnegative, and e+p=d>=1. A source is primitive
when gcd(D,E,U,d)=1. In this note the stronger gcd(D,U)=1 will hold in the
family used to prove the main obstruction. Thus the obstruction does not
come from preservation of a common divisor.

The ten complete two-nontrivial-source readouts are N,M,A,B,C,T,J0,J1,J2,J3,
with strictly positive integer matrices P_F in the frozen dependency. Their
physical source domains are all satisfied by nontrivial inputs. The C3
singleton boundaries are not operations obtained by extending two supplied
nontrivial sources, so they do not enter the extension menus below.

## 2. A three-parameter family of primitive actual sources

Starting at t, take the source word

    O^i J O^j Z J O^k,    i>=0, j>=1, k>=1.

Write z=2^(i+1)+1, y=2^j, q=2^k. Its response is

    D=(3y+1)z-2,
    E=0,
    U=q((3y+1)z-1)=q(D+1),
    d=(y+1)z.                                      (1)

Proof by successive actual operations:

    after O^i J:       (z,1,z,1),
    after O^j:         (z,0,yz,1),
    after Z:           ((y+1)z,1,2yz-1,1),
    after J:           ((3y+1)z-2,(y+1)z-1,(3y+1)z-1,(y+1)z),
    after O^k:         the response in (1).

D is odd, q is a power of two, and gcd(D,D+1)=1. Consequently gcd(D,U)=1.
Every member is therefore primitive, with an explicit physical graph of
order i+j+k+5.

This family is polynomially dense in the three-dimensional weak-source
coordinate plane E=0. Precisely: if a polynomial in D,U,d vanishes on every
member of (1), it is the zero polynomial. To prove this, substitute (1).
The three parameters z,y,q each range independently over an infinite set,
so a polynomial vanishing on the product of these sets is identically zero,
by induction on the number of variables and the elementary one-variable
root bound. The substitution is injective on polynomials because it has a
rational inverse on a nonempty open set:

    z=(3d-D-2)/2,  y=d/z-1,  q=U/(D+1).             (2)

Equivalently, the primitive physical source family is Zariski-dense in E=0.
The elementary polynomial argument, not algebraic-geometric terminology, is
what is used below. It does not say that every integer point of that plane is
physical.

## 3. No finite affine-child menu works on every primitive pair

Fix an input bilinear readout F(x,w)=x^T P w. Fix any finite list of m rules.
Rule r consists of fixed source words L_r,R_r, a chosen output bilinear family
P_r, an integer multiplier b_r, and a nonzero integer offset c_r. Interchanging
the two input sources is allowed. The proposed rule produces

    G_r(x,w)=(L_r x)^T P_r(R_r w)

or its source-swapped version, and asks for

    G_r(x,w)=b_r F(x,w)+c_r.                         (3)

The words may have arbitrary finite lengths, and branch selection may inspect
the entire supplied source pair. They are simply a fixed finite menu.

**Theorem 1.** There is a pair of primitive actual sources from (1) for which
none of the m equalities (3) holds.

**Proof.** On the product of the two coordinate planes E=0, the difference in
(3) is a polynomial

    H_r(x,w)=x^T Q_r w-c_r,

where Q_r is a fixed rational (in fact integer) matrix. This polynomial is
not zero: at x=w=0 its value is -c_r. The product of the m nonzero polynomials
is nonzero because a polynomial ring over the rationals has no zero
divisors. By the density statement of Section 2, applied to both independent
physical source families, this product cannot vanish on every physical pair.
At a pair where it is nonzero, no rule holds. All chosen sources have
gcd(D,U)=1. QED.

There is an explicit finite witness bound. Each source coordinate in (1) is
of degree at most one in each of z,y,q separately. Thus each pulled-back H_r
has degree at most one in each of the six parameters, and their product has
degree at most m in each parameter. Test the finite grid

    i in {0,...,m}, j,k in {1,...,m+1}

independently for each source. It contains m+1 distinct values of each
parameter. Repeated univariate interpolation shows that a nonzero polynomial
with these degree bounds cannot vanish throughout that grid. Therefore a
witness is guaranteed with both source orders at most 3m+7. This is a
finite constructive bound, not a claim that exhaustive grid search is fast.

The same reasoning permits fixed rational homogeneous linear source maps
in place of words, provided a proposed branch is used only where it is
physically valid. The obstruction concerns polynomial equality before that
additional validity restriction. It does not cover arbitrary affine source
maps with constant terms or freely chosen new representations.

## 4. For binary children, all extension lengths reduce to 71 tests

For this section the parent value is n=N(x,w), so n>=13. Consider every output
family among the ten bilinear readouts and arbitrary source-word extensions
on both sides. A prospective binary child is 2n-1 or 2n+1.

The following exact inequalities follow entrywise from the positive matrices
in the inherited proof. The verifier reconstructs the matrices from the
original response laws before checking the inequalities.

1. Every N or M source letter multiplies that readout by at least 4/3.
2. P_M >= (11/6)P_N.
3. For A,B,C,T, P_F-2P_N is nonnegative. Its submatrix on coordinates e,p
   has minimum entry respectively 3,4,41,4. Since e+p>=1 in both inputs,
   every such F is at least 2n+3.
4. For J0,J1,J2,J3, P_F >= (3/2)P_N. Every O or Z letter on either source
   multiplies F by at least 3/2. Every J letter is nondecreasing.
5. Three J letters on either source multiply each of J0,J1,J2,J3 by at least
   two. All ten readouts are nondecreasing under every source letter.

These give a complete, all-word reduction:

* For N, at least three letters in total give (64/27)n>2n+1. Hence only
  |L|+|R|<=2 can work: 34 ordered word pairs.
* For M, any nonempty extension gives at least (22/9)n>2n+1. Only the empty
  pair can work: one test.
* For A,B,C,T, even the empty extensions already exceed 2n+1. No tests remain.
* For J0,J1,J2,J3, any O or Z gives at least (9/4)n>2n+1. Thus both words
  consist only of J. Three J letters on one side give at least 3n>2n+1.
  Only J^a,J^b with 0<=a,b<=2 remain: nine tests per family.

In total there are 34+1+4*9=71 tests, some with identical matrices. No word
length cutoff is guessed: the growth inequalities exclude every other word.
Source swaps add no candidates: the ten matrices are closed under transpose
(B and T exchange, J0 and J3 exchange, J1 and J2 exchange, and the others are
symmetric), and the two ordered source-word slots are already both tested.

## 5. A primitive bad representation, and physical rebasing repairs

Take

    v=(5,1,5,1),  t=(1,0,2,1).

Both are primitive actual sources: v=JO(t), where the operations are read
as functions, so its source word is O then J. Their parent value is

    N(v,t)=41.

The 71 complete tests in Section 4 yield exactly these 45 distinct values:

    41,69,73,75,77,79,85,93,99,115,123,125,129,131,133,
    135,137,141,143,147,149,153,155,167,225,249,257,263,
    267,273,283,287,289,303,435,441,467,469,479,499,501,
    819,855,885,897.

Neither 81 nor 83 appears. By the all-word reduction, no extension of these
supplied source profiles, of any length, can produce either binary child
through any of the ten bilinear readouts.

This is not a scalar gap. Let h=(3,1,3,1) be the triangle source. Then

    v=Z(h),    v' = J(h)=(5,2,6,3),
    N(v,t)=N(v',t)=N(h,h)=41.                       (4)

Thus one can rebase the same scalar 41 to different actual source pairs.
The children are physically realized by

    81=N(O(v'),t),
    83=J1(h,h).                                    (5)

All sources in (4)-(5) are actual, and direct graph counts verify the stated
Hamilton-path values and cycle classes. Z(h) and J(h) use the same unrooted
four-vertex tournament with different roots; the other equality in (4)
replaces the source pair itself. These changes are exactly the additional
information a representation-local extension strategy had omitted.

The simplest twin rebasing has an exact, narrow physical domain. For an
actual source x=(a,f,e,p),

    N(Zx,t)-N(Jx,t)=4a+6f-4e.                       (6)

This vanishes exactly when x=t or x is a physically realizable positive
multiple of the triangle vector (1,0,1,0). In particular the only primitive
possibilities are t and the triangle.

Here is a last-letter proof, rather than a finite pattern extrapolation.
The seed t has zero defect. An O output has e=0 and f>0, so its defect is
positive. For a Z output, substituting its predecessor (a,f,e,p) into (6)
gives 8a+10f. Equality forces predecessor a=f=0, and the output is
(e+p,0,e+p,0), a positive triangle multiple. For a J output the defect is
4(f-a-e). A predecessor with e>0 has f<=a: its final letter is Z or J,
and the two formulas give this inequality directly. Such a predecessor
cannot give zero defect. If e=0, equality forces f=a, and the output is
(a+p,0,a+p,0), again a positive triangle multiple. Conversely every such
actual multiple has zero defect by (6). This proves the stated domain; it
does not assert that every positive integer multiple is a physical source.

Thus the H=41 repair is a real but very thin rebasing mechanism, rather than
a uniform procedure already known to work at arbitrary primitive inputs.

## 6. Scope and the precise remaining positive problem

The finite-menu theorem is not a noncoverage theorem. It quantifies over
every supplied primitive pair, whereas a universal scalar construction may
choose a special presentation for each target. The concrete H=41 example
shows why these quantifiers differ: its first pair is stuck, while two
alternative pairs supply its two children.

The positive problem is now more specific: construct a target-preserving
physical rebasing system that always reaches a presentation with the needed
inverse or child branch. Finite checks of rebasing success, modular density,
and the existence of a decreasing physical membership inverse do not prove
that this system covers every odd target. Nothing here claims that they do.

## 7. Reproduction

Run `python3 verify.py` and `python3 -O verify.py`. Both use explicit failure
checks. They verify the primitive family identities and small literal graphs,
the rational parameter inverse, all matrix bounds underlying the 71-case
reduction, the complete bad-pair output list, and literal graph witnesses for
the rebased parent and children. Finite tests support the displayed proofs;
they are not used to infer density or all-length exclusion.
