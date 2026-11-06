# The exact append-only signed-section orbit

**PROVED, independently audited in the stated scope.** For each supplied actual source, this
note classifies every signed target that can reach the physical 6/7 section
by appending physical letters. It gives a sharp minimum and an infinite
family of already represented positive targets for which append-only repair
of a supplied signed witness is impossible.

**OPEN:** a target-dependent replacement or deletion-and-reconstruction
that supplies a section witness at each desired target. Signed surjectivity
and the classification below do not prove that existence statement.

## 1. Domain and the inherited physical section

Use the actual nontrivial source orbit S of t=(1,0,2,1), with chronological
operations

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D-E+U,d,2U-d,d),
    J(D,E,U,d)=(D-E+U,D-E,D+U,D).

The nonnegative physical coordinates are

    (a,f,e,p)=(D-E-d,U-D+E-d,E,d-E).

Thus D=a+2e+p, U=a+f+2e+2p and d=e+p>=1. Both D and d are
nondecreasing along physical words. The f coordinate is preserved by Z,
is set to zero by J, and becomes strictly positive under O.

Write L(q)=6U_q-7D_q and let T be the actual physical slice

    E_q=2(U_q-D_q), d_q=3(U_q-D_q).

The independently audited dependency SIGNED-SECTION.md proves the following
facts used here. Every q in T has positive target n=L(q). It is either:

* canonical: q=r.JJJ with an actual balanced r=(m,n,m,n), m>=3n;
  every presentation of q has these last three J letters;
* exceptional: q=r.J with actual balanced r=(3n,n,3n,n).

The canonical statement forces all three predecessor letters before filling
in their deletion coordinates. Its intermediate predecessors q_1 and q_2
both have f=0. The exceptional predecessor r is balanced and hence has f=0.
Also, on every actual source v,

    W(v)=v.OJJJJ belongs to T, L(W(v))=D_v,
    U(W(v))=6D_v+14U_v.                         (1)

For a fixed actual source s, let Desc(s) include s and all its physical word
extensions. Define

    Dtail(s)={D_v:v in Desc(s)},
    A(s)={L(q):q in Desc(s) intersect T},
    C(s)={L(s.J^j):0<=j<=2 and s.J^j belongs to T}.

All these sets retain the actual physical source s; no cone point is
substituted for one. Membership in C(s) is a direct test of three actual
short extensions.

## 2. Balanced descendants: one exact scalar invariant

**Lemma 1.** The endpoint counts of balanced descendants of s are exactly

    Dtail(s) union ({d_s} if f_s=0, otherwise the empty set).       (2)

**Proof.** Every v in Desc(s) gives the balanced descendant v.OJ with
endpoint count D_v. If f_s=0, s.Z is balanced with endpoint count d_s.

Conversely, let r be a balanced descendant. If its supplied extension word
contains J, take the last J. Its output deletion count is D_v for the
physical prefix v immediately before that J. Subsequent O and Z preserve
deletion, and balance gives E_r=d_r=D_v. If the extension contains no J,
then it contains no O: an O makes f positive and subsequent Z cannot
remove it. Hence f_s=0 and its unchanged deletion is d_s; balance again
gives E_r=d_s. This includes an empty extension. QED.

## 3. Complete target classification and a sharp minimum

**Theorem 2.** For every actual source s,

    A(s)=Dtail(s) union ({d_s} if f_s=0) union C(s).                 (3)

**Proof.** Equation (1) realizes all of Dtail(s). For any balanced
descendant r, direct iteration gives r.JJJ in T with target E_r. Lemma1
therefore realizes the additional d_s when applicable. The three short
completions are actual by definition.

For the reverse inclusion take q=s.gamma in T. If q is exceptional and
gamma is nonempty, its last letter is J and the preceding profile is a
balanced descendant of s. Its endpoint is the target n. Apply Lemma1.
If gamma is empty, q is recorded in C(s).

If q is canonical and gamma has length at least three, its forced last
three J letters leave a balanced descendant of s, again with endpoint n.
Apply Lemma1. If gamma has length at most two, every letter of gamma is J,
so it is exactly one of the short completions in C(s). QED.

This is an exact append-only criterion, not a claim that Dtail(s) is all
large odd integers.

**Corollary 3 (sharp lower frontier).** If f_s>0, then

    A(s)=Dtail(s), min A(s)=D_s.                                  (4)

If f_s=0, then

    min A(s)=min({d_s} union C(s)).                              (5)

In f=0 coordinates (a,e,p), the short-completion tests simplify to

    j=0: e=2p,      target p-a;
    j=1: e=a+p,    target a;
    j=2: a=0,      target p.                                    (6)

Every listed target is positive when the tested profile is actual.

**Proof.** First suppose f_s>0. The j=0 case is impossible because T has
f=0. For j=1, the predecessor of either type of slice point has f=0, so
that case is impossible. For j=2, a canonical slice point has f=0 at its
second predecessor, again impossible. If s.JJ is exceptional, s.J is
balanced. The identity U(Js)-D(Js)=E_s forces E_s=0, and its endpoint is
D_s-E_s=D_s. Thus C(s) is contained in {D_s}. Monotonicity of D and (1)
prove (4), including sharpness.

If f_s=0, D_s>=d_s and all descendant D-values are at least D_s.
Theorem2 gives (5). The raw slice equations after J^j are, respectively,

    j=0: 2D+E-2U=0 and 3D-3U+d=0;
    j=1: D-3E=0;
    j=2: -2D+2E+U=0.

Their signed readouts are L, -D+7E-U and 5D-6E-2U.
Substitution of f=0 coordinates gives exactly (6). The inherited physical
section theorem supplies positivity; the corresponding hyperplanes in the
whole integer lattice do not. QED.

## 4. An infinite hostile family at already represented targets

For r>=0 set lambda=3^r. The actual source

    s_r=t.Z.(OJ)^r.OZ=lambda*(9,1,11,1)

has f=2lambda>0 and L(s_r)=3lambda<D_s=9lambda. Indeed t.Z=b and
b.OJ=3b, so the formula follows by iteration. Corollary3 shows that no
physical suffix gamma can make s_r.gamma lie in T while retaining the
signed target 3lambda.

These targets are already represented: t.Z.(OJ)^r is balanced with
D=3lambda, and its W extension is a section witness with the same target.
The smallest hostile is

    s_0=t.ZOZ=(9,1,11,1), L(s_0)=3.

Thus the obstruction is not the known missing target7, nor any question
about whether3 itself is representable. It is a failure of append-only
repair from a supplied signed witness. A successful universal construction
must have permission to replace part of that witness, delete and rebuild,
or choose a materially different signed representative.

## 5. A finite exact decision bound for a fixed supplied source

For completeness, (3) can be tested without searching arbitrarily long
extensions. For a prescribed positive integer n, Dtail(s) membership is
false if n<D_s; true if n=D_s; and for n>D_s it suffices to enumerate
physical descendants through U<=max(U_s,2n-3).

To prove the last claim, take any extension having D=n and strip only its
terminal O letters. Because n>D_s, at least one new Z/J remains. The
resulting positive-endpoint core is still a descendant of s and has D=n.
A Z output satisfies 2D-U>=3d>=3. For a J output,

    2D_new-U_new=D_old-2E_old+U_old>=3:

if U_old>=3 this follows from D_old>=2E_old; the sole U_old=2 source is t,
where the expression is3. Therefore every such core has U<=2n-3.
The converse follows from keeping actual words.

In particular every target n<=N in A(s) has a section descendant of cover
height at most

    B(s,N)=max(34N-42,14U_s+6N), N>=1.                           (7)

For Dtail(s), use W(s) when n=D_s, and use the core above otherwise.
The additional d_s uses s.Z.JJJ, of height14U_s-8d_s. The short completions
have height at most4U_s. These bounds also cover N=1 via the second term.
Equation (7) is a complete finite decision bound for a supplied prefix,
not an existence theorem for any unsupplied target.

## 6. Verification and the remaining constructive obligation

The verifier independently evaluates both sides of (3) for declared source
seeds and target cutoffs, using (7) to make the forward section search
complete. It checks the short faces on complete physical banks, the sharp
minimum in both f chambers, the infinite hostile formulas over a declared
range, and formal-basis matrix identities. Normal and optimized Python
runs must agree.

The remaining constructive issue can now be stated more narrowly: how can
one select or alter a globally surjective signed witness so that its target
belongs to the exact tail set (3)? Merely appending a sufficiently long
physical correction word cannot solve this, even for infinitely many
already represented positive targets. No all-target coverage is claimed.
