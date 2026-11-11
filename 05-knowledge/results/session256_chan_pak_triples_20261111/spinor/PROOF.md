# A Berggren substitution lift into actual apex tournaments

Session256. New construction/proof, with small tests of these new claims only.
No claim of literature priority. This does not prove Chan–Pak or coverage by
this restricted family. The S251 five-core formula and its binary g convention
are inherited inputs, not results re-audited here.

## 1. The graph carrier and the exact spinor dictionary

Use Euclidean spinors m>n>0, gcd(m,n)=1, m−n odd, with triple
(m²−n²,2mn,m²+n²). The branches are

 A:(m,n)→(2m−n,m), B:(m,n)→(2m+n,m), C:(m,n)→(m+2n,n).

Set x=m−n,y=n. Their matrices on (x,y)^t are

 A=[[1,0],[1,1]], B=[[1,2],[1,1]], C=[[1,2],[0,1]].

A is the elementary left shear L; C=R² for R=[[1,1],[0,1]];
B=RSR for the coordinate swap S. The negative entry in the usual A spinor
matrix disappears in these positive coordinates, while B retains determinant
−1. All three preserve gcd and odd x, starting at (x,y)=(1,1).

For any word w in O,I, define T(w): the base vertices are totally ordered,
all base edges point forward, and the apex points to precisely the O vertices.
Every I vertex points to the apex. This specifies every arc.

Define substitutions

 a(O)=OI,    a(I)=I;
 b(O)=OI,    b(I)=OIO;
 c(O)=O,     c(I)=OOI.

Their incidence matrices, with rows #O,#I, are exactly A,B,C above.
Starting with w=OI and (m,n)=(2,1), apply the same sequence of substitutions
and spinor branches. Then exactly

 |w|=m, #I(w)=n, #O(w)=m−n.

Proof: the incidence matrices agree, including their multiplication order.
Induction proves all three equations. The standard spinor-to-triple formula
then produces a primitive Pythagorean triple from each graph's cut counts.

This is a literal graph operation: replace each base vertex by the ordered
block prescribed by its letter, orient all edges between blocks in the old
base order, and reassign apex edges by the letters of that block. It is NOT
a standard module substitution, because apex adjacency can vary within a
new block. Calling these blocks modules would lose precisely the needed
coordinate. The base order and apex are part of the decorated graph state.

Thus the whole Berggren tree has an explicit lift into decorated tournaments,
not merely into a numerically labeled ternary tree. The graph has m+1 vertices.
Recovering its spinor requires only counting base vertices and apex indegree.
Unmarked tournament isomorphism classes are not asserted to remain distinct.

## 2. A four-coordinate linear transfer for Hamiltonian paths

Let H(w) be the number of directed Hamiltonian paths in T(w). Let P(w) be
the integer obtained by reading O as binary 1 and I as binary 0; let e(w)
be the last bit, taking e(empty)=0. Then starting at (H,P,e,1)=(1,0,0,1),

 append O: (H,P,e,1)→(H,2P+1,1,1),
 append I: (H,P,e,1)→(H+P+e,2P,0,1).

These are nonnegative integral linear maps on four coordinates. In the
ordered coordinates H,P,e,1 their matrices are

 U=[[1,0,0,0],[0,2,0,1],[0,0,0,1],[0,0,0,1]],
 V=[[1,1,1,0],[0,2,0,0],[0,0,0,0],[0,0,0,1]].

Proof by a closed graph-count formula: a Hamiltonian path consists of an
increasing subset S, the apex, and the increasing complement. Empty S
contributes iff the first letter is O; empty complement contributes iff the
last letter is I. Otherwise put i=max S and j=min(complement). If i<j,
necessarily j=i+1 and the adjacent letters must be IO. If j<i, the letters
at j,i must be O,I, and the i−j−1 intermediate vertices can be distributed
freely between the two segments. Therefore for a nonempty word

 H(w)=1[first O]+1[last I]+#(adjacent IO)
       +sum over j<i, w_j=O,w_i=I of 2^(i−j−1).

Appending O creates an adjacent IO exactly when it removes the old last-I
endpoint contribution, so H is unchanged. Appending I adds P(w) to the
sum, plus 1−1[last I]=e(w) to the endpoint contribution. The empty case
also obeys the recurrence. This proves the transfer without a black-box
Hamiltonian counting identity.

The recurrence gives H(w)=g(P(w)) in S251's convention. In particular this
new spinor lift lands inside the inherited one-apex/transitive-base grammar;
it does not claim a new count-coverage record. What is new here is the
explicit branch-compatible graph family and its noncommutative compiler.

## 3. A compressed noncommutative compiler, with address order explicit

For a word w=w1...wt write M(w)=M(wt)...M(w1), with M(O)=U,M(I)=V.
Keep a PAIR (X,Y) of 4×4 matrices. The three polynomial updates are

 A: (X,Y)→(YX,Y),
 B: (X,Y)→(YX,XYX),
 C: (X,Y)→(X,YXX).

If a chronological Berggren address is d1...dk, initialize (X,Y)=(U,V)
and read the address in REVERSE, dk...d1, using these updates. The final
matrix for the tournament word is YX. Apply it to (1,0,0,1)^t; its first
coordinate is H and its second is the binary cut integer P.

Proof: each update substitutes the relevant O/I image into an existing
representation, and matrix evaluation reverses concatenation. The map
w↦M(w) is an anti-homomorphism. Reverse processing computes
M(phi_dk(...phi_d1(OI)...)), exactly the chronological graph address.
This order is essential: AB and BA have different spinors and counts.

The computation needs O(k) fixed-size integer matrix multiplications, without
materializing the word of m letters. This is an arithmetic-operation bound,
NOT a claim of polynomial bit complexity in k: entries and H can have
exponentially many bits. Modulo a fixed q the same compiler is genuinely
O(k) fixed-cost work and O(1) state space.

## 4. Exact regular congruence languages on the ternary tree

Fix q≥2. Reduce X,Y modulo q. The compiler has at most q^32 states; each
A/B/C update is a well-defined deterministic state transition. Its output
is the first coordinate of YX(1,0,0,1)^t modulo q.

For any residue r, reversed Berggren addresses whose lifted tournament
counts satisfy H≡r (mod q) form a regular language. By closure under
word reversal, the same is true for ordinary chronological addresses.
The latter assertion need not have the same deterministic state bound:
reversal first gives an NFA, and determinization may enlarge it.

This creates a rigorous finite-state laboratory for congruence obstructions
and residue distribution in a genuine graph family. Combining it with
spinor residues modulo s uses the SAME branch address. One must not combine
independent branch choices in the two registers, in agreement with S255's
shared-branch CRT warning. No equidistribution or reachability is assumed.

## 5. The scalar-H obstruction is visible at H=9

Three decorated tournaments have H=9:

 address B:  word OIOIO, (m,n)=(5,2), A/B/C-child H=(107,1615,129);
 address C:  word OOOI,  (m,n)=(4,1), A/B/C-child H=(73,117,33);
 address AA: word OIII,  (m,n)=(4,3), A/B/C-child H=(17,401,545).

Therefore no functions f_A,f_B,f_C of H alone implement these graph branches.
This is an actual collision certificate, not a generic warning about lost
information. It does not rule out a different graph grammar with scalar
branches 3H+2 and 3H+4. It identifies the additional state needed for THIS
Berggren lift.

Two rays have simple exact counts:

 A^k: w=O I^(k+1), H=2^(k+1)+1;
 C^k: w=O^(2k+1) I, H=2^(2k+1)+1.

Both follow immediately from the closed formula. Thus H(C^k)=H(A^(2k))
for every k, although their spinors differ for k>0. Count collisions are
structural and infinite, not just a small-order accident.

## 6. How this connects to the five-core route rather than replacing it

The compiler returns BOTH the ordinary one-apex count g(P) and its full
cut integer P. These are exactly distinct inputs needed by the S251
five-core flank formula: that formula retains g(A), g(B), AB, cut parity,
and a finite core signature. Hence two independently chosen spinor words
can furnish actual flank cuts around any of the five cores, and the
inherited formula evaluates their actual child counts. The full integer
cut registers cannot be replaced by g(A),g(B), as the collision examples
already demonstrate.

This yields an explicit two-spinor/fixed-core graph grammar: select two
Berggren addresses and one finite core cut; compile the words to flank
codes, using the complemented/reversed convention on the right flank
specified in S251; construct that graph and apply the existing formula.
No coverage is claimed for this restricted set of codes. A plausible next
question is whether carefully chosen continued-fraction or substitution
families of flank codes force useful affine count identities after a core
signature cancels the exponential terms. That cancellation problem is now
an exact matrix identity question, not a visual analogy.

## 7. Reproducible new-claim-only experiments

Run `python test_spinor.py` in this directory (standard library only).
The committed results.json records:

- all 1,093 chronological addresses through depth 6: incidence/spinor
  agreement, direct word recurrence versus compressed matrix compiler;
- residues modulo 3,5,8,15 for every such address;
- independent subset Hamiltonian-path DP on all 511 arbitrary cut words
  of lengths 0 through 8, checking the newly proved transfer formula;
- independent graph DP on 32 morphic words having at most 13 base vertices.

These are tests of this new construction, not a verification campaign of
older corpus claims. Proofs above establish the infinite statements.

## 8. Sorting the cut gives a smaller polynomial lift and a permanent obstruction

There is a complementary graph carrier with just two arithmetic registers.
For a spinor set x=m−n (odd), y=n, and SORT its cut to O^x I^y. The
resulting graph is the literal cyclic blowup C3[TT_x,TT_y,K1]: apex→O
block→I block→apex. Its Hamiltonian count, by §2, is

 H_sort=2+(2^x−1)(2^y−1).

Write p=2^x−1,q=2^y−1. The positive spinor matrices give polynomial maps

 A: (p,q)→(p,p+q+pq),
 B: (p,q)→(p+2q+q²+2pq+pq²,p+q+pq),
 C: (p,q)→(p+2q+q²+2pq+pq²,q),

starting at (1,1). Indeed multiplication of 2^x and 2^y becomes
(p+1)(q+1), and 2^(x+2y)=(p+1)(q+1)². Every displayed coefficient is
nonnegative. Thus this is a smaller monotone polynomial graph lift.
Its graph construction resets the cut order after each branch, rather
than retaining the morphic order of §§1–3. It is not the same operation
on tournaments, even though it carries exactly the same spinor.

There is a complete terminating inverse for this sorted IMAGE, independent
of any Chan–Pak conjecture. For a target N≥3, enumerate positive ordered
factor pairs p q=N−2. Retain exactly the pairs for which p+1=2^x,
q+1=2^y with x,y≥1, gcd(x,y)=1, and x odd. Each retained pair outputs
a primitive triple with (m,n)=(x+y,y), and the displayed cyclic tournament
of count N. Conversely every sorted-spinor graph supplies such a pair.
This characterizes the image over all primitive spinors directly. Identifying
it with the image of all Berggren addresses additionally uses the inherited
completeness of that spinor tree; no completeness claim is needed for the
factorization test or explicit graph witnesses.

A permanent congruence obstruction shows why keeping cut order matters.
Since x is odd, p=2^x−1≡1 (mod3). If y is even then q≡0, and if y is
odd then q≡1. Consequently H_sort≡2 or0 (mod3), NEVER 1. In contrast,
the morphic address CA gives H=73≡1 (mod3). For this same address,
(m,n)=(7,4), (x,y)=(3,4), the sorted graph has H=2+7·15=107.
Thus ordering changes the attainable congruence classes even within a
fixed-spinor graph construction. This is a theorem-level distinction,
not evidence for unrestricted count coverage.

The sorted family is a transparent connection to safe cyclic inflation:
its two transitive blocks and singleton are retained as actual modules.
It is the singleton-third-block case of any three-block cyclic/Gram
count framework. No additional Gram identity is needed here; the direct
Hamilton-path partition proof already establishes the formula.

`test_sorted.py` and `sorted_results.json` give separate receipts, leaving
the original experiment receipt unchanged: all 1,093 addresses through
depth6 check the polynomial registers, graph formula and mod3 exclusion;
small sorted graphs are independently counted; a factor inverse is checked
against direct exponent enumeration through target511.

## 9. Integration with the Pareto compiler

The companion [compiler proof, §9](../compiler/PROOF.md) applies its
componentwise-monotone polynomial/Dickson optimization to the pair (X,Y)
of §3. All matrix-word updates and the output H=(YX)[0,0]+(YX)[0,3]
have nonnegative coefficients, so its hypotheses hold without a signed
coordinate workaround. Consequently this morphic graph image also has a
terminating all-height exact inverse, now by finite-residue Pareto search
rather than the sorted family's elementary factorization inverse.

For example, its modulo16 minimum for residue13 is3181 and for residue15
is1615. Therefore targets13 and15 are absent from this morphic family at
EVERY depth. These obstructions are specific to this graph lift, not to
tournament Hamiltonian counts in general. See the companion script and
receipt for the exact finite certificates; the spinor experiments above
do not independently claim these minima from a bounded-depth census.
