# A source-marked elliptic lift, and what the torus quotient loses

2026-10-16. **PROVED:** an exact elliptic-coordinate lift of every supplied
native Collatz word, faithful decoding of its marked affine frame, strict ROOT
substitution, and two obstructions to forgetting that frame. **CITED:** the
classical Nagell–Lutz integrality theorem and Chang's cyclic analytic quotient.
**FINITE-EXACT:** the controls below. **OPEN:** a new uniform signed-energy
budget, a new source-covering rule, and universal Collatz convergence. The lift
does not claim any of these.

## 1. Inheritance and connection contract

The current closest mechanism is
[refinement debt](collatz_refinement_debt_20261016a.md), which transports a
Gram remainder through a **common signed native map**. Its essential data are
the original source coefficients, full integer clocks, and shared integer
terminal labels. The one-reset `7X` failure at `k=23` and the missing angular
term are retained in [paired energy](collatz_paired_energy_20261015g.md).
Neither finite energy observations nor a changed quotient close its open
uniform budgets. The less-used relevant sidecar here is the full affine
translation, decoded in [translation phase](translation_phase_decoder_20261005.md).

The [golden order bridge](golden_order_bridge_20261014j.md) already separates
a cyclic transfer matrix from a uniform poset-order law. The new
[recursive interfaces](recursive_interfaces_20261016b.md) constructs infinitely
many genuine elliptic mixed contacts but retains geometric compatibility as
a separate obligation. This packet applies the same discipline to Collatz.

Board: fixed integer source; native valuation word; affine carry; marked
elliptic coordinates; shared signed terminal; cyclic versus dense quotient.
The source is a supplied legal word. The target is an isomorphic elliptic
model with a marked rational point. The map preserves the exact endpoint and
clock. Forgetting the frame destroys the native rule; forgetting the integer
terminal destroys cancellation. The cheapest controls are the actual seam
`3 -> 5 <- 5`, and the words `(1,2)` and `(2,1)` with equal slope but distinct
carries and source cylinders. No tournament orientation is imposed on ties.

## 2. A faithful guarded lift

Use the fully odd map `U(n)=(3n+1)/2^v2(3n+1)` on positive odd integers,
stopped at ROOT 1. For a positive valuation word `w=(a_1,...,a_m)`, write

\[
 F_w(x)=\frac{Px+B}{Q},\qquad P=3^m,\quad Q=2^A,\quad A=\sum_i a_i.
 \tag{1}
\]

Appending `a` changes `(P,Q,B)` to `(3P,2^a Q,3B+Q)`. The empty word has
`(1,1,0)`. Exact native valuations hold precisely on the odd cylinder

\[
 n\equiv (Q-B)P^{-1}\pmod {2Q}.                         \tag{2}
\]

Indeed a final odd integer can be reversed through each affine step; oddness
at each reverse step follows modulo 2. At the last letter, divisibility of
`Qz-B` by `P=3^m` implies `3 | (2^a z-1)`, because the preceding power of 2
is invertible modulo 3. Repeating proves integral reverse intermediates.
Equivalently, the usual inductive cylinder calculation gives (2).
Positive source values make the forward formal intermediates
positive. The additional **stopped** guard is that no proper endpoint is 1;
this is checked separately, not inferred from the cylinder.

A direct prefix check makes the exact valuation explicit: if a nonempty
tail follows a prefix of cost `A_i`, then
`B=3^(tail length) B_i+2^A_i B_tail`, with `B_tail` odd. Divisibility of
`Pn+B` by `2^(A_i+1)` therefore forces
`P_i n+B_i = 2^A_i mod 2^(A_i+1)`. The final prefix is odd by (2).

Fix an odd integer source `s>=3`. Set

\[
 A_s=s(1-s),\quad E_s:y^2=x^3+A_sx,\quad P_s=(s,s).
 \tag{3}
\]

The discriminant is `-64 A_s^3 != 0`, so this is an elliptic curve. Given
`u=P/Q` and `v=B/Q`, define

\[
 E_{s;u,v}:Y^2=u(X-v)^3+A_su^3(X-v),\qquad
 \iota_{u,v}(x,y)=(ux+v,u^2y).                         \tag{4}
\]

Substitution proves (4); its inverse is `((X-v)/u,Y/u^2)`. It extends to an
isomorphism of the smooth projective cubics preserving their point at
infinity. The affine plane transformation also preserves lines and vertical
reflection, so the chord-and-tangent group law is transported. In particular
the marked point is

\[
 P_{s;w}=(F_w(s),u^2s).                               \tag{5}
\]

If the next affine word is `cx+d`, the new frame is
`(cu,cv+d)` and the point changes by `(X,Y)->(cX+d,c^2Y)`. Thus composition
is chronological and source-preserving. On an authenticated native word its
first coordinate is the actual odd integer endpoint, and `A` is the number
of shortcut steps for `T(even)=n/2`, `T(odd)=(3n+1)/2`. It is not the number
of steps for the unaccelerated odd rule `3n+1`, whose clock is `A+m`.

**Faithfulness.** The reduced slope determines both `m` and `A` by unique
factorization. The full carry determines the word: put `N=B mod 2^A`.
For `m>1`,

\[
 Z=(N-3^{m-1})\bmod 2^A,
 \quad a_1=v_2(Z),
 \quad (N,A,m)\leftarrow(Z/2^{a_1},A-a_1,m-1).          \tag{6}
\]

The last numerator must be 1. This is the inherited fixed-length real-phase
decoder: the first summand is `3^(m-1)` and the remaining numerator is
`2^a_1` times an odd tail. Re-encoding the full `B`, not only its residue,
rejects extraneous integer translations. This proves that the **marked**
elliptic lift is a lossless representation of the native receipt. It does
not provide a cheaper word or a rule-selection oracle.

The change of elliptic coordinates (4) works for any nonzero rational
affine slope. Its specifically Collatz content is supplied by the restricted
frame, the native cylinder, the literal clock, and the stopped-source check;
the curve equation does not impose those guards by itself.

Payment still has its original integer meaning:

\[
 F_w(s)<s\quad\Longleftrightarrow\quad (P-Q)s+B<0.       \tag{7}
\]

`paid` checks (7) only after strict source authentication. A smaller child
is a dependency, not a ROOT proof. `discharge` requires a supplied first-hit
ROOT word for that exact child; the words are then concatenated and checked.
An already terminal child permits only the empty suffix.

## 3. Elliptic invariants do not become a Collatz rank

Every model (4) for one fixed source is rationally isomorphic to (3).
Consequently its `j=1728`, rational group rank, and torsion subgroup are
unchanged by the lifted Collatz word. These invariants cannot strictly rank
this representation. This does not assert that all curves for different
sources are rationally isomorphic; they merely share `j=1728`.

There is a stronger explicit torsion hostile. On the integral model (3),
the tangent slope at `P_s` is `s+1/2`, hence

\[
 x(2P_s)=s^2-s+\tfrac14.                              \tag{8}
\]

This is nonintegral for every integer `s>=2`. The classical Nagell–Lutz
integrality theorem for rational torsion on an integral short Weierstrass
equation therefore implies that `P_s` is nontorsion: otherwise its finite
double would be torsion with integral coordinates. See
[Sutherland, 18.783, Problem Set 3, Problem 2](https://math.mit.edu/classes/18.783/2025/ProblemSet3.pdf).
All points (5) are nontorsion under the isomorphisms.

For example the supplied strict word `(1,4)` takes `3 -> 5 -> 1` and gives

\[
 E_3:y^2=x^3-6x,\qquad
 (P,Q,B)=(9,32,5),\qquad P_{3;(1,4)}=(1,243/1024).
\]

Its first coordinate is ROOT while the point remains nontorsion. More
generally, for every `k>=2`, the already elementary ROOT ray
`s_k=(4^k-1)/3` has actual word `(2k)` and endpoint point
`(1,9s_k/16^k)`, again nontorsion. This is an infinite exact control, not
new Collatz coverage or a new classical rank-family claim.

Even keeping one fixed model cannot repair this by the most direct affine
ansatz. Suppose `(x,y)->(ux+v,cy)` maps `E_s` to itself with positive rational
`u` and nonzero rational `c`. Coefficient comparison gives `v=0`,
`c^2=u^3`, and `c^2 A_s=A_su`. Since `A_s!=0`, `u=1` and `c=+-1`.
No nontrivial positive Collatz slope occurs in this ansatz. This is not a
classification of arbitrary nonlinear maps or other elliptic constructions.

## 4. Where the current native energy lives

The compiler exports the projection

\[
 (s,w,E_{s;u,v},P_{s;w})\longmapsto
             (-1)^A e_{F_w(s)}.                        \tag{9}
\]

Extension changes the clock by the supplied native cost and uses the actual
shared integer endpoint, so (9) commutes with every common native extension
used in the current energy account. It retains the original source mark;
there is no population reweighting by curve points or isomorphism classes.

For `3 -> 5` with word `(1)` and the empty receipt at 5, the target columns
are `-e_5` and `+e_5`. Their curve points are `(5,27/4)` and `(5,5)`, on
different source-marked models. The cancellation occurs after the exact
projection (9); equality of models or marked points would miss it. Conversely
merging by `j` alone would identify every source, including unrelated empty
receipts at 3 and 7. Thus elliptic data can decorate the common observed
graph but cannot replace its shared integer terminal/sign interface.

There is an equally small slope hostile:

\[
 F_{(1,2)}(x)=(9x+5)/8,\quad n=11\pmod {16};\qquad
 F_{(2,1)}(x)=(9x+7)/8,\quad n=9\pmod {16}.              \tag{10}
\]

The slope and total clock agree; carry and actual source guards do not.
The word order is therefore an indispensable part of any proposed
arithmetic or elliptic selection mechanism.

## 5. A genuine golden elliptic quotient, and a non-Hausdorff replacement

**CITED.** Chang constructs an analytic elliptic curve by quotienting the
punctured sphere by the cyclic Möbius map `gamma(z)=1/(z-1)`, whose fixed
points are `phi` and `phi'=1-phi`; see
[Chang, The Fibonacci sequence and an elliptic curve, §§2,4](https://jtnb.centre-mersenne.org/item/10.5802/jtnb.1210.pdf).
Only that quotient construction is used here.

Direct calculation with `phi^2=phi+1` gives

\[
 h(z)=\frac{z-\phi'}{z-\phi},\quad q=\phi-2=-\phi^{-2},
 \qquad h(\gamma z)=q h(z).                             \tag{11}
\]

Thus the quotient is `C*/q^Z`. In logarithmic coordinates its periods are
`2 pi i Z + (log|q|+pi i) Z`, a discrete rank-two lattice because the second
generator has nonzero real part. This supplies an actual elliptic bridge
for the golden Möbius recurrence, not a numerical match of constants.

The tempting replacement by the group generated by both Collatz slopes
does not have this property. The shortcut slopes `1/2,3/2` generate
`<2,3>={2^a3^b:a,b in Z}`. The ratio `log 2/log 3` is irrational, since
a rational ratio would contradict unique prime factorization. The additive
group `Z log 2+Z log 3` has arbitrarily small nonzero positive elements by
the pigeonhole principle; integer multiples then approximate every real
number. It is dense in `R`. Therefore `<2,3>` is dense in `R_{>0}` but
countable, so the orbit of 1 is not closed in `C*`. The orbit-space quotient
`C*/<2,3>` is not T1, hence not Hausdorff and not an analytic elliptic curve.

More precisely, any continuous map from `C*` to a Hausdorff space invariant
under multiplication by both 2 and 3 is constant on each positive radial
ray: approximate a positive scale by group elements and use continuity.
Its restriction to the positive integer sources is constant. This proves
a specific loss, not a prohibition on all elliptic approaches to Collatz.
Moreover the quotient ignored the additive carry in (10), and parity guards
mean the actual dynamics was never a global scaling-group action.

## 6. API, exact controls, and remaining obligation

The script is
[collatz_curve_bridge_20261016c.py](../../04-computation/experiments/collatz_curve_bridge_20261016c.py);
the deterministic report is
[collatz_curve_bridge_20261016c.json](collatz_curve_bridge_20261016c.json).
`decode_frame` is a total membership reader. `compile_lift`, `audit_lift`,
`compose`, `paid`, and `discharge` retain the unchanged source. No production
function discovers a ROOT word. ROOT 1 itself is a terminal symbol, not the
singular specialization `s=1` of (3).

The finite universe contains every valuation word of length 0 through 5
with letters 1 through 4 (1,365 frames), and three first cylinder lifts per
word when the source is at least 3. Strict ROOT checking retains 4,073 native
lifts and verifies all 23,108 word splits. Independent literal valuations,
inverse coordinates, and direct tangent doubling check the lift. The report
also checks `2<=k<=64` on the infinite ROOT-ray theorem, 420 exact conjugacy
identities in `Q(phi)`, and 22 malformed/native/padded-ROOT hostiles.

Reproduction, from the repository root:

```text
python -B 04-computation/experiments/collatz_curve_bridge_20261016c.py
python -B -O 04-computation/experiments/collatz_curve_bridge_20261016c.py
```

Both runs reproduce **72,680 checks**, with byte-identical LF-normalized
JSON, SHA256
`bd63c6bb96beefbd1ef1584198aa22bf9a6dd42f1c6105833818d0d102476f5d`.

Independent proof/API audit and optimized-versus-saved replay: **PASS**.
The audit checked the transfer and its scope; it did not re-audit the cited
classical literature.

The useful next target is now precise: an elliptic observable must provide
a new **source-bound** inequality or a common-terminal signed cancellation
estimate beyond isomorphism-invariant data, and its pullback must preserve
(1), (2), and (9). The present compiler makes that transfer independently
testable. It supplies no such new inequality, and no ungrounded child is
declared complete.
