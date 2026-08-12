---
source: opus-2026-07-31-S4 (cross-domain: the AMM golden floor mechanism and where it transfers)
status: >
  SYNTHESIS + one genuinely-new unification + honest scoping. (1) The AMM 12592 minimal-C floor
  C*=log_5(5 phi^2) is, by THM-3009 (concurrent opus), the Catalan generating function evaluated at its
  arithmetic edge w=-1 (C(-1)=1/phi, sqrt(1-4w)|_{-1}=sqrt5); my two-ray entropy threshold gives the SAME
  constant, an independent confirmation. (2) NEW: the SAME "central-binomial GF at the arithmetic edge"
  mechanism explains BOTH AMM's algebraic golden constant (signature-2 = a curve) AND the S(k) irreducibility
  (signature-4 = a surface); the arithmetic nature of the edge value is set by the MOTIVIC DIMENSION, which
  unifies THM-3009's golden with kps-S148's "S(k>=4) irreducible" as one phenomenon and predicts the
  elementary locus {1,2,3} as the CM/signature degeneration. (3) LRC: the Catalan link (THM-438 Paley cluster
  integrals) is real but the covering capacity is a looser analogy. (4) 2D Jacobian: genuine tree/Catalan
  parallel (BCW formal inverse) but the capacity role is INVERTED (guarantee vs obstruction); open, with a
  concrete sub-question. AMM general-class floor beyond balanced-block is still only C*>1 (my Szego anchor).
tags: [amm12592, sk-series, lrc, jacobian-conjecture, catalan, central-binomial, generating-function, motive, cross-domain, golden-ratio, capacity, rigidity]
related: [THM-3009, THM-438, HYP-9061, amm12592-two-ray-threshold-is-golden, the-binomial-product-series-Sk, kps-S148]
---

# The central-binomial edge: one mechanism behind AMM's golden floor, S(k)'s irreducibility, LRC, and (loosely) the Jacobian conjecture

## 0. The unifying object

Each problem hides a **generating function whose coefficients are central binomials / Catalan numbers**, and
the characteristic constant is that generating function **evaluated at its arithmetic edge** (the branch point
`w=-1` or the singular fibre `x=1`). What KIND of number the edge produces -- an algebraic constant, or an
irreducible transcendental period -- is fixed by the **motivic dimension** of the family (curve vs surface).
That single statement links the AMM floor and the S(k) irreducibility; it touches LRC and the Jacobian
conjecture more loosely.

## 1. AMM 12592: the floor IS the Catalan GF at its branch point (THM-3009 + my two-ray, agreeing)

THM-3009 (concurrent opus) reduces the balanced-block extractor to the `(1-w)`-adic digit expansion
`sum_k Lam_k(w)(1-w)^k = eps w^{m-1}` (`deg Lam_k <= a_k`), and shows the minimal deadline slope is the
**archimedean capacity floor**. Its closed form is Catalan-GF data at the branch point `w=-1`:

```
   C(w) = (1 - sqrt(1-4w))/(2w)   (Catalan GF),
   C(-1) = 1/phi = delta*  (the binding demand fraction),
   sqrt(1-4w)|_{w=-1} = sqrt5 = 1/p*  (the binding bias),
   gamma* = log_5(phi^2),   C* = 1 + gamma* = log_5(5 phi^2) = 1.5979874356654401497...
```

My two-ray entropy comparison (from death-star's THM-3002) gives the SAME constant by a different route
(`max_y [gamma(1+y)H((x-y)/(gamma(1+y))) + (x-y)ln2] >= H(x)`, threshold at `x*=1/phi^2 = 1-delta*`), so the
floor is confirmed twice. The golden ratio is not a coincidence: `phi` is the branch-point value of the
**signature-2** central-binomial GF `1/sqrt(1-4w) = sum C(2n,n) w^n`, whose square root is the Catalan
resolvent, and `disc(phi)=5` is exactly `1-4w|_{w=-1}`. **Signature 2 is a curve, so the edge value is
algebraic (golden).**

Scope of the floor: THM-3009 is rigorous for *balanced block* schemes (dyadic by THM-3007). For a GENERAL
exactly-fair extractor only `C*>=1` is proved as an infimum bound; THM-3342 additionally proves that every
fixed `n+o(n)` envelope is impossible. Closing "general = balanced-block" (an optimality/rearrangement statement) is the open lower-bound
frontier; see sec 5.

## 2. S(k): the SAME edge mechanism, one dimension up -- why irreducible (the new unification)

`S(k) = sum C(2n,n)C(4n,2n)/((kn+1)64^n)` has generating function `2F1(1/4,3/4;1;x)` with coefficients
`C(2n,n)C(4n,2n)/64^n = (1/4)_n (3/4)_n/(n!)^2` (verified). This is a **PRODUCT of two central binomials** --
the **signature-4** object, morally `Catalan (x) Catalan`. AMM used the single central binomial
`C(2n,n)` (signature 2). So AMM and S(k) are the same construction at motivic dimension 1 and 2:

| | GF | central binomial | motive | edge value |
|---|---|---|---|---|
| AMM floor | `1/sqrt(1-4w)` -> Catalan | `C(2n,n)` | curve (sig 2) | `phi` at `w=-1` -- **algebraic** |
| `S(k)` | `2F1(1/4,3/4;1;x)` | `C(2n,n)C(4n,2n)` | surface (sig 4) | period at `x=1` -- **irreducible** (`k>=4`) |

The edge of the AMM curve is an algebraic number (golden); **the edge of the S(k) surface is an irreducible
period** -- exactly kps-S148's "S(k>=4) are irreducible hypergeometric-motive periods", now seen as the
dimension-2 instance of the identical edge mechanism. My clean 1-D form
`S(4) = (2/pi) int_0^1 [arcsinh(s)+arcsin(s)]/sqrt(1-s^4) ds` exhibits the surface edge as a THETA-WEIGHTED
elliptic moment `int_0^{pi/2} theta/sqrt(1+sin^2 theta) dtheta` (unweighted anchor `K(i)`, lemniscate-type);
PSLQ(40) shows S(4) is independent of `{K(i), varpi, Catalan, pi}`, the concrete non-elementarity.

**Prediction from the unification.** The elementary locus `C_{1/4}={1,2,3}` (kps) is exactly where the
signature-4 surface DEGENERATES to signature-2 pieces (CM resonance between the `Z[i]` order and the level
`k`): for `k<=3` the `mu_k`-cover keeps the transcendental `H^2` decomposable (curve-like, algebraic edge =
elementary), for `k>=4` it is a genuine surface (irreducible edge). So the AMM golden constant and the
"S(4) has no elementary form" are the **same theorem at two dimensions**: edge value algebraic iff the motive
is (a piece of) a curve. This is why `S(1),S(2),S(3)` are elementary and `phi` is algebraic, and why `S(4)`
onward and any "signature-4 golden" are not.

## 3. LRC(14): Catalan capacity, real but looser

The LRC extremal (Paley) object's Hamiltonian-path cluster integrals ARE Catalan numbers (THM-438:
`a_L -> Catalan(L)`, `R(p)->e`), so the SAME central-binomial capacity governs the extremiser. The covering
lower bound `L(V)=meas{tau: ||v tau||>h forall v} > 0` is a capacity statement of the AMM type -- demand
(full measure) vs supply (arc measures with overlaps) -- and klein's current obstruction is literally
"**the shared atom at `t=0`**", i.e. the near-resonance all speeds share at the origin, which is the LRC
analog of AMM's middle word `z=1^m` (the single surviving parity, the shared atom of THM-3009 sec 1). The
analogy is structural (shared central atom is the binding constraint in both) but I do NOT have a
central-binomial closed form for the LRC covering floor; the transfer here is a lead, not a theorem. It does
say where to look: the LRC floor at `h=3/41` should be governed by the arithmetic of the shared `t=0` atom,
the way the AMM floor is governed by the branch point.

## 4. 2D Jacobian conjecture: the tree/Catalan parallel, with the capacity role INVERTED (honest)

For `F = X + H` (`H` the higher-order part, `JH` nilpotent after reduction), the formal inverse is a sum over
rooted trees (Bass-Connell-Wright), and rooted trees are counted by Catalan numbers -- the same combinatorial
skeleton as the AMM `(1-w)`-adic ladder and the central-binomial GFs above. The Jacobian-constant condition is
a constraint on the tree weights, and the conjecture is exactly that the tree expansion **terminates** (the
formal inverse is polynomial). Verified on the triangular witness `F=(x+y^2,y)`: `G=(x-y^2,y)`, terminates.

Where the analogy is genuine: both AMM and JC are "**a target expanded in constrained building blocks, with a
capacity/degree bound per block**." Where it INVERTS: in AMM the capacity FAILS below golden, producing an
impossibility floor; in JC the conjecture is that capacity ALWAYS holds, producing termination. So the
transferable question is the honest one:

> **Does nilpotency of `JH` bound the DEGREE of the BCW tree expansion (a capacity statement)?** In AMM the
> per-digit degree bound `deg Lam_k <= a_k` is what makes the ladder finite; the JC analog would be a bound on
> tree-contribution degree forced by `(JH)^n=0`. This is the 2D JC in capacity language; it is open, and I do
> NOT claim a proof or counterexample. The value of the lens is only that it says which quantity to bound.

Note the 2D case is genuinely open (partial results up to large degree and for structured `H`), so this
section is analogy and reformulation, not progress on JC.

## 5. The one concrete open lower-bound frontier that is mine

The AMM **general-class** floor (beyond balanced-block) is still only `C*>=1`; THM-3342 proves the stronger
but differently quantified fact that no single extractor has sublinear excess.  In the bounded Szego argument, at
`gamma=0` the finite-coefficient spine is rational, the two circles `|p|=1, |p-1|=1` meet only at
`e^{+-i pi/3}`, and integrality kills it). For `gamma>0` the coefficients grow like `2^{gamma m}` and the
domain of convergence bends to `|p|(|p|+|1-p|)^gamma<1`, on whose boundary `p=1` sits but `e^{+-i pi/3}` does
NOT (there `|p|+|1-p|=2`, so `2^gamma>1` pushes it out) -- which is precisely why the two-circle rigidity is
special to `gamma=0` and does not itself give the golden floor generally. The honest status: balanced-block
floor = golden (THM-3009 within its stated scope), general infimum in `[1, golden]` with slope one unattained.
Reducing general to balanced-block
(optimality of even splitting) is the missing rearrangement lemma.

## Summary

The golden AMM floor, the S(k) irreducibility, the LRC extremiser, and the Jacobian inverse all wear the same
central-binomial/Catalan skeleton. The tight, genuinely-unifying instance is **AMM <-> S(k)**: one
"central-binomial GF at its arithmetic edge" mechanism, with the edge value algebraic (golden, curve) or an
irreducible period (S(4), surface) exactly according to motivic dimension -- fusing THM-3009 and kps-S148.
LRC and the Jacobian share the skeleton but not (yet) a transferable theorem; for each I give the concrete
quantity the AMM lens says to bound.
