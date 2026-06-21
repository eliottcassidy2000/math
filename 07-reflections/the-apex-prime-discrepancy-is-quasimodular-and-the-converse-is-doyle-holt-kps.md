# The apex-prime discrepancy is quasimodular (E_2 on three legs), and the converse is Doyle-Holt

**kind-pasteur-2026-06-21.** The owner handed four inspirations -- Delsarte LP, Tanner/weakly-regular
graphs, the Doyle-Holt half-arc-transitive graph, and the algebraic modular-forms picture (branched
covers over {0,1,infty}, Gamma(2)) -- and asked to continue toward the LRC proof. Each landed
concretely. The surprise is how cleanly they place the LRC objects in known mathematical homes, and
how honestly those homes explain the wall the proof keeps hitting.

## Modular: the discrepancy is quasimodular E_2, not a holomorphic form

The L7 sharp cell-discrepancy (HYP-2739) is residue-only mod the apex prime 7. Pushing it: it is
residue-only mod P for **every** prime (HYP-2743, verified P=3..43), with one uniform closed form
(HYP-2745):
```
   D_P(p,q) = G_P(||p||,||q||) / (P p q),   G_P = [2 A B (P-A)(P-B) + 2 C (P-C)] / P,
   A=||p||_P, B=||q||_P, C=||pq||_P,   ||x||=min(x mod P, P - x mod P).
```
Each "leg" `g(t)=2t(P-t)` equals `-2P^2 B_2(t/P) + P^2/3 = P * R_eff(0,t)` on the cycle graph `C_P`.
So **the discrepancy is built from cycle-graph effective resistances = second-Bernoulli (quasimodular
E_2 / weight-2 Eisenstein) values**, on three Markoff-like legs `||p||, ||q||, ||pq||`.

Is it "modular"? The honest answer is the sharp one. There is **no** nontrivial holomorphic weight-2
form; `B_2`/`E_2` is exactly the **quasimodular** generator, and the discrepancy is its **absolute /
L1** avatar. The classical **signed Dedekind sum** `s(slope,P)` -- which *does* carry the SL_2(Z)
reciprocity and the eta multiplier -- is only the **degree-1 Fourier shadow** of `G_P`. The symmetry
group is the Klein four-group `<z->-z, z->1/z>` (PROVED for all `P>=5`): the order-2 (hyperelliptic)
part survives, the order-3 (cubic/QR/doubling) part washes out, uniformly in `P`. And `G_P` is a
function not of the slope but of the **pair** `(p,q) mod P` modulo `<+-, swap>` -- the third leg
`||pq||` is a multiplicative coordinate the slope cannot see (this corrects my own HYP-2742).

This is the precise modular fact, and it is the same wall in new clothes: the LRC bound needs the
**absolute** discrepancy `|E_2|`-type object, not the holomorphic/signed one with clean transformation.
That the signed reciprocity object is the *shadow* and the absolute one is the *carrier* is exactly
MISTAKE-082's conditional-convergence lesson, now read on the modular curve.

## Doyle-Holt: half-arc-transitivity IS the converse Z_2

The owner's analogy was exact (HYP-2748). By Tutte, a half-arc-transitive graph carries an
Aut-invariant orientation `D` -- a partial tournament -- with arc-orbits `{D},{D^op}` and **no
automorphism reversing it**. That is literally the tournament converse `T <-> T^op` (the Z_2 whose
fundamental domain is the half-tiling, THM-549/550) one categorical level up: half-arc-transitive =
"the converse Z_2 is not realized by the symmetry group." The genuine tournament incarnation is a
**vertex-transitive non-self-converse tournament**, smallest at `n=21` on the Frobenius group `F_21`
-- exactly the canon THM-052/MISTAKE-013 fact. Circulant/Paley tournaments are *always*
self-converse (inversion `i->-i` realizes the converse) -- the digraph shadow of Chen-Quimpo's "no
half-arc Cayley graph on an abelian group" -- which is *why* the SC spine of the metagraph is the
circulant/Paley locus. The owner's "edge maps to edge, but only one of two ways" is the orientation
Z_2; the half-tiling was its fundamental domain all along.

## Tanner/Delsarte: honest negatives that confirm the structure

The Tanner graph of the relation code `Lambda(E)` gives **no** usable expansion bound (girth is
degenerate `=4` for every set; the spectral gap runs the wrong way -- the hard AP code is *denser* and
mixes *better*; the absolute kernel enumerator diverges). What survives is the **weight distribution**:
AP = anti-MDS (`d_min=2`, low-weight-codeword-rich), Sidon/arc = MDS (`d_min=3`, sparse) -- the same
extremal dichotomy as HYP-2723/2724. And the cover bound is a Delsarte LP (HYP-2726) whose
consec-saturation is irreducibly aggregate (HYP-2738): no nonnegative certificate certifies it,
because to certify consec-max you must *subtract* a consec-max quantity, forcing the alternating
Bonferroni signs. The Delsarte/Tanner/weakly-regular lenses sharpen *which* objects matter; they do
not supply the missing inequality.

## The meta-point

The apex prime 7 runs through all of it: the discrepancy (`E_2` on three legs over `F_7`), the
converse Z_2 (Doyle-Holt / `F_21 = 7:3`), the Delsarte dual roots (residue classes mod 7). The **three
legs** `||p||, ||q||, ||pq||` echo the project's oldest theme -- the staircase triangle's three sides.
The structure is gorgeous and now well-placed: quasimodular, hyperelliptic, effective-resistance,
Markoff. But the proof's bottleneck (HYP-2602: consec/AP minimizes `mu_{1/7}` for unbounded spread)
is, by HYP-2738, *aggregate* -- a signed cancellation with no clean certificate. The modular home tells
us why the absolute object is the carrier, and the effective-resistance bound `R_eff <= (P-1)/4` is a
concrete lead for sharpening the L7 window. But the wall is the wall: the signed-vs-absolute split,
seen now as `|E_2|` vs the holomorphic form. The mathematics keeps handing us the same coin from new
angles -- which is itself the strongest evidence that there is exactly one coin.
-> HYP-2739/2742/2743/2745/2748, HYP-2602/2738, HYP-2726, THM-052/549/550, MISTAKE-082, OPEN-Q-108,
the triangle foundation.
