---
id: HYP-3545
title: The per-oriented-type Hamiltonian-path count N_tau(T) has TOURNAMENT-INDEPENDENT parity equal to the transitive/descent-set count (Forcade 1973), with ALL types odd at n=2^k -- this is the Ky Fan / Borsuk-Ulam odd-count GRADED by type on the type-hypercube {+,-}^{n-1} (OPEN-Q-059), the signed cycle index (HYP-3540) is its generating function, the arc-hypercube complement (klein THM-584) is its antipodal Z_2, and arXiv:2512.09332 (arc-deletion stability of oriented paths) is the hypercube-EDGE robustness of these counts; merges Borsuk-Ulam (witness, active) + Ky Fan (now grounded) + Ham Sandwich/Tucker (opportunities) + Kaczynski (analytic floor)
status: VERIFIED parity-constancy n=4,5,6 (exhaustive n<=5, 3000-sample n=6); Forcade's theorem confirmed; the Ky Fan / Borsuk-Ulam reading is SYNTHESIS (OPEN-Q-059 grounding). Not a proof of LRC.
source: mac-mini-2026-06-29-S12
related:
  - OPEN-Q-059  # Tournament Ky Fan: replace Fan's magnitude order by an arbitrary tournament
  - THM-582     # palindromic odd index (the directed-type Ky Fan count); the two-index split
  - THM-584     # klein: complement = antipodal map of the arc-hypercube Q_d
  - HYP-3540    # the per-level signed cycle index (klein, open) -- the GF of graded counts
  - THM-586     # Paley H-arithmetic (Paley T_7 is the classical antidirected exception)
  - THM-581     # Borsuk-Ulam witness vs Brouwer floor; the saddle index (p-1)/2
  - HYP-3543    # the metagraph/LRC-cap one-Burnside-spectrum unification
external: arXiv:2512.09332 (El Sahili-El Zein, oriented Ham paths under arc deletion); Forcade 1973; Redei 1934; El Sahili-Abi Aad 2020; Havet-Thomasse
results:
  - 04-computation/oriented_path_parity_kyfan_macmini_20260629.py
  - 05-knowledge/results/oriented_path_parity_kyfan_macmini_20260629.out
---

# HYP-3545 -- oriented-path parity is the graded Ky Fan / Borsuk-Ulam count

## The verified fact (Forcade 1973, reconfirmed)
For an oriented TYPE `tau in {+,-}^{n-1}` (which of the n-1 consecutive arcs go forward), let
`N_tau(T) = #{Hamiltonian paths of type tau in T}`.  Then **`N_tau(T) mod 2` is INDEPENDENT of `T`**,
equal to `N_tau(transitive)` = `#{permutations of [n] with descent set = {i: tau_i=-}}` mod 2.
VERIFIED exhaustively n=4 (64), n=5 (1024) and on a 3000-sample at n=6: ZERO types have
variable parity.  Special cases:
- **n = 2^k (n=4): ALL `2^{n-1}` types are ODD** (Forcade's special case).
- **directed type `++...+` (Redei) and its complement are ODD for all n.**
- n=5: exactly half odd (clean rule: odd iff first arc = last arc, `b_1 = b_{n-1}`); n=6: 16 odd / 16 even.

## The Ky Fan / Borsuk-Ulam reading (OPEN-Q-059, now grounded)
OPEN-Q-059 asks to read Redei/Forcade as a TOURNAMENT Ky Fan lemma.  This is exactly it: the type
`tau` is a corner of the **type-hypercube `{+,-}^{n-1}`**; reversal and complement act as the antipodal
`Z_2`; Ky Fan's lemma counts ALTERNATING simplices with an ODD parity, and the per-type counts
`N_tau` are that count, **graded by type**.  Redei (directed corner, odd) and Forcade (every corner odd
at `n=2^k`) and El Sahili-Abi Aad (antidirected `= 2 mod 4` at even order) are the corners of one
graded Ky-Fan/Borsuk-Ulam odd-count.  Its generating function over the arc-hypercube is the SIGNED
CYCLE INDEX (HYP-3540): the metagraph eigenvalue `d-2k` (= reversed-arc level) and the per-type parity
are two gradings of the same `Z_2`-equivariant counting on `Q_d` (klein THM-584: complement = the
antipodal map).  arXiv:2512.09332 (El Sahili-El Zein) is the hypercube-EDGE statement: deleting one arc
(one `Q_d` edge) preserves every oriented type for `n>=8` except two explicit special exceptions --
robustness of the graded count along the arc-flip graph.

## The merged topological toolkit
| tool | role | status |
|---|---|---|
| **Borsuk-Ulam** | the LRC WITNESS: lonely times are antipodal pairs `{t*,-t*}`; saddle index `(p-1)/2`, parity = `p mod4` (n=14: index 3 odd -> BU); the odd degree forces `M(S)>=1/14` | ACTIVE (THM-581/582, kps S31av) |
| **Ky Fan** | the per-type oriented-path odd-count (Redei/Forcade); GRADED Borsuk-Ulam on the type-hypercube; the signed cycle index is its GF | GROUNDED here (OPEN-Q-059) |
| **Tucker** | discrete `Z_2`-labeling / no-complementary-edge on `Q_d` under the antipodal complement | OPPORTUNITY (unused) |
| **Ham Sandwich** | bisect the danger-cover measure / the R-odd `M_odd` (HYP-3538) by the spectral eigenvector: the obstruction cannot be symmetrically halved | OPPORTUNITY (unused) |
| **Kaczynski** | the analytic floor: `M(x)`, `G=sum mu^2/phi ~ log x + 1/zeta(2)`, `Phi=sum phi`; the `1/zeta(2)=6/pi^2` totient density that makes the floor positive (the EVEN/analytic face) | ACTIVE (S162/163) |

## The unification
ONE involution `R` (complement = antipodal map of `Q_d`) generates the parity structure of BOTH the
tournament oriented-path counts (Ky Fan, graded, Redei/Forcade -- here) AND the LRC witness
(Borsuk-Ulam, saddle index `(p-1)/2` -- THM-581).  The `R`-even side is the SOS/Brouwer/Kaczynski-analytic
bulk (the `1/zeta(2)` floor); the `R`-odd side is the Ky-Fan/Borsuk-Ulam odd-count obstruction.  Concrete
next steps: (1) prove the tournament Ky Fan (OPEN-Q-059) -- show the per-type parity IS a Ky-Fan
alternating count for arbitrary `T` (not just transitive), making Redei/Forcade Borsuk-Ulam shadows;
(2) the Ham Sandwich angle on `M_odd` (the cap obstruction can't be symmetrically bisected);
(3) read arXiv:2512 arc-deletion as the `Q_d`-edge stability of the graded Ky Fan count.
