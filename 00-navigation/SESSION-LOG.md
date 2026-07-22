## kind-pasteur-2026-07-21-S128c148 -- THE ADDITIVE ROUTE BYPASSES THM-1550: its algebraic core formalized kernel-pure (corrects HYP-8975)

**Owner:** work the additive route (THM-2101) that bypasses THM-1550 entirely; get creative; explore
past work.

**CORRECTION to my HYP-8975.** Last session I concluded "all DvdK routes bottom out at THM-1550" --
WRONG (death-star-S113 cited it). THM-2101 sec6 is explicit: the product/log/Wiener-Hopf bottleneck
(THM-1550) is created by INTEGRATING the observable (Abel A(G-1)=int(G-1)ds/s). The additive route
stays BEFORE integration (eq 15) and closes by ADDITIVE ORBIT INCIDENCE -- a purely algebraic
group-action identity, no product/Hensel/log/Wiener-Hopf.

**FORMALIZED (GMC2AdditiveOrbitSum.lean, kernel-pure [propext,Classical.choice,Quot.sound], pushed):**
`additive_orbit_contradiction` -- a finite group G acting TRANSITIVELY on a finite set Omega over a
char-0 field K, weight w:Omega->K; if hB (full-root Lagrange sum ∑_Omega w = 0) and hA (∀σ, ∑_{α∈S}
w(σ•α) = 1), then False (∑_σ ∑_{α∈S} w(σα) = |G| by hA; = ∑_α |Stab|·(∑_Omega w) = 0 by hB+transitivity;
so |G| = 0 in char 0, absurd). This is the ADDITIVE ANALOG of thm2067_contradiction (THM-2067's
multiplicative closer) and needs NO THM-1550. It was NOT in the tree (the residue files are a different
mod-p route); freshly formalized. Helpers sum_smul_eq, orbitSum_const (transitivity), orbitSum_zero.

**REFRAMED FRONTIER -- TWO independent routes to DvdK1, NOT sharing a crux:**
| route | closer | remaining input | THM-1550? |
| multiplicative (THM-2067) | thm2067_contradiction (done) | small-root PRODUCT Π=c·t (THM-1550, DEEP) | YES |
| additive (THM-2101) | additive_orbit_contradiction (done THIS session) | hA + hB (ALGEBRAIC) | NO |
hB = elementary Lagrange (∑ α^{M-1}/Φ_u(α)=0, since M-1<=deg-2, N>=1); hA via THM-2101 proof-3 t-adic
formal partial fractions ("needs neither continuation nor specialization") => algebraic.

**HONEST.** Did NOT finish DvdK1 on the additive route: hA and hB still need concrete formalization on
Φ's rootSet, and the wrapper instantiated at Φ.Gal ↷ Φ.rootSet (transitivity already
GMC2GalRootAction.isPretransitive_rootAction; declined the concrete instantiation to avoid the
CharZero/instance-diamond rabbit hole death-star hit in S113). Closed the additive route's one
STRUCTURAL gap kernel-pure + corrected the fleet's belief that THM-1550 is unavoidable. Files
GMC2AdditiveOrbitSum.lean + reflection the-additive-route-bypasses-thm1550-formalized-kps-S128c148.
HYP-8980. NAMED NEXT: formalize hB (Lagrange), hA (t-adic partial fractions), concrete instantiation.

## kind-pasteur-2026-07-21-S128c147 -- DvdK1 elementary reductions (zero-charge + both-signs + gcd) + the honest bypass verdict
## boxeph-2026-07-22-S238 -- concrete root-packet lemma: the additive-route core, kernel-pure (HYP-8980)

**Owner:** work the additive route bypassing THM-1550; explore additive work; aim at completing.

**DELIVERED kernel-pure (GMC2RootPacketConcrete):** codex had the abstract additive machinery (GMC2LaurentShiftCheckA); I built the concrete instantiation at the Galois action on roots.
- weight_equivariant: w(a)=a^k/Phi'(a) Galois-equivariant (aeval_algHom_apply + GMC2GalRootAction.coe_smul).
- root_packet_eq_zero: THE root-packet lemma (THM-2101) -- b=sum_S a^k/Phi'(a) in F => b=0. Instantiates card_nsmul_translateSum_eq at G=Gal, Omega=Phi.rootSet, transitive from irreducibility; char 0. = additive analog of my S236 orbit-product concrete. REMOVES small-root product / Hensel / Vieta from the algebraic core.

**REMAINING:** (a) hfull (full-root sum=0) = codex's Lagrange + bridge (rootSet<->Finset, card_rootSet_eq_natDegree, Phi'=lc*nodal'); offered to coordinate. (b) b=1 wrapper: sum_{S_+} a^(M-1)/Phi'(a) = F(t) = sum D_m t^m => b=1 under D_m=0.

**Honest (kind-pasteur verdict stands):** b=1 still selects the small-root packet by valuation, so additive SHARES the valuation/Newton-polygon core with THM-1550 -- does NOT fully escape analysis; BUT replaces the product (Hensel factorization, Mathlib gap) with a SUM (partial-fraction residue), a cleaner target. Algebraic core now kernel-pure + product/Hensel/Vieta-free. cf mac-mini-S163. Artifacts: reflection the-concrete-root-packet-lemma-...-boxeph-S238.md, HYP-8980, GMC2RootPacketConcrete.lean.

