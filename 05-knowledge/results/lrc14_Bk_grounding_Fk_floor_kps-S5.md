# B(k) grounding — the iid ceiling F(k) and the bounded-spread μ-floor (kind-pasteur-2026-06-18-S5)

Target: prove `inf_E μ(E) > 0` (then ∩ G_P), `μ(E)=meas{x: maxgap{frac(e_i x)}>2/7}`, `0∈E`, `k=|E|≤13`.

## The iid equidistribution CEILING F(k) = P(k uniform pts have maxgap>2/7), exact
`Σ_{j≥1}(−1)^{j+1}C(k,j)(1−2j/7)_+^{k−1}`:
- k=3: 1; k=4: 342/343=.9971; k=5: 2325/2401=.9683; k=6: 15125/16807=.8999; k=7: 13443/16807=.7998;
- k=8: .6846; k=9: .5689; k=10: .4621; k=11: .3688; k=12: .2904; **k=13: 3132376013/13841287201 = .2263**.

F(k) is the LARGE-SPREAD limit (orbit → iid on the full subtorus when the offsets are relation-free).

## Consecutive μ_k (grid Q=2520)
k=3:1, k=7:.394, k=8:.332, k=9:.282, k=10:.252, k=11:.215, k=12:.202, **k=13:.179** (≈ mac-mini 829/4620).

## Bounded-spread min μ (k-subsets of {0..k+4} ∋ 0, grid)
k=7: .368 at {0,2,3,4,5,6,8}={0..8}\{1,7}; k=8: .319; k=9: .221; **k=10: .173** at {0,1,3,5,6,8,10,11,12,14}.
The minimizers are **perforated near-APs** (confirms HYP-2585: consecutive is NOT extremal for k≥7).
Floor `inf_k μ_min(k)` ≈ μ_min(13) ~ 0.10–0.12 (pure); with G_P, mac-mini's min ρ* ~ 1/90 ≈ 0.011. **Positive.**

## Key realizations for B(k) (the proof structure)
1. **The 0-anchored gap μ_0(E) = meas{(1/7,3/7) empty} is the SMALL-spread mechanism, NOT uniform.**
   For consecutive it ≈ near-0 window 1/84 (= mac-mini's floor); for large spread it → 0 like 1/maxE.
   So neither the four-window (THM-528) nor the 0-anchored bound gives the uniform floor — both shrink
   with spread. The uniform floor MUST come from the BULK / equidistribution mass.
2. **Large spread ↛ iid automatically.** The orbit `{(e_i x)}` equidistributes on the subtorus
   `T(E)=Λ(E)^⊥`, `Λ(E)={m: ⟨m,e⟩=0}`. Relation-free E → full torus → μ→F(k). But structured high-spread
   E (e.g. `{0,1,N,N+1}` has `e_1+e_2=e_3`) → lower-dim subtorus → μ = a different subtorus integral.
   So B(k) is NOT "large spread → F(k)"; it is "min over the FINITELY many subtorus types (k≤13), each
   positive". The short-relation lattice structure is the real object.
3. **μ>0 always** (near-0 window, L2), but on a positive-dim subtorus the near-0 is measure-0, so the
   subtorus-integral floor is the genuine question — and it is what the bounded-spread finite check +
   the equidistribution-rate (Erdős–Turán) must jointly pin.
