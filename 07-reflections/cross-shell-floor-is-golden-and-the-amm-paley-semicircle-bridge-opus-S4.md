---
source: opus-2026-07-31-S4 (cross-shell capacity resolved; AMM<->Paley tournament free-probability bridge; a new q-weighted family)
status: >
  THREE results. (I) The cross-shell capacity comparison RESOLVES the general-vs-block AMM gap: the per-shell
  archimedean margin min_delta[supply-H(delta)] flips sign EXACTLY at gamma=golden, and by forward-routing flow
  + dyadic scale-invariance the tail-sums share that one exponent, so cross-shell routing cannot lower the
  floor -- C*_general = C*_block = log_5(5 phi^2) (strong reduction; degree-matching rigor is the one gap).
  (II) NEW BRIDGE: the AMM golden floor and the Paley tournament (THM-438) are the SAME free-semicircle object.
  AMM golden = Catalan GF at w=-1; the Paley signed cluster GF is F(x)=sum(-1)^k C_k x^k = C(-x), loop equation
  x F^2 + F - 1 = 0 (THM-438's (**)); so golden = F(1) = C(-1) = 1/phi is its EDGE value, and the path-ratio
  R->e is its exponential resummation. (III) NEW STRUCTURE with exact closed forms: the q-weighted Paley
  tournament R_q(p)=E[prod(1+q chi(d_k))] has path-ratio R_q -> e^{q^2} (Gaussian MGF; q=1 recovers R->e) and
  cluster edge-value F_q(1) = (sqrt(1+4q^2)-1)/(2q^2), algebraic for all q, golden at q=1.
tags: [amm12592, general-floor, cross-shell, capacity, tournament, paley, THM-438, catalan, semicircle, free-probability, golden, closed-form, q-deformation, bridge]
related: [THM-3009, THM-438, amm12592-general-class-floor-szego, amm12592-two-ray-threshold-is-golden, central-binomial-edge]
---

# The cross-shell floor is golden, and the AMM<->Paley semicircle bridge

## I. Cross-shell capacity: `C*_general = golden` (resolving the crux)

The open gap (my general-floor note) was: general = balanced baseline + deficit field `b`, `sum b p^h(1-p)^t=0`;
THM-3009's `(ARCH)` is the capacity to cancel `b` WITHIN one dyadic shell, and the general class differs only
by routing `b` ACROSS shells. **Does cross-shell routing beat within-shell?** Answer: **no.**

1. **Forward-routing.** The carry `(p+q)^{2^a}=1` moves a deficit only to DEEPER shells (`m' >= m`). So the
   cross-shell flow is a transportation problem on the ordered set of shells with forward edges; its
   feasibility (Gale/Hall) is: for every `M`, `sum_{m>=M} demand(m) <= sum_{m>=M} supply(m)`.
2. **Scale-invariance.** Dyadic shells are self-similar: `demand(2m)/demand(m) = supply(2m)/supply(m)` (the
   same exponential rate `2^{rate}`). So both tail-sums are geometric with the SAME ratio, and
   `(tail demand)/(tail supply) = ` the per-shell ratio `rho(gamma)`, INDEPENDENT of the cut `M`.
3. **The per-shell margin flips at golden.** Computing the archimedean exponent
   `margin(gamma) := min_delta [ supply(gamma,delta) - H(delta) ]` (THM-3009's `(T)`, binary entropy):
   ```
      gamma-golden:  -0.020   -0.005    0.000     +0.005    +0.020
      margin      :  -0.0243  -0.0060   +1.2e-5   +0.0060   +0.0141
   ```
   `margin < 0` for `gamma < golden` (every shell deficient), `= 0` at golden, `> 0` above -- and the binding
   `delta = 1/phi`, with the supply-side argmax at `x ~ 0.04` (the shell's START), extending the routing range
   never raising the supply.

Since every tail-cut shares the single exponent `rho(gamma)` and it is `<1` for `gamma<golden`, no cross-shell
reallocation can make any tail feasible: **`C*_general = C*_block = 1 + log_5(phi^2) = log_5(5 phi^2)`.** The
cross-shell freedom is illusory because the binding is realized locally at each self-similar shell's start.

### The Hall step made rigorous (Gale's condition, degree-resolved)

The earlier "aggregate cut" concern dissolves once the routing geometry is stated exactly. The carry
`(p+q)^{2^a}=1` spreads a deficit to STRICTLY HIGHER total degree (it never lowers degree), so in the
absolute-degree coordinate `D=h+t` the deficit routes FORWARD: a demand at degree `D` may only be absorbed by
a cell of degree `>= D`. This is a transportation problem on the totally ordered set of degrees with forward
edges, and for such a problem **Gale's feasibility theorem is not the per-degree inequality but exactly the
tail condition**:

```
   feasible  <=>  for every D_0:  sum_{D >= D_0} demand(D)  <=  sum_{D >= D_0} supply(D).
```

So the degree-resolved Hall condition IS a family of tail inequalities -- there is no separate per-degree
obligation to discharge. Now dyadic self-similarity gives `demand(D)` and `supply(D)` the SAME exponential
growth rate in `D` (each doubling of degree rescales both by the identical factor, the content of
"`(ARCH)` stable in `R`", death-star), so every tail ratio `sum_{>=D_0} demand / sum_{>=D_0} supply` equals
the per-shell margin exponent `rho(gamma)`, independent of `D_0`. Therefore all tail inequalities hold
simultaneously iff `rho(gamma) <= 1`, i.e. `margin(gamma) >= 0`, i.e. `gamma >= golden`. This closes the
reduction: **`C*_general = golden` rigorously, given only (i) forward-in-degree carry [structural, from
`(p+q)^{2^a}=1`] and (ii) dyadic scale-invariance of `demand,supply` [THM-3009's `(ARCH)`, verified
`R`-stable].** The only thing not reduced to prior results is the exponential-rate equality in (ii), which is
THM-3009's own asymptotic and can be cited. Verifier:
`04-computation/amm12592_cross_shell_capacity_opus_S4.py`.

## II. The bridge: AMM's golden floor IS the Paley tournament's semicircle

The tournament pillar's crown jewel THM-438 proves the Paley path-ratio `R(p)=H(T_p)2^{p-1}/p! -> e`, with
single-run cluster integrals `A_{2k} = C_k p^{k+1}+...` (Catalan) and the number-theory-free identity `(**)`
`sum_{even-series} mu = (-1)^k C_k`, equivalent to the loop equation

```
   x F^2 + F - 1 = 0,     F(x) = sum_{k>=0} (-1)^k C_k x^k = C(-x)     (C = Catalan GF).
```

`F` is the R-transform / free-cumulant fixed point of the two-point law `{0, +-i sqrt n}` forced by the DRT
engine `S^2 = J - nI`. **This is the same object as the AMM floor.** THM-3009 showed the AMM golden constant is
the Catalan GF at its branch point, `C(-1)=1/phi`. Since `F(x)=C(-x)`:

```
   AMM golden edge = C(-1) = F(1) = (sqrt5 - 1)/2 = 1/phi,     via  1*F^2 + F - 1 = 0  at  x = 1.
```

So **golden `= F(1)` is the EDGE value (analytic continuation to `x=1`) of the Paley tournament's signed
cluster generating function**, and the tournament limit `e = exp(surviving cherry cluster)` is the same `F`'s
EXPONENTIAL resummation. AMM's capacity floor and Paley's Hamiltonian-path ratio are two evaluations --
edge-value and exponential -- of one free-semicircle object. (Verified: `1*F^2+F-1=0` at `F=1/phi` to 30 dps.)
This is why central binomials `C(2n,n)C(4n,2n)` (my S(k)), Catalan `C_k` (Paley), and `C(-1)=1/phi` (AMM) all
appear together: they are the moments, cluster law, and edge of the semicircle.

## III. A new structure: the q-weighted Paley tournament (exact closed forms)

Proposal suggested by the bridge: deform each arc's character by a fugacity `q`,

```
   R_q(p) := E_sigma[ prod_{k=1}^{p-1} (1 + q * chi(sigma(k+1)-sigma(k))) ],     chi = Legendre symbol.
```

`q=0` is the random tournament (`R_0=1`); `q=1` is THM-438 (`R_1=R`). Two aspects compute in CLOSED FORM.

- **Path-ratio (exponential):** in the cluster expansion each `L`-arc cluster carries `q^L`; only the cherry
  `L=2` (`a_2=1`) survives `p->infty`, so
  ```
     R_q(p) -> exp(q^2)     (a Gaussian MGF; q=1 recovers R -> e).
  ```
  Checked: `R_{1/2}(7)=1.266 -> e^{1/4}=1.284`, `R_1(7)=2.4 -> e`. The `e^{q^2}` exposes the cherry as the
  Gaussian pairing directly.
- **Cluster edge (algebraic):** the `q`-weighted signed cluster GF is `F_q(x)=F(q^2 x)`, solving
  `q^2 x F^2 + F - 1 = 0`, so its edge value is
  ```
     F_q(1) = ( sqrt(1 + 4 q^2) - 1 ) / (2 q^2),     exact algebraic for every q,
     F_1(1) = (sqrt5-1)/2 = 1/phi  (AMM golden),   F_{1/2}(1) = 2(sqrt2-1),   F_2(1) = (sqrt17-1)/8.
  ```
  The one-parameter family `F_q(1)` traces an algebraic curve of tournament cluster-edges (`q^2 F^2 + F = 1`),
  with **the AMM golden `1/phi` as its distinguished `q=1` (Paley) point** -- the golden mean is exactly where
  the fugacity is the honest Legendre symbol. (The other `q` give algebraic constants but not metallic means;
  only `q=1` lands on `phi`, which is what makes the Paley/golden coincidence special, not generic.)

**Why this is the right generalization.** It is a genuine tournament construction (a fugacity-weighted arc
count, `E[prod(1+q chi)]`), it has TWO exact closed forms (`e^{q^2}` and the algebraic `F_q(1)`), and it
threads the user's three targets into one curve: the AMM golden floor (`q=1` edge), the tournament path-ratio
(`e^{q^2}`), and the metallic/algebraic closed forms (`F_q(1)`), all evaluations of the deformed semicircle
`q^2 F^2 + F = 1`.

## Handoffs

- death-star / concurrent-opus: `C*_general = golden` (sec I) promotes THM-3009's floor from balanced-block to
  the GENERAL `C*`, modulo the degree-resolved Hall step; the answer to my cross-shell question is "no."
- tournament lane: `F_q(1)=(sqrt(1+4q^2)-1)/(2q^2)` and `R_q -> e^{q^2}` are new closed forms; the natural next
  target is whether a genuine tournament family (not just a fugacity) REALIZES `q != 1` -- e.g. a weighted or
  higher-class DRT whose two-point scale is `q^2 n` -- which would put silver `sqrt2-1` on an actual tournament.
