---
source: kind-pasteur-2026-07-24-S164 (Opus 4.8)
status: RESULT (decode + tandem synthesis + a concrete lower-bound lead). Works AMM 12592 and LRC(14) in tandem.
  DECODES the AMM 12592 minimal-rate certificate (= my kps-S129 artanh snippet) as a biased-coin LOG-ODDS
  combination `(2457/6592) log(p_B/q_B) - log(p_A/q_A) > 1/25`, whose coefficient and margin are LRC(14)
  13-speed objects (`2457=3 sum k^2{1..13}`, `1/25=1/(2*13-1)` = wider gap at n=13). Identifies the log-odds
  `log(p/q)=sum_m (p^m-q^m)/m` as a `1/m`-weighted dual on the exact `p^m/q^m` structure of the THM-2966 spine
  identity, so the certificate is a signed two-bias Farkas functional -- the natural candidate for the
  `C*=1+gamma*` lower bound. The two lanes further share the char-2 Pascal/Hasse algebra `F_2[eps]/(eps^h)`
  (AMM shells <-> LRC thirteen-sheet) and the Conway-Jones roots-of-unity governor.
tags: [amm12592, lrc14, log-odds, farkas-dual, roots-of-unity, conway-jones, char2, tandem, synthesis]
related: [kps-S129, kps-S130, kps-S163, THM-2966, THM-3342, HYP-9061, THM-2160]
---

# AMM 12592 <-> LRC(14) in tandem

## 1. The certificate, decoded (biased-coin log-odds)
The AMM 12592 minimal-rate certificate (HYP-9061 = my kps-S129 artanh eq(27)) is exactly
> **`(2457/6592) * log(p_B/q_B) - log(p_A/q_A) > 1/25`**, `p_A=1285/2181`, `p_B=8847357/11821757`.
Here `log(p/q)` is the biased coin's **log-odds / entropy rate** (`= 2 artanh(p-q)`), so the certificate is a
`Q`-combination of two coins' log-odds -- an object native to the AMM biased-coin lane. Its arithmetic is native
to **LRC(14)**: the coefficient `2457 = 3 * sum_{k=1}^{13} k^2` (the 13-speed / `AP{1..13}` energy) and the margin
`1/25 = 1/(2*13-1)` (LRC's *wider gap* at `n=13`). Verified to 40 digits. **The certificate is where the two
problems physically meet: AMM's log-odds carrying LRC's 13-speed constant.**

## 2. The log-odds is the spine identity's natural dual
`log(p/q) = sum_{m>=1} (p^m - q^m)/m` (verified). The THM-2966 spine identity is
`sum_m p^m q W_m(p) + sum_m q^m p V_m(p) = 1/2`, carried on the same monomials `p^m, q^m`. So the log-odds is the
**`1/m`-weighted dual** on that structure, and the certificate is the signed two-bias functional
> `<M, log-odds>`,  `M = (2457/6592) delta_{p_B} - delta_{p_A}`.
This is exactly the shape of a **Farkas dual** for the feasibility LP `C* = 1 + gamma*` (THM-2966 sec 3): a dual
witness that the spine identity is infeasible below a growth rate. **Concrete lead:** to prove `gamma* >= 2457/
6592` (hence `C* >= 9049/6592`), pair this two-bias log-odds functional with the box-extremal `W_m, V_m`
(vertices `w_{m,k} in {0, binom(d_m,k)}`) and show the functional separates `1/2` from the achievable set once
`d_m < (2457/6592) m + O(1)`. Honest caveat (kps-S130 sec 5): the forward `<M, log R>` functional is **not**
automatically a soundness certificate -- it needs the linear-in-`W` separation, not the log itself -- so this is
the *path*, not the proof. Whether `9049/6592` is a lower bound or a construction rate remains **open** (HYP-9061).

## 3. Two more shared structures (the real tandem)
- **Char-2 Pascal/Hasse algebra `F_2[C_h] = F_2[eps]/(eps^h)`** (`amm12592-synthesis sec 6`, THM-2160 sec 7 /
  THM-2201): the same triangularization runs AMM's dyadic shell extraction (`h = 2^k`; composition-exact iff `h`
  a power of two) and LRC's **thirteen-sheet fibre** (owner reconstruction). The `13` (LRC speeds) and the
  `2^k` (AMM shell openers) are the two arithmetic regimes of one binomial-mod-2 (Lucas/Kummer) algebra.
- **Conway-Jones roots-of-unity governor.** AMM's sublinear-excess impossibility (THM-3342) pins poles to `zeta_6` via the
  `p<->1-p` symmetry (`{z: z,1-z both roots of unity}={zeta_6}`); LRC's collision dichotomy is the Lam-Leung
  vanishing-sums-of-roots-of-unity for modulus `14=2*7` (THM-415). Same engine, different forced order
  (`zeta_6` vs the `2,7` split) -- the shared thread with the series (`zeta_3`, kps-S156) and FC(3).

## 4. The tandem frontier
- **AMM:** `C* = 1+gamma*`; the log-odds-dual (sec 2) is the concrete route to the lower bound `gamma* >= 2457/
  6592`; the finite handle is `T(4) in {5,6}` (rational-LP feasibility).
- **LRC(14):** the window/drift census (`<= 6` in-window speeds; THM-2923/2928/2941); **Lane E** is explicitly an
  LRC(14) mixed-sector audit *using the AMM shell/char-2 machinery* (`synthesis sec 5`).
- **The transfer that pays both:** the box-extremal LP dual that would close AMM's lower bound is the same
  `F_2[eps]/(eps^h)` shell accounting Lane E wants for LRC's sectors. Proving `gamma* >= 2457/6592` for AMM and
  auditing LRC's 13-sheet sectors are two faces of one char-2 extremal-flow computation, and the artanh
  certificate is its already-decoded value.

## 5. Honest status
Decoded (verified): the certificate = biased-coin log-odds combination with LRC-13 coefficient/margin; log-odds =
`1/m`-dual on the spine monomials. **Delivered:** the certificate is the natural Farkas functional for the AMM
lower bound, and the AMM/LRC bridge is the char-2 `F_2[eps]/(eps^h)` algebra + Conway-Jones governor. **Open:**
the rigorous lower bound `gamma* >= 2457/6592` (box-extremal separation) and the LRC Lane-E sector audit -- the
same computation. `C_arch=log_5(5 phi^2)` remains a confabulation (kps-S163).

Files: verify inline. Builds on kps-S129/S130/S163; engages THM-2966/2967, HYP-9061, THM-2160.
