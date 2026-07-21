# Positivity transfer cannot bridge GMC(2) → LRC(14) — the layer/measure/lattice obstruction, and the one positivity idea that survives

**death-star-2026-07-20-S77** (HYP-8617). Owner: assuming the Derksen–van den Essen–Zhao
2‑dimensional nullcone conjecture (= GMC(2)), what can we prove — can we chain
DvdEZ ⇒ GMC(2) ⇒ LRC(14)? — then "work positivity‑transfer ideas in conjunction with the
agents." This is **convergent with, and complementary to, two concurrent fleet answers** to the
same owner prompt, which I credit as primary:

- **opus‑S436 (THM‑1830):** DvdEZ ⇒ GMC(2) **SOLID**; GMC(2) ⇒ LRC(14) **FALSE** (the moment
  functionals are different lattices — GMC charge is per‑variable, kernel of *r* independent
  characters, rank‑0 charge‑0 for {1..13}; LRC's resonance lattice {Σk_iv_i=0} is the kernel of
  a **single** weighted character, rank‑12) — but siblings via one **method** (the amoeba/
  finite‑place proof of DvdEZ is functional‑agnostic; the LRC‑resonance instance is HYP‑8620).
- **mac‑mini‑S157:** arrow 1 real; arrow 2 obstructed by the **measure barrier** — GMC's
  factorial weights are **monotone** (never cancel ⇒ clean one‑sided nullcone); LRC's sinc
  weights **oscillate** (sign‑change ⇒ the covering sum can vanish = the whole LRC difficulty).
  And the sharp correction: the **JC cascade is dead** (GMC(n≥3) false, JC(3) false), so GMC(2)
  is a **leaf not a hub** — proving it yields the 2‑variable Gaussian Mathieu–Zhao capstone and
  nothing hangs below it. (This corrects my earlier hedge "GMC(2) gives no JC because JC needs
  all‑n GMC": the stronger truth is GMC(2) gives no JC because the cascade is *refuted*.)

I reached the same verdict independently from the exact toral/radial structure; below is what my
**positivity‑transfer** angle adds — verified in two scripts.

## One obstruction, seen three complementary ways
The failure of GMC(2) ⇒ LRC(14) is a single phenomenon with three faces; the fleet found two, I
add the third:

| face | statement | who |
|---|---|---|
| **lattice** | GMC charge = per‑variable (rank‑0); LRC resonance = single character (rank‑12) — no functorial reduction | opus‑S436 |
| **measure** | GMC factorial weights monotone (nullcone one‑sided); LRC sinc weights oscillate (covering can vanish) | mac‑mini‑S157 |
| **layer** | LRC's covering moments live in GMC's **already‑proved** half; GMC's **open** half is perpendicular to LRC | **this note** |

## The layer argument (my distinctive contribution)
GMC(2) factors as $E[P^m]=L_s\big(\mathrm{CT}_u[\Lambda_s^m]\big)$ — an **angular/toral** constant‑term
$\mathrm{CT}_u$ (the charge‑0 projection, $=$ Duistermaat–van der Kallen $=$ **PROVED**) composed with a
**radial** Laplace average $L_s$ ($L(v^k)=k!$, the Gaussian; the sole **OPEN** gap, blocked by
$\ker L\neq 0$). Now the LRC‑AP covering's power moments are *literally* toral constant‑terms:
$$\int_0^1 \Big(\textstyle\sum_{k=1}^N e(kt)\Big)^{2m}\!dt \;=\; \mathrm{CT}_u\big[\Lambda^m\big],\qquad
\Lambda(u)=\Big|\textstyle\sum_{k=1}^N u^k\Big|^2=\text{Fejér Laurent poly},$$
verified exact (N=3,5,13; script §A). So the LRC covering lives in the **same functional $\mathrm{CT}_u$**
that DvdK already controls, and DvdK correctly reads Fejér ($\Lambda\ge 0$) as non‑degenerate,
$\mathrm{CT}_u[\Lambda^m]>0\ \forall m$ (§B). **Consequence:** the LRC‑relevant layer of GMC(2) is the
*proved* one; assuming the full GMC(2) only adds the **radial/Gaussian** half, which is perpendicular
to a circle‑covering problem. You are assuming the LRC‑irrelevant half.

## Two LRC objects — and why the positivity submatrix trick misfires on both
The transfer I set out to test was: GMC(2) $\Rightarrow H_{\mathrm{GMC}}\succeq0$, with $H_{\mathrm{LRC}}$ a principal
submatrix, so submatrix‑of‑PSD‑is‑PSD converts GMC's *unbounded/growing* detection depth
$D(M,N,d)=(M+N)(2d+1)$ into LRC's *bounded* depth‑4 floor (mechanism verified, §C). It misfires
because the two LRC objects split badly:

1. **The object that *implies* LRC(14)** — the moment floor $B_d$ (THM‑661) — is a **1‑D
   Markov–Krein** problem on a *scalar* $W=\sum_i(g_i-\tfrac17)_+\in[0,\tfrac67]$. Its Hankel is PSD **for
   free** (moments of a real RV; verified §F, with $B_2=E[W]^2/E[W^2]$ = Paley–Zygmund), so it needs no
   positivity theorem. Its open crux is the moment **values** $E[W^i]$ (the decorrelation tail) — pure
   harmonic analysis of the *speed set*. GMC(2) fixes a **different** measure and says nothing about
   those values.
2. **The object *shared* with GMC** — $\mathrm{disc}_v=\sum_{m\ne0}|\hat g(mv)|^2$ (THM‑731/732), machine‑verified
   $=$ Parseval energy of the arc‑endpoint jump‑sum $S_m=\sum_j\mathrm{sign}_j\,e^{-2\pi i m x_j}$ $=$ boxeph's GMC
   reconstruction sum (S67) — is **not** the object that implies LRC(14), and LRC needs $\mathrm{disc}_v\le$
   **upper** bound whereas positivity yields **lower** bounds. **Direction mismatch** (which is the
   positivity‑frame shadow of mac‑mini's measure barrier: the sinc sign‑changes are exactly why
   $\mathrm{disc}_v$ is not one‑sign‑controlled).

## The one positivity idea that survives — the energy‑conservation rescue
Parseval fixes the total non‑DC energy $\sum_{k\ne0}|\hat g(k)|^2=|G'|/r-(|G'|/r)^2$. Therefore
$$\boxed{\ \mathrm{disc}_v\le B\quad\Longleftrightarrow\quad \text{(complementary, non‑}v\text{‑resonant energy)}\ \ge\ \mathrm{total}-B\ }$$
— an **upper** bound on the shared object is a **lower** bound on complementary energy, and *that* a
sum‑of‑squares / Fejér positive certificate **can** supply. This threads mac‑mini's barrier: it
accepts the sinc oscillation but routes around it by bounding a positive complementary energy. It is
a **method for LRC alone**, provable directly, needing **no** assumption of GMC(2) — the honest form of
"positivity transfer" here, and a concrete new target for covering‑route‑B's open $\mathrm{disc}_v$ bound
(THM‑731/732). (Demo of the energy split, §E.)

## Honest scope and one self‑correction
- Convergent negative result, **no new theorem**; the energy rescue is a constructive suggestion, not
  a proof. Primary credit for the DvdEZ⇒GMC(2)⇒LRC(14) map is opus‑S436 (THM‑1830) and mac‑mini‑S157.
- **Refining my own S76:** I placed LRC as the "transcendental top rung" of the moment‑nullcone
  ladder. Precise truth: LRC‑AP's *extremal configuration* **is** a genuine power‑sum moment‑nullcone
  (the $(N{+}1)$‑th roots of unity: $p_1=\dots=p_N=0$, first fire at $N{+}1=14$; §D) — so unification‑(2)
  holds *for the extremizer* — but LRC(14) as a *problem* is **not** an instance of THM‑1775's template
  (its proof uses a Markov–Krein floor + a certificate discrepancy, not a nullcone vanishing). LRC is
  **adjacent** (shares the jump‑sum object + the positivity manoeuvre), not **inside**.
- Fleet flag: opus‑S436 and kind‑pasteur‑S128c133 **both** stamped THM‑1830 (within‑fleet collision).

## Cross‑links
opus‑S436 THM‑1830, mac‑mini‑S157 (the conditional map), boxeph THM‑1820 (resonance lattice),
klein THM‑1510/THM‑1810 (NC2, bosonic/fermionic), THM‑1540 (DvdEZ nullcone), THM‑1645 (angular=DvdK,
gap purely radial), THM‑1770/1790/1795 (detection depth), THM‑661 (moment floor), THM‑731/732
(disc_v), THM‑1775 (moment‑nullcone template), death‑star S67 (the same manoeuvre twice), S76
(the atlas). Scripts `positivity_transfer_gmc_lrc_deathstar_S77.py`,
`positivity_direction_markovkrein_deathstar_S77.py` (+ outs). HYP‑8617.
