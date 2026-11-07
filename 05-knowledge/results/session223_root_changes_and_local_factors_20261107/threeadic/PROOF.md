# The mod-3 derivation for fixed-length two-root-change words

Status: PROVED for the fixed-length ensemble below. The all-precision extension is proved in [ALL_K.md](ALL_K.md). The subsequent ALL_K.md extends this result to every 3-adic precision and gives the complete fixed-length valuation factor. The height-truncated local factor in HYP-9494(E) remains OPEN. No fit to the census means enters the proof.

## 1. Definitions and inheritance

Use the physical operators, initial state, and F from THM-4904, `05-knowledge/results/root_change_carmichael_chamber_20261107/PROOF.md`, section 0:

O(D,E,U,d)=(D,0,2U,d), Z(D,E,U,d)=(D-E+U,d,2U-d,d),
J(D,E,U,d)=(D-E+U,D-E,D+U,D), t=(1,0,2,1), F=D-E+U.

The parent packet's independent `audit/AUDIT.md` fixes the source-prefix/final-close convention. Every J here is in the source prefix. The final close is not part of the sampled word. This is word multiplicity, not distinct-state multiplicity.

Closest mechanism: finite physical-state transfer plus the exact O/Z dynamics. Hostile: the height restriction is not a length restriction. Corrected near miss: the packet's mod-3 class minima prediction failed, so a class-mean theorem must not be read as a minimum theorem. Sidecar retained: length and J count. The mod-3 quotient discards height completely.

## 2. The theorem

Choose a word uniformly from the binomial(L,2)*2^(L-2) words of length L containing exactly two J's. As L tends to infinity,

Pr(F=0 mod 3)=34/81+O(1/L),
Pr(F=1 mod 3)=47/162+O(1/L),
Pr(F=2 mod 3)=47/162+O(1/L).

Thus the residue local factors relative to uniform residue mass are

(34/27, 47/54, 47/54).

For exactly one J the corresponding limiting probabilities are

(41/108, 8/27, 35/108).

For no J the limiting distribution is uniform. In particular, the 1-versus-2 asymmetry at one J vanishes in the leading fixed-length two-J law. Finite-length asymmetry need not vanish identically.

## 3. O/Z Markov projection, proved directly

All coordinates in this section lie in F_3. Let P be the transition matrix choosing O or Z independently with probability 1/2. Let Pi be its Cesaro stationary projection.

Define pi_a, for a=1,2, to be uniform on the eighteen states

(D,E,U,d)=(D,e*a,U,a), with D,U in F_3 and e in {0,1}.

Define eta to be uniform on the six states (D,0,U,0) with U nonzero, and delta_c the point mass on (c,0,0,0).

For every state s=(D,E,U,d), its projection is:

- if d is nonzero, s Pi=pi_d;
- if d=0 and U is nonzero, s Pi=eta;
- if d=U=0, s Pi=(delta_D+delta_(D-E))/2.

Proof for d nonzero: one step puts E in {0,d}. Every state in the eighteen-state set has exactly two incoming O/Z edges, so uniform measure is invariant. The E=0 slice is strongly connected: O acts there by S(D,U)=(D,-U), and the word ZO acts by T(D,U)=(D+U,U+d). The permutations T^{-1} and STS differ by the nonzero D-translation d. Thus their generated group contains all D-translations and, using T, all U-translations. Every permutation inverse is realized by a positive power on this finite set. Applying Z reaches all states of the E=d slice. Hence the chain is irreducible. At E=U=0, O is a self-loop, so the chain is aperiodic. Any starting E enters this chain after one step.

For d=0,U nonzero, one step sets E=0. Each step negates U. On the six-state set, two-step words OO and ZO respectively fix D and translate D by U. Each parity class is therefore irreducible and aperiodic for P^2. The whole chain has period two and Cesaro projection eta. For d=U=0, the first O or Z leads to the two stated absorbing points. This proves the formula, including transient states.

## 4. Exact rational calculation

Let R be deterministic application of J. Directly from the projection formula,

pi_a R Pi = (pi_1+pi_2)/3 + (2/9) eta + (1/12) delta_0 + (1/36) delta_(-a).

Indeed D=0 and U=0 has probability 1/9. In that event the first O/Z after J lands at delta_0 with probability 3/4 and delta_(-a) with probability 1/4. D=0,U nonzero has probability 2/9 and projects to eta; D=1 or 2 each has probability 1/3 and projects to pi_D.

Also,

eta R Pi=(pi_1+pi_2+eta)/3,
delta_0 R Pi=delta_0,
delta_c R Pi=pi_c for c nonzero.

F is uniform under each pi_a and eta. Under delta_c, F=c. Since t Pi=pi_1, the one-J probability of residue zero is 8/27+1/12=41/108, residue one is 8/27, and residue two is 8/27+1/36=35/108.

After the second J and projection, the absorbing mass at zero is

(2/3)*(1/12)+1/12=5/36,

and the absorbing mass at each nonzero residue is 1/108. The remaining mass 91/108 has uniform F. Therefore the zero probability is 91/324+5/36=34/81, and each nonzero probability is 91/324+1/108=47/162.

## 5. Why this is the fixed-total-length limit, not only independent long blocks

Put n=L-2. Every word has the unique decomposition A J B J C with O/Z block lengths a,b,c satisfying a+b+c=n. Each triple has exactly 2^n words, so the length triples are uniform among binomial(n+2,2) possibilities. The modular law is exactly

[binomial(n+2,2)]^(-1) sum_(a+b+c=n) t P^a R P^b R P^c.

Here t denotes a point-mass row vector. Its generating function before normalization is

t (I-zP)^(-1) R (I-zP)^(-1) R (I-zP)^(-1).

The initial d=1 chain is aperiodic and mixing, so the first row resolvent is pi_1/(1-z) plus a function analytic on a disk of radius greater than one. For the full finite chain, the only peripheral eigenvalues are 1 and -1, as the preceding classification proves, and they are semisimple. Thus at z=1 the coefficient of the order-three pole is pi_1 R Pi R Pi. All other terms have poles of order at most two on the unit circle, including z=-1 because the first row resolvent has no pole there. Their coefficients are O(n); poles strictly outside the unit circle give exponentially decaying terms times fixed polynomials. The order-three pole has coefficient binomial(n+2,2). Division proves the claimed O(1/L) error. This argument also proves the one-J statement with the two-resolvent expression.

The possible period-two component must not simply be ignored: the initial mixing block is precisely what prevents a competing order-three pole at -1.

## 6. Exact transfer at every 3-adic precision

For any q=3^k define integer matrices A_q=O_q+Z_q and R_q=J_q on (Z/qZ)^4; transitions are interpreted with multiplicity. If e_t is the initial row vector and h_r selects F=r, then

sum_(a+b+c=n) e_t A_q^a R_q A_q^b R_q A_q^c h_r

is exactly the number of length n+2 words with exactly two J and F=r mod q. Equivalently it is the coefficient of z^n in

e_t (I-z A_q)^(-1) R_q (I-z A_q)^(-1) R_q (I-z A_q)^(-1) h_r.

This is a rigorous modular counting identity at every precision. ALL_K.md subsequently proves the higher-adic asymptotic by classifying the O/Z communicating classes; no k>1 mixing assertion is needed for the shorter mod-3 proof here.

## 7. Precise unresolved gap and strongest legitimate consequence

M_2(N)=[N=1]+2 B_0(N)+2 B_1(N)+B_2(N). The theorem concerns B_2 words at fixed length, not the sum of M_2(N) over a height window. The event X<=F<2X reweights both O/Z lengths and the real normalized state within each block. Neither Pi nor the modulo-q transfer records this real state. A proof that this height weighting has the required conditional residue equidistribution, or an explicit joint real/modular renewal theorem, is missing.

Consequently the theorem does not derive the observed height-truncated factors in HYP-9494(E), does not imply coverage or a multiplicity lower bound, and does not derive height-truncated valuation-by-valuation factors. ALL_K.md does derive the valuation-by-valuation factors for the fixed-length ensemble. It does rigorously identify a nonuniform local law in a specified natural ensemble and isolates the singular mechanism: J can send D=U=0 mod 3 to the absorbing zero state, retaining mass through a subsequent J.

## 8. Reproduction and controls

Run `python3 verify.py` in this directory. `VERIFY.txt` records exact rational projection calculations, invariant-distribution and strong-connectivity checks, self-loop checks, a fully literal word enumeration control through length nine, and independent exact integer transfer counts through length 1000. Large-length decimal readouts summarize exact integers; they are diagnostic, not proof inputs. Exploratory floating-point counts were not used by the proof or controls and are not included as evidence.

Independent reviews are in `../audit/SCOPE.md`, `../audit/THREEADIC_AUDIT.md`, and `../audit/HIGHER_AUDIT.md`.
