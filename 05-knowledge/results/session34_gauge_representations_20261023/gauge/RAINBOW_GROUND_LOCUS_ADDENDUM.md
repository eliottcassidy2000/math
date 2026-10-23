# Full rainbow action: exact ground locus and its nongauge code symmetry

Status: PROVED by finite linear-algebra certificate; FINITE-EXACT partition histogram. Addendum to SIMPLICIAL_WILSON_FANO.md, preserving that main file for concurrent audit.

Retain its fixed Fano lines, symmetric link signs, Paley baseline, vertex-cut space B, Fano line space L, and line-code C=ker(syndrome). Embed C into edge space via full-line flips, denoting its image L_C. Let R(u) be the sum of21 normalized rainbow four-cycle Wilson loops, now on ALL2^21 edge fields.

**Theorem.** The span W of the21 rainbow edge-cycle vectors has dimension11, and

W^perp = B direct-sum L_C.                                    (A1)

Therefore R(u)<=21, with equality exactly when u/u^0 belongs to B+L_C. The full model's1024 ground fields consist of16 ordinary vertex-gauge orbits, represented precisely by the zero-syndrome Fano orientations. The particular energy -R has extra nongauge symmetries L_C on the entire configuration space, not merely on the Fano slice. This does NOT identify L_C as the physical vertex gauge group or assert every observable is invariant under it.

Proof: every cycle annihilates cuts. Every rainbow cycle restricts on L to one of the seven dual simplex-code characters, so it annihilates L_C. Earlier L intersect B=0 proves dim(B+L_C)=6+4=10, hence rank(W)<=11. The following11 original rainbow masks, with edges ordered lexicographically as combinations(range(7),2) and least-significant bit edge(0,1), have independent XOR echelon pivots:

8273,8707,524577,525321,34826,36870,1114148,1179668,134208,147648,295552.

The elimination certificate's (pivot index,row mask) pairs are
(20,1114148),(19,524577),(18,295552),(17,196656),(16,68720),
(15,34826),(14,19584),(13,8273),(12,6156),(10,1320),(9,594).

Distinct leading pivots prove rank>=11. The script independently enumerates geometric rainbow cycles, verifies these rows and computes the reduction, so this is an exact finite certificate rather than a float rank. Dimension equality proves(A1). A sum of21 signs is21 iff every sign is+1, exactly the annihilator condition. QED.

The complete unnormalized counting partition function is

Z_R(beta)=1024 e^(21 beta)+21504 e^(13 beta)+143360 e^(9 beta)
          +322560 e^(5 beta)+688128 e^beta
          +637952 e^(-3 beta)+282624 e^(-7 beta).               (A2)

These counts sum to2097152=2^21. Enumeration fixes all six root edges, checks all32768 rooted fields and multiplies every count by64. Background normalization means one may enumerate relative fields directly; translation by the Paley baseline is a bijection. A second normal/-O run is identical. This exact full-field census is NOT an enumeration of only8192 constrained fields from the main note.

The quartic action therefore extends the Fano code symmetry to a larger ambient system, while distinguishing it from ordinary gauge symmetry. On the selected128 Fano representatives only R=21 and R=-3 occur. The full model also has -7,1,5,9,13; discarding those values is a genuine model restriction, not gauge fixing.

Reproduction: `python3 full_rainbow.py > full_rainbow.json`; rerun with `-O`. The main note's proposed first next experiment is now completed by this addendum. No finite mass eigenvalue or continuum claim is inferred from a finite energy landscape.
