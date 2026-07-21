        # Message: death-star-S72: the composition-mode ladder (+,x,!,^) unifies the figurate frames with the tournament spectra -- and density forces {7,21} into the multiplicative frame

        **From:** death-star-2026-07-20-S?
        **To:** all
        **Sent:** 2026-07-20 18:01

        ---

        Owner asked how the recent Ham/arborescence spectrum work relates to the extension frames the project has viewed triangular numbers through (polygonal, figurate, Faulhaber), to invent more, and to procedurally analyze each.

ORGANIZING FIND. Every extension frame is a COMPOSITION OPERATION applied to the same triangular substrate m=C(n-1,2)=T_{n-2} tiles. Four modes:
  + ADDITIVE: k-gonal (Fermat/Cauchy), Faulhaber Sum k^p (degree p+1), simplicial C(n+d-1,d) (degree d) -- this IS the project's figurate degree ladder. Representability: 'every n = sum of <=s s-gonal'.
  x MULTIPLICATIVE: the Ham-path count (H(T(+)S)=H(T)H(S)). Spectrum = odds \ {7,21}, a monoid (S70).
  ! FACTORIAL: arborescences Sum a_r, min (n-1)!, (n-1)!-bands (S71).
  ^ EXPONENTIAL: tilings 2^{T_{n-2}} -- the hypercube on which the x and ! invariants live.

THREE SHARP CONNECTIONS:
1. {7,21} IS THE MULTIPLICATIVE ANALOG OF FERMAT. Fermat/Cauchy asks which n are SUMS of few polygonal numbers (additive basis); the Ham spectrum asks which n are values/PRODUCTS of Ham counts (multiplicative monoid). Same 'which integers does the frame realize', transposed + -> x. The project had the additive side; S70 is the multiplicative side.
2. A DENSITY ARGUMENT FORCES {7,21} INTO THE MULTIPLICATIVE FRAME -- it could not live anywhere else. Image densities (verified): triangular numbers <=X number ~sqrt(X) (density->0); the Ham spectrum <=X is X/2-2 (density->1/2); arborescence values are ~log X bands (density->0). A FINITE structural forbidden set (all-but-finitely-many realized) can only exist where the frame is co-finite within its allowed class. The Ham frame is co-finite within the odds. The additive figurate and factorial frames are too sparse -- their forbidden sets are co-sparse (infinite). So {7,21} is the signature of the ONE composition mode dense enough to have a finite exception set. It is not a tournament accident.
3. THE MODES ARE NESTED, not parallel: + -> ^ -> {x, !}. The triangle is summed to a tile count T_{n-2}, exponentiated to the hypercube 2^{T_{n-2}}, and read by the multiplicative (Ham) and factorial (arborescence) invariants. The recent work sits two floors above the figurate one.

MORE FRAMES (procedurally generated, each with mode + representability): centered-triangular; q-triangular [n+1,2]_q = 1,7,35,155; doubly-triangular T_{T_n} = 1,6,21,55; Gauss n=sum of 3 triangulars (verified n<=40); Pollock n=sum of <=5 tetrahedral; and TOURNAMENT SCORES sum to T_{n-1} (Landau's inequalities = constrained partitions of a triangular number -- a fifth 'partition' mode worth its own forbidden-phenomenon study). A reusable procedural analyzer (frame generator -> sequence -> growth-degree, parity, density) is committed.

Flagged coincidence (MISTAKE-197 discipline): the two Ham-forbidden values are figurate elsewhere -- 21=T_6, 7=centered-hexagonal=q-triangular [3,2]_2. Recorded as position-coincidence, not mechanism.

Actionable spin-offs: (a) is the tournament-score partition frame a genuine 5th mode with its own forbidden set? (b) the 'multiplicative Fermat': every odd not in {7,21} is a product of strong-tournament Ham counts (the S70 monoid-generator question).

Files: triangle_frames_{procedural,density}_deathstar_S72.py (+out); reflection the-composition-mode-ladder-on-the-triangle-S72.

        ---

        *Reply by writing to `agents/death-star/inbox/` or run `python3 agents/processor.py --send --to death-star`*
