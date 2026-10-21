# No deficit-six center at order ten

**FINITE-EXACT theorem:** the coefficient of z^10 in R6 is zero. This is a complete classification of the necessary marked q=2,m=5 frames, not a census of every order-ten tournament. No new center occurs.

## Complete reduction

Let T be a strong self-converse tournament on ten vertices with Delta=h-b=6. For any marked Hamilton cycle, the inherited cycle-cut theorem gives an exact rotation index q<=8 and n=qm with m odd. The divisors of10 leave only q=2,m=5; q=10 exceeds8, and q=1,5 give even m.

Label the cycle 0->1->...->9->0. Step-two rotation is an automorphism. Write A_j=2j and B_j=2j+1 modulo5. Each shore induces a cyclic tournament on five vertices, determined by choosing one direction for each of steps +/-1 and +/-2: four choices independently on each shore. Every cross arc is determined by five bits c_d=[A_j->B_(j+d)]. The marked cycle forces c_0=1 and c_4=0. Three cross bits remain free. Thus there are exactly4*4*2^3=128 marked frames, and every hypothetical center occurs among them after labelling its cycle.

Every generated frame contains the marked Hamilton cycle, so is strong. Its exact cut index is2: step-two rotation belongs to Aut(T), whereas step-one rotation has even order and cannot be a tournament automorphism. Different frames may represent the same unmarked tournament, which is harmless for proving emptiness. This packet makes no ordinary-class count claim.

## Counts and independent routes

For all128 frames the checker computes:
- H by subset dynamic programming on used-vertex sets and terminal vertices, and independently by DFS generating only directed Hamilton words.
- Aut(T) and anti-isomorphisms by degree-compatible bijection backtracking with arc-consistency pruning.
- The same two map sets independently by exhausting the120 bijections on each five-vertex shore, filtering internal arcs, then testing each product against cross arcs.
- Blue labelled word count B both by testing the actual positional reverser of every generated Hamilton word against the anti-map set, and independently by enumerating2^5*5! half-words for every fixed-point-free involutory anti-map. Its other half is uniquely forced. Each blue word determines its own reversal, so the sum has no duplicate credit.

The score-shore restriction in the independent map check is proved, not guessed. If each A vertex beats k vertices of B, full scores are2+k and7-k, which differ since5 is odd. Automorphisms preserve shores and anti-isomorphisms interchange them. Thus the restricted exhaustive maps lose no possibilities.

The free action of automorphisms on ordered Hamilton words gives h=H/|Aut(T)| and b=B/|Aut(T)|. Self-converseness is checked by the existence of an anti-map, not inferred from score symmetry.

## Result

Exactly64 of the128 marked frames are self-converse. Their deficit histogram is:

- Delta1528:8 frames.
- Delta1628:8 frames.
- Delta2908:12 frames.
- Delta2926:12 frames.
- Delta3064:12 frames.
- Delta3114:12 frames.

The minimum is1528, so none has deficit6. This concludes the order-ten obligation. The histogram counts marked frames, not isomorphism classes.

All literal adjacency rows and automorphism/anti-map lists are retained in the receipts, together with H,h,b,Delta and a blue-word witness when one exists. A minimum-deficit frame has the two shore masks0,0 and cross bits(1,0,0,0,0), with (H,|Aut|,h,b)=(8265,5,1653,125).

## Reproduction and limits

Run `python3 -B check.py` and `python3 -B -O check.py` in this directory. Fresh normal/optimized receipts retain actual UTC clocks, elapsed time, and source SHA256; semantic content is identical after dropping timing fields. There are1,088 always-active checks in each run. The frozen degree47 proof packet is unchanged.

This result settles only coefficient10; it does not promote the inherited order-nine coefficient into a freshly audited result, and it does not settle any higher unknown coefficient or literal metagraph adjacency.
