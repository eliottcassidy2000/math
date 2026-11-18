# Universal two-line structural lift, and why heights must be chosen before lifting

The [rim exchange](INVERSE_METAGRAPH.md) gives an exact two-line repair when the wrong closing arc lies between two adjacent twins. Every count-changing one-arc edit can be put into this geometry by the same physical template. This is a structural lift, not universal arithmetic placement.

Suppose T,T' differ at a single arc with H(T)≠H(T'). One of the orientations has a Hamiltonian frame P using that arc. (If neither does, the two sets of Hamilton paths coincide.) Insert the tournament in the middle block of the cyclic substitution C3[TT_w,T,singleton c], w≥2. Let a,b be the first two vertices of TT_w. A frame of the whole host is
```
a, the other TT vertices except b in increasing order, P, c, b.
```
It contains the edited arc. Its endpoint arc a→b closes the backwards rim after frame reversal. Reverse a second frame going backwards around that rim and omitting the original edited arc. This flips just that arc and a→b. Swapping the adjacent twins a,b removes the latter change and gives exactly C3[TT_w,T',c]. These are two genuine merged metagraph lines, up to relabeling and converse, and can be traversed either way if the selected frame originally belongs to T'. The same template lifts a whole sequence of count-changing single arcs, one edge at a time.

Write h=H(T), B=f_T(2), and a_w=2^w−1. The ordinary count of the lifted tournament is
```
L_w = 2h + a_w(B−h).
```
This is the balanced three-cycle substitution with one singleton: the two main blocks have one or two runs, and their weights give 2h+(B_TT−1)(B−h). The replay independently counts the lifted graphs on both ends and the intermediate frames, e.g. 549,457→444,295→440,143 at w=2 from the 10,007/7,505 one-spoke pair. No monotonicity of the intermediate or of a lifted multi-edit path is assumed.

## The change of height is essential

Suppose the original endpoints h,h' satisfy an exact signed shortcut 2^v h'−3h=s with s=±1. Then for fixed w
```
2^v L'_w−3L_w = a_w(2^v B'−3B−s)+2s.
```
This cannot equal ±1 at the **same strip v**. To see it, a_w≥3 would have to divide the difference of ±1 and 2s, which has absolute value 1 or 3. The only possible case is a_w=3 and desired sign −s, forcing 2^vB'=3B. But B'≡(−1)^{|T'|}h' mod 3 and 2^v h'≡s mod3, so its left side is nonzero modulo 3. This rules out using generic padding as if it preserved an arithmetic shortcut. The all-width beta=59 example works because the twin vertices are already built into the exact count construction, not appended by this generic lift.

A target-first use of this geometry must therefore prescribe the **lifted** count and find a core pair that meets the two height equations. For each proposed core h and width w, the required B is fixed by L_w=2h+a_w(B−h); a merely congruent answer cannot replace it. The structural part (two actual lines) is solved once a physical pair is supplied. General supply of the needed core profiles remains open.

## The short-width boundary of the exact beta=59 lane

The core-built beta=59 family in the inverse atlas used a pair of twins and hence width w≥2. At width one the count request is 187→35. The separate order-eight source (4,17,34,7,109,203,143,31), with frame (0,2,5,3,1,4,6,7), has a literal black line with those exact endpoint counts. Thus the beta=59 shortcut has a genuine metagraph route at **every positive width**: one different scaffold for w=1, the two-line rim exchange for w≥2. The w≥32 exact distance statement remains specific to that latter pair of graph classes.

Replay: `python3 rim_replay.py`.
