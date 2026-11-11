# Extension: single-letter binary returns fail at every seed length

**PROVED, 2026-11-11.** This extends the finite scan in `PROOF.md`; no general two-apex count obstruction is claimed.

Fix core a_1->a_2. Let a physical seed word w have count h>=3. There do not exist two distinct letters p,q in {0,1,2,3} such that every continuation over these two letters follows the exact digit rules h->2h+1 and h->2h-1, with the signs assigned once and for all to p,q.

## 1. Three live segments defeat a binary count cylinder

Suppose a fixed continuation prefix can initialize all three spine segments for some core order pi, and a fixed continuation suffix can make the two constrained last endpoints legal. Then between those fixed blocks insert m arbitrary copies of any allowed letter. Each can be assigned independently to any of the three initialized segments. All earlier first endpoints and later last endpoints remain valid. Consequently the resulting tournament has at least3^m Hamilton paths.

Its continuation length is m+C, where C is fixed independently of m. In a proposed binary cylinder every such word has at most

    2^(m+C)(h+1)-1

paths, by iterating h->2h+1. This is incompatible with3^m for sufficiently large m. Every object is physical: the assignments count different literal Hamilton paths of a single chosen word, not an average over words.

## 2. Eliminate the six unordered alphabets

If letter3 is present, it appends a universal sink, since both core vertices and all earlier spine vertices point to it. It leaves H unchanged. For h>=3 this is neither2h+1 nor2h-1. Thus all pairs containing3 fail immediately.

For {1,2}, both bits take both values. Choose core order(a_1,a_2). After the seed, initialize the first segment with any new letter, the middle segment with1 (a_1 points to it), and the last segment with2 (a_2 points to it). After the free block, close the first segment with2 (it points to a_1) and the middle segment with1 (it points to a_2). Assign all old seed vertices to the first segment. This activates the three-segment argument and rules out the alphabet.

For {0,1}, bit1 (the a_2 bit) is always0.

* If the seed contains a letter with bit1=1, use core order(a_2,a_1). Put all seed vertices preceding the selected vertex in segment0 and the selected vertex and all later seed vertices in segment1. Its first endpoint is legal. Start segment2 with a new1. Initialize segment0 if necessary, then after the free block close segment0 with0 and segment1 with0. These are legal endpoints into a_2 and a_1 respectively. Thus three segments are live, a contradiction.
* Otherwise a_2 is a universal sink in the entire seed and all its continuations. Remove it. The remainder is one apex a_1 over the spine. Appending1 is a universal sink in this remainder, so H again stays fixed. This contradicts either required digit rule at h>=3.

For {0,2}, bit0 (the a_1 bit) is always0.

* If the seed contains bit0=1, use core order(a_1,a_2), start the middle segment at a selected such seed vertex, assign earlier seed vertices to segment0 and later ones to segment1, and initialize segment2 with2. Close segment0 with0 and segment1 with0 after the free block. This again gives three live segments.
* Otherwise all spine vertices point to a_1. Deleting a_2 gives the transitive chain consisting of the spine followed by a_1. Thus the tournament is the one-apex word s0, where s records the a_2 bits of the old spine. Appending0 or2 to the original spine changes this one-apex word to s00 or s10, respectively. This is a genuine coordinate identification, not a scalar analogy.

## 3. The residual one-apex calculation

For a binary word s, let A be its binary integer, p its last bit (p=0 for the empty word), and g=H(s). The inherited elementary one-apex append laws are

    H(s1)=g,       H(s0)=g+A+p.

They can alternatively be checked directly from the k=1 compiler. Hence

    H(s00)=g+3A+p,
    H(s10)=g+2A+2,
    H(s00)-H(s10)=A+p-2.

For the two outputs to be2h+1 and2h-1, their difference must be plus or minus2. Thus A+p is0 or4. Since p=A mod2, the only possible A are0,3,4.

Now append the original letter2. This replaces s by s1, so A becomes2A+1 and p becomes1. The necessary value A+p becomes2A+2, namely2,8,10 in the three cases. None is0 or4. Therefore the two required digit rules already fail at the following step. This eliminates the last alphabet, completing the proof.

## 4. Scope and next representation

The obstruction is stronger than the depth8 finite seed scan: it holds at all seed lengths and does not assume a particular invariant-line representation. It does **not** eliminate:

* multi-letter digit blocks;
* state-dependent digits or changes of core order/mark between digits;
* constructions that replace a spine segment rather than append;
* a finite union of different physical presentation types;
* two-apex count completeness by a non-binary synthesis algorithm.

The mechanism is a resource/spectral mismatch: a stationary free alphabet either opens a third independent path slot, forces an inert append, or leaves a one-apex arithmetic coordinate that prevents a return. Keeping physical segment resources is what turns the finite failed scan into an all-size theorem.
