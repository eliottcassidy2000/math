# Four consecutive odd counts with matched continuation data

A selective port can read a true +2 without keeping its full continuation state, as in [OUTLET_WIDTH](OUTLET_WIDTH.md). Working backward from a window that can survive every common continuation gives a stronger positive result: four literal seeds realizing consecutive odd counts with the **same** right interface. This construction uses several orders and different selective masks rather than forcing a common core, common initial port, or identical tail width.

## The exact target condition

Retain the [four-coordinate cut state](../chanpak_selective_target_20261114/CUT_STATE.md) (h,N,S,d): ordinary spine count h, marked-apex count N, weighted marked starts S, and its odd-minus-even part d. For an appended ordinal singleton with mark b, the transfer is
```
(h,N,S,d) → (h,N+(1−b)(S+d),2S+bh,bh).
```
Thus to obtain a width-r unit-step window preserved by every common transitive suffix, solve for r different N in steps of two with exactly the same (h,S,d). Every common step then adds the same displacement to all N. For a shared ordinal suffix block of Hamilton count a, the earlier cut transfer instead multiplies the gaps by a; keeping the unit step specifically calls for a=1. This is an actual absolute interface, not marginal agreement modulo a finite modulus.

Here are four witnesses. Rows are outgoing bitmasks (bit i corresponds to vertex i), and the mask describes the initial vertices the apex beats. The appended ordinal words precede one final common 0.

| Spine rows | Mask | Initial (h,N,S,d) | Tail before common 0 | State before common 0 |
|---|---:|---|---|---|
| (48,29,25,1,8,30) | 25 | (15,27,53,11) | 1 | (15,27,121,15) |
| (48,29,25,1,8,30) | 16 | (15,29,19,5) | 11 | (15,29,121,15) |
| (18,20,17,7,8) | 1 | (15,31,19,3) | 11 | (15,31,121,15) |
| (48,29,25,1,8,30) | 1 | (15,33,19,1) | 11 | (15,33,121,15) |

The replay computes the initial signed weights directly. The initial d values do not enter the displayed masked transfers because the alignment tails contain only ones. The common 0 adds 136 to each N and doubles S, resetting d=0. Both underlying spines are strong, the masks nonempty, and the common final 0 supplies a return arc to the apex from the end. Hence all four initialized glues are strong. The shared right state is (h,S,d)=(15,242,0), with N=163,165,167,169.

Write t for the canonical integer address of any subsequent transitive suffix (including 0 for the empty suffix), and p(t) for the exact [binary suffix correction](../chanpak_collatz_exact_addresses_20261114/BINARY_SUFFIX_DECODER.md). The common suffix increments every count by 242t+15p(t). Thus the four consecutive odd values
```
163+242t+15p(t)+{0,2,4,6},       t≥0
```
are literal Hamilton counts for **every** such t. For a nonempty suffix use the convention ending in 0 to retain strength. No low-part size restriction is imposed. This is one window translated along infinitely many addresses; it is not a proof that those translations cover every location or can be extended by one more odd count.

A smaller pair uses the single spine (6,8,10,1) with h=5, masks 8 and 6, and respective alignment words 11 and 1. Their states before the common 0 are (5,9,35,5) and (5,11,35,5). After that 0 the common right state is (5,70,0) and N=49,51, so both values 49+70t+5p(t) and 51+70t+5p(t) are realized for every t.

## Keep the preparation height free: the backward signature is oddpart(S+h)

An all-one tail of length w does not change N; it sends S+h to 2^w(S+h) and sets d=h if w≥1. Thus at fixed h a whole list of initial counts can be synchronized this way **if and only if** the painted bases share the odd part of S+h. To work backward from an anticipated consecutive window, compute that odd-part list for each base count and intersect them; then choose the physical mask and the neutral all-one width. This is an exact compiler for a specified base bank, not a search over suffix words.

In the four witnesses above, h=15 and the initial S+h values are 68,34,34,34. Their common odd part is 17. At every j≥3 independently choose the all-one widths so the matched S+h is 17·2^j, always with at least one appended 1 on every arm. A common 0 then adds 17·2^j to all four N values, resets d=0, and leaves S=2(17·2^j−15). Thus for every j≥3 and t≥0 the whole physical window is
```
27+17·2^j(2t+1)−30t+15p(t)+{0,2,4,6}.
```
The previous formula is the j=3 section. The exact two-parameter family is the advantage of fixing the *whole* right interface: further common right growth need not consist solely of zeros or of a gap-separated pattern.

## Why the reverse alignment search is finite at fixed input cores

For h fixed and any colored transitive word of length l and usual binary value b, the start readout is S'=2^l S+hb, 0≤b<2^l. Suppose two states with the same odd h reach equal S values by possibly different suffixes. Split the longer word at its excess length r. After that prefix, let its start value be S*. At equal remnant lengths l, start equality gives 2^l(S*−S_other)=h(b_other−b_long). Since h is odd, h divides S*−S_other; since |b_other−b_long|<2^l, their difference in absolute value is less than h. Thus S*=S_other, followed by the **same** remaining word. The alignment prefix, if one exists, has a uniquely forced value (S_other−2^r S)/h. Its length is bounded by 2^r S≤S_other when S>0. If S=0, the marked set is empty, (N,d)=(h,0), and one can strip arbitrary neutral initial zeros; the first nonzero digit supplies the finite bound instead.

Once S aligns, a common first 1 resets d with no change in the gap; a common first 0 resets d while adding the old d difference to the gap. If d already agrees, the empty remaining word works. This gives an exact all-length finite test for right matching from any fixed painted base bank. It also explains why the construction can tolerate changing the widths: the longer arm aligns before any common continuation must start.

Replay: independently compute the painted starts, biases and full Hamilton counts, including several literal common tails, with `python replay_matched_windows.py`.
