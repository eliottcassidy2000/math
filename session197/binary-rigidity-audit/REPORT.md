# Independent binary-rigidity review

**PASS after one boundary clarification, with no formula or verifier correction.** Bound to `../../session197-orbit-wildcard/binary-rigidity/PROOF.md`, corrected SHA256 `65102468958b198f8538022288238c40d4e3ca33684cf0d419bc86b0993681a9`; verifier `9cb8898f86656c60d55460b66b817611af6bde4db64f9e06772f4db76070e92c`; receipt `2741f8b773e725348be5e0b3867df5d96ca22de0a11628e1ba0c5b48dac58b81`.

## All-order argument

The root-fixed anti-involution pairs O with I and makes the crossing matrix symmetric. Matching rank one really does force a single diagonal crossing edge. The exterior-cycle proof then has a unique backward arc, fixing each cycle's cut and eliminating rotational multiplicity. All odd external cycles decompose into one I-side path, an increasing clone subset, and one O-side path. Remaining odd-cycle collections independently give the two endpoint responses, so the square formula h^2+q^2*2^m follows for every positive m.

The original wording asserted that every path/subset choice produced a cycle. The singleton/singleton/empty-clone choice instead asks for the absent joining arc and would have length2. The corrected proof now asserts the converse only for odd-total-length choices and explicitly identifies this filtered-out boundary. This does not change the parity count or any formula.

The general rank-one7p exclusion is valid for weak and strong sources alike: modulo7, -2^m is a nonsquare, so a multiple of7 forces both h and q to be multiples of7. The resulting lower bound rules out49 as well. Unit offset then forces H(S)=1, hence transitive S, and q=2^j. The unique single-backward-arc core, its middle root, equal transitive padding and the converse construction all follow exactly. Eventual identities become exact by finite exponential independence. The root-fixed and positive-clone hypotheses remain essential.

## Independent computational controls

The independent checker constructs each source as the ordered join S^op -> root -> S and reverses just the distinguished cross edge. Its endpoint response uses a reverse subset recurrence, not the source's literal path-prefix enumeration.

-44 side sources, including nontransitive side graphs through order5.
-88 physical clone graphs counted using literal odd directed cycles and disjoint packings, rather than the source's Hamilton endpoint DP.
-29 complete scalar shift operators checked by the independently developed odd-cycle skeleton compiler. The positive bases are exactly1:h^2 and2:q^2.
-1,099 symmetric crossing matrices through size4 checked by direct disjoint-edge tests.
-15 complete transitive core/padding controls and294 residue controls.

The scalar operator may retain a base-zero term, already present for C3 as -f(0). The independent comparison therefore discards only that term on the asserted nonempty polynomial ideal. This confirms why the source correctly makes no m=0 deletion claim.

Source normal and optimized runs reproduce the frozen receipt. The independent checker also agrees normally and with `-O`. Its SHA256 is `b914e2fa35e87e13ba401d0fc19173096f8d7595e1e39f3589e832a90b7851d6`; its receipt SHA256 is `57a0ed3b3f94dc965977dc6273d15311c134eecee37ceea63c3ff1a67a6efe32`.

The source files were not modified by this audit. The author's single explanatory correction is explicitly bound by the corrected proof hash above. No claim about unmarked self-converse roots, arbitrary nontransitive root modules, or the unrestricted Hamilton spectrum is added.
