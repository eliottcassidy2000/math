# Independent review: fixed-prefix synchronous nonregularity

**PASS on the corrected frozen snapshot.** The complete binary function graphs are nonregular, including eventual graphs and nonempty complete tagged finite banks. The theorem does not imply a nonregular range, a scalar omission, or failure of a selected-witness strategy.

## Final source binding and correction

Source directory: `session197-prime-atoms/fixed-prefix-nonregularity/`.

- PROOF.md: `fe0c822872555e54db26f313967fd97cdf0939737f39bd379c5a432dbaec8058`.
- verify.py: `7bd4525faf311fa5de6d1e65b51a85f69fcca82daa6ce6b21c614fd161f35116`.
- Both source receipts and copied normal/optimized replays: `3c5de3d665ef2b2bc049c7dcf45c4d0f57690ccf413ea6af48545e9558da22df`.

The first review identified a real formal normalization gap: after projecting away a longer auxiliary number track, canonical intersection alone does not remove the common leading-zero columns left on the retained tracks. It can discard the desired pair. The author repaired both compositions by explicitly taking the regular left quotient by(0,0)* before intersecting with canonical encodings. The final theorem also states integer coefficients explicitly. A forced-padding control was added. The final proof and verifier above incorporate these repairs and pass. No frozen author file was edited by the reviewer.

The recurrence and physical final-gadget law are inherited from the already reviewed ARITHMETIC_INVERSE.md, SHA256 `4ac1727eb27c85f8a39f5598eb9f0e1b020bc4020cbca0a94e6288cb0676c3b7`. The nonregularity argument itself is an all-order language proof, not a finite automaton experiment.

## Balanced-language argument

Writing f(n)=g(n+1), the append-one step leaves f unchanged, while appending a zero adds the current prefix plus its parity. Starting from the all-one word1^a, only the first appended zero has an odd current prefix. This yields exactly

    f((2^a-1)2^b)=2^(a+b)-2^a-2^b+3.

For a=b=k>=3, its binary form is1^(k-1)0^(k-1)11. For unequal a,b, setting s=min(a,b),t=max(a,b) gives1^(s-1)01^(t-s)0^(s-2)11. Summing the positional values verifies both forms. The constraints s>=3 and t>s make every separating run positive, so the unequal form cannot pass the equal-case output filter. Both numbers have the same a+b digit length on this slice.

Thus intersection of a hypothetical regular f graph with the two regular track filters, followed by input projection, gives precisely{1^k0^k:k>=3}. The pairwise distinguishing suffixes0^i for prefixes1^i prove that language nonregular. This establishes the contradiction at all lengths, independently of numerical tests.

## Arithmetic relations, padding and eventual graphs

Addition by one is regular using a finite carry from the least-significant end. Reversal transfers it to the stated most-significant-first convention. For fixed integer coefficients, the equation y+e epsilon=hz+RL with epsilon=L mod2 is also regular. One explicit least-significant-first carry update is

    c_next=(c+y_i-h z_i-R L_i+e epsilon*[i=0])/2,

allowed only when the numerator is even. Start and finish at carry zero. Its carry magnitude stays bounded by1+h+R+e, so no unbounded register is concealed here. A finite state also distinguishes the first digit. Reversing this regular relation gives the required reading direction.

When auxiliary tracks are eliminated, allow common zero padding first, project, take the regular left quotient by common-zero columns, and then impose canonical pair encoding. This restores exactly the intended numerical graph. For example L=1,z=1 and the(13,2,18) prefix give auxiliary y=29: projection leaves four common zero columns before the pair(1,1). Canonical intersection without the quotient would lose it. The repaired proof handles this case and its arbitrary-length analogue correctly.

Consequently a regular g graph would yield a regular f graph through L=n+1. A regular fixed-prefix F graph would yield the g graph by the fixed-coefficient equation, because h>0 makes z unique and the assumed nonnegative F supplies a legal y. Both contradict the balanced slice. Removing a finite initial domain interval cannot repair the contradiction, since its finitely many canonical pairs form a regular language and could be restored. Fixing a tag in a nonempty complete tagged bank recovers one forbidden component graph.

The scope exclusions are essential and correctly stated: projection of a nonregular relation can be regular; an untagged union need not reveal a component; a proper witness subrelation need not encode every pair. Finite-dimensional arithmetic with unbounded values and the target-dependent square-root inverse are not synchronous finite-state graphs.

## Independent finite controls

The standalone `check.py` imports no source code. Its f computation accumulates contributions from zero digits while scanning a binary word. For small L, a separate physical one-apex count assigns each transitive-deletion vertex before or after the root and checks the two boundary incidences directly.

The independent normal and optimized runs agree on:

- 2,048 literal physical gadget counts;
- 65,025 exact large-integer two-run bit forms, with255 balanced filter hits;
- 129,024 fixed-coefficient carry checks, including false-output and false-recovered-value rejection;
- 13,543 projected pairs that actually require common-zero trimming;
- 8,192 explicit addition-by-one carry controls.

The physical and carry controls corroborate different parts of the algebra and normalization. The infinite distinguishability proof, not the sample size, establishes nonregularity.

Hashes:
- check.py: `add9211e493c12459a62ec28ea0a15a1033c5b16002b49bd18886154233eba05`.
- independent-normal.json and independent-optimized.json: `e9fc7b0fd5bd5fa3f98140c81ff2d42578fcfc9e05177c384736d4f864cfb959`.

Run `python3 check.py` and `python3 -O check.py`; neither writes files. The copied source verifier likewise writes no files, and both source modes reproduce the final frozen receipt exactly.
