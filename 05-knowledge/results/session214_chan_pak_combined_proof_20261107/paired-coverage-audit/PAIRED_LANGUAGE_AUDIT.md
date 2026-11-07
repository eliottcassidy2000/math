# Session214 paired-carrier exact-language audit

**Verdict: PASS.** For the four fixed carrier pairs677,759,1547,821, the synchronized common-word language is exactly {O,Z}*. Any first J permanently separates adjusted bases under every common suffix. No mathematical correction was required.

## Frozen author bytes

Staging packet: `session214-paired-coverage/session214/paired-carrier-language`.

- README.md: `a2a8e238274fae507df1a8674d305fba5150590418f0a6e109a5b9924a16bf3e`
- SHA256SUMS: `54fd3ad6f530d883056c5cf8e15bf2a97ecbed69bf5cf76614d51f4d20d88a53`
- PROOF.md: `cd8505a54bc7f30dae3ff88ca29cda33acb05988ea8fd7e17672969d03b1c042`
- verify.py: `6f957916e0cb7af137741a8a9d1602635dcce1fad8cc98cd4736ac1765fd5c1f`
- RESULT.json: `9f12776c3c9aca0f012d2d7897110e44129b8bee9cb35a17f3054e174c6ab881`

The full author SHA256SUMS manifest was successfully rehashed. The frozen rooted-closure and paired-carrier sources are explicitly inherited, not claimed as new results.

## All-order proof check

The linear change to cover coordinates(L,C-2L,T-b,b) gives precisely the three displayed nonnegative matrices. Every new L contains the old L, so strict sign of L persists within either closed orthant. The separate ordinary transformation(a,f,e,p) was independently expanded from the inherited h,e,R,d law and agrees term by term.

Pairs759 and1547 have identically zero ordinary differences. Every O/Z prefix keeps L=b=0 and both other cover differences negative. The first J makes L negative and enters the negative orthant, which remains invariant.

For677 the prefix identities L=-1-u and b=-v remain exact under O and Z. The strict inequalities v<u<0 and C,T<0 persist. The first J makes every cover coordinate strictly negative, and the ordinary difference becomes(v-u,0,u,0). Its h is negative, and under Z/J its change R-e is nonpositive in that cone. Thus h+L+1 remains strictly negative; the constant offset has not been dropped.

For821, the first J gives ordinary difference twice the singleton physical profile(h,e,R,d)=(1,1,1,1). The ordinary laws do hold at singleton inputs; the Z singleton exception concerns only s, which this argument never uses. Every suffix therefore retains positive h. Its cover difference enters the strictly positive cone, giving positive h+L+1 forever.

The adjusted-base offsets follow from the actual fixed s and j values: zero for759/1547 and one for677/821. All source graphs are nontrivial and their start values stay fixed under the three physical operations. Empty words and arbitrary O/Z prefixes are covered. The theorem concerns common words on these fixed pairs, not asynchronous words, different starts, new sources or a new operation alphabet.

## Independent computation

`paired_language.py` directly enumerates all directed simple vertex words for the eight literal seeds and their depth<=2 children. Their path and ordered-cover profiles are assembled from those actual words. This is a separate implementation from the author's subset recurrence.

- 104 independently recounted literal graph profiles through order10.
- Every seed difference and both adjusted-base offsets verified.
- Explicit integer matrices exhaust all common words through length11:265,720 per pair,1,062,880 total.
- Exactly4,095 O/Z words per pair synchronize; every other word has the proved strict sign.

Reproduce:

    python3 paired_language.py > /tmp/session214-paired-language-audit.json
    cmp paired_language.json /tmp/session214-paired-language-audit.json

The all-order conclusion rests on the invariant proof, not the finite depth check. This audit does not cover the separate exploratory search9 program or any claimed coverage extension. Chan–Pak completeness and F(8,137) remain open.

The final README records the completed independent PASS while preserving the earlier candidate-status line in the unchanged proof. Its scope and replay counts agree with this audit, and the regenerated manifest rehashed successfully.
