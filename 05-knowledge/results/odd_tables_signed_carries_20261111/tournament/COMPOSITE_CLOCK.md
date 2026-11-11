# A341-phase clock gives an all-prefix H25-to-C3 paired reset

**PROVED IN NOTE, 2026-11-11:** an actual simultaneous source rebase
driven by a composite residue clock, a seven-bit arbitrary-prefix
preparation, and uniform improvement after changing the retained root.
**FINITE-EXACT:** the stated small controls. **OPEN:** completion of an
arbitrary prescribed target. This extends [PROOF.md](PROOF.md).

## 1. The actual sources determine the arithmetic channel

Use the inherited H25 cores and apex cuts

    C-=(50,8,3,5,46,14),     cut38;
    C+=(50,40,35,5,46,8),    cut35.

At left address0 and odd right address B=2n+1, their counts are

    P(n)=25g(n)+244n+186,       P(n)+2.             (1)

Their shared apex-root state has E=d=25 and

    U=147·2^bitlength(B)-25rev(B).

The C3 pair from [PROOF.md, equation1](PROOF.md) has counts
`N(C)=3g(C)+12C+10`, N(C)+2, using cuts3/1 and right addresses2C+1/4C+3.
Choose C=19n. The two actual target differences are the same function

    D19(n)=P(n)-N(19n)
          =25g(n)-3g(19n)+16n+176.                (2)

Define the shifted phase

    F19(n)=D19(n)+22(n mod2)-18.

The coefficient25 is the H25 deletion count, not a formal rescaling
of the earlier H13 formula. The named graphs and cuts in (1) remain
part of the construction throughout.

## 2. A synchronized composite clock adds four to the defect

**Lemma1.** For every x>=0,

    D19(1024x+3)=F19(x),
    F19(1024x+3)=F19(x)+4.                        (3)

Proof. Ten-bit concatenation gives

    g(1024x+3)=g(x)+255x+(x mod2),
    g(19(1024x+3))=g(19x)+7581x+(x mod2)+22.

Substitution proves the first identity. The appended number is odd,
so F19 at that address is D19+4, proving the second. QED.

The block is binary0000000011. Its periodic fraction is

    3/(2^10-1)=1/341,       341=11·31.

Multiplication by19 takes it to the ten-bit block0000111001, with
periodic fraction19/341. The orders of2 modulo11 and31 are10 and5,
so their common CRT clock has period10. The two component phases are
synchronized by the same digit word; two independently selected cycles
would not supply (3).

The separate g-drift rates per block are86/341 and262/341. Their
weighted difference is

    (25·86-3·262)/341=4.

This explains why a composite modulus supplies a useful operation:
the exact weighted defect survives the residue projection and is
changed by a known integer. It is not an inference from modular
reachability to exact integer coverage.

## 3. Seven preparation bits suffice for every prefix

**Lemma2.** For any p>=0, the two prepared prefixes q0=128p+27 and
q1=128p+28 satisfy

    F19(q0)=25g(p)-3g(19p+4)-1068p+22(p mod2)-44,
    F19(q1)=25g(p)-3g(19p+4)-1576p+22(p mod2)-182. (4)

Both are strictly negative, and exactly one is divisible by4.

Proof. The equalities are seven-bit expansions of the binary recurrence.
Both values are even. Their difference is508p+138, which is2 modulo4.
Using g(p)<=p+1 and g(19p+4)>=1 bounds the first expression by
`-1043p+22(p mod2)-22<0`, and the second is still smaller. QED.

Select the q whose value is divisible by4 and put m=-F19(q)/4>=1.
Define

    n=q·1024^(m+1)+3(1024^(m+1)-1)/1023.           (5)

**Theorem1 (arbitrary-prefix H25 rebase).** The integer n in (5)
satisfies D19(n)=0. Thus both actual H25 targets P(n),P(n)+2 have
the C3 witnesses at C=19n. For p>0 the word of n is

    prefix(p) · (0011011 or0011100) · (0000000011)^(m+1).

Proof. Apply the first identity of (3) to the terminal block and the
second identity to the preceding m blocks. The resulting defect is
F19(q)+4m=0. Equation(2) supplies both target equalities and the named
graphs supply the actual realizations. QED.

The empty-prefix case p0 selects q28, with F19(q)=-172 and m43.
The resulting address is large. This is an exact completion algorithm,
not a claim of efficiency in the prefix bit length. As before, the
completed target differs from the target at the unpadded prefix.

## 4. Every equality in this channel has a beneficial root choice

The H25 odd-address pair satisfies

    U-P(n)>=58,       U-(P(n)+2)>=56.              (6)

One way to retain the exact quantifiers is to start before the right
word: the two states have H73, U147, d25 and E9/7. Their invariant
`U-H-d+E` is58/56. O increases this invariant by U-E>=0 and Z preserves
it. A nonempty canonical right word ends in Z, making E=d and proving
(6). This is the same physical recurrence used for the H13 bound.

Whenever `P(n)=N(C)` with C>=2, reroot each C3 witness at its last
right-flank vertex as in the main note. Write C0 for C with its leading1
removed, and epsilon for its second leading bit. The resulting states
are

    (P(n),   N(C0)-3epsilon,   P(n)+3epsilon,   N(C0)),
    (P(n)+2, N(C0)+2-3epsilon, P(n)+2+3epsilon, N(C0)+2).    (7)

**Theorem2.** The replacement in (7) preserves the two exact current
targets and improves each old H25 q state componentwise, where
`q=(H,U,U-d,H+U-E)`. The savings in its last three coordinates are at
least55,43,40 on the lower branch and53,43,40 on the upper branch.

Proof. With delta=U-H and the new deletion count d', old minus new is

    (0,delta-3epsilon,
       delta+d'-25-3epsilon,
       delta+d'-25-6epsilon).

Use (6), d'>=13 on the lower branch, and d'>=15 on the upper branch.
QED. This comparison holds for any equality P(n)=N(C); C=19n is one
new constructive channel producing such equalities.

Every common future O/Z/J program or fixed right-root amalgamation
preserves weak H dominance by the inherited order theorem. As in the
H13 rebase, the new pair difference is `2(1,1,1,1)`, not pure H. Its
general common-address gap is2g(address), and an all-Z continuation
preserves the unit gap2. Neither common deletion ports nor arbitrary
unit-gap continuation is asserted after this root change.

## 5. A moderate literal meeting and the remaining gate

The small input n890 has

    B=1781,       C=16910,
    P(n)=N(C)=231421.

The old H25 states are

    (231421,25,265981,25),
    (231423,25,265981,25).

The rerooted C3 states are

    (231421,7213,231421,7213),
    (231423,7215,231423,7215).

The q savings are `(0,34560,41748,41748)` and
`(0,34558,41748,41748)`. This is a concrete source change with a large
future-response improvement and the same two current targets.

The theorem still does not say that a prescribed arbitrary odd target
lies in this completed language. It supplies a synchronized actual
family, a short residue preparation, and a root choice with a uniform
future comparison. The missing all-target inverse remains the same
construction obligation stated in the main note.

Reproduce:

    python3 -B 05-knowledge/results/odd_tables_signed_carries_20261111/tournament/composite_clock.py

The receipt records128 exact identity/preparation inputs and four full
completed prefixes. The n890 graph certificates are counted by literal
endpoint subset DP on18,18,20,21 vertices. Old U is also counted from
root-start covers; the new U uses the proved one-out-neighbour identity
after independently counting H,E,d. No coverage range is rescanned.
