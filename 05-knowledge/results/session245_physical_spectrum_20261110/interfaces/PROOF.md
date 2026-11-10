# Reversal bridges cancel physical response errors

**PROVED conditional composition law with explicit unconditional physical instances. FINITE-EXACT discovery. Arbitrary anchored interval repair and Chan–Pak coverage remain OPEN.**

## Inheritance and concept board

The inherited source is Session244 `composition/INTERFACES.md`, especially its seven-coordinate ordinal product and raw JJJJ/JOZJ module. Session233 THM-5000 already composes coprime physical moving cores into long intervals, with exponential baselines. Session199 already uses converse sources in `source-norm-capacity/dependencies/C3_TWO_SOURCE.md`; merely using converse is not a novelty claim. The scoped contribution here is exact cancellation of the raw response errors using an opposite-facing module, with a physical witness. The Session244 fixed-chain two-adic obstruction applies to similarly oriented E=1 modules and does not cover this architecture.

Board: (1) exact anchored repair, still open; (2) full-interface rather than scalar composition; (3) converse exchanges response directions; (4) separator kills unwanted endpoint interaction; (5) coprime deletion does not itself create intervals. Broad targeted searches in canon and Session199/Session233 found the inherited mechanisms above, but no matching numerical witnesses. This is not a claim of literature novelty.

## 1. Interface and converse

Write an interface as x=(h,q,A,G,T,S,E), with conventions of Session244. Its converse is

    x*=(h,q,G,A,T,E,S).

Reverse each path in an ordered cover and exchange the two path positions. This preserves T and exchanges A,G and S,E, proving the identity for every physical interface. The selected cut is retained, not complemented.

For nonempty modules the inherited ordinal law is

    xy=(hj,qs,As+hB,Fq+jG,AF+Tj+hL+EV,Sj,hW)

for y=(j,s,B,F,L,V,W). A transitive chain of t>=1 unselected vertices has interface (1,2^t,0,0,0,0,0). The restriction t>=1 matters: the empty graph is not a genuine identity for the endpoint-coordinate formulas. No zero-length separator is used in the theorem below.

The protected closure adds an unselected sink and then a root beating exactly the selected cut. Its full rooted state (count, ending-at-root, two-path-start response, deletion count) is

    (h+2A-T, h, h+2A, h).

Indeed the sink forces the deletion Hamilton path to end outside the cut, the extra sink doubles A and leaves T unchanged, and the inherited one-root formulas give the displayed coordinates.

## 2. Raw module data

Let source i have D_i,E_i,U_i,d_i and higher response W_i, as defined in Session244. Delete the terminal root in the words JJJJ and JOZJ. These yield the same underlying deletion graph, and two singleton cuts with interface change

    delta x_i=(0,0,E_i,-4W_i,0,0,-U_i).

Their common h_i=3D_i-E_i+3U_i and q_i=12Q_i+14W_i-6V_i are inherited identities. Denote the lower-state interface coordinates by A_i,G_i,S_i,F_i (F_i means interface E, to avoid confusion with source E_i). Both raw T values are zero because the same singleton cannot lie on two disjoint paths.

## 3. Exact positive cancellation bridge

Place raw module 1, then t>=1 unselected transitive vertices, then the converse of raw module 2. Compare the simultaneous choice 00 with 11. Suppose

    2^t E_1 q_2 = 4 h_1 W_2.                   (C)

Then their h,q,A,S,E coordinates coincide, while

    delta G = -4 h_2 W_1 + 2^t q_1 E_2,
    delta T = 2^t (A_1 E_2 + E_1 A_2 + E_1 E_2).   (D)

Here A_2 is the original raw module's lower-state A, hence the converse module's G. In particular, if E_1,E_2>0 then delta T>0, and protected closure gives a pure count decrement delta T, with all three other rooted coordinates unchanged.

**Proof.** After the separator, module 1 has A=2^t A_1, G=G_1, T=0, E=0 and S=S_1. The converse second module has A=G_2,G=A_2,T=0,S=F_2,E=S_2. Thus product A=2^t A_1 q_2+h_1 G_2, G=2^t q_1 A_2+h_2 G_1, and T=2^t A_1 A_2. Its S=S_1 h_2 and E=h_1 S_2 are unchanged under either raw switch. Taking simultaneous differences gives (C) and (D), including the product E_1 E_2 term. Protected closure then gives the asserted state equality. All statements concern actual graphs and cuts, not freely assigned integer states.

This cancellation is a useful positive escape from the earlier similarly oriented-chain obstruction. It does not imply an arbitrary prescribed value of delta T or its baseline.

## 4. Literal witness, fully recounted at order19

Source1: outgoing bit rows [8,1,3,6], root0, with (D,E,U,W,h,q)=(5,1,5,17,29,580).
Source2: rows [24,1,3,6,14], root2, with (11,3,17,72,81,2088).
Set t=2. Equation(C) reads 4*2088=4*29*72.

Product interfaces are

    (2349,4844160,300440,240077,13328,81,58),
    (2349,4844160,300440,241529,14140,81,58).

Their protected closures are strong tournaments on19 vertices with rooted states

    (589901,2349,603229,2349),
    (589089,2349,603229,2349).

The gap812 is not divisible by3. The executable constructs both actual adjacency lists, checks strong connectivity, and independently recounts both full Hamilton-path counts and endpoint counts with subset DP. This is not merely verification of the product formula.

## 5. Coprime-deletion witness

Source1 rows [44,21,48,54,33,2], root2 has (D,E,U,W,h,q)=(31,7,45,236,221,6880).
Source2 rows [58,20,25,34,40,6], root5 has (33,7,41,224,215,7072).
Again t=2, and 4*7*7072=4*221*224.

The product interfaces are

    (47515,194621440,7947160,7446740,272384,860,884),
    (47515,194621440,7947160,7436420,287196,860,884).

Their protected strong closures have counts15669451 and15654639, gap14812, deletion47515 (coprime to3), common endpoint47515 and common response15941835. Each graph has22 vertices. The verifier literally constructs these graphs and checks strongness, and recounts each constituent interface. The order22 total is proved by the exact composition law; no independent order22 full subset-DP recount is claimed.

## 6. A smaller cut-cell and composition consequence

A direct cut search found the fixed deletion core C=[18,8,3,5,14]. Independent interface enumeration verifies

    cut3: (13,104,32,41,6,5,6),
    cut6: (13,104,32,44,8,4,6).

Protected closures are strong counts71 and69, sharing rooted sidecars (13,77,13). This unit gap2 is stronger than the reversal witnesses for immediate repairs, but its baseline is still fixed. It is not a full-interface pure-T switch: delta G=3 and delta S=-1 remain.

Append an unselected sink to EACH deletion core before ordinal composition. Every resulting module has h13,q208,A64,E0, and switch delta=(0,0,0,3,2,-1,0). In an arbitrary fixed chain the exact scalar coefficients are

    -h_suffix(i) [3 A_prefix(i)+2 h_prefix(i)].

The S discrepancy does not enter because every preceding product has E=0. There are no mixed switch terms since A is unchanged. Two copies give steps -26 and -218, whose gcd is2. Their four values do not constitute a consecutive interval; gcd is not coverage. A chain of k copies has h-prefix=13^j and A-prefix=64*13^(j-1)*(16^j-1)/15 for j>=1. In particular the coefficients have a growing common divisor at least2*13^(k-2) for k>=2, so repeated identical cores do not provide arbitrary-length unit intervals.

## Scope and receipt

Run `python3 -I -B verify.py` and `python3 -O -I -B verify.py`; compare outputs. `--fail` must terminate nonzero. The verifier is self-contained. Exploratory source searches motivated the explicit witnesses, but are not required to reproduce any proved statement. The positive certificates do not assert a census beyond their stated universes.

An early exploratory t=0 ranking briefly used an invalid endpoint identity and produced an odd gap105. This was rejected by tournament count parity, repaired to the true gap54, and is not used by any theorem or final verification. The t>=1 witness and proof were unaffected.

Remaining bridge: realize arbitrary target N in one of these anchored physical families, or embed a sufficiently wide interval of switch increments inside a controlled near-target baseline. Reversal cancellation and coprime deletion alone prove neither. Chan–Pak remains open.

## 7. An h37 cell removes the growing divisor at every length

The cut search subsequently supplied a second physical core C37=[32,17,3,7,13,30], with cut17 and cut33. Independent recomputation gives respectively

    (37,352,96,102,12,12,10),
    (37,352,96,96,32,10,10).

After an unselected sink the common h,q,A,E are37,704,192,0, and the switch changes G by−6,T by20,S by−2. Compose one protected h13 core followed by k−1 protected h37 cores. Let h_j,A_j be the common product values for the first j modules. Every switch is affine, because A stays fixed and E=0 eliminates S interactions. The first scalar coefficient is −2*37^(k−1); later coefficients, at positions i=2,...,k, are

    c_i=37^(k−i)(6 A_(i−1)−20 h_(i−1)).

For k=1 the sole coefficient is−2. For k=2 the last coefficient is124, not divisible by37. For k>=3, h_(k−1)=13*37^(k−2) is zero modulo37. Meanwhile A_1=64, A_2=64*704+13*192=47552≡7 (mod37). All subsequent updates A_(j+1)=704 A_j+192 h_j preserve the residue7 because704≡1 (mod37) and37 divides h_j for j>=2. Therefore c_k≡42≡5 (mod37). All coefficients are even. Since their gcd divides the first coefficient2*37^(k−1) but is not divisible by37, their gcd is exactly2 for every k.

This is a genuine all-length escape from a growing scalar divisor in a fixed physical affine cube. It does not show that the binary subset sums cover an interval. The first boundary contribution13*192 is essential in the modular calculation; discarding it gives an incorrect intermediate residue even though the conclusion happens to survive.

The supplementary `verify_chain.py` independently checks both explicit cells, exhausts all binary choices through k=12, verifies the exact affine formulas and gcd2, and checks both h13 seven-vertex closures by direct permutation enumeration independent of subset DP.

## 8. Exact exponential spacing despite primitive coefficient gcd

**PROVED.** For every k>=2, the mixed h13/h37 chain of Section7 has exactly2^k distinct scalar counts, and the minimum difference between two distinct counts is

    50*37^(k−2).

Thus its coefficient gcd2 does not even produce a unit odd-count pair once k>=2.

**Proof.** Put S=37^(k−2). The first two coefficients are c_1=−74S and c_2=124S. Their four subset sums are−74S,0,50S,124S, so their minimum spacing is50S.

For a prefix of length j let b_j=6A_j−20h_j. We have b_1=124, h_j>0, and the recurrence from Section7 gives

    b_(j+1)=704b_j+14492h_j.

Consequently b_j>0 and, for successive positive coefficients of the fixed k-chain,

    c_(i+1)/c_i = b_i/(37b_(i−1)) >704/37>19,  i>=2.

Let L_i=sum_(j<=i)|c_j|. Initially L_2=198S<2c_2. By induction, whenever i>=3,

    c_i−L_(i−1)>19c_(i−1)−2c_(i−1)
                 =17c_(i−1)>=2108S>50S,

and L_i< c_i+2c_(i−1)<2c_i. For any two binary choices whose largest differing index is i>=3, the magnitude of their scalar difference is at least c_i−L_(i−1)>50S. If the largest differing index is at most2, the minimum is exactly50S by the four explicit values. Hence no collision occurs and the minimum is attained by choices differing only in the first two positions. This proves both conclusions for all k.

The positive consequence remains exact primitive coefficient gcd2 for every length. The hostile consequence is equally exact exponential separation in this particular fixed family. Neither statement limits other module sequences or proves an obstruction for arbitrary tournaments. `verify_chain.py` exhausts this family's2^k counts for each k<=12, checks that all are distinct, sorts them, and verifies their minimum adjacent difference equals the formula.
