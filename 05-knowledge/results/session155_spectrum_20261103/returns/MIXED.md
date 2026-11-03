# Mixed-return obstruction and a descendant boundary

**Status.** A finite-state proof for the seven named macros below, of arbitrary word length. It applies to uniformly 3-divisible bulk quotients and the indicated seed residue class. It does not settle mixed returns from the full 65-symbol bank or all possible target states. A separate exact integer witness shows the original-seed obstruction must not be extended to all descendants.

## Setup

Use row states and processing order. Put

    t=(2,-1,0),
    beta(h,x,y)=(y,h+2x-y)=vB,
    B=[[0,1],[0,2],[1,-1]].

For a return macro `X`, write `P_X B=B Q_X`. The seven macros and their exact 2-by-2 bulk matrices are:

    macro  processing word             Q_X
    L2     r24 a a a                   [[240,339],[90,127]]
    L6     r24 q4                      [[720,1101],[270,413]]
    L45    r19 r32 r48 q4              [[1020432,1560769],[366879,561148]]
    L13    r19 r19 a a r18             [[211965,378780],[76220,136205]]
    b      b                            [[5,8],[2,3]]
    q2     q2                           [[11,18],[4,7]]
    Lp2    r24 a a r21                 [[4773,8628],[1789,3234]]

Each word has a literal tournament construction in the source atom bank and returns the error line `Zt`. The script `mixed_audit.py` contains the ten atomic 3-by-3 matrices, independently multiplies the words, derives these quotients, and checks each return-line identity. The three original target states are

    U=(31,14,25), U2=(95,41,75), U4=(97,40,75),
    beta(U)=(25,34), beta(U2)=beta(U4)=(75,102).

The first bulk reduces to `(1,0)` modulo 2 and to a nonzero multiple of `(0,1)` modulo 5, as do the other two.

## The mixed-word theorem

For every word `W` over those seven macros, let `Q_W` be the ordered product of their bulk matrices. Suppose all four entries of `Q_W` are divisible by 3. Then for **every** integer target state `v` such that

    beta(v) = (1,0) mod 2,
    beta(v) is a nonzero multiple of (0,1) mod 5,

there do not exist `w ∈ Z³, m ∈ Z` satisfying

    w P_W = 3v + m t.                                      (1)

In particular this excludes exact uniform-elevator resets to each of U, U2 and U4, even with arbitrary formal integer predecessors and arbitrary integral `t` carries.

**Finite closed certificate.** Starting with the 2-by-2 identity and repeatedly multiplying on the right by each `Q_X` modulo 30 yields exactly 6,400 matrices. The included verifier checks closure of that set under all seven generators. Thus it is the entire generated monoid, not a bounded-length sample. Exactly 184 members vanish modulo 3. Among those 184, the numbers for which the row space contains `(1,0)` modulo 2 and `(0,1)` modulo 5 are, respectively:

    mod2 fails, mod5 passes: 12
    mod2 passes, mod5 fails: 111
    mod2 fails, mod5 fails: 61
    both pass: 0.

**Proof from the certificate.** Write `Q_W=3R_W`. Applying beta to (1) and using `beta(t)=0` gives `beta(w) R_W=beta(v)`. Since multiplication by 3 is a unit modulo 2 and modulo 5, `R_W` and `Q_W` have the same row space over each of those fields. The seed residue conditions would require both row-space gates to pass. The closed finite certificate proves this never occurs for any `Q_W=0 mod 3`. No condition on `det R_W` is required.

Both prime gates are needed. For `W=Lp2 L2`, the normalized bulk matrix is `[[640680,904601],[240140,339063]]`, determinant `-1300`, which passes the mod-5 gate but fails mod 2. For `W=Lp2 Lp2`, it is `[[12739007,23028132],[4774841,8631416]]`, determinant `16900`, passing mod 2 but failing mod 5. Both determinants are prime to 3, so each is itself a uniform exact 3-adic level elevator.

## Mixed descendant escape at the formal level

Let `D=L2 b L2` and consider the mixed word `W=b D`. Direct 3-by-3 matrix multiplication gives

    P_D=[[336359,155087,268050],[672726,310170,536100],[560922,258618,447000]],
    P_W=[[2803655,1292663,2234250],[5607302,2585330,4468500],[4373654,2016542,3485400]],
    v=U P_D=(33868343,15615527,26989950),
    w=(69,0,-21),
    w P_W = 3v + 216t.                                     (2)

Thus the original-seed mod-30 obstruction cannot be generalized to all D descendants. Equation (2) is only a **formal** reset: `beta(w)=(-21,90)` and `y=-21`, incompatible with a physical tournament observer. In fact the normalized bulk matrix of `W` is `[[1906550,2692023],[744750,1051577]]`, with determinant 100, so the predecessor bulk solving (1) is unique over the rationals. Its first coordinate is `-21`; therefore **no physical predecessor exists for this particular W and v**, regardless of the integer carry. Other mixed words, targets and re-markings are not excluded.

## Reproduction

Run `python mixed_audit.py` and `python -O mixed_audit.py`; both use explicit checks and their output is identical. The finite monoid count, exhaustive closure and row-space tests, named macros and separate-gate examples, and all entries of (2) are checked. Input atoms and source construction are inherited from `05-knowledge/results/return_channels_20261103/INCOMING.md` and `.../channels.py` at commit `aa7514b988295fe0caac7fc06359a28b330dccf3`; the graph and source-label validity beyond matrix identities are not independently reproved here.
