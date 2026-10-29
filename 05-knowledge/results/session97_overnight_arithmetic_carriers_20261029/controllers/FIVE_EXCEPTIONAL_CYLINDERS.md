# A complete five-cylinder classification of expanding first returns

**Status: PROVED.** Exact scope: the cheaper ternary first-stop map F of `INTEGER_TERMINATION.md`, on every positive odd integer. This is not a classification of divergent Collatz integers.

## Theorem

The controller F is everywhere defined on the positive odd integers, and

    F(n)>n

if and only if n belongs to one of the following five pairwise disjoint source cylinders:

| Exact valuation word | Source guard | First-return affine map |
|---|---|---|
| (1) | n≡3 (mod4) | (3n+1)/2 |
| (3,1,1,1) | n≡125 (mod128) | (81n+179)/64 |
| (3,3,1,1,1,1,1) | n≡1357 (mod4096) | (2187n+16177)/2048 |
| (3,1,3,1,1,1,1) | n≡2269 (mod4096) | (2187n+12289)/2048 |
| (3,1,1,3,1,1,1) | n≡2877 (mod4096) | (2187n+9697)/2048 |

F(n)=n holds exactly at n=1. Every other positive odd integer satisfies F(n)<n. Thus a **12-bit arithmetic classifier** decides the first-return comparison globally, even though computing the entire return still retains its exact source, ports, carry, and integer clocks.

## Proof: successively eliminate expensive letters

Use P=(n+1)(3/2)^Q and its branch contraction from the preceding note. At termination Q changes from1 to0, so one must account for a factor3/2 when comparing n_final+1 with n_initial+1.

1. If any leaf letter has a≥2, its P ratio is at most2/3. That pays the factor3/2. If there is any branch, its ratio is strictly less than1, so the endpoint is strictly smaller. If there is no branch, the word has length1: a≥2 gives strict descent except n=1,a=2, which is equality.
2. If any branch has a≥6, its input is at least13 and

       P'/P ≤ (9/4)(3n+65)/(64(n+1)) ≤117/448.

   Since (3/2)(117/448)<1, this alone pays the terminal factor and gives strict descent.

Therefore every still-unresolved nontrivial word uses only a=3 at branches and a=1 at leaves.

## A finite arithmetic exclusion improves every remaining branch

For a=3, positive inputs satisfy n≡13 mod16. The seven such inputs below125 cannot occur at a branch of a completed word over{1,3}. The following short native prefixes show a forbidden letter before the subtree can close:

| Input | Prefix ending at first forbidden letter |
|---:|---|
|13|(3,4)|
|29|(3,1,2)|
|45|(3,2)|
|61|(3,1,1,5)|
|77|(3,3,1,2)|
|93|(3,1,5)|
|109|(3,2)|

In each row, starting Q=1, the proper prefix queue stays positive. If the occurrence has a larger queue, exactly the same prefix stays active. Thus every branching input of a completed{1,3}-only word is at least125. This is a finite integer certificate, not an empirical inference from testing a range.

At such an input,

    P'/P=(9/4)(3n+9)/(8(n+1)) ≤6/7,

with equality at125. If there are m≥3 branches,

    (n_final+1)/(n_initial+1) ≤(3/2)(6/7)^m
                                 ≤324/343 <1.

So every unresolved nontrivial expansion has m=1 or2. There is one ordered ternary tree with1 branch and three with2 branches. Their chronological words are exactly the four nontrivial rows of the theorem. Each has positive carry and slope strictly greater than1, so every positive source realizing its exact word expands. The source guards come from the exact affine cylinder formula, proving both directions and excluding additional cases. The zero-branch expansion is exactly the leaf a=1, n≡3 mod4.

## Exact density and source-preserving retry interface

Normalized among odd integers, the expanding set has density

    1/2 +1/64 +3/2048 =1059/2048.

The density of strictly descending inputs is989/2048; the singleton equality has density0. These are first-controller statistics, not orbit-independence statements. In particular an iterated source can revisit the expanding cylinders with a highly nonuniform distribution.

The five rows are an exact finite dispatcher. On any expanding row, keep its original congruence, affine carry and clocks, then pull the next row's guard back through the displayed map. On a descending row, the controller certifies descent relative to its current input. A stronger dispatcher seeking descent relative to an earlier source must retain that earlier source and compose its full affine carry. Neither a resampled Haar law nor a remembered12-bit output alone proves that stronger statement.

The finite comparison classifier is substantially smaller than the faithful carrier: it decides only F(n)−n, not arbitrary future guarded observations. This is an explicit instance of a short certificate coexisting with an information-rich exact representation.

## Unbounded clock hostile and next target

For n=2^t−1, t≥2, the next t−1 native steps all have a=1 and strictly increase. Hence no fixed number of repeated cheap first returns can guarantee a descending step uniformly. An exact compressed run keeps t=v_2(n+1), u=(n+1)/2^t, giving

    n_after(t−1)=2·3^{t−1}u−1.

The next valuation is1+c where c=v_2(3^t u−1)≥1, and the next odd endpoint is(3^t u−1)/2^c. This transition retains an unbounded integer clock and the source guard on u; a bounded-state substitution that forgets t is not claimed faithful.

Useful next research: compose these five guarded rays into source-marked return languages, paying the exact carry and stopping at a drop below the **original** source. The present theorem gives the complete failure alphabet for this controller, but does not prove global termination of that retry process.

## Reproduction

`python3 -B check_exceptional_cylinders.py`, repeated with `-O`, produces `exception_checks.json`. It independently checks all500001 odd sources through1000001 against the five-cylinder iff, verifies the equality locus, and tests initial expanding runs for2≤t≤300. The proof uses only the seven finite obstruction rows and the exact inequalities above; it does not extrapolate the million-scale check.
