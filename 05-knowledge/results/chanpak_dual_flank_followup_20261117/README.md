# Dual marked flanks: balanced exact sevenfold pairs, and their counting limit

Follow-up to [the single-marker bank](../chanpak_bank_fusion_20261117/readers/PROOF.md) and the [fixed arbitrary-module sevenfold reader](../chanpak_changed_contacts_20261117/PLUS_RAY.md). The two-flank witnesses and the small pump came from the background pass; this packet replays their signatures, adds an exactly balanced prime-to-sevenfold pair, and proves an obstruction to turning a fixed finite marker bank into a universal selector.

## Construction and exact signatures

Take the binary transitive word of u followed by a marked block E, and the word of w followed by a marked block J. They are the two flanks X and Y. Between X and the reversed, mark-complemented Y put the center with outgoing bitmasks (54,32,2,119,6,20,55), mark mask 26. Add an apex beating just the marked vertices. All joins between individual blocks before the apex are oriented forward. The second construction exchanges X,Y. These are strong for positive u,w: the apex beats the first vertex of the first transitive word and is beaten by the last vertex of the reflected other word, and each of these spans the intervening blocks by the forward joins. Inside either outer word the same reachability follows from its transitive order.

Use the terminal profile (h,p,Q,L,e), where h=H, p=P₂, Q=P₃, L is the number of nonempty paths starting marked with complement weighted by P₂, and e counts full paths ending marked. As derived in the single-marker packet, the difference (swapped minus original) in H is 2(h_J e_E−h_E e_J), and P₂ is unchanged exactly when

    sum_j (A w_j − B u_j)(3/2)^j = C,
    A=p_E Q_J, B=p_J Q_E, C=2(p_J L_E−p_E L_J).

Bits here are read least-significant first. The carry starts at C and advances by c' = 2(c−A w_j+B u_j)/3. Complete binary solutions have integer carries (the yet-unprocessed sum has denominator a power of two, while each proposed division has denominator dividing three).

* **Gap six at equal length.** Use E=(22,12,56,49,34,3), mask 57, and J=(10,20,25,18,1), mask 21. Their profiles are (43,424,2259,1230,30) and (13,104,471,222,9). The integer reduction is (A,B,C)=(8321,9789,2816). For (u,w)=(6,7), the carry visits 2816, −3670, −1468, 0. The hosts read (570968953,3203502193724) and (570968959,3203502193724) in (H,P₂). Add two ordinal singleton pads and apply the published 128P₂+H reader and its sevenfold companion. Both lower counts, 1640193694155641 and 1640193694155647, are prime by exhaustive trial division. Their exact sevenfold values are 11481355859089487 and 11481355859089529.
* **A pump with equal lengths.** For E=TT₃, mask 6, and J=(10,12,9,0), mask 1, the profiles are (1,8,27,10,1) and (3,24,99,21,0). The reduced system is (11,9,2). The least-significant column pairs (1,0),(0,1) loop 2→−6→2 and (1,1) exits 2→0. Hence u=(5·4^k−2)/3, w=(4·4^k−1)/3, k≥0, yield the exact gap six and equal P₂ at all heights. For k=0 the counts 93859 and 93865, with common P₂ 58530716, also follow from a literal order-17 subset DP, independent of the signed transfer.
* **Exact balance and prime quotients on both sides.** On both sides use (22,60,56,49,32,1), with marks 24 and 36. Both profiles start (19,256,1611,404), but their full end counts are 10 and 7. Thus every positive common address u=w yields H and H+114 at equal P₂. At u=48 both raw counts 39155020907 and 39155021021 are prime. For fusion with the sevenfold reader, instead choose u=w=39: H=25924835857, P₂=2843189815014472. One ordinal singleton doubles P₂; the physical paired reader reads the two primes **727856618568540689** and **727856618568540803** and their exact sevenfolds **5094996329979784823** and **5094996329979785621**. The replay certifies both primes by Lucas's criterion with complete factorizations of N−1.

## Why fixed marker banks cannot select arbitrary targets

For fixed E,J and positive address integers of bit lengths m,n, the unsigned P₂ seam equality is

    (Q_E D_m(u)/2^m + L_E)/p_E = (Q_J D_n(w)/2^n + L_J)/p_J,
    D_m(u) = sum_i b_i 2^i 3^(m−1−i),

where these bits are read most-significant first. D_m is odd, as are Q_E,Q_J by the slot packing identity at q=3. Since m,n are positive and the L terms are integral, the valuations of the two normalized reads force n−m=v₂(p_E)−v₂(p_J). For each u there is **at most one w**. Indeed the equation fixes D_n(w), and its digits are decoded from its parity, then by subtracting the leading contribution and dividing by two. In particular, even if both words are nominally free, exact two-cover compatibility has already spent one full address.

For the raw count, take just those odd cycles through the apex whose other vertices lie in the two transitive words, with at least one vertex in each. One such cycle is determined by its ascending subpath from a marked vertex in the left word and the corresponding subpath ending unmarked in the reflected right word. There are u and w choices of the two subpaths. Signed by the size of each subpath the sums have magnitude at most one: after fixing the first vertex, choices later in the transitive order cancel in pairs unless that vertex is last. Hence imposing even total subpath size leaves at least (uw−1)/2 choices. The leftover tournament has a Hamiltonian path, so the odd-cycle packing formula at q=1 gives H≥uw−1. As the bit lengths differ by a fixed amount, each fixed pair of markers supplies only O(sqrt(X)) matched raw counts up to X, even allowing all positive addresses and both reflections.

For the 128P₂+H reader the same O(sqrt(X)) bound holds with **arbitrary ordinal padding**. With a pad of j vertices, P₂ of the input is at least 2^(j+m+n+1) p_E p_J p_center: count the packings entirely inside the pre-apex blocks, with the apex a singleton. Thus n−m fixed gives O(sqrt(X)·2^(−j/2)) addresses for fixed j. Sum the geometric series. The companion sevenfold reader only makes the bound stronger. Any fixed finite bank of terminal pairs therefore has zero density in both the raw spectrum and this paired-reader spectrum. This does not bound a growing bank, or a construction which relinquishes exact P₂ equality and compensates in the reader.

## Next target

Use the forced length difference as a design equation rather than a restriction to work around afterwards. To get two independently selectable addresses, either let the terminal profile vary with the target (so that at least one of the carry coefficients actually grows), or solve the compensated relation α ΔP₂+β ΔH=small instead of ΔP₂=0. The balanced example proves that the sevenfold scaling and actual prime quotients survive the first change of port; the counting bound says why a fixed collection of such port changes cannot finish Chan–Pak.

Run `python replay.py`. It checks local signatures, both readers, an independent literal host, all cited prime certificates, and the initial levels of the exact pump.
