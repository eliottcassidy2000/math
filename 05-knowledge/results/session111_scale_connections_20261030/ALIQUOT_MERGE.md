# The doubled289 shell rejoins the1001 aliquot branch

**PROVED elementary transport and complete finite fibre; CITED perfect-power classification.**

The earlier square-sum examples include289=8²+15². Gaussian doubling gives
578=7²+23²=2·17².
The aliquot sum then gives
578 →343 →57,
because s(2p²)=p²+3p+3 for odd prime p, and at p=17 this is7³. Also343+57=400=20².

This joins the older branch1001 →343 →57. The first arrow289→578 is doubling, not an aliquot edge; indeed s(289)=18. Keeping edge types prevents a false orbit claim.

## Complete fibre of343

For s(n)=a>1, n is composite. Its largest proper divisor is at least√n, so a≥1+√n and n≤(a−1)². Thus the343 fibre is fully contained in1≤n≤116964.

Two independent enumerations, a divisor sieve and direct divisor-pair sums, agree on exactly:
578,1001,1685,3641,9077,9641,16697,17141,19241,19637,
20777,24341,24617,25397,25877,27641,28757,28841,29177.

Only578 is even. The other exceptional factorization is1001=7·11·13. The remaining17 sources are products of two distinct primes whose sum is342. This is a finite Goldbach-pair classification of this fibre, not a proof of Goldbach in general.

An independent parity guard reduces the even search further: σ(n) is odd iff n is a square or twice a square, since every odd prime exponent must be even. For even n with odd aliquot sum, that criterion applies. No assumption that every such n is2p² is needed for completeness.

## Classical perfect-power restriction on the2p² subfamily

For an odd prime p,
s(2p²)=(p+1)²+(p+1)+1=Φ3(p+1).
The classical theorem of Nagell says x²+x+1=y^e, x,y positive and e≥2, has the unique solution (x,y,e)=(18,7,3). A primary research-paper statement and attribution are in [Terai–Takakuwa, Lemma2](https://www.m-hikari.com/imf/imf-2013/17-20-2013/teraiIMF17-20-2013.pdf), citing Nagell's1921 paper.

Consequently, **using that cited theorem**, 578=2·17² is the unique member of the odd-prime2p² subfamily whose aliquot sum is a nontrivial perfect power. This classification is not proved anew here and says nothing comparable about all composite sources. The full19-element fibre above is established independently of that citation.

[aliquot_fibre.js](aliquot_fibre.js) gives the complete finite certificate, factors every source, checks the17 semiprime cases and verifies both merging branches and the square-total edge.
