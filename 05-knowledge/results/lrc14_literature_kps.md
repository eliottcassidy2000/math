# LRC(14) literature dossier (kind-pasteur, 2026-06-16)

## INDEX CONVENTION (load-bearing)
Three conventions in the literature; OUR problem = 13 distinct speeds, gap 1/14.
- Tao / Wills "n runners, gap 1/(n+1)" with delta_n: OUR case = **n=13**, delta_13 = 1/14.
- Goddyn-Wong "n-1 speeds (runners), gap 1/n", LR(n): OUR case = **n=14**, LR(14)=1/14.
- Recent computer-proof papers "LRC(k): k speeds, gap 1/(k+1)": OUR case = **k=13**.
All agree: 13 speeds, gap 1/14.

## PROVEN FRONTIER (as of 2026)
- Survey (Perarnau-Serra, arXiv:2409.20160, Sep 2024): proved for n<=6 (=7 runners) via Barajas-Serra 2008; n+1=7 prime exploited.
- Rosenfeld arXiv:2509.14111: 8 runners (k=7). Trakulthongchai arXiv:2511.22427: 9,10 runners (k=8,9).
- Sungkawichai-Trakulthongchai arXiv:2604.23906 (Apr 2026): LRC(k) for k in {10,11,12}.
- => Largest proven: k=12 (gap 1/13). **OUR k=13 (gap 1/14) is the FIRST OPEN CASE.**

## TIGHT INSTANCES (Goddyn-Wong, Integers 6 (2006) A38; survey Sec 4)
- Tight = LR(v)=1/n exactly (equality in conjecture).
- AP [n-1]={1,...,n-1} always tight. Three sporadics (Wills/Flor): {1,3,4,7}(n=5), {1,3,4,5,9}(n=6), {1,2,3,4,5,7,12}(n=8). For n<=6 these are ALL tight instances up to dilation (Cusick/Chen/Bohman-Holzman-Kleitman).
- Goddyn-Wong Thm 2.3: [n-1] with r->mr is tight IFF n=2, or (n,r,m)=(3,1,4), or gcd(r,b)>1 for all b in {n-r,...,m(n-r)-1}.
  (Survey form, kappa=1/(n+1): replace r by mr in [n] tight iff r shares factor with all of [(n+1)-r, m(n+1-r)-1].)
- Family: {1,...,n-2, n, 2(n-1)} tight for every n=6t+1 (survey eq 12). For OUR case n=13=6*2+1: gives {1,..,11,13,24} (=Goddyn-Wong T5, n=14 their convention). MATCHES our sporadic.
- Goddyn-Wong Thm 2.4 (FINITENESS, single accel): "For any fixed speed r there are only finitely many pairs (n,r') s.t. [n-1]_{r->r'} is tight." Tight forces n<12r and r'<12r.
- Corollary 3.3 + remark: infinitely many tight vectors EXIST but n GROWS with #accelerated runners; two accel runners (subseq (2,3)) needs n >= 1.2e14. Infinitude is across n->infty, NOT at fixed n.
- Pomerance bound (survey 8.2): exists c>0, if n < v_n < 2n - c log^2 n then V NOT tight. Bounds max speed of tight configs near the AP regime.
- Survey verdict: "complete characterization of tight instances is still widely open." Converse of multi-accel theorem FALSE in general.

## BOUNDED-SPEED REDUCTION (makes compactness rigorous)
- Tao arXiv:1701.02048 Thm 1.3 / Cor 1.4: LRC for n<=n0 EQUIVALENT to checking only tuples with speeds of size n^{O(n^2)}; decidable in finite time. Convention: Wills n runners.
- Sharper explicit (Malikiosis-Santos-Schymura, via 2604.23906 Lemma 2.6): if LRC(k-1) holds, any counterexample to LRC(k) has product u_1...u_k < B_k = ( C(k+1,2)^{k-1}/k )^k.
- Recent proofs reduce to gcd=1 coprime tuples + congruence (u_i ≡ i mod p) sieve.

## GAP LOWER BOUNDS (NOT measure-of-lonely-set bounds)
- Trivial kappa(n) >= 1/(2n) via union bound = EXACTLY our L(S)>0 union-bound (Pr(union of danger sets)<1).
- Chen: kappa(n) >= 1/(2n-1-1/2^{n-3}). Chen-Cusick: kappa(n)>=a/p for primes p>=max{2a(n-1)-1,2(a-1)n+2}.
- Tao Thm 1.2: delta_n >= 1/(2n) + c log n/(n^2 (log log n)^2).
- Bedert arXiv:2511.16636 Thm 1.3 (Riesz products): ML(V) >= 1/(2n) + 1/n^{5/3+o(1)}, all n. Riesz-product measure mu with controlled Fourier coeffs; Lemma 2.1 integral>=1.
- NONE of these lower-bound the MEASURE of the lonely set (our L). They bound the max gap (sup_t min_v ||tv||). No literature result directly bounds inf L below.
