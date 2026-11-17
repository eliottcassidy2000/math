# A congruence and ratio-band obstruction for pairing by a normalized receiver change

The [8-to-1 receiver catalog](../chanpak_prime_dichotomy_20261117/RECEIVER_GATE.md)
leaves just two prime-capable normalized reads from a common source B. In one orientation they are

    p=1482P₃−128P₂−31h,
    q=28626P₃−5264P₂−403h.

The same receiver on both sides just gives the diagonal; swapping receiver order exchanges p,q. The mixed pair gives exact identities

    q−13p=720(13P₃−5P₂)>0,
    367p−19q=120(442P₂−31h)>0.

This already requires 13<q/p<367/19. For primes p,q>3 one obtains the stricter **q≡3613p (mod 8640)**. Indeed p≡P₂−h modulo 3, so |B| is odd and 3∤h. A singleton cannot give prime p by this pair of receivers (q=74947=149·503). For order ≥2 we have P₂≡0 mod4; the odd-pack formula gives P₃≡(−1)^|B|h mod4 since slot variable 3≡−1. Thus t=P₃≡p mod4 and t≡0 mod3. Also P₂≡−p mod3. Set K=(q−13p)/720; then K≡p modulo4 and K≡−p modulo3, i.e. K≡5p modulo12. Substitute back for the claimed congruence.

For a supplied pair passing these phases, the entire exact integer inverse is (t=P₃) with

    P₂=(13t−K)/5,     h=(5746t+128K−5p)/155.

The candidate is a finite arithmetic progression after integrality, divisibilities and positivity, not a guarantee of a physical profile. For instance 13t≡K mod20, t≡14(128K−5p) modulo155, and t≡0 modulo3.

The width of the possible ratio band tends to zero along large sources. In fact P_{k+1}≥(n+k)P_k/k for n=|B|: split a chosen slot at any of its cuts and put the suffix in a new slot. There are n+k choices for any k-cover and at most k preimages per output, by choosing the slot to restore. Hence P₂/P₃≤2/(n+2), and the second determinant identity makes 367/19−q/p=O(1/n) because p≥1323P₃. This rules out using only the normalized receiver switch to accommodate two arbitrary independent prime inputs. Deliberately moving the defect, in contrast, gives the new [signed pairing](SIGNED_PAIR.md).
