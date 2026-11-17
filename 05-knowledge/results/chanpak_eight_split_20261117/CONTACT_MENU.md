# Two decompositions of 8n inside the same recontacting family

Use the second scaffold in [LITERAL_LIFT](LITERAL_LIFT.md), with any tournament B in its variable slot and n=22P₂(B)+5H(B). Keep that same physical n-path tournament R and the same two newborn vertices q,d. Vary only the mask Q of old vertices q beats; d still beats only q. The following is the complete menu of integer scale factors that are independent of B when q loses to the whole slot:

| Q (mask) | H(R+q+d) | Strong? |
|---:|---:|:---|
| 0 | n | no |
| 32 | 3n | no |
| 4 | 5n | yes |
| 34 | 7n | yes |
| 38 | 9n | yes |

In particular **8n=n+7n=3n+5n** has two decompositions into actual tournament counts extending the same labelled source. No single tournament realizes the even number 8n. There is a one-contact continuation 3n→7n (q begins beating old vertex 1) and then 7n→9n (q begins beating old vertex 2); no old arc is changed. The last recontact stays strong throughout. At B=TT₃, n=181, so the three strong companions have counts 905,1267,1629. This is an arbitrary-B construction, not a proof that one of its n values equals every admissible prime.

For verification, each of the 128 contact masks not touching slot 7 keeps the old unique entrance into B; thus all words over exterior vertices and B are classified by whether they visit B in one or two runs. The read of any mask is a linear function of H(B) and P₂(B), as in [LITERAL_LIFT](LITERAL_LIFT.md). The actual substitutions TT₁ and TT₂ distinguish these two coefficients. The exact replay tests all 128 masks for a common integer multiplier on these two evaluations, and checks the strength flags on the contracted graph. For comparison, the first scaffold's complete analogous mask bank has only the ×1 (mask 0) and ×7 (mask 3) invariant multipliers. No claim is made for tags in the variable slot.
