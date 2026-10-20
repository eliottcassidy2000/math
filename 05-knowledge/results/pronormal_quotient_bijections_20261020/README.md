# Pronormal bases and prime-local normalizer quotients

Status: proof and exact controls staged for independent review.

Read `PRONORMAL_QUOTIENT_THEOREM.md` for the theorem, scope, and failure boundary. The new mixed-prime survivor is: if K is odd and pronormal in S_n and q is odd, all symmetry predicates for K<=H<=N(K) with H/K a q-group are preserved simultaneously by one bijection on the K-admitting tournament/sign-even-graph types.

The fixed-cube result without pronormality is relative to N(K), and is not silently substituted for a global type count.

Reproduce:

```sh
python3 n6_c3_certificate.py > n6_c3_certificate.json
python3 -O n6_c3_certificate.py > n6_c3_certificate.optimized.json
cmp n6_c3_certificate.json n6_c3_certificate.optimized.json
python3 c7_normalizer_quotient.py > c7_normalizer_quotient.json
python3 -O c7_normalizer_quotient.py > c7_normalizer_quotient.optimized.json
cmp c7_normalizer_quotient.json c7_normalizer_quotient.optimized.json
```

Both scripts are standard-library-only. The first independently checks full S6 actions. The second checks the entire 24-by-24 Walsh matrix on the C7-fixed cube at n=10, including the whole quotient Sylow action and all its C3 predicates. Its global interpretation uses the proved ambient-Sylow descent.

Result highlights:

- C3 at n=6: single cycle 11=11, double cycle 3=3, intersection 1=1.
- C15 at n=8 and 9 has a suitable ambient Sylow C5 base; at n=10 that base fails pronormality and the known K_(5,5) obstruction appears.
- C7 at n=10: the F21, C7 x C3, diagonal F21, and F21 x C3 counts are 12,4,2,2 on each side, preserved by one quotient-prime matching.
- C5 wr C3 at n=15: a prime-power normal base yields the unique tournament C3[C5,C5,C5] paired in this sector with the empty graph.

No repo or Git mutation was performed by this worker. Publication is centralized.
