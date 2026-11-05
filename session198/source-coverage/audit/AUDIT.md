# Independent audit: five-form C3 inverse and sharp cutoff

Verdict: PASS on final frozen bytes.

PROOF.md SHA256: cf7294fd2759b1708d1850377d6c1b46019cea637b73f1585c9620b6697edc4e
coverage.py SHA256: fb1144fbea71e2f19c9a2583914a33bea6475c408ab1890bd444eee52b76b562
verify.py SHA256: 7c7ad84e6d9ffc649a6e68a71908a3d26b7f0852a8d3ce76b3c8478533b7f5a9

All new arguments were independently reviewed: deletion gives E<=d; the two
cover injections give U>=D+d-E; the invertible y coordinates are nonnegative;
all product bounds are entrywise weighted matrix bounds; every singleton
boundary is handled outside that cone; and the scalar transpose identities
preserve the admitted five-form domains. The uniform U<=floor(2H/15) cutoff
is correct. The p_z family with z=2^k+1, k>=1, proves representation-wise
sharpness along infinitely many physical examples. No stronger assertion about
minimum source size over all presentations is made.

Separately recomputed all six distinct matrices directly from inherited scalar
rows, without importing coverage.py. Their exact weighted minima are7,7,14,
7,6,29/4 for A,B,C,T,J0,J1. Results are in
INDEPENDENT_MATRIX_AUDIT.json. Confirmed the sharp J0 formula is
12z*2^n+3z and the k>=1 boundary restriction is present.

Replayed the full frozen verifier in an isolated copy, normal and optimized
modes; outputs match byte-for-byte. All release manifest hashes still match
after replay. The old and new C3 supports agree through5,001; the new search
uses6,339 sources and94,591 candidate evaluations rather than24,975 sources
and758,909 evaluations. Every126,946 source-bank recipe was replayed, and18
selected actual graphs were counted independently, up to order14.

Finite support remains finite: C1/C2 miss7,21,67 through5,001, and adding C3
leaves7,21. No argument here proves universal or cofinite odd coverage.
