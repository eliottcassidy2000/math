# Hollow-frame second-moment supplement

Read `HOLLOW_FRAME_SECOND_MOMENT.md` for the reviewed all-k formula and proof.
The main maximum-census packet remains unchanged.

With R(z)=(1+z^2/4)/(1-z^2/4), the formal coefficient identity is

    E[h^2]/k! = [z^k] R(z) sum_(m>=0) m! (z R(z))^m.

Subtract (k!)^2 for the variance. Exact values through k=7 are recorded in
`second_moment_checks.json`; only k<=5 is supported by a complete frame
census. The k=6,7 entries are formula and independent overlap-sum values.

From the published `second_moment` directory:

    python3 verify_second_moment.py --census-dir ..

Without the finite census files:

    python3 verify_second_moment.py

Only the Python standard library is needed. Checks remain active under
`python3 -O`. See `SECOND_MOMENT_REVIEW.md` for review and scope boundaries.
