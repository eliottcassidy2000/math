# Physical source fibers and canonical target descent

The main result is the canonical physical descent in PROOF.md, Section 7.
Every nonbalanced represented C1 target has an actual balanced parent below
half its size, together with a required odd slope. The pair of parent count
and slope must come from one physical source. The target 443 proves that the
one-step version of this descent is insufficient.

The earlier sections derive a physical twin operation, its exact dyadic
continuation defect, a boundary restriction on hypothetical singleton source
fibers, and a shifted Mersenne-product subfamily.

Reproduce the fast complete census and graph controls:

    python3 verify.py --cap 10001
    python3 -O verify.py --cap 10001

Reproduce the full reported census, serially:

    python3 verify.py --cap 50001
    python3 -O verify.py --cap 50001

No third-party Python packages are needed. The full run keeps 3,816,297
source profiles with recipes in memory. Normal and optimized JSON outputs
must agree. The enumeration has an explicit U cutoff and is complete for
the stated F cutoff by F>=U+1.

The singleton-fiber statement beyond this finite cutoff is a conjecture.
The all-odd {7,21} coverage statement remains open. Profile multiplicity is
not word multiplicity or graph-isomorphism multiplicity.
