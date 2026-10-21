# Independent support arithmetic review

Date: 2026-10-21 UTC. Verdict: PASS as a necessary support restriction, conditional on explicitly inherited low-order data.

The even q cases give the result directly: q=2 permits only m=1,3 after the insertion exclusion, and order2 is not strong; q=4 permits orders4,12,20,28; q=6 permits6,18; q=8 permits8,24. Their union is {4,6,8,12,18,20,24,28}. Odd q cannot yield even n because m is odd. The q=1 lane leaves every odd order from11 through47 potentially unresolved.

Removing orders through9 using the explicitly inherited low-order census leaves exactly the nineteen odd orders11,13,...,47 and five even orders12,18,20,24,28, hence24 unresolved coefficients. The term 'exactly' applies to this parameter-filter output, not to existence of tournaments at any listed order. No individual unknown coefficient is computed or proved positive.

Copied support_check.py runs normally and under -O and uses always-active checks. This is arithmetic bookkeeping, not a fresh tournament census. In particular the review does not reauthenticate order9 or the inherited coefficients at5,6,7. The proposed decomposition z^5+z^6+z^7 plus unknown nonnegative coefficients at the24 listed orders is valid under the expressly inherited low-order input.
