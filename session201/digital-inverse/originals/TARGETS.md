# Targeted arithmetic hostiles, not a new census

Status: author-verified exact certificates; pending independent review.

The separate divisor-gate theorem shows that a prime target n cannot use
JJZJ if3n-2 is also prime, or JJZZJ if5n-4 is also prime. We therefore chose
three individual prime triples, rather than extending the earlier interval
census. They exclude those two templates at every source height.

The digital search nevertheless found actual witnesses in the remaining menu:

    100000837       = Psi_JOZJJ(3,1890994),
    10000001533     = Psi_empty(2,7840792990),
    1000000001911   = Psi_JOZJ(8,22798770134).

The last witness was normalized from the full-menu word JOZOJ with left
index4 by the already reviewed boundary absorption. All three displayed
centers belong to the actual twelve-word menu. Their compiled graphs have
29,35,43 vertices respectively and exactly one Hamilton cycle by the
reviewed physical compiler. The positive path counts follow from full
source-profile replay, not a claimed literal enumeration of these large
graphs. `TARGETED-WITNESSES.json` records the complete graph words/profiles.

`PRIME-CERTIFICATES.json` supplies26 recursive Lucas certificates, including
all nine primes n,3n-2,5n-4. Verification requires only exact modular powers,
gcds and complete factorizations of p-1 into previously certified primes.
For completeness, the criterion is elementary: if a^(p-1)=1 modulo p and
gcd(a^((p-1)/q)-1,p)=1 for every prime q dividing p-1, then for any prime
r dividing p, the order of a modulo r contains every prime-power factor of
p-1. Thus p-1 divides r-1, forcing r>=p and hence p prime. Two is the base
case. No probabilistic primality test is trusted by the checker.

Run `python check_targets.py` and `python -O check_targets.py`. The SymPy
factorization used to find certificates is not needed for their verification.

The long first-witness search receipts retain their author node statistics.
Checking a positive index certificate does not require replaying that search.
The exhaustive result at2134117 is different: its complete length-five
inverse has no witness, agrees with the older finite proof, and uses97131
search nodes over25400 fixed rows in the author implementation.

None of these prime triples is a counterexample to occupancy. Neither the
positive probes nor the algorithm proves the uniform twelve-form lemma.
