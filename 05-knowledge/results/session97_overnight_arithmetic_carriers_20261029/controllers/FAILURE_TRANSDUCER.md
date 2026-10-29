# Exact failure charts versus their rank-one measure quotient

**Status: PROVED finite transducer and source-measure identities.** No every-integer exit theorem is inferred.

The five expanding words W from `FIVE_EXCEPTIONAL_CYLINDERS.md` form a prefix code of complete cheaper-ternary returns. For each ordered pair u,v in W, concatenate the native valuation words without changing their order. The resulting25 charts are generated in `transducer_checks.json`; each retains `(J,D,C,residue,modulus)` and its intermediate port.

Composition is

    J_uv=J_u+J_v,
    D_uv=D_u+D_v,
    C_uv=3^(J_v) C_u+2^(D_u) C_v,
    n0≡(2^(D_uv)−C_uv)3^(−J_uv) mod2^(D_uv+1).

This single guard is equivalent to the first guard and the second guard pulled back through the first affine map. Every positive representative has two strictly increasing first returns. The script checks this through direct native execution, independently of the guard formula, on five representatives of all25 charts. In arbitrary length, the same recurrence supplies5^k pairwise disjoint guarded histories. This is a finite symbolic alphabet, not a claim that the cumulative arithmetic sidecar has only finitely many states.

## Why the measure quotient is rank one

A complete word w maps its odd source cylinder bijectively onto the odd2-adics. If the cylinder is parametrized n=r+2^(D+1)k, its endpoint is

    y=y0+2·3^J k.

The coefficient3^J is odd. Consequently normalized Haar measure conditioned on that particular original cylinder pushes to uniform odd Haar measure at the endpoint. This proves the product law for future finite valuation words while retaining their original-source pullbacks; it does not authorize resampling a fixed integer.

The five failure weights are

    p=(1/2,1/64,1/2048,1/2048,1/2048),
    rho=sum p=1059/2048.

The measure-level failure transition matrix has all five rows equal to p, so it has rank1. Its k-step retained source mass is exactly rho^k. Thus the first syntactic exit from this five-word language is Haar-almost-sure, with mean number of attempts2048/989. These statements use the exact prefix code, not an independence assumption about a fixed integer orbit.

The accumulated clocks of the failed words before exit have means

    E[J_failure]=1173/989,
    E[D_failure]=1249/989,
    E[J_failure+D_failure]=2422/989.

Indeed the one-attempt unnormalized weighted sums are1173/2048 and1249/2048, and summing their geometric repetitions divides by1−rho. This does not include the terminal nonfailure word's clock.

## Loss ledger and arithmetic boundary

The rank-one quotient preserves source mass and these averaged clocks. It discards the actual source residue, intermediate integer, cumulative carry, and unbounded integer clock. Those discarded data determine whether a future local decrease crosses the **original** source. The25 exact charts retain them.

An infinite word over the failure alphabet has Haar mass0 as an event, but may still encode an exceptional individual source. Source2^t−1 realizes an arbitrarily long prefix of the first failure word. Hence a finite clock bound cannot be recovered from the geometric measure tail. The universally terminating auxiliary queue theorem applies to each first-return word individually; it does not prove that the failure-word retry language has no positive integer infinite path.

This makes the distinction precise: an exact small transition matrix can calculate aggregate laws while faithful composition still needs growing arithmetic memory. The proof is specific to dyadic affine word cylinders and is not an identification with a tournament centralizer or an assertion that a generic finite-state quotient preserves its ports.

Reproduce with `python3 -B check_failure_transducer.py` and `python3 -O -B`; outputs are byte-identical. All25 rows and125 native source witnesses are recorded/checked.
