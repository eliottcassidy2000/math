# Reviewed weak start-one source oracle

The [proof](WEAK_START_ONE_ORACLE.md) gives exact all-order membership and every possible root-deletion count for a weak start-one response triple. It uniquely normalizes the transitive sink tail, then existentially recovers the erased endpoint in the strong core. The [independent review](audit/REPORT.md) passes with the arithmetic-versus-bit complexity clarification incorporated.

The literal weak triple(9,0,18) has deletion counts1 and3. Thus the oracle returns a finite set, with one physical graph recipe for each deletion count, rather than falsely treating that response as a function of the triple. It does not establish cofinite scalar coverage.

Run the source normally or with -O:

    python3 verify_weak_start_one.py
    python3 -O verify_weak_start_one.py

The verifier rewrites its certificate identically; use a copy to preserve a pristine checkout. The previously reviewed strong inverse and its proof dependencies are included as pinned inputs. Its earlier independent audit is published at the sibling ../backward-profile-inverse/ checkpoint.

From audit/, run the independent forward-quadruple and literal-permutation controls:

    python3 independent.py
    python3 -O independent.py

The exact copied source directory required by that checker is included. The independent tests cover82,416 weak numerical candidates, all1,695 start-one roots through order5, and exact tail exponents up to200. All code uses standard libraries. FILES.sha256 binds every released byte. No navigation or index is changed.
