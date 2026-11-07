# Connect tournament recursions

**Status: PROVED elementary statements, independently audited; exact finite verification is recorded separately.**

The two inherited recursions act on vertex count, not Hamilton-path count:

- Source/sink padding: n to n+2, preserving H.
- Skew-Sylvester doubling: n to 2n, retaining a reflected copy of every Hamilton path but generally creating many others.

The [proof](PROOF.md) connects them to Session210's gated cycle carrier and coherent permutohedron:

1. A single apex-arc flip turns padding into the canonical strong gate and creates exactly P_T(2) new permutohedron sinks.
2. At core order n, the gate's first two H values are 1+2^n and 3+3*2^(n-1). Equality means a transitive core or exactly one C3 component, respectively. A triangle-resolved lower bound gives the exact small gated H spectrum {3,5,9,15,17,27} through 31.
3. The doubled reflected word has energy binom(n,2)+2E and forms an affine permutohedron section; an original adjacent swap lifts through a square face.
4. The two one-switch sheet sectors have exact sizes (H^2+sum e_v^2)/2 and (H^2-sum s_v^2)/2.
5. Even when the arithmetic vertex counts agree, the operations fail to commute, witnessed by sink or gate degrees.

These are structural statements and scoped gate-spectrum results. They do not prove unrestricted spectrum completeness beyond the known omissions 7 and 21. The familiar scalar identities, operations, reflected path injection, and 333/189 hostile pair are inherited, with sources pinned in [SOURCES.json](SOURCES.json).

## Reproduce

Run from the repository root:

    python3 session211/recursion-bridge/verify.py
    python3 -O session211/recursion-bridge/verify.py

Both outputs match [RESULT.json](RESULT.json). The exhaustive universe comprises all 1,099 labelled tournaments through order five and all 124,469 permutation words of those inputs. An independent product-of-path-lists test checks both sheet sectors through order four. Separate independent audit results are in `audit/`.
