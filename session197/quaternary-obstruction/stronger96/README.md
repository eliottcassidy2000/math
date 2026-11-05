# Quaternary unit-clone amplitude is at least96

This reviewed Session197 strengthening proves: if an actual rooted finite tournament has an eventual transitive-clone identity H_m=A*4^m+epsilon, with A>0 and epsilon either-1 or+1, then A>=96. There is no symmetry assumption. This does not assert existence at amplitude96, nor exclude other graph-changing constructions.

- [Proof and exact order cutoff](PROOF.md)
- [Independent review](audit/REPORT.md)
- [Complete8740-class strong bank](strong-bank.jsonl)
- [Standalone reproduction driver](verify.py)
- [Independent uncompressed gap counter](audit/check.cpp)

The general unit-clock reduction removes arbitrary weak singleton padding and forces6 to divide A. Below96 the strong core has H<=361;361^3<5^11 forces order at most11. The included generator enumerates every mixed one-vertex strong extension and uses exact degree-preserving canonicalization. The whole bank is regenerated and byte-hash-checked. All12263 eligible roots fail the required second/third clone equations. A separate gap counter and a third independent uncompressed counter agree.

## Reproduce the complete result

Python3 and a C++17 compiler are the only requirements. No network or external data is needed.

    python3 verify.py
    python3 -O verify.py

Both outputs match receipt.json and receipt-optimized.json. The driver compiles its programs in a temporary directory, regenerates the complete bank from the directed triangle, then executes both root counters. The bank is not merely accepted as an unverified external census.

For the independent counter, convert the included bank to its CSV input:

    python3 -c 'import json; from pathlib import Path; a=[json.loads(s) for s in Path("strong-bank.jsonl").read_text().splitlines()]; Path("bank.csv").write_text("".join(",".join(map(str,[x["n"],x["H"]]+x["rows"]))+"\n" for x in a))'
    g++ -O3 -std=c++17 audit/check.cpp -o independent
    ./independent bank.csv

The output matches audit/normal.txt and audit/debug.txt. An assertion-enabled unoptimized build gives the same output. The independent review also compared the cap217 sub-bank to the previously reviewed smaller bank; that extra historical consistency check is not required to run the main counter or to regenerate this complete bank.

## Reviewed dependency and source names

[dependencies/UNIT_CLOCK_REDUCTION.md](dependencies/UNIT_CLOCK_REDUCTION.md) is an unchanged copy of the earlier reviewed proof containing the general eventual-to-exact, factorial-divisibility and SCC reduction. Its earlier amplitude60 finite census is superseded by the new amplitude96 bound here, not rewritten. Its [published checkpoint](https://github.com/eliottcassidy2000/math/blob/f633072db124b398e2563b9e3ba7b12e4fc2c7ec/session197/quaternary-obstruction/README.md) contains the independent review and original controls.

The frozen review retains its original session-local source names. Source files are at this release root; its independent checker is audit/check.cpp and receipts are audit/normal.txt and audit/debug.txt. All copied proofs, sources, reports and receipts retain their reviewed bytes. RELEASE_FILES.sha256 binds every other release file.
