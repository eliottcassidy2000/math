# Session208: Higher atoms and the global efficiency boundary

Explicit strong tournaments through order12 improve the global bound to lim m(15^k)/k<=8/3. An order9 witness refutes any global extension of the preceding order-eight dual prices. The earlier finite-dictionary optimality remains valid. An independently proved permutation-orbit bound supplies a universal sublinear lower bound; no exact global efficiency is claimed.

Run `python3 verify.py` with assertions enabled. This reconstructs nine witnesses, verifies strong connectivity and reverse-endpoint path counts, and checks all four construction identities. It rewrites WITNESSES.json and RESULT.json deterministically; stdout should match verify.out.

Optional discovery replay, using a C++17 compiler:

    g++ -O3 -std=c++17 search9_10.cpp -o /tmp/s208-small
    /tmp/s208-small 9 100000 > /tmp/search9.out
    /tmp/s208-small 10 200000 > /tmp/search10.out
    g++ -O3 -std=c++17 search11.cpp -o /tmp/s208-eleven
    /tmp/s208-eleven 11 100000 > /tmp/search11.out
    g++ -O3 -std=c++17 search.cpp -o /tmp/s208-twelve
    /tmp/s208-twelve 12 100000 > /tmp/search12.out

These searches are not exhaustive censuses. The frozen witnesses, not discovery probabilities, establish the results. See PROOF.md for precise scope, code convention and the inherited dominance-join mechanism.

Independent replay: `python3 independent_check.py` should match `independent.out`. `REVIEW.md` records the independent mathematical audit. `DEPENDENCIES.json` pins the prior published packet. Check frozen payloads with `sha256sum -c SHA256SUMS`.
