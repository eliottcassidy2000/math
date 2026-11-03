# Order-eight projected certificate

See PROOF.md for the exact theorem, dependence on the order-seven lifting proof, scope, and improved all-order bound. This is a quotient certificate, not a full 13028-by-13028 rank computation: the previously proved order-seven lifting supplies the other 10795 dimensions. The two ordered context parts concatenate to a fixed 2233-row certificate. To replay:

    cat contexts8.part00 contexts8.part01 > contexts8.txt
    sha256sum contexts8.txt
    g++ -O3 -std=c++17 -o verify8 verify8.cpp
    ./verify8 8 contexts8.txt reps8.txt

Expected SHA-256 is 94df48ba7e434651698928b156d81ef253c39cf76f9df8c20a4c87b6adcad57b. The verifier independently generates the full-support signature representatives by a path exact-cover method, checks them against `reps8.txt`, validates each exterior context, and computes exact prime-field ranks. Expected ranks and determinants are 2233/70 modulo 101 and 2233/94 modulo 103. The duplicate-row control has rank 2232. Optional deletion control: `./verify8 8 contexts8.txt reps8.txt --delete-last` reports rank 2232 over each prime.

Each line of the concatenated certificate is `N k u1 v1 ... uk vk`; old vertices are 0 through 7, and every listed edge has at least one new endpoint. Contexts have total vertex counts from 8 to 28 and 0 to 19 exterior edges. There are 2233 unique contexts. The earlier order-seven packet contains additional direct Hamilton, exact fixed-cut Walsh-extraction, and hostile-parser controls.

`generate_quotient8.cpp` records the discovery procedure from the first 2000 rows of the saved bank; verification does not depend on running it. To replay the continuation, in a temporary directory create `start.txt` from the first 2000 lines of `contexts8.txt`, compile the generator with `g++ -O3 -std=c++17`, and execute it. It writes `reps8.generated.txt` and `contexts8.generated.txt`. The random shuffle implementation may depend on the standard library, so byte identity of generated output is not a cross-platform correctness criterion.
