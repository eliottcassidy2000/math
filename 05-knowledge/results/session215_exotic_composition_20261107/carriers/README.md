# Positive paired carriers from an endpoint reset

## Main result

`PROOF.md` gives two explicit physical carrier pairs at adjusted bases 4225
and 4357. Each has bounded starts with j in [0,8], and every common O/Z word
preserves equality of the two adjusted bases. In particular, for every k >= 1,

- X = 8715 * 2^k - 4295;
- X = 9585 * 2^k - 4977.

Each X supplies actual tournaments with Hamilton-path counts 3X+2 and 3X+4.
The positive mechanism is a final O operation resetting unequal endpoint
counts after new J-containing source constructions. The two initial bases
are outside the four inherited O/Z scalar orbits. Their descendant four-states
stay outside the old four-state orbits because d differs.

No all-odd coverage theorem or new scalar Hamilton-count value is claimed.
The two scalar rays are not asserted to avoid every old scalar base.

`NEW_CARRIERS.json` contains the four explicit rooted graphs, complete
nine-coordinate profiles, compiled graphs, and counts. The independent literal
checks in `verify_new_carriers.py` are the main reproducibility entry point.

## Supporting bounded searches

`new_bank.cpp` searches the supplied order-at-most-eight bank through three
O/Z/J moves and records 20 selected synchronized four-state keys. Completeness
of that bank is not needed for the two explicit constructions.

The old eight-source searches are separate exploratory evidence:
`search.py` records 122 selected asynchronous scalar collisions through depth 9;
`families.py` and `jrays.py` test selected exponential-ray signatures through
depth 10. These finite probes have the representative-selection limitations
spelled out in `PROOF.md`. They are not global impossibility or coverage claims.
`WITNESSES.json` independently verifies two selected finite collisions.

## Reproduce from this directory

Check the delivered files first:

```sh
sha256sum -c SHA256SUMS
```

Reproduce the positive constructions without modifying the saved receipts:

```sh
python3 -B verify_new_carriers.py > /tmp/carriers-normal.json
python3 -B -O verify_new_carriers.py > /tmp/carriers-optimized.json
cmp NEW_CARRIERS.json /tmp/carriers-normal.json
cmp NEW_CARRIERS_OPTIMIZED.json /tmp/carriers-optimized.json
```

Reproduce the discovery catalogue and selected asynchronous witnesses:

```sh
g++ -O2 -std=c++17 new_bank.cpp -o /tmp/carriers-new-bank
/tmp/carriers-new-bank BANK.txt > /tmp/carriers-new-bank.txt
cmp NEW_BANK.txt /tmp/carriers-new-bank.txt
python3 -B verify_witnesses.py > /tmp/carriers-witnesses.json
python3 -B -O verify_witnesses.py > /tmp/carriers-witnesses-optimized.json
cmp WITNESSES.json /tmp/carriers-witnesses.json
cmp WITNESSES_OPTIMIZED.json /tmp/carriers-witnesses-optimized.json
```

The optional old-bank searches rewrite their matching JSON receipts:

```sh
python3 -B search.py > /tmp/carriers-search.txt
python3 -B families.py > /tmp/carriers-families.txt
python3 -B jrays.py > /tmp/carriers-jrays.txt
cmp SEARCH.txt /tmp/carriers-search.txt
cmp FAMILIES.txt /tmp/carriers-families.txt
cmp JRAYS.txt /tmp/carriers-jrays.txt
sha256sum -c SHA256SUMS
```

## Provenance and publication boundary

`PROVENANCE.json` pins the supplied bank, frozen transfer and paired-carrier
proofs, and imported verification definitions to commit
`0f8750677dd20b594a0e032a3b56aa23f58a950c`, with exact repository paths, Git
blob IDs, and SHA-256 hashes. The local copies make reproduction self-contained.
`inherited.py` is frozen; importing it does not run its main program.

`PUBLICATION_FILES.txt` is the explicit publication list. `SHA256SUMS` covers
every listed file except itself, including the publication list. Binaries,
`__pycache__`, and the unused exploratory `old_probe.cpp` are excluded.
