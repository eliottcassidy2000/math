# Session182 independent closure-candidate audit

Status: FINITE-EXACT / VERIFIED. Neither witness proves repeatable closure.

The pinned `source.json` is Session181's 17-name state with H=512425.
`audit.py` recomputes endpoint subset dynamic programming from candidate rows,
checks the tournament, strong connectivity, physical inports, transitive
prefix and prefix-to-suffix orientation, old-name retention, and actual changed
pairs. It derives h,F,J and E2 independently and verifies H=h-J+2^m F.

`native_audit.py` independently finds a strong guarded arc-flip schedule. Each
flip is realized by a temporary helper birth, two Hamiltonian-path complement
operations, and deletion of only that helper. It checks both path predicates,
actual resulting rows, every completed state's strong connectivity, and exact
counts using the independent C++ endpoint subset DP in `count.cpp`.
Temporary helper names start at100, disjoint from all retained names.

From this directory run:

```sh
python3 audit.py 2
python3 audit.py 4
python3 native_audit.py 2
python3 native_audit.py 4
```

Repeat with `python3 -O` for optimization-mode validation. No checks rely on
Python assertions. `g++` is required and no software is installed by these
scripts. Candidate files are sibling `../closure/candidate2.json` and
`../closure/candidate4.json`.

Results:

- +2: H=1537277=3*512425+2, (h,F,J)=(80637,182080,0),
  E2=1537277; root39 has sole inport21. Three old-arc flips and birth55;
  six native complements; no old-name deletion.
- +4: H=1537279=3*512425+4, (h,F,J)=(58063,190302,43200),
  E2=1580479; root39 has inports33,40. Six old-arc flips and birth55;
  twelve native complements; no old-name deletion.

The counts along completed native steps are in `native_audit2.json` and
`native_audit4.json`. This audit does not establish a phase-transfer theorem,
preservation of another designated frame, or uniform continuation for arbitrary
input states. Joint phase recomputation is separate.
