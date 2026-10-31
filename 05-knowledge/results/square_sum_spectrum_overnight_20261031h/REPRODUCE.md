# Reproduction and audit scope

Run from the repository root. The packet directory below is abbreviated only in the shell variable. Python uses the standard library. C++ requires C++17. The independent representative provider is installed nauty `gentourng`2.9.3; its resolved path, executable digest and generated-input digests are frozen in `external_receipts.json`.

```bash
p=05-knowledge/results/square_sum_spectrum_overnight_20261031h
c++ -std=c++17 -O2 "$p/strong_profiles.cpp" -o /tmp/square-strong-profiles
/tmp/square-strong-profiles 7 31 > /tmp/strong_profiles.json
/tmp/square-strong-profiles 8 63 > /tmp/strong_profiles128.json
python3 "$p/strong_reps.py" 9 > /tmp/square-strong-reps9.txt
python3 "$p/strong_reps.py" 8 > /tmp/square-strong-reps8.txt
c++ -std=c++17 -O2 "$p/low_strong_ears.cpp" -o /tmp/square-low-strong-ears
/tmp/square-low-strong-ears 127 9 < /tmp/square-strong-reps9.txt > /tmp/strong_profiles256.json
/tmp/square-low-strong-ears 127 8 < /tmp/square-strong-reps8.txt > /tmp/strong_profiles256.ears.json
python3 "$p/bounded_sources.py"
python3 "$p/cost_census.py" 128
python3 "$p/cost_census.py" 256
python3 "$p/audit_sources.py"
python3 "$p/audit_costs.py"
python3 "$p/pythagorean_glue.py"
python3 "$p/norm_class_glue.py"
python3 "$p/port_synchronization.py"
python3 "$p/verify_packet.py"
```

The consumers intentionally read the frozen packet banks. Compare regenerated `/tmp` files with their corresponding packet files before using them as replacements. Their JSON schemas retain the universe, bounds and witnesses. The two large strong-ear outputs use different label gauges and have different rows/bytes; their profile sets, rather than their presentation lists, are equal.

Every Python consumer also runs under `python3 -O`, with identical output bytes. `strong_profiles.cpp` has an O0/O2 equality control through order7. The full order8 census and low-H strong-ear cover use O2. Absence is justified by the complete declared universes and the all-order reduction; no sampled graph bank supplies an absence claim.

Independent computational paths:

- fixed-Hamilton-path enumeration versus nauty representatives plus strong-ear augmentation;
- direct order9 representatives versus two strong-ear stages from order8;
- subset-cover profile production versus literal permutation covers for all9499 admissible profiles;
- formal block actions versus actual adjacency path counts on all737 small states and all48 positive dyadic states of the larger window;
- edge-set constructions versus degree/connected-traversal consumers and explicit restoration of both cut input edges.

The two strong-ear routes share the same C++ profile/count implementation; they are not two fully independent implementations. Literal permutation profiles and adjacency DP supply the separate arithmetic audits. This packet has same-agent mathematical review, not an independent subagent audit.
