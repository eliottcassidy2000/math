# Session137 complete portable publication replay

**PASS — all requested verification is terminal.** This supersedes the earlier reconciliation snapshot's statements that the coverage and complete independent suites had not yet been rerun. The 233-file publication inventory is unchanged.

## Results

- **65 scientific executions pass**, each byte-for-byte equal to its frozen receipt: 14 coverage Python executions, 42 independent Python executions, 8 independent C++ O2/UBSan executions, and one complete coverage k=7 O3 census.
- **9/9 C++ compilations pass** with `-Wall -Wextra -Werror`.
- **12/12 regenerated JSON sidecar comparisons pass**: three phase-two atlas files and three full-star/deck files, each regenerated in normal and optimized mode.
- **No stderr or UBSan diagnostics** on any scientific run or build.
- **All 233 frozen original files remain unchanged.** All 9 manifests, containing 242 entries, pass both in the original tree and in the portable tree after replay.
- The complete additional execution ledger has **74 records** (65 scientific + 9 builds). The earlier 26 bounded replays are separate; they were not silently added to the 65 count.

## Full k=7 result and scope firewall

Command: `g++ -O3 -std=c++17 -Wall -Wextra -Werror zero_escape.cpp -o zero_escape`, then `./zero_escape 7 quotients7.txt`, from the copied coverage directory.

Elapsed runtime: **1038.263 seconds** (~17 min 18 s); peak resident memory **345,696 KiB** (~337.6 MiB). Output: **301,724 bytes**, SHA256 `b0fa02026c74e5f55b88560d80627764d1632a064efaa253539c198c43e70777`. Exit status 0; stderr empty.

The complete terminal line is:

```
PASS classes 353 rows 2471 profiles 40484864 represented_cuts 2591031296 zero_escape_rows 0
```

This is a **portable replay of the producer's same source code**, not an independently implemented full census. It confirms that the completed producer census relocates and reproduces exactly; it does not remove the producer-only algorithmic scope label, promote scalar zero-rigidity to an all-order theorem, or replace the existing independently implemented selected k=6 audit. The compressed 40,484,864 profiles represent 2,591,031,296 cuts; they are not falsely described as billions of separately recounted literal children.

## Execution controls and terminal records

Coverage ran in dependency order, normal then `-O`; all 21 independent Python checkers ran in both modes. Four independent C++ programs ran under O2 and UBSan. The mathematical producer sources were unchanged. The timing harness used Python `resource.getrusage` for peak child RSS; the terminal ledger excludes harness-only preflight failures.

Toolchain: Python 3.12.14; SymPy 1.14.0 for the symbolic witness checker; g++ 14.2.0 / C++17. Other primary Python checkers use the standard library.

- [Terminal summary](full-replay-summary.json): exact totals and original-file preservation.
- [Execution ledger](full-replay-results.json): commands, terminal status, exit codes, timing, RSS, output hashes, and sidecar comparisons. Expected-receipt paths are relative to the packet root.
- [Postflight](full-replay-postflight.json): both original and portable manifest checks and diagnostics.
- [Source-integrity reconciliation](RECONCILIATION.md): frozen inventory and inherited/reviewed pins.

This publication copy omits transient execution paths and harness bookkeeping; scientific measurements and terminal results are unchanged. Compiled binaries and temporary run products are not included.

## Per-command measurements

All statuses below are PASS. Hashes and literal commands are in the JSON ledger; the tables summarize runtime and memory. Concurrent timings are observed, not a platform-independent performance promise.

| Scientific command tag | Seconds | Peak RSS KiB | Output bytes | SHA256 |
|---|---:|---:|---:|---|
| `coverage-audit-normal` | 14.383 | 18,868 | 1,927 | `63073f4afd3893c311ef7f3f1d957255fa0fbbc743cd8e99b97fe7678c2a92d6` |
| `coverage-audit-O` | 14.213 | 22,252 | 1,927 | `63073f4afd3893c311ef7f3f1d957255fa0fbbc743cd8e99b97fe7678c2a92d6` |
| `coverage-phase2-normal` | 3.610 | 15,708 | 2,481 | `18e7a0722604b08ac076fbc1285ecace44485f4dd8bf2b5caaa545001e851d63` |
| `coverage-phase2-O` | 3.355 | 18,248 | 2,481 | `18e7a0722604b08ac076fbc1285ecace44485f4dd8bf2b5caaa545001e851d63` |
| `coverage-fullstar-normal` | 0.501 | 10,624 | 2,269 | `c2f8c0cf52ceb84026bb91de228b6aeaa2a2b6dedbd633a2c8c12ca3072fcfe6` |
| `coverage-fullstar-O` | 0.519 | 12,972 | 2,269 | `c2f8c0cf52ceb84026bb91de228b6aeaa2a2b6dedbd633a2c8c12ca3072fcfe6` |
| `coverage-deck_leaf_audit-normal` | 2.519 | 9,472 | 1,004 | `bbb7241e6fbe7d2087c82c5d10cdc09f7f9a1d9d07ff71dfae9b9b086ec28641` |
| `coverage-deck_leaf_audit-O` | 2.370 | 12,928 | 1,004 | `bbb7241e6fbe7d2087c82c5d10cdc09f7f9a1d9d07ff71dfae9b9b086ec28641` |
| `coverage-symmetry_edit-normal` | 34.444 | 18,868 | 4,945 | `76eabd8eeb47059ef7201bd9f092f36fa17aee8eeabcce8adc8f83c45c2c6f8e` |
| `coverage-symmetry_edit-O` | 32.452 | 21,864 | 4,945 | `76eabd8eeb47059ef7201bd9f092f36fa17aee8eeabcce8adc8f83c45c2c6f8e` |
| `coverage-zero_closure_controls-normal` | 9.452 | 10,112 | 703 | `b34195a60ef05462c049ccaa2d51e1d6c02ead4b2b5cc9f5280096b09a85cff1` |
| `coverage-zero_closure_controls-O` | 9.588 | 13,152 | 703 | `b34195a60ef05462c049ccaa2d51e1d6c02ead4b2b5cc9f5280096b09a85cff1` |
| `coverage-card_rigidity_audit-normal` | 0.434 | 9,472 | 268 | `8dca41e7598bd5c4977f2ca2131642ad118634493e77820c62dbddc75ef44292` |
| `coverage-card_rigidity_audit-O` | 0.466 | 13,020 | 268 | `8dca41e7598bd5c4977f2ca2131642ad118634493e77820c62dbddc75ef44292` |
| `independent-bounded_interior-normal` | 4.382 | 10,496 | 548 | `4a5138ed1396bcf551b6cf123a27ac8098ceb0b638c14f6470936ba703074ee5` |
| `independent-bounded_interior-O` | 4.209 | 12,740 | 548 | `4a5138ed1396bcf551b6cf123a27ac8098ceb0b638c14f6470936ba703074ee5` |
| `independent-bridge-normal` | 3.517 | 11,900 | 790 | `5bea34390da9814ac9bf449672ea3928a611d4e1fa3296ae1e2b28533e4b2282` |
| `independent-bridge-O` | 3.351 | 14,680 | 790 | `5bea34390da9814ac9bf449672ea3928a611d4e1fa3296ae1e2b28533e4b2282` |
| `independent-budget_and_faces-normal` | 0.092 | 9,984 | 7,245 | `dd1756c37be0bbe920664f624e9423da9093c539941b94dcf021c60cb1be5170` |
| `independent-budget_and_faces-O` | 0.210 | 13,184 | 7,245 | `dd1756c37be0bbe920664f624e9423da9093c539941b94dcf021c60cb1be5170` |
| `independent-card_rigidity-normal` | 3.573 | 9,088 | 407 | `cd49b7174bb4f505e0a712775998ca4b5668d256e04dacf242f45cf383d3b53a` |
| `independent-card_rigidity-O` | 4.047 | 13,156 | 407 | `cd49b7174bb4f505e0a712775998ca4b5668d256e04dacf242f45cf383d3b53a` |
| `independent-constructor_cone-normal` | 0.599 | 9,088 | 738 | `fb53b437e200cf4cf9814747fd7b997a4881188bac129e4ea46f90c295d5232d` |
| `independent-constructor_cone-O` | 0.487 | 12,736 | 738 | `fb53b437e200cf4cf9814747fd7b997a4881188bac129e4ea46f90c295d5232d` |
| `independent-counting_certificates-normal` | 0.675 | 11,904 | 784 | `50a0c0c34ad163d55a208940e0096ec7dcc8c1513a84d1619f078286944d7088` |
| `independent-counting_certificates-O` | 0.830 | 13,484 | 784 | `50a0c0c34ad163d55a208940e0096ec7dcc8c1513a84d1619f078286944d7088` |
| `independent-cycles-normal` | 1.479 | 10,112 | 10,424 | `11cf0305e4f21d6b90176295714dd60b9b7f699413370f9dc336747bf6991543` |
| `independent-cycles-O` | 1.698 | 12,704 | 10,424 | `11cf0305e4f21d6b90176295714dd60b9b7f699413370f9dc336747bf6991543` |
| `independent-cycles15-normal` | 64.894 | 14,208 | 5,685 | `f35ec89ec403f816212b36adf4b82b4c323e6dac4f6e4a954b41ac71b738b52f` |
| `independent-cycles15-O` | 65.727 | 17,392 | 5,685 | `f35ec89ec403f816212b36adf4b82b4c323e6dac4f6e4a954b41ac71b738b52f` |
| `independent-discriminant-normal` | 0.486 | 45,268 | 1,783 | `dfae09faf7f40581f6fb6e5b463b3e93df25cac633530c64b52a331ddf060104` |
| `independent-discriminant-O` | 2.253 | 52,376 | 1,783 | `dfae09faf7f40581f6fb6e5b463b3e93df25cac633530c64b52a331ddf060104` |
| `independent-feedback2-normal` | 0.714 | 9,728 | 1,163 | `dffb2e7de5702d97026cddc7352be2b8465ddc1b9574cb99b9b726f9c34e9272` |
| `independent-feedback2-O` | 0.941 | 12,908 | 1,163 | `dffb2e7de5702d97026cddc7352be2b8465ddc1b9574cb99b9b726f9c34e9272` |
| `independent-formal_capacity-normal` | 0.415 | 21,924 | 508 | `5a7cb8f9b028f387f2e3d38c71b39312fcb42666b71a58445df7687e1c44f7f6` |
| `independent-formal_capacity-O` | 0.570 | 23,920 | 508 | `5a7cb8f9b028f387f2e3d38c71b39312fcb42666b71a58445df7687e1c44f7f6` |
| `independent-fullstar-normal` | 2.479 | 10,496 | 3,160 | `42d8343d99297aa11a8a58e76270cdbd69940423f90b4d378c7dfe6bbc2c3ba1` |
| `independent-fullstar-O` | 2.497 | 13,048 | 3,160 | `42d8343d99297aa11a8a58e76270cdbd69940423f90b4d378c7dfe6bbc2c3ba1` |
| `independent-harmonic-normal` | 5.590 | 9,088 | 976 | `b61b6180cb536dcef664c9f03c5608c85a0cb66f9f59eb8d0b5d2fac2d80bf00` |
| `independent-harmonic-O` | 6.035 | 12,992 | 976 | `b61b6180cb536dcef664c9f03c5608c85a0cb66f9f59eb8d0b5d2fac2d80bf00` |
| `independent-harmonic_tails-normal` | 2.717 | 9,088 | 252 | `b3716715eea22593a1909eec7f4b91e4383469194d198ed244ce547997168742` |
| `independent-harmonic_tails-O` | 2.815 | 12,880 | 252 | `b3716715eea22593a1909eec7f4b91e4383469194d198ed244ce547997168742` |
| `independent-tournaments-normal` | 8.084 | 13,156 | 1,883 | `bdc3cfb885360c5341d356374e00b637f417f76389eed6fe4fabdca13b9ba506` |
| `independent-tournaments-O` | 7.712 | 14,764 | 1,883 | `bdc3cfb885360c5341d356374e00b637f417f76389eed6fe4fabdca13b9ba506` |
| `independent-weighted_gap-normal` | 2.221 | 10,880 | 510 | `74c80b7bedbbdfb2363dd83600fa70d694af5e592ae1900ede43f5297e30554c` |
| `independent-weighted_gap-O` | 2.118 | 13,184 | 510 | `74c80b7bedbbdfb2363dd83600fa70d694af5e592ae1900ede43f5297e30554c` |
| `independent-weighted_polynomial-normal` | 1.152 | 10,752 | 769 | `5a953fb038f72f19b256e0962706907efdc48f15f4003b243b267c896293fa7b` |
| `independent-weighted_polynomial-O` | 1.310 | 13,188 | 769 | `5a953fb038f72f19b256e0962706907efdc48f15f4003b243b267c896293fa7b` |
| `independent-weighted_stability-normal` | 1.490 | 10,624 | 584 | `d179728dde319d9ed7b8bd200e89dc4449e772ae90ed82d17b6055fc17166ca3` |
| `independent-weighted_stability-O` | 1.633 | 13,404 | 584 | `d179728dde319d9ed7b8bd200e89dc4449e772ae90ed82d17b6055fc17166ca3` |
| `independent-width-normal` | 0.834 | 9,088 | 3,079 | `1861122050abf2ca05a55f6c5b6277d98a9e901a9eb0175cbbe501f925706308` |
| `independent-width-O` | 0.925 | 12,732 | 3,079 | `1861122050abf2ca05a55f6c5b6277d98a9e901a9eb0175cbbe501f925706308` |
| `independent-width_states-normal` | 1.231 | 11,264 | 2,061 | `6937bca782f26c5f3bc0cbbded43ccb725d011a3ee62b087369b531d1bc07ad6` |
| `independent-width_states-O` | 1.289 | 13,304 | 2,061 | `6937bca782f26c5f3bc0cbbded43ccb725d011a3ee62b087369b531d1bc07ad6` |
| `independent-width_universality-normal` | 0.787 | 11,520 | 701 | `523f754eeb2ed25a23918d7c064fec66c12eeedf78f531bc3a649768c0ea7987` |
| `independent-width_universality-O` | 0.939 | 14,876 | 701 | `523f754eeb2ed25a23918d7c064fec66c12eeedf78f531bc3a649768c0ea7987` |
| `independent-global_gap-O2` | 2.542 | 9,088 | 459 | `7ea06c6757380f62310fb5768bd87c0b906b0a20aa3d82072b079b4e7711b6cc` |
| `independent-global_gap-UBSan` | 5.146 | 9,088 | 459 | `7ea06c6757380f62310fb5768bd87c0b906b0a20aa3d82072b079b4e7711b6cc` |
| `independent-audit_defects-O2` | 0.405 | 9,088 | 322 | `c207ce558bbda4231f0d69e35af5224624596ee0f5847fa46157c3799dc99369` |
| `independent-audit_defects-UBSan` | 0.696 | 9,088 | 322 | `c207ce558bbda4231f0d69e35af5224624596ee0f5847fa46157c3799dc99369` |
| `independent-feedback2_full-O2` | 0.627 | 9,088 | 774 | `b73ddebe5ae896a72b6c30fb93bbc523e4ffd82c6ff93e6a934e1f06978e4634` |
| `independent-feedback2_full-UBSan` | 1.569 | 9,088 | 774 | `b73ddebe5ae896a72b6c30fb93bbc523e4ffd82c6ff93e6a934e1f06978e4634` |
| `independent-audit_large_blowup-O2` | 0.377 | 74,204 | 388 | `577ccc83c5ef273e2f5c248582da749660fd366bc4358986bad8f4bf302629e2` |
| `independent-audit_large_blowup-UBSan` | 0.659 | 74,228 | 388 | `577ccc83c5ef273e2f5c248582da749660fd366bc4358986bad8f4bf302629e2` |
| `coverage-zero_escape-k7` | 1038.263 | 345,696 | 301,724 | `b0fa02026c74e5f55b88560d80627764d1632a064efaa253539c198c43e70777` |

### Builds

| Build tag | Seconds | Peak RSS KiB |
|---|---:|---:|
| `coverage-zero_escape-k7-compile` | 1.099 | 103,040 |
| `independent-global_gap-O2-compile` | 0.918 | 91,040 |
| `independent-global_gap-UBSan-compile` | 1.430 | 129,892 |
| `independent-audit_defects-O2-compile` | 1.130 | 94,820 |
| `independent-audit_defects-UBSan-compile` | 1.894 | 141,816 |
| `independent-feedback2_full-O2-compile` | 0.688 | 73,624 |
| `independent-feedback2_full-UBSan-compile` | 0.998 | 97,548 |
| `independent-audit_large_blowup-O2-compile` | 0.707 | 74,976 |
| `independent-audit_large_blowup-UBSan-compile` | 1.053 | 97,800 |

Overall two-lane elapsed time: 1039.381 seconds. No active replay process remains, no verification blocker remains, and no Git mutation, fetch, commit or push was performed.
