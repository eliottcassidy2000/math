# Bounded probes and large-target scales

Seed: 27120261111. Each growth row uses 1,200 iid words; standard deviations are empirical, not confidence guarantees. No atlas was regenerated.

| Law | m | r | mean H/r^m | median H/r^m | sd(log H) |
|---|---:|---:|---:|---:|---:|
| uniform | 16 | 3 | 0.122653 | 0.0805449 | 1.28184 |
| uniform | 32 | 3 | 0.122026 | 0.082547 | 1.27743 |
| uniform | 64 | 3 | 0.124188 | 0.07857 | 1.35146 |
| uniform | 128 | 3 | 0.126995 | 0.0765739 | 1.35784 |
| uniform | 256 | 3 | 0.127007 | 0.0797839 | 1.38395 |
| uniform | 512 | 3 | 0.127039 | 0.0855916 | 1.38333 |
| skewed | 16 | 3 | 0.0521863 | 0.0104175 | 3.06992 |
| skewed | 32 | 3 | 0.0515649 | 0.00969817 | 3.4278 |
| skewed | 64 | 3 | 0.053573 | 0.00976704 | 3.47077 |
| skewed | 128 | 3 | 0.0512261 | 0.00942834 | 3.46595 |
| skewed | 256 | 3 | 0.0523773 | 0.00902196 | 3.46672 |
| skewed | 512 | 3 | 0.0520735 | 0.00965094 | 3.45543 |
| binary03 | 16 | 3 | 0.124279 | 0.0469035 | 2.37972 |
| binary03 | 32 | 3 | 0.118837 | 0.038289 | 2.38777 |
| binary03 | 64 | 3 | 0.117294 | 0.0437514 | 2.40947 |
| binary03 | 128 | 3 | 0.122496 | 0.0429049 | 2.39771 |
| binary03 | 256 | 3 | 0.128462 | 0.0417132 | 2.45283 |
| binary03 | 512 | 3 | 0.118477 | 0.0401048 | 2.3851 |
| binary01 | 16 | 2 | 0.243357 | 0.179871 | 1.38473 |
| binary01 | 32 | 2 | 0.258845 | 0.192151 | 1.42357 |
| binary01 | 64 | 2 | 0.248115 | 0.192094 | 1.40967 |
| binary01 | 128 | 2 | 0.252406 | 0.198057 | 1.40827 |
| binary01 | 256 | 2 | 0.247181 | 0.182831 | 1.41084 |
| binary01 | 512 | 2 | 0.256798 | 0.19025 | 1.40409 |
| constant2 | 16 | 2 | 1.00002 | 1.00002 | 1.69477e-21 |
| constant2 | 32 | 2 | 1 | 1 | 0 |
| constant2 | 64 | 2 | 1 | 1 | 0 |
| constant2 | 128 | 2 | 1 | 1 | 0 |
| constant2 | 256 | 2 | 1 | 1 | 0 |
| constant2 | 512 | 2 | 1 | 1 | 0 |

## Target order scales, not existence certificates

For H=10^d (even: deliberately only a numerical scale), use a nearby eligible odd target if seeking a witness. The lower bound holds for every word. The proxy is the order where the ensemble mean reaches H; it is neither a minimum nor an upper bound.

| d | rigorous language-specific order lower bound | ensemble mean-scale order proxy |
|---:|---:|---:|
| 6 | 14 | 17 |
| 12 | 27 | 30 |
| 30 | 65 | 67 |
| 100 | 211 | 214 |
| 1000 | 2098 | 2100 |
| 1000000 | 2095905 | 2095908 |

## Complete support classification for the two-core iid carrier

| Positive-probability letters | exponential base r |
|---|---:|
| [0] | 1 |
| [1] | 1 |
| [0, 1] | 2 |
| [2] | 2 |
| [0, 2] | 2 |
| [1, 2] | 3 |
| [0, 1, 2] | 3 |
| [3] | 1 |
| [0, 3] | 3 |
| [1, 3] | 2 |
| [0, 1, 3] | 3 |
| [2, 3] | 2 |
| [0, 2, 3] | 3 |
| [1, 2, 3] | 3 |
| [0, 1, 2, 3] | 3 |
