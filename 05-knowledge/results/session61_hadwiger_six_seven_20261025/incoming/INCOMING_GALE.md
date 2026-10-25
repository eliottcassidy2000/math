# Incoming canonical balance and Gale return

The concurrent checkpoint `c05e56c19147903a108d2d06b7bb5532508f8ca3` was retained before Session 61's first publication. The parent reviewed its [canonical-balance/Gale proof](https://github.com/eliottcassidy2000/math/blob/c05e56c19147903a108d2d06b7bb5532508f8ca3/05-knowledge/results/six_seven_gale_cube_20261025/BUDGET_AND_GALE.md). Its smooth complete-intersection formula is `K=(rc-n)H`, with `c=n-d-1`; for square-root surfaces the balance compares independent relations and parameters. Equal dimensions do not by themselves identify objects.

The exact six-line Gale self-return retains the minus-one arithmetic twist in the double-plane presentation. This is consistent with, but does not replace, the Session 60/61 marked branch and metric sidecars. No new surface classification or HN bound is inferred here.

The pinned checker `04-computation/six_seven_gale_cube_20261025/verify_budget_gale.py` was materialized directly from that Git commit and replayed in normal and optimized modes. Both outputs match the stored intake receipts. Source SHA256: `975f86c34a7958ed340c1681bcf5101c086d899a0d84d391d29aa0415464a3d8`. `INCOMING_REPLAY.json` also records the later [cube-geometry intake](INCOMING_CUBE.md). Upstream checker sources are not duplicated in this packet.

To reproduce an intake, extract its exact source with `git show <source_commit>:<source_path>` using the fields of `INCOMING_REPLAY.json`, verify the listed SHA256, then run the extracted script normally and with `python3 -O`. These two checkers do not require local producer inputs. The initial source-read claims and independently replayed computations are distinguished from a new independent audit of every upstream theorem.
