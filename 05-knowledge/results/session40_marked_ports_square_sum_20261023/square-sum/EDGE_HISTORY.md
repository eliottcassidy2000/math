# Exact old-edge support through32

F = every path; O = some but not every path; U = unused; – = edge not yet born. Only orders with Hamiltonian paths are columns. No extrapolation beyond32.

|edge|15|16|17|23|25|26|27|28|29|30|31|32|
|---|---|---|---|---|---|---|---|---|---|---|---|---|
|(1, 3)|U|U|U|F|O|U|O|O|O|O|O|O|
|(1, 8)|F|F|F|F|O|O|O|O|O|O|O|O|
|(1, 15)|F|F|F|U|U|O|O|O|O|O|O|O|
|(1, 24)|–|–|–|–|O|O|O|O|O|O|O|O|
|(2, 7)|F|F|F|O|O|O|O|F|U|U|O|O|
|(2, 14)|F|F|F|O|O|F|O|O|F|F|O|O|
|(2, 23)|–|–|–|F|F|O|O|O|F|F|O|O|
|(3, 6)|F|F|F|U|O|O|O|O|O|O|O|O|
|(3, 13)|F|F|F|U|O|O|O|O|O|O|O|O|
|(3, 22)|–|–|–|F|O|O|O|O|O|O|O|O|
|(4, 5)|F|F|F|U|O|O|O|O|F|F|O|O|
|(4, 12)|F|F|F|F|O|O|O|O|O|O|O|O|
|(4, 21)|–|–|–|F|F|F|O|O|O|O|O|O|
|(5, 11)|F|F|F|F|O|O|O|O|F|F|O|O|
|(5, 20)|–|–|–|F|F|F|F|F|U|U|O|O|
|(6, 10)|F|F|F|F|O|O|O|O|O|O|O|O|
|(6, 19)|–|–|–|F|F|F|O|O|O|O|O|O|
|(7, 9)|F|F|F|O|O|O|O|U|U|U|U|O|
|(7, 18)|–|–|–|F|F|F|F|F|F|F|O|O|
|(8, 17)|–|–|F|F|F|F|F|O|O|O|O|O|
|(9, 16)|–|F|F|F|F|F|O|F|F|F|F|O|
|(10, 15)|F|F|F|F|F|O|O|O|O|O|O|O|
|(11, 14)|F|F|F|F|O|O|O|O|O|O|O|O|
|(11, 25)|–|–|–|–|F|F|O|O|O|O|O|O|
|(12, 13)|F|F|F|F|O|O|O|O|F|F|O|O|
|(12, 24)|–|–|–|–|O|O|O|O|O|O|O|O|
|(13, 23)|–|–|–|F|O|O|O|O|O|O|O|O|
|(14, 22)|–|–|–|O|F|O|O|O|O|O|O|O|
|(15, 21)|–|–|–|F|F|F|O|O|O|O|O|O|
|(16, 20)|–|–|–|F|F|F|F|F|F|F|O|O|
|(17, 19)|–|–|–|F|F|F|F|F|F|F|F|O|
|(24, 25)|–|–|–|–|F|F|F|O|O|O|O|O|

## Exact minimum-bank witnesses

Actual unordered endpoints are supplied. These are globally fixed edge-bit banks, not adaptive queries or arbitrary invariants.

|n|minimum|one certified bank|
|---|---:|---|
|25|0|[]|
|26|1|[[3, 6]]|
|27|3|[[1, 15], [2, 14], [4, 12]]|
|28|4|[[2, 14], [1, 15], [3, 22], [3, 6]]|
|29|1|[[4, 12]]|
|30|1|[[4, 12]]|
|31|5|[[3, 6], [2, 14], [4, 5], [1, 15], [12, 13]]|
|32|6|[[4, 21], [2, 14], [3, 6], [5, 11], [1, 8], [12, 24]]|
