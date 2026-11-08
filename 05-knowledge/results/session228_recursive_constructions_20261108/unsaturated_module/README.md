# Exponential strong gcd-one physical fibers

Read `PROOF.md` for the all-order theorem and its limits. `verify.py` is a self-contained standard-library checker, requiring no other source files or repository imports. It binds `PROOF.md` by hash in its receipt.

Reproduce:

    python verify.py > CHECK.json
    python -O verify.py > CHECK.optimized.json
    cmp CHECK.json CHECK.optimized.json

The exact category-basis and annihilating-polynomial certificates establish the all-order results. Literal whole-tournament subset DP checks every binary choice for t=1,2. Transfer checks exhaust choices only through t=10 and are regression controls, not the all-order argument. The order-four negative control exhausts all64 labeled tournaments and applies only to H5,Q32 blocks. Every validation uses explicit exceptions, so optimization does not remove checks.

Freeze/deliver only PROOF.md, verify.py, README.md, CHECK.json, CHECK.optimized.json. Other scratch files are exploratory, not dependencies or canonical verification inputs.
