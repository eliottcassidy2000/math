# Common-flag periodic tape calculus

Read PROOF.md for the all-k annihilator, generalized-eigenspace argument, minimal mixed-period example, and physically sharp k=1,2 bounds. General-k sharpness and all-target Chan–Pak coverage are not claimed.

    python3 -I -B verify.py > replay.json
    cmp RECEIPT.json replay.json
    python3 -I -B -O verify.py > replay_opt.json
    cmp RECEIPT.json replay_opt.json
    python3 -I -B verify.py --fail

The last command must fail; observed exit code 1. The verifier is standalone stdlib Python, with exact integers and fractions, no external repository imports, no randomness and no optimization-sensitive assertions. Receipt checks are structural, not a numerical spectrum expansion.
