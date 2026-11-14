Reproduction (needs nauty's `gentourng`; run inside this directory):

    cc -O2 -o hl hl.c && cc -O2 -o hubcut hubcut.c
    for n in 3 4 5 6 7 8 9; do gentourng $n 2>/dev/null | ./hl > hl$n.txt; done   # ~1 min, ~25 MB, not committed
    python3 staircase.py            # Busch floors, staircase steps, strong runs (section 1)
    python3 cover_cuts.py 400000    # shallow Fermat-block cuts: no miss in [69, 398129]; first miss 398131 (section 2)
    python3 gen_onecentre.py 10 | ./hl > oc10.txt && python3 hubL.py   # spine-port transfer vs brute force (section 3)
    ./hubcut 200000 20              # hub-cut family coverage (section 3)
    python3 transfer_check.py       # R_1 transfer at the three new roots, 300 random checks (section 4)

`hl` reads gentourng lines (upper-triangle bits, 1 = i->j for i<j) and prints, for every root v:
H(A-v), L(A,v), end_v(A), start_v(A); L(A,v) = ordered pairs (P,Q) of disjoint paths covering A-v,
empty allowed, with Q empty or v->first(Q).
