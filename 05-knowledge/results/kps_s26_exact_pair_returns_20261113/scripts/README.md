# Replay of the finite-exact claims (THM-5151)

All sources are plain text. Compile the C files with `gcc -O2`. Run every step in one scratch directory, so that
generated data never enters the repository. The quick standard-library checks live one level up in
`verify_s26.py`.

1. **Complete class lists, orders 1-9.** Build `classgen` from
   `../../first_creators_genealogy_20261107/classgen.c`. Write `cls_n.txt` for n = 1..9, one class code per line.
   Codes are bit words over the pairs i < j in order; bit set means i -> j.
2. **Profiles.** Run `fprof cls_n.txt n > prof_n.txt`. Each line is
   `n code strong nsink nsource H C F_1 .. F_n`. Order 9 takes about a minute.
3. **Proposition 4 (census).** Run `python census.py 4 5 6 7 8 9`. It prints the profile groups that mix strong and
   non-strong classes.
4. **Proposition 5 (C = 1 slice).**
   - Seed with `echo "3 2 4 1" > c1_3.txt`.
   - For n = 3..11 run `c1ext < c1_n.txt > c1_(n+1).txt`. Counts go to stderr: 1, 3, 8, 21, 55, 144, 377, 987,
     2584.
   - For each order run `rprof < c1_n.txt > c1p_n.txt`.
   - Run `python match.py`. It needs `prof_3,5,7,9.txt` from step 2 for the odd strong blocks.
5. **Proposition 3 (no small odd merge).** Run `python refill2.py`.
6. **Pre-registered searches S1 and S2.** Run `ysearch 11 5 targets.txt`, where each target line is
   `c3 F_1 .. F_11`.
   - S1 targets are Q_X / t over the C = 1 order-12 classes.
   - S2 targets are Q(Y5 => Y7) / t over the odd strong pairs.
   - The PROOF section 7 notes give the construction. Validation: `ysearch 7 4` and `ysearch 9 4` on the order-8 and
     order-10 targets return exactly the known partners (H = 95, 317).

`trade.py` (odd-cycle census and the shared `decode`) is used by the census and by the one-off structure checks.
