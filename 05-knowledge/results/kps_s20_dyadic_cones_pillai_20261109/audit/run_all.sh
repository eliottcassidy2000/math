#!/bin/sh
# S20 audit: rebuild and rerun everything (Git Bash / POSIX sh; gcc and python3 required).
# usage: sh run_all.sh CLASS_DIR SCRATCH_DIR AUTHOR_GADGET_DIR
#   CLASS_DIR         : directory with cls_1.txt .. cls_8.txt (S19 classgen output; certified complete by a02)
#   SCRATCH_DIR       : any writable directory for the regenerable gadget tables and value lists (kept out of the repo)
#   AUTHOR_GADGET_DIR : directory with the author's gadgets5.txt and gadgets7c.txt (only for the a09 comparison)
set -e
CLS="$1"; S="$2"; AUTH="$3"; mkdir -p "$S"
gcc -O2 -o a_paths.exe a_paths.c
gcc -O2 -o a_classes.exe a_classes.c
gcc -O2 -o a06_abstract.exe a06_abstract.c
gcc -O2 -o a09_gadgets.exe a09_gadgets.c
gcc -O2 -o a10_closure.exe a10_closure.c
python3 a01_validate_and_tableB.py > a01_validate_and_tableB.out
(for n in 1 2 3 4 5 6 7 8; do ./a_classes.exe $n "$CLS/cls_$n.txt"; done; echo "A000568(1..8) = 1 1 2 4 12 56 456 6880") > a02_class_lists_certified.out
python3 a03_census.py "$CLS" > a03_census.out            # also writes a03_triples_le8.txt (sha256 139d73af...deb6)
python3 a04_cones.py > a04_cones.out
python3 a05_pillai.py > a05_pillai.out
(./a06_abstract.exe 9 40 1; ./a06_abstract.exe 9 40 0; ./a06_abstract.exe 12 60 0; ./a06_abstract.exe 14 60 1) > a06_abstract.out
python3 a07_realised.py > a07_realised.out
python3 a08_tiling.py > a08_tiling.out
./a09_gadgets.exe labelled 5 > "$S/agad5_labelled.txt"
./a09_gadgets.exe classes 7 "$CLS" > "$S/agad7_classes.txt"
python3 a09_compare_gadgets.py "$S" "$AUTH" > a09_gadgets.out
sh a10_run_closures.sh "$S" > a10_closures.out 2>&1
python3 a11_cone_reland_anywhere.py 4 > a11_cone_reland_anywhere.out
python3 a12_reroot_chains.py 8 4 4 15 "$CLS" > a12_reroot_chains.out
python3 a13_misc.py "$CLS" > a13_misc.out
python3 a14_session236_seeds.py > a14_session236_seeds.out
python3 a15_scope_witnesses.py > a15_scope_witnesses.out
sh a16_stress.sh "$S" > a16_stress.out 2>&1
gcc -O2 -o a10b_closure.exe a10b_closure.c
sh a17_tail.sh "$S" > a17_tail.out 2>&1
