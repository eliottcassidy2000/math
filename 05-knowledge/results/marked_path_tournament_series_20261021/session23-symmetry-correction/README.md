# Scoped THM-260 denominator repair

Prepared outside the repository; not applied. The frozen research and independent-review packets are unchanged.

`THM260-correction.patch` changes only:
1. THM-260's general mixed-cycle denominator, with a dated explanation and the one-line rooted Burnside proof.
2. A dated unnumbered MISTAKES entry immediately below the maintained purpose paragraph.

The existing correct `(5,3)` table entry, all numerical prefixes, original file headers, and unrelated graded-ring interpretation are untouched. The research packet explicitly does not inherit the latter interpretation. No scarce ID is reserved.

`correction-receipt.json` pins original files, proposed patched files, and the exact patch. `dry-run.txt` records successful nonmutating patch checks against both files. `prepare.py` only reads the chosen repository and writes patch artifacts beside itself; it does not apply the patch or run Git.

Reproduce preparation with `python3 prepare.py /path/to/math`. Publisher should recheck concurrent changes before applying the patch. If MISTAKES has gained another entry, preserve its current heading and purpose paragraph and insert the dated entry beneath them rather than replacing the file.
