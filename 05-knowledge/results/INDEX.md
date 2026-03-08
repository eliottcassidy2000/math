# Script Results Index

Maps every script in `04-computation/` to its stored output.
Results are stored in `05-knowledge/results/` with matching filenames (`.out` extension).

**Convention:** When you run a script, save its output:
```bash
python3 04-computation/SCRIPT.py 2>&1 | tee 05-knowledge/results/SCRIPT.out
```

---

## Results catalog

Scripts with stored results are marked [STORED]. Scripts without results are marked [MISSING].

### High-priority results to capture
These scripts have been run but results exist only in /tmp or session transcripts:

| Script | Status | Key finding | Result file |
|--------|--------|------------|-------------|
| `dc_symmetry_path.py` | /tmp only | DC verified; M(T\e) NOT symmetric | MISSING |
| `deletion_contraction_M.py` | /tmp only | M submatrix != M(T-v) | MISSING |
| `M_n7_structure.py` | /tmp only | Diagonal all odd at n=7 | MISSING |
| `reversal_proof_attempt.py` | /tmp only | M(T) != M(T^op), M(T) != M(T^op)^T | MISSING |

### How to bulk-capture results

Run this to capture results for any script:
```bash
for f in 04-computation/*.py; do
    base=$(basename "$f" .py)
    if [ ! -f "05-knowledge/results/${base}.out" ]; then
        echo "Running $f..."
        timeout 300 python3 "$f" > "05-knowledge/results/${base}.out" 2>&1
    fi
done
```

Note: Some scripts take >5 minutes. Use `timeout` appropriately.

---

## Convention for new scripts

1. Save the script in `04-computation/`
2. Run it and pipe output to `05-knowledge/results/SCRIPT.out`
3. Add a row to this INDEX
4. If the result confirms/refutes a hypothesis, update `05-knowledge/hypotheses/INDEX.md`
