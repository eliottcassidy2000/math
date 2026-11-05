# Reproduce the binary-module rebasing packet

The theorem uses deletion followed by a birth to transfer a vertex between two modules while preserving H. It is not an induced extension or a universal coverage theorem. The root and four retained old ports are specified in PROOF.md. All mathematical and reviewer files preserve their frozen bytes.

Run in this directory:

```sh
sha256sum -c RELEASE_FILES.sha256
python verify.py > /tmp/rebase-normal.json
python -O verify.py > /tmp/rebase-optimized.json
cmp /tmp/rebase-normal.json VALIDATION_NORMAL.json
cmp /tmp/rebase-optimized.json VALIDATION_OPTIMIZED.json
```

The verifier rewrites VALIDATION.json with the same receipt. The graph counter is bundled. The two named inherited mathematical dependencies are copied under dependencies/.

Run the independent audit from audit/:

```sh
g++ -std=c++17 -O3 check.cpp -o /tmp/check-rebase-audit
/tmp/check-rebase-audit > /tmp/rebase-audit.json
cmp /tmp/rebase-audit.json independent-normal.json
python check_tensor.py > /tmp/rebase-tensor.json
cmp /tmp/rebase-tensor.json tensor-normal.json
```

The independent checker has no source imports. The tensor checker uses only the bundled independent-normal.json and Python standard-library exact fractions. REPORT.md gives the debug/UBSan and optimized Python commands. No large data download or graph bank is required.
