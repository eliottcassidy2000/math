# Reviewed low-cycle arithmetic inverse and finite coverage

This is a self-contained frozen reproduction bundle. The mathematical proof is [ARITHMETIC_INVERSE.md](source/session197-lowcycle/ARITHMETIC_INVERSE.md); the independent full audit is [REPORT.md](REPORT.md), and the separate graph-count review is [C2_REPAIRS_REVIEW.md](C2_REPAIRS_REVIEW.md).

The all-order results are a necessary physical response cone and an exact final-gadget inverse for a fixed reachable prefix. The finite conclusion constructs every odd count from3 through10,000,001 except7,21,67 in a strong tournament with at most two Hamilton cycles. This does not prove cofinite or all-order coverage. The three exclusions rely on the separately audited complete two-port compiler; the positive constructions here are independently verified.

The source folder preserves the reviewed local relative layout so that frozen scripts run unchanged. Its helper copies are dependencies, not new theorem promotions. Original audit paths beginning session197-lowcycle or session197-hybrid-constructions are relative to source/ here. The original two-port module lies at source/session197/lowcycle-two-port/.

Run from this directory (Python3 and a C++17 compiler, standard libraries only):

    python3 source/session197-lowcycle/verify_arithmetic_inverse.py
    python3 -O source/session197-lowcycle/verify_arithmetic_inverse.py
    g++ -O3 -std=c++17 independent_coverage.cpp -o /tmp/s197-arithmetic-coverage
    /tmp/s197-arithmetic-coverage image-normal 10000001
    python3 check_inverse.py
    python3 check_hard_witnesses.py
    python3 source/session197-hybrid-constructions/audits/c2-repairs/check.py
    g++ -O3 -std=c++17 source/session197-hybrid-constructions/audits/c2-repairs/walk.cpp -o /tmp/s197-c2-walk
    /tmp/s197-c2-walk < source/session197-hybrid-constructions/audits/c2-repairs/walk-input.txt

The source verifier writes its certificate beside itself; the graph checker rewrites its graph/word receipt and walk input with identical bytes. Run a copy if preserving a pristine checkout matters. The independent coverage executable regenerates packed image files and the two prefix CSV inputs. Packed images are reproducible outputs and are not included; their exact hashes are in REPORT.md and the source certificate. The supplied prefix CSVs allow the independent Python inverse check to run immediately.

The source normal/optimized replay receipts, independent optimized/debug results, both audit reports, all26 C2 graph recipes, all four larger C1 recipes, and their required code dependencies are included. FILES.sha256 binds the bundled bytes. No navigation file is changed by this checkpoint.
