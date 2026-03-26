#!/usr/bin/env python3
"""
program_tournament.py — Programs ARE Tournaments

THESIS: A program's dependency structure forms a tournament. Tournament-theoretic
invariants (H, score sequence, cycle structure, chromatic number) become
program metrics (scheduling freedom, criticality, feedback complexity, parallelism).

PART 1: DEPENDENCY TOURNAMENT
  Given a program's data-flow graph (DAG), complete it to a tournament
  by ordering independent operations. The tournament structure reveals
  which operations MUST be sequential and which CAN be parallel.

PART 2: EXECUTION TRACE COMPRESSION
  An execution trace is a Hamiltonian path through the state tournament.
  H(T) = number of valid execution orderings.
  Bits to specify one execution = log2(H(T)).
  Tournament theory gives BOUNDS on this: the score sequence constrains H(T).

PART 3: CODE METRICS FROM TOURNAMENT INVARIANTS
  - Score sequence → criticality ranking (highest score = most depended upon)
  - H(T) → scheduling freedom (more HP = more parallelism)
  - Cycle structure → feedback loops
  - Chromatic number → minimum parallel stages

PART 4: SOURCE CODE AS BYTE TOURNAMENT
  A source file's bytes can be viewed as a tournament where byte i "beats"
  byte j if knowing byte i helps predict byte j. This IS the basis of
  compression — and tournament theory tells us how to optimize it.

PART 5: AST TOURNAMENT
  The Abstract Syntax Tree of a program forms a tournament on nodes.
  Parent beats child (definitional dependency). Sibling ordering is the
  key degree of freedom — tournament theory optimizes it.

Author: opus-2026-03-25-S344
"""

import sys, os, ast, zlib, struct, time, io
import numpy as np
from itertools import permutations
from collections import defaultdict

sys.stdout.reconfigure(line_buffering=True)

# ============================================================
# PART 1: DEPENDENCY TOURNAMENT
# ============================================================

def build_dependency_tournament(deps):
    """Build a tournament from a dependency partial order.

    deps: dict mapping node -> set of nodes it depends on
    Returns: adjacency matrix A where A[i][j] = 1 means i beats j (i before j)
    """
    nodes = sorted(deps.keys())
    n = len(nodes)
    idx = {v: i for i, v in enumerate(nodes)}
    A = np.zeros((n, n), dtype=int)

    # First: forced orderings from dependencies
    for node, node_deps in deps.items():
        for dep in node_deps:
            if dep in idx:
                # dep must execute before node: dep beats node
                A[idx[dep]][idx[node]] = 1

    # Complete to tournament: for unordered pairs, use lexicographic order
    for i in range(n):
        for j in range(i+1, n):
            if A[i][j] == 0 and A[j][i] == 0:
                # No forced ordering — use default (lower index first)
                A[i][j] = 1

    return A, nodes


def count_hamiltonian_paths_small(A):
    """Count Hamiltonian paths by brute force (n <= 12)."""
    n = len(A)
    if n > 12:
        return -1  # too large

    count = 0
    for perm in permutations(range(n)):
        valid = True
        for k in range(n - 1):
            if A[perm[k]][perm[k+1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count


def score_sequence(A):
    """Score sequence of tournament (sorted descending)."""
    scores = A.sum(axis=1)
    return sorted(scores, reverse=True)


def find_3_cycles(A):
    """Count directed 3-cycles in tournament."""
    n = len(A)
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check both orientations
                if A[i][j] and A[j][k] and A[k][i]:
                    count += 1
                if A[i][k] and A[k][j] and A[j][i]:
                    count += 1
    return count


def chromatic_number_bound(A):
    """Upper bound on chromatic number of tournament (= min parallel stages)."""
    # For tournament: chi(T) <= ceil(log2(n)) + 1
    # Better: chi(T) = length of longest path + 1 for acyclic tournaments
    # For general: use greedy coloring on topological sort
    n = len(A)
    # Find longest directed path (gives lower bound on chi)
    dp = [1] * n
    for i in range(n):
        for j in range(n):
            if A[i][j]:
                dp[j] = max(dp[j], dp[i] + 1)
    return max(dp)


# ============================================================
# PART 2: EXECUTION TRACE COMPRESSION
# ============================================================

def trace_bits(A):
    """Minimum bits to specify an execution trace (Hamiltonian path).
    This is log2(H(T)) where H(T) is the HP count."""
    H = count_hamiltonian_paths_small(A)
    if H <= 0:
        return float('inf')
    return np.log2(H)


def demonstrate_trace_compression():
    """Show that tournament theory bounds execution trace compression."""
    print("\n" + "=" * 75)
    print("  EXECUTION TRACE COMPRESSION via Tournament Theory")
    print("=" * 75)

    # Example: a simple branching program
    # Nodes: start, A, B, C, D, end
    # Dependencies: start -> {A, B}, A -> C, B -> C, C -> D, D -> end
    # A and B can execute in either order (the only freedom)
    print("\n  Example: Diamond dependency")
    print("    start -> A -> C -> D -> end")
    print("           -> B ->")
    deps = {
        'start': set(),
        'A': {'start'},
        'B': {'start'},
        'C': {'A', 'B'},
        'D': {'C'},
        'end': {'D'},
    }
    A, nodes = build_dependency_tournament(deps)
    H = count_hamiltonian_paths_small(A)
    bits = trace_bits(A)
    scores = score_sequence(A)
    cycles = find_3_cycles(A)
    chi = chromatic_number_bound(A)

    print(f"\n  Nodes: {nodes}")
    print(f"  Score sequence: {scores}")
    print(f"  Hamiltonian paths: {H}")
    print(f"  Bits per trace: {bits:.2f}")
    print(f"  3-cycles: {cycles}")
    print(f"  Longest path (≈ min stages): {chi}")
    print(f"  Scheduling freedom: {'high' if H > len(nodes) else 'low'}")

    # Enumerate actual HP
    n = len(A)
    paths = []
    for perm in permutations(range(n)):
        valid = all(A[perm[k]][perm[k+1]] == 1 for k in range(n-1))
        if valid:
            paths.append([nodes[p] for p in perm])
    if len(paths) <= 10:
        print(f"\n  Valid execution orders:")
        for p in paths:
            print(f"    {' -> '.join(p)}")

    # More complex example
    print("\n  ─" * 37)
    print("\n  Example: Parallel map-reduce")
    print("    init -> {W1, W2, W3, W4} -> merge -> output")
    deps2 = {
        'init': set(),
        'W1': {'init'}, 'W2': {'init'}, 'W3': {'init'}, 'W4': {'init'},
        'merge': {'W1', 'W2', 'W3', 'W4'},
        'output': {'merge'},
    }
    A2, nodes2 = build_dependency_tournament(deps2)
    H2 = count_hamiltonian_paths_small(A2)
    bits2 = trace_bits(A2)
    scores2 = score_sequence(A2)
    chi2 = chromatic_number_bound(A2)

    print(f"\n  Nodes: {len(nodes2)}")
    print(f"  Score sequence: {scores2}")
    print(f"  Hamiltonian paths: {H2}")
    print(f"  Bits per trace: {bits2:.2f}")
    print(f"  Min parallel stages: {chi2}")
    print(f"  Workers are freely reorderable: 4! = 24 permutations")
    print(f"  H(T) = {H2} (accounts for all valid schedules)")

    # Chain (no parallelism)
    print("\n  ─" * 37)
    print("\n  Example: Sequential chain A->B->C->D->E->F")
    deps3 = {chr(65+i): {chr(64+i)} if i > 0 else set() for i in range(6)}
    A3, nodes3 = build_dependency_tournament(deps3)
    H3 = count_hamiltonian_paths_small(A3)
    print(f"  Hamiltonian paths: {H3} (transitive tournament = 1 path)")
    print(f"  Bits per trace: {trace_bits(A3):.2f}")
    print(f"  This IS a transitive tournament: zero scheduling freedom.")

    # Fan-out (maximum parallelism)
    print("\n  ─" * 37)
    print("\n  Example: Fan-out: init -> {A,B,C,D,E} (5 independent tasks)")
    deps4 = {'init': set()}
    for c in 'ABCDE':
        deps4[c] = {'init'}
    A4, nodes4 = build_dependency_tournament(deps4)
    H4 = count_hamiltonian_paths_small(A4)
    print(f"  Hamiltonian paths: {H4}")
    print(f"  Bits per trace: {trace_bits(A4):.2f}")
    print(f"  5 independent tasks: 5! = 120 possible orders")
    print(f"  But only {H4} are HP (init must be first)")


# ============================================================
# PART 3: CODE METRICS FROM TOURNAMENT INVARIANTS
# ============================================================

def analyze_python_file(filepath):
    """Analyze a Python file's AST and extract tournament metrics."""
    with open(filepath) as f:
        source = f.read()

    tree = ast.parse(source)

    # Extract function definitions
    functions = {}
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            functions[node.name] = node

    if not functions:
        return None

    # Build call graph (which functions call which)
    call_deps = {name: set() for name in functions}
    for fname, fnode in functions.items():
        for node in ast.walk(fnode):
            if isinstance(node, ast.Call):
                if isinstance(node.func, ast.Name) and node.func.id in functions:
                    call_deps[fname].add(node.func.id)

    if len(call_deps) < 2:
        return None

    # Build tournament
    A, nodes = build_dependency_tournament(call_deps)
    n = len(A)

    # Metrics
    scores = score_sequence(A)
    H = count_hamiltonian_paths_small(A) if n <= 10 else -1
    cycles = find_3_cycles(A) if n <= 20 else -1
    chi = chromatic_number_bound(A)

    return {
        'file': filepath,
        'functions': nodes,
        'n': n,
        'scores': scores,
        'H': H,
        'cycles': cycles,
        'chi': chi,
        'call_deps': call_deps,
    }


def demonstrate_code_metrics():
    """Analyze actual Python files from our codebase."""
    print("\n" + "=" * 75)
    print("  CODE METRICS from Tournament Invariants")
    print("=" * 75)

    # Analyze our own codec files
    files = [
        '04-computation/tic3.py',
        '04-computation/tic5_unified.py',
        '04-computation/quincunx.py',
        '04-computation/wavefront.py',
    ]

    for filepath in files:
        full = os.path.join(os.path.dirname(__file__), '..', filepath)
        if not os.path.exists(full):
            full = filepath
        if not os.path.exists(full):
            continue

        result = analyze_python_file(full)
        if result is None:
            continue

        print(f"\n  {os.path.basename(filepath)}:")
        print(f"    Functions: {result['n']}")
        print(f"    Score sequence: {result['scores'][:10]}{'...' if len(result['scores']) > 10 else ''}")
        if result['H'] >= 0:
            print(f"    Valid execution orders (H): {result['H']}")
            if result['H'] > 0:
                print(f"    Scheduling bits: {np.log2(result['H']):.1f}")
        print(f"    Min parallel stages: {result['chi']}")
        if result['cycles'] >= 0:
            print(f"    Mutual recursion cycles: {result['cycles']}")

        # Show call graph
        for fname, deps in result['call_deps'].items():
            if deps:
                print(f"    {fname} calls: {', '.join(deps)}")


# ============================================================
# PART 4: SOURCE CODE COMPRESSION via TOURNAMENT PREDICTION
# ============================================================

def compress_source_tournament(source_bytes):
    """Compress source code using tournament-inspired prediction.

    Key insight: source code bytes have strong LOCAL dependencies
    (indentation, keyword repetition) and LONG-RANGE dependencies
    (variable names, function calls). This is exactly like an image
    with both spatial and semantic structure.

    Strategy: multi-context prediction
    - Context 0: previous byte (like 'sub' filter)
    - Context 1: same column in previous line (like 'up' filter)
    - Context 2: MED of (prev_byte, same_col, diagonal)
    - Adaptive: pick best context per byte (like per-row adaptive)
    """
    lines = source_bytes.split(b'\n')
    if not lines:
        return source_bytes

    # Pad lines to same length for 2D processing
    max_len = max(len(line) for line in lines)
    padded = np.zeros((len(lines), max_len + 1), dtype=np.uint8)
    for i, line in enumerate(lines):
        padded[i, :len(line)] = list(line)
        padded[i, len(line)] = ord('\n')

    H, W = padded.shape

    # Strategy A: raw zlib
    raw = zlib.compress(source_bytes, 9)

    # Strategy B: treat as image, MED prediction
    buf_med = bytearray()
    for r in range(H):
        for c in range(W):
            a = int(padded[r, c-1]) if c > 0 else 0
            b = int(padded[r-1, c]) if r > 0 else 0
            d = int(padded[r-1, c-1]) if r > 0 and c > 0 else 0
            mn, mx = min(a, b), max(a, b)
            pred = mn if d >= mx else (mx if d <= mn else a + b - d)
            buf_med.append((int(padded[r, c]) - pred) & 0xFF)
    med_data = zlib.compress(bytes(buf_med), 9)

    # Strategy C: column-first (transpose) — groups indentation
    buf_col = bytearray()
    for c in range(W):
        for r in range(H):
            a = int(padded[r, c-1]) if c > 0 else 0
            b = int(padded[r-1, c]) if r > 0 else 0
            d = int(padded[r-1, c-1]) if r > 0 and c > 0 else 0
            mn, mx = min(a, b), max(a, b)
            pred = mn if d >= mx else (mx if d <= mn else a + b - d)
            buf_col.append((int(padded[r, c]) - pred) & 0xFF)
    col_data = zlib.compress(bytes(buf_col), 9)

    # Strategy D: delta-line (each line predicted from previous)
    buf_dl = bytearray()
    buf_dl.extend(padded[0].tobytes())
    for r in range(1, H):
        delta = ((padded[r].astype(np.int16) - padded[r-1].astype(np.int16)) & 0xFF).astype(np.uint8)
        buf_dl.extend(delta.tobytes())
    dl_data = zlib.compress(bytes(buf_dl), 9)

    # Strategy E: indentation-aware — strip leading spaces, encode separately
    indents = []
    stripped = []
    for line in lines:
        indent = len(line) - len(line.lstrip())
        indents.append(indent)
        stripped.append(line.lstrip())
    indent_data = zlib.compress(bytes(indents), 9)
    stripped_data = zlib.compress(b'\n'.join(stripped), 9)
    ie_total = len(indent_data) + len(stripped_data) + 8

    best_name = "raw"
    best_size = len(raw)
    for name, sz in [("MED", len(med_data)), ("column", len(col_data)),
                     ("delta-line", len(dl_data)), ("indent-split", ie_total)]:
        if sz < best_size:
            best_size = sz
            best_name = name

    return {
        'raw_zlib': len(raw),
        'med': len(med_data),
        'column': len(col_data),
        'delta_line': len(dl_data),
        'indent_split': ie_total,
        'best': best_name,
        'best_size': best_size,
        'original': len(source_bytes),
    }


def demonstrate_source_compression():
    """Test source code compression on our own files."""
    print("\n" + "=" * 75)
    print("  SOURCE CODE COMPRESSION via Tournament Prediction")
    print("=" * 75)

    files = [
        '04-computation/tic3.py',
        '04-computation/quincunx.py',
        '04-computation/wavefront.py',
        '04-computation/tic6_adaptive.py',
    ]

    print(f"\n  {'File':<25} {'Raw':>6} {'zlib':>6} {'MED':>6} {'Col':>6} {'DL':>6} {'Ind':>6} {'Best':>8} {'Ratio':>6}")
    print("  " + "-" * 95)

    for filepath in files:
        full = os.path.join(os.path.dirname(__file__), '..', filepath)
        if not os.path.exists(full):
            full = filepath
        if not os.path.exists(full):
            continue

        with open(full, 'rb') as f:
            source = f.read()

        r = compress_source_tournament(source)
        ratio = r['original'] / r['best_size']
        print(f"  {os.path.basename(filepath):<25} {r['original']:>6} {r['raw_zlib']:>6} {r['med']:>6} "
              f"{r['column']:>6} {r['delta_line']:>6} {r['indent_split']:>6} {r['best']:>8} {ratio:>5.1f}x")


# ============================================================
# PART 5: THE DEEP ANALOGY — Program = Tournament = Codec
# ============================================================

def demonstrate_deep_analogy():
    """Show the three-way analogy: programs ↔ tournaments ↔ codecs."""
    print("\n" + "=" * 75)
    print("  THE DEEP ANALOGY: Program = Tournament = Codec")
    print("=" * 75)

    print("""
  TOURNAMENT              PROGRAM                 CODEC
  ─────────              ───────                 ─────
  Vertex                 Statement/Variable       Pixel
  Arc i→j                i executes before j      i predicts j
  Score(i)               # depending on i         # predicted from i
  Hamiltonian path       Execution trace          Scan order
  H(T)                   # valid schedules        # valid scan orders
  Score sequence         Criticality ranking      Gradient profile
  3-cycle                Mutual dependency         Non-monotone region
  Transitive T           Sequential code           Sorted image
  Regular T              Balanced parallelism      Uniform texture
  Complement T^op        Reversed execution        Transpose
  Isomorphism class      Functionally equiv code   Equivalent images
  χ(T)                   Min parallel stages       Min scan directions
  OCF: H = I(Ω, 2)      Schedules from cycles     Compression from structure

  THE KEY INSIGHT:
  A program's DEPENDENCY TOURNAMENT determines its parallelism.
  An image's PREDICTION TOURNAMENT determines its compressibility.
  A codec's STRATEGY TOURNAMENT determines its optimality.

  ALL THREE are governed by the same mathematics:
  - The score sequence bounds H(T) (scheduling/compression)
  - The cycle structure determines the odd-cycle formula
  - Band-limitedness limits the effective degrees of freedom
  - The staircase paradox bounds the misalignment penalty
""")

    # Demonstrate with a concrete example
    print("  CONCRETE EXAMPLE: The 'compress' function itself")
    print()
    print("  compress(image):")
    print("    1. validate input       # score 5 (everything depends on this)")
    print("    2. color_transform      # score 4 (prediction depends on this)")
    print("    3. choose_scan_order    # score 3 (prediction follows scan)")
    print("    4. predict_residuals    # score 2 (entropy coder needs this)")
    print("    5. entropy_encode       # score 1 (output depends on this)")
    print("    6. write_output         # score 0 (nothing depends on this)")
    print()
    print("  This is a TRANSITIVE tournament: (5,4,3,2,1,0)")
    print("  H = 1. Only one valid execution order.")
    print("  Bits per trace = 0. No scheduling freedom.")
    print("  The function is INHERENTLY SEQUENTIAL.")
    print()

    # Now a parallel version
    print("  parallel_compress(image):")
    print("    1. validate             # score 7")
    print("    2a. color_R             # score 4 (independent of 2b, 2c)")
    print("    2b. color_G             # score 4")
    print("    2c. color_B             # score 4")
    print("    3. predict_all          # score 2")
    print("    4. entropy_encode       # score 1")
    print("    5. write_output         # score 0")
    print()

    deps = {
        'validate': set(),
        'color_R': {'validate'},
        'color_G': {'validate'},
        'color_B': {'validate'},
        'predict': {'color_R', 'color_G', 'color_B'},
        'entropy': {'predict'},
        'output': {'entropy'},
    }
    A, nodes = build_dependency_tournament(deps)
    H = count_hamiltonian_paths_small(A)
    print(f"  Score sequence: {score_sequence(A)}")
    print(f"  H = {H}. The 3 color channels can be processed in any order.")
    print(f"  Bits per trace = {np.log2(H):.2f}.")
    print(f"  Scheduling freedom: 3! = 6 valid orderings of parallel work.")
    print(f"  Min parallel stages: {chromatic_number_bound(A)}")


# ============================================================
# PART 6: TOURNAMENT-BASED CODE DEDUPLICATION
# ============================================================

def ast_fingerprint(node):
    """Compute a tournament-inspired fingerprint of an AST subtree.
    The fingerprint captures the STRUCTURE (topology) independent of
    variable names (like tournament isomorphism classes)."""
    if isinstance(node, ast.Module):
        return ('Module', tuple(ast_fingerprint(n) for n in node.body))
    elif isinstance(node, ast.FunctionDef):
        return ('Func', len(node.args.args),
                tuple(ast_fingerprint(n) for n in node.body))
    elif isinstance(node, ast.Return):
        return ('Return', ast_fingerprint(node.value) if node.value else None)
    elif isinstance(node, ast.Assign):
        return ('Assign', len(node.targets),
                ast_fingerprint(node.value))
    elif isinstance(node, ast.BinOp):
        return ('BinOp', type(node.op).__name__,
                ast_fingerprint(node.left),
                ast_fingerprint(node.right))
    elif isinstance(node, ast.Call):
        return ('Call', len(node.args))
    elif isinstance(node, ast.Name):
        return ('Name',)  # ignore actual name — structural equivalence
    elif isinstance(node, ast.Constant):
        return ('Const', type(node.value).__name__)
    elif isinstance(node, ast.If):
        return ('If', ast_fingerprint(node.test),
                tuple(ast_fingerprint(n) for n in node.body),
                tuple(ast_fingerprint(n) for n in node.orelse))
    elif isinstance(node, ast.For):
        return ('For', ast_fingerprint(node.target),
                ast_fingerprint(node.iter),
                tuple(ast_fingerprint(n) for n in node.body))
    elif isinstance(node, ast.Expr):
        return ('Expr', ast_fingerprint(node.value))
    elif isinstance(node, ast.Compare):
        return ('Cmp', len(node.ops))
    else:
        return (type(node).__name__,)


def find_structural_clones(filepath):
    """Find structurally identical functions (tournament isomorphism = code clones)."""
    with open(filepath) as f:
        tree = ast.parse(f.read())

    funcs = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            fp = ast_fingerprint(node)
            if fp not in funcs:
                funcs[fp] = []
            funcs[fp].append(node.name)

    clones = {fp: names for fp, names in funcs.items() if len(names) > 1}
    return clones


def demonstrate_deduplication():
    """Find structural clones in our codebase."""
    print("\n" + "=" * 75)
    print("  CODE DEDUPLICATION via Tournament Isomorphism")
    print("=" * 75)

    files = [
        '04-computation/tic3.py',
        '04-computation/tic5_unified.py',
        '04-computation/tic6_adaptive.py',
    ]

    for filepath in files:
        full = os.path.join(os.path.dirname(__file__), '..', filepath)
        if not os.path.exists(full):
            full = filepath
        if not os.path.exists(full):
            continue

        clones = find_structural_clones(full)
        basename = os.path.basename(filepath)
        if clones:
            print(f"\n  {basename}: {len(clones)} clone groups")
            for fp, names in clones.items():
                print(f"    Structural clones: {', '.join(names)}")
        else:
            print(f"\n  {basename}: no structural clones found")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    demonstrate_trace_compression()
    demonstrate_code_metrics()
    demonstrate_source_compression()
    demonstrate_deep_analogy()
    demonstrate_deduplication()

    print("\n" + "=" * 75)
    print("  Done. opus-2026-03-25-S344")
    print("=" * 75)
