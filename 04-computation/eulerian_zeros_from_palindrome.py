"""
eulerian_zeros_from_palindrome.py
kind-pasteur-2026-03-07-S37

KEY QUESTION: Which F_k are FORCED to be 0 mod 3 for ALL tournaments,
purely from palindrome + c_0=c_1=c_2=0 mod 3 (THM-085)?

METHOD: Over F_3, F(T,x) = sum_{j>=3} c_j (x-1)^j. Palindrome F_k = F_{n-1-k}
imposes linear constraints on c_3,...,c_{n-1}. These constraints force certain
c_j = 0 mod 3, which in turn forces certain F_k = 0 mod 3.

We compare the palindrome-forced zeros with the actual Eulerian zeros A(n,k) = 0 mod 3.
"""

import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
import numpy as np


def eulerian_numbers(n):
    """Compute Eulerian numbers A(n,k) for k=0,...,n-1."""
    if n == 1:
        return [1]
    prev = eulerian_numbers(n - 1)
    A = [0] * n
    for k in range(n):
        A[k] = (k + 1) * prev[k] if k < len(prev) else 0
        if k > 0:
            A[k] += (n - k) * prev[k - 1]
    return A


def analyze_palindrome_constraints_mod3(n):
    """
    Given c_0=c_1=c_2=0 mod 3 and the palindrome F_k=F_{n-1-k},
    determine which c_j and F_k are forced to be 0 mod 3.

    Work over F_3.
    """
    d = n - 1  # degree of F(T,x)

    # Compute (x-1)^j mod 3 for j=3,...,d
    # Store as coefficient arrays of length n (degrees 0 to d)
    powers = {}
    # (x-1)^0 = 1
    p = [0] * n
    p[0] = 1
    for j in range(1, d + 1):
        # Multiply by (x-1)
        new_p = [0] * n
        for k in range(n):
            # x * p[k] -> new_p[k+1]
            if k + 1 < n:
                new_p[k + 1] = (new_p[k + 1] + p[k]) % 3
            # -1 * p[k] -> new_p[k]
            new_p[k] = (new_p[k] - p[k]) % 3
        p = new_p
        if j >= 3:
            powers[j] = list(p)

    # F_k = sum_{j=3}^{d} c_j * [(x-1)^j]_k mod 3
    # This gives n equations (k=0,...,d) in (d-2) unknowns (c_3,...,c_d)

    num_vars = d - 2  # c_3, c_4, ..., c_d
    # Build the matrix M where M[k][j-3] = coefficient of c_j in F_k (mod 3)
    M = [[0] * num_vars for _ in range(n)]
    for k in range(n):
        for j in range(3, d + 1):
            M[k][j - 3] = powers[j][k]

    # Palindrome constraints: F_k = F_{d-k} for k=0,...,floor((d-1)/2)
    # This gives: sum_j c_j * ([(x-1)^j]_k - [(x-1)^j]_{d-k}) = 0 mod 3
    num_constraints = d // 2  # number of palindrome pairs
    constraints = []
    for k in range(num_constraints):
        row = [0] * num_vars
        for j in range(3, d + 1):
            row[j - 3] = (M[k][j - 3] - M[d - k][j - 3]) % 3
        constraints.append(row)

    # Solve the constraints over F_3 using Gaussian elimination
    # Find which variables (c_j) are forced to 0

    # Augmented matrix (constraints | 0)
    mat = [row[:] for row in constraints]
    nrows = len(mat)
    ncols = num_vars

    # Gaussian elimination over F_3
    pivot_cols = []
    row_idx = 0
    for col in range(ncols):
        # Find pivot
        found = -1
        for r in range(row_idx, nrows):
            if mat[r][col] % 3 != 0:
                found = r
                break
        if found == -1:
            continue
        # Swap
        mat[row_idx], mat[found] = mat[found], mat[row_idx]
        pivot_cols.append(col)
        # Scale pivot row
        inv = pow(mat[row_idx][col], -1, 3)  # modular inverse mod 3
        for c in range(ncols):
            mat[row_idx][c] = (mat[row_idx][c] * inv) % 3
        # Eliminate
        for r in range(nrows):
            if r != row_idx and mat[r][col] % 3 != 0:
                factor = mat[r][col]
                for c in range(ncols):
                    mat[r][c] = (mat[r][c] - factor * mat[row_idx][c]) % 3
        row_idx += 1

    # Identify free variables (not pivot columns)
    free_vars = [col for col in range(ncols) if col not in pivot_cols]

    # For each pivot column, the variable c_{col+3} is determined by free vars
    # A variable c_j is FORCED to 0 iff it's a pivot variable whose expression
    # in terms of free variables is identically 0.
    # Actually, c_j is forced to 0 iff for ALL assignments of free variables,
    # c_j = 0. This means c_j has no dependence on any free variable (row is all 0
    # except the pivot).

    forced_zero_cj = []
    for i, col in enumerate(pivot_cols):
        # Check if row i of reduced matrix has any nonzero entries in free columns
        depends_on_free = False
        for fc in free_vars:
            if mat[i][fc] % 3 != 0:
                depends_on_free = True
                break
        if not depends_on_free:
            forced_zero_cj.append(col + 3)  # c_{col+3} is forced to 0

    # Now determine which F_k are forced to 0 mod 3
    # F_k = sum_{j in free} c_j * M[k][j-3] + sum_{j in pivot, forced 0} 0 + ...
    # F_k is forced 0 iff its expression in free variables is identically 0

    # Express F_k in terms of free variables only
    # First, express each c_j (pivot) in terms of free variables
    # From reduced matrix: c_{pivot_col+3} = -sum_{free fc} mat[pivot_row][fc] * c_{fc+3}
    # (in our case, forced_zero means that expression is just 0)

    # Build expression: for each pivot variable, its coefficients in terms of free vars
    pivot_expr = {}  # pivot_var_col -> {free_col: coeff}
    for i, col in enumerate(pivot_cols):
        expr = {}
        for fc in free_vars:
            if mat[i][fc] % 3 != 0:
                expr[fc] = (3 - mat[i][fc]) % 3  # negate
        pivot_expr[col] = expr

    # Now F_k = sum_{all j=3..d} c_j * M[k][j-3]
    # = sum_{free fc} c_{fc+3} * M[k][fc] + sum_{pivot col} c_{col+3} * M[k][col]
    # = sum_{free fc} c_{fc+3} * (M[k][fc] + sum_{pivot col depending on fc} coeff * M[k][col])

    # Coefficient of free variable c_{fc+3} in F_k:
    def fk_free_coeff(k, fc):
        total = M[k][fc]  # direct contribution
        for col, expr in pivot_expr.items():
            if fc in expr:
                total = (total + expr[fc] * M[k][col]) % 3
        return total % 3

    forced_zero_fk = []
    fk_depends = {}
    for k in range(n):
        coeffs = {}
        all_zero = True
        for fc in free_vars:
            c = fk_free_coeff(k, fc)
            if c != 0:
                all_zero = False
                coeffs[fc + 3] = c
        if all_zero:
            forced_zero_fk.append(k)
        fk_depends[k] = coeffs

    return {
        'forced_zero_cj': forced_zero_cj,
        'forced_zero_fk': forced_zero_fk,
        'free_vars': [fc + 3 for fc in free_vars],
        'fk_depends': fk_depends,
        'num_palindrome_constraints': num_constraints,
        'num_vars': num_vars,
        'rank': len(pivot_cols),
    }


# Main analysis
print("=" * 70)
print("Which F_k are FORCED to 0 mod 3 by palindrome + c_0=c_1=c_2=0?")
print("=" * 70)

for n in range(5, 20):
    result = analyze_palindrome_constraints_mod3(n)
    A = eulerian_numbers(n)
    eulerian_zeros = [k for k in range(n) if A[k] % 3 == 0]
    forced = result['forced_zero_fk']

    match = set(forced) == set(eulerian_zeros)
    extra_forced = set(forced) - set(eulerian_zeros)
    missing = set(eulerian_zeros) - set(forced)

    status = "EXACT MATCH" if match else f"MISMATCH (extra={sorted(extra_forced)}, missing={sorted(missing)})"

    print(f"\n  n={n:2d}: Eulerian zeros = {eulerian_zeros}")
    print(f"        Forced zeros  = {sorted(forced)}")
    print(f"        Status: {status}")
    print(f"        Forced c_j=0: {result['forced_zero_cj']}, "
          f"Free c_j: {result['free_vars']}")

    if not match and missing:
        print(f"        ** MISSING zeros need c_j=0 for j in {result['free_vars']} **")
        # Show what each missing F_k depends on
        for k in sorted(missing):
            deps = result['fk_depends'][k]
            terms = [f"{v}*c_{j}" for j, v in sorted(deps.items())]
            print(f"           F_{k} = {' + '.join(terms)} mod 3")


# Summary table
print("\n" + "=" * 70)
print("Summary: palindrome-forced zeros vs Eulerian zeros")
print("=" * 70)
print(f"{'n':>3s} | {'#Eul zeros':>10s} | {'#Forced':>7s} | {'Match?':>6s} | {'#Missing':>8s}")
print("-" * 50)

mismatches = []
for n in range(5, 25):
    result = analyze_palindrome_constraints_mod3(n)
    A = eulerian_numbers(n)
    eulerian_zeros = set(k for k in range(n) if A[k] % 3 == 0)
    forced = set(result['forced_zero_fk'])
    missing = eulerian_zeros - forced
    match = "YES" if not missing else "NO"
    if missing:
        mismatches.append(n)
    print(f"{n:3d} | {len(eulerian_zeros):10d} | {len(forced):7d} | {match:>6s} | {len(missing):8d}")

if mismatches:
    print(f"\nMismatches at n = {mismatches}")
    print("These require ADDITIONAL c_j divisibility (beyond palindrome + THM-085).")
else:
    print("\nAll match! Palindrome + THM-085 explains all Eulerian zeros mod 3.")


print("\nDONE")
