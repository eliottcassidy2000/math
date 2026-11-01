# Fixed source: no calls, arrays, recursion, or sequence tables.
def count(n):
    residue = n % 3
    k = n // 3 + 2
    if residue == 0:
        return 1
    if residue == 1:
        return k
    return (k * (k - 1)) // 2
