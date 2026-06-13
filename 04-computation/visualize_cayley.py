#!/usr/bin/env python3
"""
visualize_cayley.py — ASCII visualization of the Cayley line and its structure.

Usage: python3 visualize_cayley.py [--primes] [--eps] [--golden] [--all]
"""

import sys
from math import sqrt, log, atanh, pi

phi = (1+sqrt(5))/2
tau = 1.8392867552141612

def cayley_line(n_max=20, width=70):
    """ASCII visualization of natural numbers on the Cayley line [0,1)."""
    print("THE CAYLEY LINE: Natural Numbers on [0, 1)")
    print("=" * width)
    print()

    # Scale: x from 0 to 1
    bar = [' '] * width
    labels = {}

    for n in range(1, n_max + 1):
        x = (n - 1) / (n + 1)
        pos = int(x * (width - 1))
        pos = min(pos, width - 1)
        if pos not in labels or n < labels[pos]:
            labels[pos] = n
        bar[pos] = '|'

    # Print the bar
    print("0" + "-" * (width - 2) + "1")
    print(''.join(bar))

    # Print labels
    label_line = [' '] * width
    for pos, n in sorted(labels.items()):
        s = str(n)
        for i, c in enumerate(s):
            if pos + i < width:
                label_line[pos + i] = c
    print(''.join(label_line))
    print()

    # Print arctanh distances
    print("Addresses and distances:")
    for n in [1, 2, 3, 5, 7, 11, 13, 17, 19]:
        if n <= n_max:
            x = (n-1)/(n+1)
            d = log(n)/2
            print(f"  n={n:2d}: x={(n-1)}/{(n+1):2d} = {x:.4f}, arctanh = {d:.4f}")

def prime_view(width=70):
    """Show primes on the Cayley line."""
    print("\nPRIMES ON THE CAYLEY LINE")
    print("=" * width)

    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    composites = [n for n in range(2, 50) if n not in primes]

    bar = ['.'] * width
    for p in primes:
        x = (p-1)/(p+1)
        pos = min(int(x * (width-1)), width-1)
        bar[pos] = 'P'

    print("0" + "-" * (width - 2) + "1")
    print(''.join(bar))
    print("P = prime, . = empty")
    print()
    print("Primes ACCUMULATE near x=1 (denser in hyperbolic metric).")
    print(f"Bertrand's postulate: gap [n, 2n] has constant width ln(2)/2 = {log(2)/2:.4f}")

def ep_view():
    """Show the exceptional points."""
    print("\nEXCEPTIONAL POINTS OF M(x)")
    print("=" * 60)
    print()
    print("Discriminant: Delta(x) = 4x(x^2 - 11x - 1)")
    print()

    x_ep1 = 8 - 5*phi
    x_ep2 = 3 + 5*phi

    print(f"  EP0: x = 0.000     d = 0       (trivial)")
    print(f"  EP1: x = {x_ep1:+.4f}   d = 1/phi = {1/phi:.4f}  (golden reciprocal)")
    print(f"  EP2: x = {x_ep2:+.4f}  d = -phi = {-phi:.4f}  (negative golden)")
    print()
    print(f"  EP equation: d^2 + d - 1 = 0 (GOLDEN RATIO)")
    print(f"  EP locations: product = -1, sum = 11, diff = 5*sqrt(5)")
    print()
    print(f"  Delta = -44 at x = -1, 1, 11 (zero, pole, Paley)")
    print()

    # ASCII diagram of Delta(x)
    print("  Delta(x) vs x:")
    xs = [i/10 for i in range(-10, 120)]
    for x in xs:
        delta = 4*x*(x**2 - 11*x - 1)
        if -50 < delta < 50:
            bar_pos = int((delta + 50) / 100 * 50)
            bar_pos = max(0, min(49, bar_pos))
            line = ' ' * bar_pos + '*'
            if abs(x - 0) < 0.05 or abs(x - x_ep1) < 0.05:
                line += f"  <- EP (x={x:.2f})"
            elif abs(x - 1) < 0.05:
                line += f"  <- x=1 (tribonacci)"

def golden_view():
    """Show the golden shadow of each number."""
    print("\nGOLDEN SHADOWS")
    print("=" * 60)
    print()
    print("n   f_n = metallic - 1    CF              Special name")
    print("-" * 60)

    names = {1: "golden reciprocal", 2: "silver - 1", 3: "bronze - 1",
             4: "phi^2 = phi+1"}

    for n in range(1, 13):
        f = (n - 2 + sqrt(n**2 + 4)) / 2
        name = names.get(n, "")
        print(f"  {n:2d}   {f:10.6f}          [{n-1}; {n},{n},{n},...]   {name}")

    print()
    print(f"  Golden ratio phi = {phi:.6f}")
    print(f"  Tribonacci tau   = {tau:.6f}")
    print(f"  Price of memory  = tau - phi = {tau-phi:.6f}")

def formula_card():
    """Print the essential formulas."""
    print("\nESSENTIAL FORMULAS")
    print("=" * 60)
    print()
    print("  Q(x) = (1+x)/(1-x) = exp(2*arctanh(x))")
    print()
    print("  Q^m = 1 + 2*sum g_k(m)*x^k")
    print("  g_k(m) = sum C(k-1,j-1)*C(m,j)*2^{j-1}")
    print()
    print("  k*g_k(m) = m*g_m(k)                    [duality]")
    print("  g_k(-m) = (-1)^k * g_k(m)              [parity]")
    print("  Q^m * Q(-x)^m = 1                      [functional eq]")
    print()
    print("  CV^2(H) = sum 2*g_k(n-2k)/(n)_{2k}")
    print("          = 2/n + 0/n^2 - 14/(3n^3) + ...")
    print()
    print("  arctanh(i) = i*pi/4                     [Wick rotation]")
    print("  Delta(x) = 4x(x^2-11x-1)               [discriminant]")
    print("  EP eigenvalues: 1/phi and -phi           [golden EPs]")
    print()
    print("  Z = (H - n!/2^{n-1}) / (n!/2^{n-1} * sqrt(2/n))")
    print("  [one-number ranking significance test]")

def main():
    args = set(sys.argv[1:])
    if not args or '--all' in args:
        args = {'--line', '--primes', '--eps', '--golden', '--formulas'}

    if '--line' in args:
        cayley_line()
    if '--primes' in args:
        prime_view()
    if '--eps' in args:
        ep_view()
    if '--golden' in args:
        golden_view()
    if '--formulas' in args:
        formula_card()

if __name__ == "__main__":
    main()
