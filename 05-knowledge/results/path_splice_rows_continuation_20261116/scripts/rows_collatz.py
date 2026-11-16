#!/usr/bin/env python3
"""rows_collatz.py: the odd numbers in rows of doubling chains L_c = {1 + c 2^s : s >= 1} (c odd), the Collatz map in
count coordinates h = 2x + 1 as the Syracuse map on row indices, the forbidden pair and the near-join exceptions by row,
and the identities 3q + 2 = 3(q+1) - 1, 3q + 4 = 3(q+1) + 1 against the Collatz and 3x-1 odd steps."""
def row(N):            # N odd > 1: (c, s) with N = 1 + c 2^s, c odd
    m = N - 1; s = 0
    while m % 2 == 0: m //= 2; s += 1
    return m, s
def T(h):              # Collatz in count coordinates (CONDITIONAL.md): h = 3 mod 4 -> 3h, h = 1 mod 4 -> (h+1)/2
    return 3*h if h % 4 == 3 else (h + 1)//2
def Tm(h):             # 3x-1 in count coordinates: x odd -> 3x-1 gives h -> 3h - 4; x even -> x/2 gives (h+1)/2
    return 3*h - 4 if h % 4 == 3 else (h + 1)//2
def syr(c):            # Syracuse map on odd c
    m = 3*c + 1
    while m % 2 == 0: m //= 2
    return m
print("rows of the exceptional sets:")
for N in [7, 21, 1, 11, 13, 47, 67]:
    print(f"  N={N}: row L_{row(N)[0]} position s={row(N)[1]}" if N > 1 else "  N=1: on no row (the root)")
print("Collatz in count coordinates descends a row to its head (s = 1) and multiplies the head by 3:")
for N in [7, 21, 13, 25, 49]:
    c, s = row(N); print(f"  N={N} = 1 + {c}*2^{s}: T(N)={T(N)} -> row L_{row(T(N))[0]} s={row(T(N))[1]}; head of L_{c} = {1+2*c}, 3*head = {3*(1+2*c)} = 1 + {row(3*(1+2*c))[0]}*2^{row(3*(1+2*c))[1]}")
print("row index map c -> oddpart(3c+1) is the Syracuse map:", [(c, syr(c)) for c in [1,3,5,7,9,11,23,33]])
print("conjugacy check (x -> h = 2x+1) for x <= 2000:",
      all(2*(3*x+1)+1 == T(2*x+1) for x in range(1, 2001, 2)) and all(2*(x//2)+1 == T(2*x+1) for x in range(2, 2001, 2)),
      all(2*(3*x-1)+1 == Tm(2*x+1) for x in range(1, 2001, 2)))
print("3q+2 = 3(q+1)-1 and 3q+4 = 3(q+1)+1 for odd q: the 3y-1 and 3y+1 images of the even number y = q+1;")
print("in count coordinates the children 3q (digit 0), 3q+2, 3q+4 and the two famous steps 3q (Collatz) and 3q-4 (3x-1) are 3q + 2*sigma with sigma = 0, 1, 2 and 0, -2;")
print("in half-counts x = (q-1)/2 the three children are x -> 3x+1 (Collatz odd step), 3x+2, 3x+3 = 3(x+1):",
      all(((3*q+2)-1)//2 == 3*((q-1)//2)+2 and ((3*q+4)-1)//2 == 3*((q-1)//2)+3 and (3*q-1)//2 == 3*((q-1)//2)+1 for q in range(1, 2001, 2)))
# which rows contain the C2 children of q: 3q+2 = 1 + (3q+1): row index oddpart(3q+1) = syr(q) when ... ; 3q+4 = 1 + (3q+3) = 1 + 3(q+1)
print("row of 3q+2 is L_{oddpart(3q+1)} = L_{syr(q)} (the Syracuse image of q!); row of 3q+4 is L_{3 * oddpart(q+1)}; row of 3q is L_{oddpart(3q-1)}:")
for q in [3, 5, 7, 9, 11, 13, 15, 17]:
    print(f"  q={q}: 3q+2={3*q+2} in L_{row(3*q+2)[0]} (syr(q)={syr(q)}), 3q+4={3*q+4} in L_{row(3*q+4)[0]}, 3q={3*q} in L_{row(3*q)[0]}")
print("so the digit-1 child of q sits on the row of the Syracuse image of q, at height v_2(3q+1).")
