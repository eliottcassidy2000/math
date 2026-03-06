#!/usr/bin/env python3
"""Test: Paley tournament vertex deletion gives maximizer at n-1."""
import sys
sys.path.insert(0, r'C:\Users\Eliott\Documents\GitHub\math\03-artifacts\code')
from tournament_lib import hamiltonian_path_count, find_odd_cycles

# Paley T_7: QR = {1, 2, 4} mod 7
T7 = [[0]*7 for _ in range(7)]
qr7 = {1, 2, 4}
for i in range(7):
    for j in range(7):
        if i != j and (j - i) % 7 in qr7:
            T7[i][j] = 1

h7 = hamiltonian_path_count(T7)
print(f'Paley T_7: H={h7}')

# Delete v=0
T6 = [[T7[i][j] for j in range(7) if j != 0] for i in range(7) if i != 0]
h6 = hamiltonian_path_count(T6)
print(f'T_7 - v=0: H={h6} (global max at n=6 = 45)')

# Claim A check
diff = h7 - h6
print(f'H(T_7) - H(T_7 - v=0) = {diff}')

cycles = find_odd_cycles(T7)
through_0 = [c for c in cycles if 0 in c]
c3 = sum(1 for c in through_0 if len(c) == 3)
c5 = sum(1 for c in through_0 if len(c) == 5)
c7 = sum(1 for c in through_0 if len(c) == 7)
print(f'Cycles through v=0: c3={c3}, c5={c5}, c7={c7}, total={c3+c5+c7}')
print(f'All mu=1 at n=7, sum mu = {c3 + c5 + c7}')
print(f'Check: 2*sum_mu = {2*(c3+c5+c7)}, diff = {diff}, match = {2*(c3+c5+c7)==diff}')

# What about Paley T_3?
T3 = [[0]*3 for _ in range(3)]
for i in range(3):
    for j in range(3):
        if i != j and (j - i) % 3 == 1:
            T3[i][j] = 1

h3 = hamiltonian_path_count(T3)
T2 = [[T3[i][j] for j in range(3) if j != 0] for i in range(3) if i != 0]
h2 = hamiltonian_path_count(T2)
print(f'\nPaley T_3: H={h3}')
print(f'T_3 - v=0: H={h2} (max at n=2 = 1)')

# Beautiful pattern: Paley T_p -> delete v -> H-maximizer at p-1
# T_3 (H=3) -> T_2 (H=1) = max at n=2
# T_7 (H=189) -> T_6 (H=45) = max at n=6

# OEIS A038375: max H(T) = 1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745, 95095
# T_7 -> T_6: 189 -> 45. Check.
# T_11 -> T_10: 95095 -> ???
# max at n=10 = 15745 (from OEIS A038375)

print(f'\nOEIS A038375 max H(T): 1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745, 95095')
print(f'T_7 -> T_6: {h7} -> {h6}. Max at n=6 = 45. Match: {h6 == 45}')
print(f'Diff / 2 = {diff//2} = number of odd cycles through v in T_7')
