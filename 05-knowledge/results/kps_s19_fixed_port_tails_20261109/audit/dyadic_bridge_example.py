"""S19 audit: brute-force the bridge-appending pure-dyadic ray with negative offset (finding M4).
Order-7 core (out-neighbour masks) [94, 84, 32, 6, 12, 27, 60], port 0, found by find_core.c as a core with base-2
state (h0, R0, D0, hu) = (71, 287, -5, 31). After one bridge the state is (353, 543, -31) with R - hu = 512 = 2^9, so
appending further bridges gives h = 353 + 512 (2^(k-1) - 1) = 2^(k+8) - 159 for the word B^k, k >= 1 (k >= 2 shown).
Usage: python dyadic_bridge_example.py TBRUTE_EXE
"""
import sys, subprocess
sys.path.insert(0, '.')
from check_transfer import extend
exe = sys.argv[1]
core = [94, 84, 32, 6, 12, 27, 60]; u = 0
items = [(extend(core, u, 'B' * k), u) for k in range(0, 9)]
inp = "".join(f"{len(o)} {u} " + " ".join(map(str, o)) + "\n" for o, u in items)
res = subprocess.run([exe], input=inp, capture_output=True, text=True).stdout.split("\n")
f = list(map(int, res[0].split()))
print(f"Example: order-7 core rows {core}, port 0: I(Omega,2) = {f[2]}, E = {f[3]}, O = {f[4]}, hu = {f[5]}; "
      f"state (h0, R0, D0, hu) = ({f[2]}, {f[3] + f[4]}, {f[3] - f[4]}, {f[5]})")
for k in range(1, 9):
    H = int(res[k].split()[0])
    tail = f";  2^{k + 8} - 159 = {2 ** (k + 8) - 159}" if k >= 2 else ""
    print(f"  core + B^{k} (order {7 + k}): H = {H}{tail}")
