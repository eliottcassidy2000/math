"""Auditor: is the rotational tournament R_q (i -> i+1..i+(q-1)/2) / its one-vertex deletion the t0-minimiser?
Compares with the complete-census minima (regular 5..11, almost-regular 4..10).  Usage: python rotational_check.py"""
import os, sys, subprocess
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tourn
HERE = os.path.dirname(os.path.abspath(__file__))
def t0(A):
    n = A.shape[0]
    (H, c, _), = tourn.forest_numbers([tourn.to_code(A)], n)
    return tourn.smallest_zero(tourn.N_from_c(c, n), n)[0], H
def canon(A):
    n = A.shape[0]
    out = subprocess.run([os.path.join(HERE, 'canon_brute.exe'), str(n)], input=tourn.to_code(A) + '\n',
                         capture_output=True, text=True, check=True).stdout.split()
    return out[0]
for q in (5, 7, 9, 11, 13, 15, 17):
    R = tourn.circulant(q, list(range(1, (q - 1) // 2 + 1)))
    D = tourn.delete(R, 0)
    print(f'R_{q}: t0 = {t0(R)[0]:.4f};  R_{q} - v (almost-regular {q - 1}): t0 = {t0(D)[0]:.4f}', flush=True)
# identify the census minimisers
for fn, n in (('ar8_classes.txt', 8), ('ar10_classes.txt', 10)):
    rows = [l.split() for l in open(os.path.join(HERE, '..', fn))]
    D = tourn.delete(tourn.circulant(n + 1, list(range(1, n // 2 + 1))), 0)
    cD = canon(D)
    mine = [l.split() for l in open(os.path.join(HERE, f'canon_ar{n}.out'))]
    idx = [i for i, m in enumerate(mine) if m[0] == cD]
    print(f'n={n}: R_{n + 1} - v is class #{idx} of the census (H = {rows[idx[0]][1] if idx else "?"})')
