#!/usr/bin/env python3
"""Span-3 literal language (pinned minimal 62-dim transfer): exact facts behind SPAN3_LITERAL.md.
usage: verify_span3_code.py span3_intrep.txt"""
import sys, itertools
tok=open(sys.argv[1] if len(sys.argv)>1 else '../../../04-computation/experiments/ribbon_spectrum_20261104/span3_intrep.txt').read().split(); pos=0
def nx():
    global pos; v=int(tok[pos]); pos+=1; return v
dim=nx(); nl=nx(); init=[nx() for _ in range(dim)]; fin=[nx() for _ in range(dim)]
M={}
for _ in range(nl):
    li=nx(); M[li]=[[nx() for _ in range(dim)] for _ in range(dim)]
def mul(A,B): return [[sum(A[i][k]*B[k][j] for k in range(dim) if A[i][k]) for j in range(dim)] for i in range(dim)]
def rank_le_1(P):
    rows=[r for r in P if any(r)]
    if not rows: return True
    r0=rows[0]
    return all(r0[i]*r[j]==r0[j]*r[i] for r in rows[1:] for i in range(dim) for j in range(i+1,dim))
P=M[0]; P2=mul(M[0],P); P3=mul(M[0],P2)
print('1. empty letter: M_0^3 has rank <= 1 over Q (exact 2x2 minors): %s; M_0^2 rank <= 1: %s'%(rank_le_1(P3),rank_le_1(P2)))
print('   so a bottleneck-free literal word never contains 000, and after a trailing 0-run (<= 2) it parses uniquely into')
print('   the 21-word prefix code {0^i a : i <= 2, a in 1..7}')
def H(word):
    v=init[:]
    for a in word: v=[sum(M[a][i][k]*v[k] for k in range(dim)) for i in range(dim)]
    return sum(f*x for f,x in zip(fin,v))
trib=[1,1,1]
for _ in range(20): trib.append(trib[-1]+trib[-2]+trib[-3])
print('2. H(1^k), k = 3..16, equals the tribonacci numbers T(k+1): %s'%all(H((1,)*k)==trib[k] for k in range(3,17)))
