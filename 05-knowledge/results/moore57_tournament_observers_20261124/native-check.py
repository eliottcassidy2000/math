#!/usr/bin/env python3
"""Exact native Gram calibration, two-marker and three-layer route verification."""
from itertools import combinations
import numpy as np
from moore_grid_pipeline import moore

for n in (5,7,9,11,19):
    S=np.array([[0 if i==j else 1 if (j-i)%n<=n//2 else -1 for j in range(n)] for i in range(n)])
    H=S@S.T
    for h in range(1,n//2+1):
        assert H[0,h]==n-4*h

for d in (2,3,7):
    A=moore(d);n=len(A)
    B=A+np.eye(n,dtype=np.int64)
    T=np.zeros((2*n+2,2*n+2),dtype=np.int64)
    for i,j in combinations(range(n),2):
        T[i,j]=T[n+i,n+j]=1
    T[:n,n:2*n]=B
    T[n:2*n,:n]=1-B.T
    T[:n,2*n:]=T[2*n:,n:2*n]=T[2*n,2*n+1]=1
    assert np.array_equal(T+T.T,1-np.eye(2*n+2,dtype=np.int64))
    twins=[(i,j) for i,j in combinations(range(2*n+2),2)
           if all(T[i,k]==T[j,k] for k in range(2*n+2) if k not in (i,j))]
    assert twins==[(2*n,2*n+1)]
    C=B@B.T-1
    np.fill_diagonal(C,0)
    squares=(C@B)*B
    assert np.array_equal(squares,d*np.eye(n,dtype=np.int64)+A)
    # From cross blocks of the three-layer tournament, count cyclic triangles
    # L_u -> M_w -> R_v -> L_u for each ordered pair (u,v).
    U=np.zeros((3*n,3*n),dtype=np.int64)
    for c in range(3):
        for i,j in combinations(range(n),2):U[c*n+i,c*n+j]=1
    U[:n,n:2*n]=B;U[n:2*n,:n]=1-B.T
    U[n:2*n,2*n:]=A;U[2*n:,n:2*n]=1-A.T
    U[2*n:,:n]=1
    assert np.array_equal(U+U.T,1-np.eye(3*n,dtype=np.int64))
    counts=U[:n,n:2*n]@U[n:2*n,2*n:]
    assert np.array_equal(counts,(d-1)*np.eye(n,dtype=np.int64)+1)
    print("degree",d,": unique twins, alignment rectangles and directed-triangle routes verified")
