"""Check: at n=9 with rank(C_R) = r, is top eigenvalue still n?"""
import numpy as np
from math import comb
import random

def tournament_from_bits(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj

def full_constraint_matrix(n):
    E = n*(n-1)//2
    edge_idx = {}
    for i in range(n):
        for j in range(i+1, n):
            edge_idx[(i,j)] = len(edge_idx)
    rows = []
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                row = np.zeros(E)
                row[edge_idx[(b,c)]] = 1
                row[edge_idx[(a,c)]] = -1
                row[edge_idx[(a,b)]] = 1
                rows.append(row)
    return np.array(rows), edge_idx

n = 9
E = n*(n-1)//2  # 36
r = (n-1)*(n-2)//2  # 28
Cn3 = comb(n,3)  # 84
C_full, edge_idx = full_constraint_matrix(n)

random.seed(42)
violations = 0
top_eig_violations = 0

for trial in range(50000):
    bits = random.randint(0, (1<<E)-1)
    adj = tournament_from_bits(n, bits)

    trans_idx = []
    cyclic_idx = []
    row = 0
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                wins = [adj[a][b]+adj[a][c], adj[b][a]+adj[b][c], adj[c][a]+adj[c][b]]
                if max(wins) == 2:
                    trans_idx.append(row)
                else:
                    cyclic_idx.append(row)
                row += 1

    c3 = len(cyclic_idx)
    if c3 < r:
        continue  # rank(C_R) ≤ c3 < r, proof works

    C_R = C_full[cyclic_idx]
    C_T = C_full[trans_idx]
    rk_R = np.linalg.matrix_rank(C_R)

    if rk_R >= r:
        violations += 1
        # Check top eigenvalue
        CTC = C_T.T @ C_T
        top = np.max(np.linalg.eigvalsh(CTC))

        if abs(top - n) > 0.01:
            top_eig_violations += 1
            print(f"  TOP EIG VIOLATION: bits=..., c₃={c3}, rank(C_R)={rk_R}, top_eig={top:.6f}")
        elif violations <= 5:
            print(f"  rank(C_R)={rk_R}=r, c₃={c3}, t₃={len(trans_idx)}, "
                  f"top_eig={top:.6f} ← STILL n!")

            # HOW? If ker(C_R) ∩ im(C^T) = {0}, where does the eigenvector come from?
            # Check: is im(C_T^T) ∩ ker(C_R) nonempty?
            cross = C_R @ C_T.T
            rk_cross = np.linalg.matrix_rank(cross)
            rk_T = np.linalg.matrix_rank(C_T)
            eff = rk_T - rk_cross
            print(f"    rank(C_T)={rk_T}, rank(C_R C_T^T)={rk_cross}, "
                  f"dim(im(C_T^T)∩ker(C_R))={eff}")

            # Also check dim(ker(C_R) ∩ im(C^T))
            gap = r - rk_R
            print(f"    dim(ker(C_R)∩im(C^T)) = r-rank(C_R) = {gap}")

            if eff > 0:
                print(f"    *** Proof via im(C_T^T) still works! ***")
            elif gap > 0:
                print(f"    *** Proof via im(C^T) still works! ***")
            else:
                print(f"    *** NEITHER proof works. Need alternative. ***")
                # But top eig IS n... so what's happening?
                # Check eigenvector
                eigs, vecs = np.linalg.eigh(CTC)
                top_vec = vecs[:, np.argmax(eigs)]
                # Is it in im(C^T)?
                U, S, Vt = np.linalg.svd(C_full, full_matrices=False)
                nz = sum(1 for s in S if s > 1e-8)
                P_basis = Vt[:nz]
                proj = P_basis.T @ (P_basis @ top_vec)
                residual = np.linalg.norm(top_vec - proj)
                print(f"      Top eigenvec in im(C^T)? residual={residual:.6f}")
                # What about C_R applied to top eigenvec?
                cr_val = np.linalg.norm(C_R @ top_vec)
                print(f"      ||C_R · top_eigvec|| = {cr_val:.6f}")

if violations == 0:
    print(f"No rank(C_R)≥r cases found in 50000 samples.")
else:
    print(f"\nSummary: {violations} cases with rank(C_R)≥r out of {50000} with c₃≥r")
    print(f"Top eigenvalue violations: {top_eig_violations}")
