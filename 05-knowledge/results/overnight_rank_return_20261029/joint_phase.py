"""Audit one explicit direct-sum Frobenius lift against native S3 transport."""
import json
import numpy as np

from complement_kernel import HERE, need, rank_mod
from terminal_selector import krylov, polynomial


def solve(A, Y):
    A, Y = np.array(A, dtype=np.int64) % 3, np.array(Y, dtype=np.int64) % 3
    if Y.ndim == 1:
        Y = Y[:, None]
    width = A.shape[1]
    Z = np.concatenate((A, Y), axis=1)
    for j in range(width):
        choices = [i for i in range(j, len(Z)) if Z[i, j]]
        need(bool(choices), 'independent chart columns')
        pivot = choices[0]
        Z[[j, pivot]] = Z[[pivot, j]]
        Z[j] = Z[j]*int(Z[j, j]) % 3
        for i in range(len(Z)):
            if i != j:
                Z[i] = (Z[i]-Z[i, j]*Z[j]) % 3
    need(not np.any(Z[width:, width:]), 'exact chart membership')
    return Z[:width, width:]


def mpow(A, n):
    out = np.eye(len(A), dtype=np.int64)
    while n:
        if n & 1:
            out = out @ A % 3
        A = A @ A % 3
        n //= 2
    return out


def xpower(n, coefficients):
    v = [1]+[0]*(len(coefficients)-2)
    for _ in range(n):
        carry = v[-1]
        v = [(x-carry*c) % 3 for x, c in zip([0]+v[:-1], coefficients)]
    return v


def run():
    r = json.loads((HERE/'terminal_selector.json').read_text())
    old = json.loads((HERE/'periodic_defect.json').read_text())
    B, betas = old['transfer'], r['betas']
    h, q = [2, 2, 2, 1, 0, 2, 1], old['ternary_recurrence_low_first']
    base = krylov(B, betas[0], 12)
    columns = np.array([polynomial(B, v, h, 3) for v in base[:6]], dtype=np.int64).T
    c = solve(columns, polynomial(B, betas[1], h, 3))[:, 0]
    need(c.tolist() == [2, 0, 1, 1, 2, 2], 'explicit complementary generator')
    correction = np.array([[a % 3 for a in v] for v in base[:6]], dtype=np.int64).T @ c % 3
    v = [(a-int(b)) % 3 for a, b in zip(betas[1], correction)]
    need(not any(polynomial(B, v, h, 3)), 'squarefree complementary cyclic module')
    K = np.array([[a % 3 for a in row] for row in base+krylov(B, v, 6)], dtype=np.int64).T
    need(rank_mod(K.tolist(), 3) == 18, 'literal common observation chart')
    T = solve(K, np.array(B, dtype=np.int64) @ K % 3)
    S = np.zeros((18, 18), dtype=np.int64)
    S[:12, :12] = np.array([xpower(487*i, q) for i in range(12)], dtype=np.int64).T
    S[12:, 12:] = np.array([xpower(487*i, h) for i in range(6)], dtype=np.int64).T
    need(np.array_equal(mpow(S, 5), np.eye(18, dtype=np.int64)), 'invertible order-five lift')
    need(not np.array_equal(S, np.eye(18, dtype=np.int64)), 'nontrivial Frobenius')
    need(np.array_equal(S @ T % 3, mpow(T, 487) @ S % 3), 'full eighteen-dimensional clock conjugacy')
    actions, defects = [], []
    for permutation, transport in zip(r['permutations'], r['transports']):
        M = solve(K, K[transport])
        need(np.array_equal(T @ M % 3, M @ T % 3), 'native names commute with bulk clock')
        defect = (S @ M-M @ S) % 3
        nonzero = np.argwhere(defect)
        actions.append(M.tolist())
        defects.append(dict(permutation=permutation, rank=rank_mod(defect.tolist(), 3),
                            witness=nonzero[0].tolist() if len(nonzero) else None,
                            value=int(defect[tuple(nonzero[0])]) if len(nonzero) else 0))
    need(defects[0]['rank'] == 0 and all(d['rank'] > 0 for d in defects[1:]),
         'this chosen lift fails every nonidentity native name commutator')
    source = np.array(old['initial'], dtype=np.int64) @ K % 3
    need(not np.array_equal(source @ S % 3, source), 'source functional requires transport')
    selectors = [solve(K, beta)[:, 0].tolist() for beta in betas]
    need(np.array_equal(S @ np.array(selectors[0]) % 3, np.array(selectors[0])),
         'identity selector fixed in this observable-dual chart')
    return dict(complement_coefficients=c.tolist(), chart=K.tolist(), clock=T.tolist(),
                frobenius=S.tolist(), native_actions=actions, commutator_defects=defects,
                source=source.tolist(), source_after= (source @ S % 3).tolist(),
                selectors=selectors)


if __name__ == '__main__':
    import sys
    receipt = json.loads(json.dumps(run()))
    path = HERE/'joint_phase.json'
    if '--write' in sys.argv:
        path.write_text(json.dumps(receipt, separators=(',', ':'))+'\n')
    else:
        need(receipt == json.loads(path.read_text()), 'unchanged joint phase certificate')
    print('PASS explicit joint Frobenius lift and native-name/source obstruction')
