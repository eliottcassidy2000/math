"""AUDIT: isomorphism tests with nauty (dreadnaut): canonical forms of digraphs.
CD_k (audit's construction) against the doubling tower D(CD_(k-1)) of THM-448 and against converses."""
import subprocess, sys, hashlib
from aud_common import cd_tournament, double_tournament, paley


def canon(A):
    n = len(A)
    lines = ["n=%d d" % n, "g"]
    for i in range(n):
        lines.append("%d : %s;" % (i, " ".join(str(j) for j in range(n) if A[i][j])))
    lines += ["c x", "b", "q"]
    out = subprocess.run(["dreadnaut"], input="\n".join(lines) + "\n", capture_output=True, text=True).stdout
    grp = [l for l in out.split("\n") if "grpsize" in l]
    # canonical graph = everything after the line containing 'grpsize' ... take adjacency block lines containing ':'
    idx = out.rfind("grpsize")
    tail = out[idx:]
    adj = "\n".join(l.strip() for l in tail.split("\n") if ":" in l and "grpsize" not in l)
    return hashlib.sha256(adj.encode()).hexdigest()[:16], (grp[-1].strip() if grp else "?")


def conv(A):
    n = len(A)
    return [[A[j][i] for j in range(n)] for i in range(n)]


print("sanity: P_7 vs CD_3:", canon(paley(7))[0] == canon(cd_tournament(3))[0], "; P_7 vs its converse:", canon(paley(7))[0] == canon(conv(paley(7)))[0])
for k in (4, 5, 6):
    C = cd_tournament(k)
    D = double_tournament(cd_tournament(k - 1))
    cC, cD, cCc, cDc = canon(C), canon(D), canon(conv(C)), canon(conv(D))
    print("k = %d (order %d): |Aut(CD_k)|: %s ; |Aut(D(CD_(k-1)))|: %s" % (k, len(C), cC[1], cD[1]))
    print("   CD_k isomorphic to D(CD_(k-1)): %s; to the converse of D(CD_(k-1)): %s; CD_k self-converse: %s; D(CD_(k-1)) self-converse: %s" % (
        cC[0] == cD[0], cC[0] == cDc[0], cC[0] == cCc[0], cD[0] == cDc[0]))
