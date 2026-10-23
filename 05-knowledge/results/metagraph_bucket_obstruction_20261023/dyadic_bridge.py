"""The saved counterfeit is invisible mod 2 but visible in half energy."""
from pathlib import Path
from hashlib import sha256
import json
import argparse

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def need(ok, why):
    if not ok:
        raise ValueError(why)


def read(path):
    raw=path.read_bytes().replace(b'\r\n',b'\n')
    return json.loads(raw),sha256(raw).hexdigest()


def lap(M):
    n=len(M)
    A=[[M[i][j]+M[j][i] if i!=j else 0 for j in range(n)] for i in range(n)]
    return [[sum(A[i]) if i==j else -A[i][j] for j in range(n)] for i in range(n)]


def main():
    source,ss=read(ROOT/'04-computation/experiments/electric_metagraph_20261021/electric_n6.json')
    proof,ps=read(HERE/'receipt.json')
    B=source['black_closed_to_open']; L=source['blue_closed_to_open']; n=len(B)
    fake=[row[:] for row in B]
    for (i,j),delta in proof['counterfeit_n6']['merged_delta']:
        fake[i][j]+=delta
    need(all(v>=0 and v%2==0 for row in fake for v in row),'literal black twins')
    G=[[B[i][j]+L[i][j] for j in range(n)] for i in range(n)]
    F=[[fake[i][j]+L[i][j] for j in range(n)] for i in range(n)]
    A=lap(G); C=lap(F); blue=lap(L)
    need(all((A[i][j]-C[i][j])%2==0 for i in range(n) for j in range(n)), 'mod-two agreement')
    need(all(A[i][i]==C[i][i] for i in range(n)), 'degree margins')
    u,v=16,12
    need(all(blue[i][u]==blue[i][v]==0 for i in range(n)), 'blue-harmonic unit potentials')
    need(A[u][v]%2==C[u][v]%2==0, 'integral half pairing')
    original=(A[u][v]//2)%2; changed=(C[u][v]//2)%2
    need(original!=changed and (C[u][v]-A[u][v])==-2, 'second-layer witness')
    result={'status':'FINITE-EXACT pinned class-labelled witness', 'source_sha256_lf':ss,
            'switch_receipt_sha256_lf':ps, 'potentials':[u,v],
            'original_half_energy_mod2':original, 'counterfeit_half_energy_mod2':changed,
            'laplacian_entry_difference':-2, 'same_laplacian_mod2':True}
    raw=(json.dumps(result,indent=2,sort_keys=True)+'\n').encode()
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    target=HERE/'dyadic_receipt.json'
    if parser.parse_args().write:target.write_bytes(raw)
    else:need(target.read_bytes().replace(b'\r\n',b'\n')==raw,'frozen bridge receipt')
    print(raw.decode())


if __name__=='__main__':main()
