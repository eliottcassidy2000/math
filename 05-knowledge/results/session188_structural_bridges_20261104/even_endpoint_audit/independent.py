#!/usr/bin/env python3
"""Integer, not bitset, reconstruction of the endpoint parity hostile."""
from pathlib import Path
import json,hashlib
R=(598,508,344,865,40,69,656,317,113,438)
PARENT=Path(__file__).resolve().parent.parent/'even_endpoint_boundary'

def need(ok,why):
    if not ok:raise RuntimeError(why)

def main():
    n=len(R);N=1<<n;G=[[0]*n for _ in range(n)];H=[];starts_by_subset=[0]*N;ends_by_subset=[0]*N
    for start in range(n):
        d=[[0]*n for _ in range(N)];d[1<<start][start]=1
        for S in range(1,N):
            if not S>>start&1:continue
            if S!=(1<<start):
                for last in range(n):
                    if not S>>last&1:continue
                    rest=S^(1<<last)
                    d[S][last]=sum(d[rest][before] for before in range(n) if rest>>before&1 and R[before]>>last&1)
            for last,z in enumerate(d[S]):
                G[start][last]+=z
                if z%2:starts_by_subset[S]^=1<<start;ends_by_subset[S]^=1<<last
        H.append(d[-1])
    gr=[sum((x%2)<<j for j,x in enumerate(row)) for row in G];hr=[sum((x%2)<<j for j,x in enumerate(row)) for row in H]
    published=json.loads((PARENT/'normal.json').read_text());need(gr==published['all_simple_path_parity_rows'] and hr==published['hamilton_path_parity_rows'],'full parity matrices')
    for matrix in (G,H):
        start=sum((sum(row)%2)<<i for i,row in enumerate(matrix));end=sum((sum(matrix[i][j] for i in range(n))%2)<<j for j in range(n));need(start==end==296,'balanced parity vectors')
    need(G[0][1]%2==1 and G[1][0]%2==0,'asymmetric physical pair');need(sum(map(sum,G))%2==1 and sum(G[i][i] for i in range(n))%2==0,'total and diagonal parity')
    hulls=[]
    for i in range(n):
        for j in range(i+1,n):
            M=(1<<i)|(1<<j)
            while True:
                old=M
                for u in range(n):
                    if not M>>u&1 and R[u]&M not in (0,M):M|=1<<u
                if M==old:break
            need(M==N-1,'pair hull proves primality');hulls.append([i,j,M])
    need((PARENT/'normal.json').read_bytes()==(PARENT/'optimized.json').read_bytes(),'producer receipts')
    need(hashlib.sha256((PARENT/'verify.py').read_bytes()).hexdigest()==published['source_sha256'],'producer source pin')
    return dict(status='PASS',rows=R,all_simple_path_integer_matrix=G,hamilton_path_integer_matrix=H,all_simple_path_parity_rows=gr,hamilton_path_parity_rows=hr,balanced_endpoint_mask=296,asymmetric_pair=[0,1],all_simple_total=sum(map(sum,G)),hamilton_count=sum(map(sum,H)),pair_hulls=hulls,source_sha256=published['source_sha256'],subset_start_parities=starts_by_subset,subset_end_parities=ends_by_subset)
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,separators=(',',':')))
