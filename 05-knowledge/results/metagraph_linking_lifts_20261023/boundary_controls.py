"""Small exact controls for what even all electrical layers can forget."""
from itertools import combinations
from pathlib import Path
from fractions import Fraction
from collections import deque
import argparse,json,math

CHECKS=0
def need(ok,message):
    global CHECKS
    CHECKS+=1
    if not ok:raise ValueError(message)

def laplacian(n,edges):
    A=[[0]*n for _ in range(n)]
    for a,b,w in edges:
        A[a][a]+=w;A[b][b]+=w;A[a][b]-=w;A[b][a]-=w
    return A

def det(A):
    B=[[Fraction(x)for x in row]for row in A];ans=Fraction(1)
    for i in range(len(B)):
        p=next((j for j in range(i,len(B))if B[j][i]),None)
        if p is None:return 0
        if p!=i:B[i],B[p]=B[p],B[i];ans=-ans
        v=B[i][i];ans*=v
        for j in range(i+1,len(B)):
            u=B[j][i]/v
            for k in range(i,len(B)):B[j][k]-=u*B[i][k]
    need(ans.denominator==1,'integer determinant')
    return ans.numerator

def tree_sum(n,edges):
    total=0
    for selected in combinations(edges,n-1):
        components=[{i}for i in range(n)]
        for a,b,w in selected:
            left=next(s for s in components if a in s);right=next(s for s in components if b in s)
            if left is not right:components.remove(left);components.remove(right);components.append(left|right)
        if len(components)==1:total+=math.prod(w for a,b,w in selected)
    return total

def describe(n,edges):
    A=laplacian(n,edges);D=det([row[:-1]for row in A[:-1]])
    need(D==tree_sum(n,edges)>0,'independent tree enumeration and cofactor determinant')
    return {'degrees':[A[i][i]for i in range(n)],'tree_count':D,'reduced_laplacian':[row[:-1]for row in A[:-1]]}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    one=describe(4,[(0,1,3),(1,2,1),(2,3,3)])
    two=describe(4,[(0,1,1),(1,2,3),(2,3,1),(0,3,2)])
    need(one['degrees']==two['degrees']==[3,4,4,3],'same labelled degrees')
    need([one['tree_count'],two['tree_count']]==[9,17],'both full primary groups odd and distinct')
    need(all((x-y)%2==0 for rx,ry in zip(one['reduced_laplacian'],two['reduced_laplacian'])for x,y in zip(rx,ry)),'same binary reduced matrix')
    need(all(row[0]%3==0 for row in one['reduced_laplacian']) and any(row[0]%3 for row in two['reduced_laplacian']),'one ternary unit potential rejects the dyadically invisible switch')
    first=[(0,1,1),(1,2,1),(2,3,1),(0,4,1),(0,5,1)]
    second=[(0,1,1),(1,3,1),(0,2,1),(2,4,1),(0,5,1)]
    t1=describe(6,first);t2=describe(6,second)
    need(t1['degrees']==t2['degrees']==[3,2,2,1,1,1],'same labelled tree degrees')
    need(t1['tree_count']==t2['tree_count']==1,'all marked discriminant groups trivial')
    def leaf_distances(edges):
        adjacency=[set()for _ in range(6)]
        for a,b,w in edges:adjacency[a].add(b);adjacency[b].add(a)
        distance={0:0};q=deque([0])
        while q:
            for j in adjacency[q.popleft()]:
                if j not in distance:distance[j]=min(distance[k]+1 for k in adjacency[j]if k in distance);q.append(j)
        return sorted(distance[j]for j in range(6)if len(adjacency[j])==1)
    d1,d2=leaf_distances(first),leaf_distances(second)
    need(d1==[1,1,3] and d2==[1,2,2],'nonisomorphic trees by unique branch vertex leaf distances')
    result={'status':'FINITE-EXACT hostile controls in weighted graphs, not claims of tournament realizability',
            'fixed_blue_path_same_degrees':[one,two],
            'same_marked_trivial_linking_pairing':[t1,t2],
            'branch_leaf_distances':[d1,d2],'checks':CHECKS}
    raw=(json.dumps(result,indent=2,sort_keys=True)+'\n').encode();out=Path(__file__).with_suffix('.json')
    if args.write:out.write_bytes(raw)
    else:need(out.read_bytes().replace(b'\r\n',b'\n')==raw,'frozen boundary controls')
    print(json.dumps({'checks':CHECKS,'odd_tree_counts':[9,17],'trivial_group_nonisomorphic_trees':True}))

if __name__=='__main__':main()
