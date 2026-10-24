"""Literal endpoint-recording suspension: a new carrier, not a guard transfer."""
from common import *
from itertools import permutations
def partner(B,p):
    E={frozenset(e) for e in zip(p,p[1:])};n=len(B)
    return [[0 if i==j else B[i][j] if frozenset((i,j)) in E else B[j][i] for j in range(n)] for i in range(n)]
def suspend(B):
    n=len(B);a=n;b=n+1;S=[[0]*(n+2) for _ in range(n+2)]
    for i in range(n):
        for j in range(n):S[i][j]=B[i][j]
        S[a][i]=S[i][b]=1
    S[a][b]=1
    return S
def worder(B,p):
    n=len(B);return sum(B[p[i]][p[j]]<<k for k,(i,j) in enumerate((i,j) for i in range(n) for j in range(i+2,n)))
def closed(w,n):
    C=[(i,j) for i in range(n) for j in range(i+2,n)];idx={p:i for i,p in enumerate(C)};M=(1<<len(C))-1
    r=sum(((w>>k)&1)<<idx[n-1-j,n-1-i] for k,(i,j) in enumerate(C));return sorted({w,r,M^w,M^r})
def main():
    nodes,owner,merged,X,Y=data();r=json.loads((HERE/'contexts570.json').read_text())['best']['indexed_internal_arc'];s=r['source']
    need(s==2393,'frozen cheapest internal-arc witness')
    B=matrix(nodes[s]['canon_a']);p=tuple(r['before']['path']);q=tuple(r['after']['path']);fixed=[i for i in range(8) if p[i]==q[i]]
    aut=[g for g in permutations(range(8)) if all(B[i][j]==B[g[i]][g[j]] for i,j in PAIRS)]
    need(aut==[tuple(range(8))],'literal source rigid; added extremes preserve its automorphism group')
    raw=[h for h in paths(B) if all(h[i]==p[i] for i in fixed)]
    S=suspend(B);out=[];ports=set()
    for h in raw:
        z=(8,)+h+(9,);K=partner(S,z);inner=partner(B,h)
        need(all(S[u][v] for u,v in zip(z,z[1:])),'actual suspended Hamiltonian path')
        need([row[:8] for row in K[:8]]==inner,'literal marked inner target restriction')
        need([i for i in range(10) if K[8][i]]==[h[0]],'outer source port records the moving inner first owner')
        need([i for i in range(10) if K[i][9]]==[h[-1]],'outer sink port records the moving inner last owner')
        need([sum(K[i]) for i in range(8)]==[8-sum(B[i]) for i in range(8)],'suspended target inner score vector independent of path')
        need(sum(K[8])==1 and sum(K[9])==8,'fixed outer target scores')
        need(partner(K,z)==S,'literal complement involution')
        C=list(map(list,zip(*S)));rev=z[::-1]
        need(partner(C,rev)==list(map(list,zip(*K))),'converse square on actual enlarged matrices')
        w=worder(S,z);O=closed(w,10);need(len(O)==4 and not ports&set(O),'distinct enlarged normalized resources')
        ports.update(O);oldw=word(B,h);t=merged[owner[FULL^oldw]]
        out.append({'path':z,'inner_target':t,'inner_word':oldw,'word':w,'closed_ports':O,'target_matrix':K})
    need(len(raw)==3 and {o['inner_target'] for o in out}=={410,570,2234},'complete fresh common-clock domain includes the third path')
    need(sum(a!=b for a,b in zip(p,q))==3,'three literal inner owner-clock changes; outer clocks fixed')
    # The hostile distinguishes a source suspension from a two-sided quotient map.
    K=partner(S,(8,)+p+(9,));naive=suspend(partner(B,p))
    need(K!=naive,'suspension does not commute with complement if the boundary is erased')
    return {'status':'PROVED literal source-record suspension with FINITE-EXACT witness; new-carrier dyadic transfer OPEN',
            'source':s,'source_matrix':B,'suspended_source':S,'old_fixed_clock_owners':[(i,p[i]) for i in fixed],
            'fixed_clock_owners':[(0,8)]+[(i+1,p[i]) for i in fixed]+[(9,9)],'owner_clock_cost':3,
            'source_aut_order':len(aut),'domain':out,'closed_resources':3,'closed_ports':len(ports),
            'projected_target_labels':{str(o['inner_target']):label(X,Y,o['inner_target']) for o in out},
            'guard_boundary':'inner deletion transports the marked original target; harmonicity in the full order10 metagraph has not been established',
            'input_revision':REV,'inputs':INPUTS}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
