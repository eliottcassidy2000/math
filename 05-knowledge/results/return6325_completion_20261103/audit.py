"""Independent explicit-matrix native replay and forward endpoint counting."""
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
def need(b,s):
    if not b:raise ValueError(s)
def matrix(T):
    n=len(T);A=[[bool(T[i]>>j&1) for j in range(n)] for i in range(n)]
    need(all(not A[i][i] for i in range(n)) and all(A[i][j]!=A[j][i] for i in range(n) for j in range(i)),'tournament relation')
    return A
def birth(T,mask):
    A=matrix(T);n=len(T)
    for i in range(n):A[i].append(bool(mask>>i&1))
    A.append([not bool(mask>>i&1) for i in range(n)]+[False])
    return rows(A)
def rows(A):return tuple(sum(1<<j for j,b in enumerate(row) if b) for row in A)
def native(T,P):
    A=matrix(T);n=len(T);need(sorted(P)==list(range(n)),'full path permutation')
    arcs=set(zip(P,P[1:]));need(all(A[a][b] for a,b in arcs),'literal legal Hamilton guard')
    return rows([[False if i==j else A[i][j] if ((i,j) in arcs or (j,i) in arcs) else not A[i][j] for j in range(n)] for i in range(n)])
def endpoint(T,root=None):
    n=len(T);D=[[0]*n for _ in range(1<<n)]
    for i in range(n):D[1<<i][i]=int(root is None or root==i)
    for s,row in enumerate(D):
        for v,c in enumerate(row):
            if not c:continue
            q=T[v]&~s
            while q:
                bit=q&-q;q-=bit;D[s|bit][bit.bit_length()-1]+=c
    return D
def hp(T):return sum(endpoint(T)[-1])
def profile(T,r):
    D=endpoint(T);R=endpoint(T,r);h=[sum(x) for x in D];h[0]=1;full=len(D)-1
    F=sum(D[m][r]*h[full^m] for m in range(1,len(D)))
    S2=sum(sum(R[m])*h[full^m] for m in range(1,len(D)));ee=D[-1][r]
    need((S2-ee)%2==0,'integral observer')
    return (h[-1],sum(R[-1]),h[full^(1<<r)],F),(h[-1],(S2-ee)//2,(S2+ee)//2)
def strong(T):
    for i in range(len(T)):
        reached=1<<i
        while True:
            nxt=reached
            for j in range(len(T)):
                if reached>>j&1:nxt|=T[j]
            if nxt==reached:break
            reached=nxt
        if reached!=(1<<len(T))-1:return False
    return True
def main():
    summary=[]
    for filename in ('return6325.json','return6325_plus2.json'):
        d=json.loads((HERE/filename).read_text());T=tuple(d['source']);oldnames=d.get('old_names',d['names']);names=d['names'];r=d['root'];spare=d['spare']
        need(oldnames[r]==31 and oldnames[spare]==6,'actual old root and spare')
        need(profile(T,r)[0]==(6325,2107,2107,6325),'inherited source identity')
        if 'permanent_birth' in d:
            b=d['permanent_birth'];T=birth(T,b['cut']);need(T==tuple(b['rows']) and hp(T)==b['H'],'permanent literal birth')
        for step in d['trace']:
            n=len(T);local=names+step['helpers'];a,b=(names.index(x) for x in step['arc']);before=T
            X=birth(T,step['cuts'][0]);P,Q=[[local.index(x) for x in p] for p in step['paths']]
            def frame(p):
                j=p.index(r);return (p[j-1] if j else None,p[j+1] if j+1<len(p) else None)
            need(frame(P)==frame(Q),'same actual root frame')
            Y=native(X,P);Z=native(Y,Q);T=tuple(z&((1<<n)-1) for z in Z[:n])
            need(T==tuple(step['rows']),'literal replay of both native commands')
            delta={(i,j) for i in range(n) for j in range(i+1,n) if ((before[i]>>j)^(T[i]>>j))&1}
            need(delta=={tuple(sorted((a,b)))} and r not in (a,b),'one old arc only; entire root star retained')
            need(strong(T) and hp(T)==step['H_after'],'strong completed repair with exact count')
        need(T==tuple(d['target']),'target identity');p,g=profile(T,r)
        need(p==tuple(d['profile']) and g==tuple(d['target_phases']['golden']),'independent complete endpoint profile')
        V=birth(birth(T,1<<r),1<<len(T));need(V==tuple(d['next_rows']),'literal two-birth reset')
        p2,g2=profile(V,d['next_root']);need(p2[0]==d['next_H']==3*6325+2*(p[3]-p[0]),'actual ternary continuation')
        need(g2==tuple(d['next_phases']['golden']) and strong(V),'next joint golden observer and strongness')
        need(set(oldnames)<=set(d['next_names']) and d['old_deletions']==0,'every old name survives')
        summary.append(dict(file=filename,native_guards=2*len(d['trace']),returned_profile=p,next_profile=p2,next_golden=g2))
        print('independent native audit',filename,p,p2[0],flush=True)
    (HERE/'audit.json').write_text(json.dumps(dict(status='FINITE_EXACT_INDEPENDENT_MATRIX_AND_ENDPOINT_REPLAY',branches=summary),indent=2)+'\n')
    completion=json.loads((HERE/'completion.json').read_text())
    for rec in completion:
        T=tuple(rec['source_rows']);n=len(T);maps=rec['states'];chosen=rec['chosen']
        def response(mask):
            out=[]
            for p in maps:
                U=tuple(sum(1<<j for j in range(n) if T[p[i]]>>p[j]&1) for i in range(n))
                out.append(hp(birth(U,mask)))
            return out
        for c in chosen:need(response(c['mask'])==c['responses'],'literal independent completion query')
        need(len({tuple(c['responses'][i] for c in chosen) for i in range(8)})==8,'joint responses identify all states')
        if n==6:
            allresponses=[response(m) for m in range(1,63)]
            need(all(len(set(row))<8 for row in allresponses),'complete independent one-query impossibility')
            need(all(sum((-1)**i.bit_count()*v for i,v in enumerate(row))==0 for row in allresponses),'complete independent parity cancellation')
    print('independent completion audit passed',flush=True)
if __name__=='__main__':main()
