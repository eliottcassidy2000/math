"""Exact finite audit of the all-order structural strong-bridge/forbidden-factor proof."""
from actual_grammar import build
from native_kernel import need, observer, sectors, det, rref_kernel
from pathlib import Path
import collections,json,hashlib

def strong(A,V):
    if len(V)<3:return False
    for rev in (False,True):
        reached={V[0]};todo=[V[0]]
        while todo:
            i=todo.pop()
            for j in V:
                if j not in reached and (A[j][i] if rev else A[i][j]):
                    reached.add(j);todo.append(j)
        if len(reached)!=len(V):return False
    return True

def bridges(A,cyc):
    n=len(A);out=[]
    for off in range(n):
        c=cyc[off:]+cyc[:off]
        for k in range(3,n-2):
            R,Q=c[:k],c[k:]
            cross=[(r,q) for r in R for q in Q if A[r][q]]
            if len(cross)==1 and strong(A,R) and strong(A,Q):out.append({'R':R,'Q':Q,'bridge':list(cross[0])})
    return out

def canonical(A,cycle):
    n=len(A)
    return min(''.join(str(A[c[i]][c[j]]) for i in range(n) for j in range(i+1,n)) for off in range(n) for c in [cycle[off:]+cycle[:off]])

def main():
    here=Path(__file__).resolve().parent;raw=(here/'pinned/tournament_results.json').read_bytes();data=json.loads(raw)
    counts=[];cuts=0;persist=0;guard_steps=0
    for nn,layer in data['layers'].items():
        n=int(nn);bad=0
        for rec in layer['records']:
            w=rec['word'];A,cycle,history=build(w)
            need(canonical(A,cycle)==rec['bits'],('chronological reconstruction mismatch',w))
            bs=bridges(A,cycle);f='be' in w or 'ee' in w
            need(bool(bs)==f,('structural equivalence mismatch',w))
            need(f==(rec['det']==0),('pinned finite observer mismatch',w))
            if f:
                last=max(i for i in range(len(w)-1) if w[i:i+2] in ('be','ee'))
                birth=last+2;c=4+birth-1;u=c-1;r=history[birth-2]['parent'];Q={r,u,c}
                for step in history[birth:]:
                    need(step['parent'] in Q and step['parent']!=c,('last forbidden guard failed',w,step))
                    Q.add(step['new_vertex']);guard_steps+=1
                need(all(A[v][c] for v in Q if v not in (r,c)),('protected sink failed',w))
            if n>4:
                top=[u for u,row in enumerate(A) if sum(row)==n-2]
                need(n-1 in top,'newest vertex not maximal')
                if len(top)==2:
                    other=next(u for u in top if u!=n-1)
                    need(A[other][n-1]==1,'newest vertex not second maximum')
            letters='ab' if sum(sum(row)==n-2 for row in A)==1 else 'cde'
            for b in bs:
                cuts+=1
                for letter in letters:
                    B,_,h=build(w+letter);p=h[-1]['parent']
                    need(p in b['Q'],'legal parent escaped high strong block')
                    Q=b['Q']+[n]
                    need(strong(B,Q) and strong(B,b['R']),'strong bridge not preserved')
                    need(sum(B[r][q] for r in b['R'] for q in Q)==1,'unique bridge not preserved')
                    persist+=1
            bad+=bool(bs)
        counts.append({'order':n,'words':len(layer['records']),'strong_bridge_words':bad})
    kernel_controls=[]
    for w in ('ee','eec','eed','eee'):
        A,cycle,_=build(w);ports,M=observer(A);P,S=sectors(ports,M)
        kernel_controls.append({'word':w,'native_determinant':det(M),'plus_nullity':len(rref_kernel(P)),'minus_nullity':len(rref_kernel(S))})
    out={'status':'FINITE-EXACT audit of separately proved combinatorial equivalence','pinned_commit':'d3132eafbfa9c43c6d9f1110566b1e9c3a3286e1','pinned_data_sha256':hashlib.sha256(raw).hexdigest(),'layers':counts,'strong_bridge_cuts_checked':cuts,'legal_cut_persistence_checks':persist,'last_forbidden_suffix_steps':guard_steps,'kernel_nonpersistence_controls':kernel_controls,'all_order_observer_equivalence':'OPEN; this checker does not prove a spectral theorem'}
    (here/'actual_bridge_grammar_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
