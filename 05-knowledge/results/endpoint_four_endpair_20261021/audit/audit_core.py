"""Independent symbolic tournament audit, using GF(2) row elimination.
The constant is bit zero; higher bits are independent forward-arc variables.
No producer module is imported. No asserts or external dependencies are used.
"""
from itertools import combinations

def need(ok, msg):
    if not ok: raise RuntimeError(msg)
def pairs(n): return combinations(range(n), 2)
def pair(u,v): return min(u,v),max(u,v)
def path_edges(q): return {pair(a,b) for a,b in zip(q,q[1:])}
def reversal(q): return dict(zip(q,reversed(q)))
def eliminate(equations):
    basis={}
    for equation in equations:
        v=equation
        while v>1:
            p=v.bit_length()-1
            if p in basis: v^=basis[p]
            else: basis[p]=v; break
        if v==1: return None
    return basis

def reduce_form(v,basis):
    answer=v&1
    v&=~1
    while v:
        p=v.bit_length()-1
        if p in basis:
            v^=basis[p]
            answer^=v&1
            v&=~1
        else:
            bit=1<<p;answer^=bit;v^=bit
    return answer

def raw_matrix(n):
    out=[[0]*n for _ in range(n)]
    for k,(u,v) in enumerate(pairs(n),1): out[u][v]=1<<k;out[v][u]=(1<<k)^1
    return out

def blue_equations(A,q):
    r=reversal(q)
    return [A[u][v]^A[r[v]][r[u]] for u,v in pairs(len(q))]
def constraints(n,A):
    return [A[i][i+1]^1 for i in range(n-1)]+[A[n-1][0]^1]+blue_equations(A,list(range(n)))
def k_move(A,q):
    n=len(A); E=path_edges(q)
    return [[0 if u==v else A[u][v]^int(pair(u,v) not in E) for v in range(n)] for u in range(n)]
def symbolic_source(n,I,guards):
    A=raw_matrix(n); eq=constraints(n,A)
    eq += [A[0][v]^int(v not in I) for v in range(1,n-1)]
    eq += [A[u][v]^bit for u,v,bit in guards]
    basis=eliminate(eq);need(basis is not None,('empty source family',n,I,guards))
    return [[0 if u==v else reduce_form(A[u][v],basis) for v in range(n)] for u in range(n)]
def replay_family(A,words):
    n=len(A);source=A; xor=set();black=[]
    for index,q in enumerate(words):
        need(sorted(q)==list(range(n)),('word not permutation',q))
        need(all(A[u][v]==1 for u,v in zip(q,q[1:])),('not uniformly directed',index,q))
        black.append(eliminate(blue_equations(A,q)) is None)
        A=k_move(A,q);xor^=path_edges(q)
    need(xor==path_edges(list(range(n))), 'wrong edge xor')
    need(A==k_move(source,list(range(n))), 'wrong literal labelled target')
    need(all(black),('blue possible',n,words,black))
    return black

def literal_rows(rows):
    n=len(rows);need(all(isinstance(x,int) and 0<=x<1<<n for x in rows),'row range')
    A=[[int(rows[u]>>v&1) for v in range(n)] for u in range(n)]
    need(all(A[u][u]==0 for u in range(n)),'loop')
    need(all(A[u][v]^A[v][u]==1 for u,v in pairs(n)),'not complete tournament')
    return A

def check_literal(A,I,words,expected=None):
    n=len(A);P=list(range(n));need(not any(blue_equations(A,P)),'source not blue')
    need(all(A[i][i+1]==1 for i in range(n-1)) and A[-1][0]==1,'source not closed path')
    need({v for v in range(1,n-1) if A[v][0]}==set(I),'literal actual ports')
    need(sum(A[0])==n-5 and sum(A[-1])==4,'literal endpoint scores')
    B=k_move(A,P);xor=set();flags=[]
    for q in words:
        need(sorted(q)==list(range(n)) and all(A[u][v] for u,v in zip(q,q[1:])), 'literal illegal word')
        flags.append(not any(blue_equations(A,q)));xor^=path_edges(q);A=k_move(A,q)
    need(A==B and xor==path_edges(P),'literal endpoint mismatch')
    if expected is not None: need(flags==expected,('wrong literal flags',flags))
    return flags

def permutation_cycles(p):
    unseen=set(p);out=[]
    while unseen:
        v=min(unseen);cy=[]
        while v in unseen: unseen.remove(v);cy.append(v);v=p[v]
        out.append(cy)
    return out
