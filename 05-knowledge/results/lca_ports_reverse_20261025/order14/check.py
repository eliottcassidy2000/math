"""Four isolated order14 symmetry exceptions and one literal repair.

Complete affine families are checked by GF(2) equations. The new repair
is independently realized from the two seven-cycle symmetry and row0.
No tournament or Hamiltonian-path census is performed.
"""
from itertools import combinations
import json

N=14
P=tuple(range(N))
F=(0,13)+tuple(range(1,13))
PAIRS=list(combinations(P,2))
INDEX={pair:i+1 for i,pair in enumerate(PAIRS)}
A_CYCLE=(0,11,9,7,5,3,1)
C_CYCLE=(2,13,12,10,8,6,4)
G={cycle[i]:cycle[(i+1)%7] for cycle in (A_CYCLE,C_CYCLE) for i in range(7)}
PORTS=((3,5,10),(3,7,10),(5,9,10),(7,9,10))
NEW=((0,1,2,13)+tuple(range(3,13)),
     (4,0,1,2,13,3)+tuple(range(5,13)),
     (4,0,1,2,3)+tuple(range(5,14)))

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def edge(a,b):return tuple(sorted((a,b)))
def edges(q):return {edge(a,b)for a,b in zip(q,q[1:])}
def reversal(q):return dict(zip(q,q[::-1]))

def raw():
    out=[[0]*N for _ in P]
    for a,b in PAIRS:out[a][b]=1<<INDEX[a,b];out[b][a]=out[a][b]^1
    return out

def blue_equations(T,q):
    r=reversal(q)
    return [T[a][b]^T[r[b]][r[a]]for a,b in PAIRS]

def K(T,q):
    E=edges(q)
    return [[0 if a==b else T[a][b]^int(edge(a,b)not in E) for b in P]for a in P]

def eliminate(equations):
    basis={}
    for row in equations:
        while row>1:
            p=row.bit_length()-1
            if p in basis:row^=basis[p]
            else:basis[p]=row;break
        if row==1:return None
    return basis

def reduce(row,basis):
    while row>1:
        p=row.bit_length()-1
        if p in basis:row^=basis[p]
        else:
            head=1<<p
            return head^reduce(row^head,basis)
    return row

def source_equations(T,I,chord=None):
    eq=[T[i][i+1]^1 for i in range(N-1)]+[T[13][0]^1]
    eq+=blue_equations(T,P)
    eq+=[T[0][v]^int(v not in I)for v in range(1,13)]
    eq+=[T[1][12]]  # closed middle
    if chord is not None:eq+=[T[2][12]^chord]
    return eq

def reduced(T,basis):
    return [[reduce(T[a][b],basis)for b in P]for a in P]

def old_words(chord):
    W=tuple(range(2,12))
    return ((12,1,13,0)+W,(13,12,1,0)+W,F) if chord else (
        (12,)+W+(1,13,0),(13,12)+W+(1,0),F)

def rows(T):
    return tuple(sum(T[a][b]<<b for b in P)for a in P)

def literal_check(T,I,words):
    need(all(T[a][a]==0 for a in P) and
         all(T[a][b]^T[b][a]==1 for a,b in PAIRS),'literal tournament')
    need(all(T[a][b]for a,b in zip(P,P[1:])) and T[13][0] and
         not any(blue_equations(T,P)),'literal closed blue source')
    need({v for v in range(1,13)if T[v][0]}==set(I),'literal port set')
    need(sum(T[0])==9 and sum(T[13])==4,'literal endpoint degree')
    target=K(T,P);Q=T;xor=set();report=[]
    for q in words:
        need(sorted(q)==list(P) and all(Q[a][b]for a,b in zip(q,q[1:])),
             'literal directed word')
        report.append({'endpoints':[q[0],q[-1]],
                       'endpoint_degrees':[sum(Q[q[0]]),sum(Q[q[-1]])],
                       'blue':not any(blue_equations(Q,q))})
        Q=K(Q,q);xor^=edges(q)
    need(Q==target and xor==edges(P),'literal target and XOR')
    return report

def symmetry_model(I):
    # G-invariance and P-anti reconstruct every arc from B(0,*).
    outgoing={1,13}|set(I)
    def power(v,k):
        for _ in range(k):v=G[v]
        return v
    def arc(a,b):
        if a==b:return 0
        if a in A_CYCLE:
            k=(-A_CYCLE.index(a))%7
            return int(power(b,k)in outgoing)
        if b in A_CYCLE:return 1-arc(b,a)
        return arc(13-b,13-a)
    B=[[arc(a,b)for b in P]for a in P]
    need(all(B[a][b]^B[b][a]==1 for a,b in PAIRS),'carrier tournament')
    need(all(B[a][b]==B[G[a]][G[b]]for a,b in PAIRS),'carrier automorphism')
    need(not any(blue_equations(B,P)) and not any(blue_equations(B,F)),
         'carrier two anti reversals')
    return K(B,P)

def main():
    rF=reversal(F)
    need(G=={v:rF[13-v]for v in P},'printed symmetry permutation')
    X={(2,3),(12,13),(2,13),(3,13)}
    Y={(3,4),(4,5),(0,4),(3,5)}
    need(edges(NEW[0])==edges(P)^X and
         edges(NEW[1])==edges(P)^X^Y and
         edges(NEW[2])==edges(P)^Y,'relocation rectangle')
    rawT=raw()
    # Exhaust every three-port set with c outside I, closed middle.
    # The x,c chord remains free in this independent carrier classification.
    classification=[]
    for I in combinations(range(2,12),3):
        eq=source_equations(rawT,I)
        need(eliminate(eq)is not None,'full cell is nonempty')
        basis=eliminate(eq+blue_equations(K(rawT,P),F))
        if basis is not None:
            T=reduced(rawT,basis)
            need(len(basis)==91,'failed symmetry must fix all arcs')
            need(all(T[a][b]in(0,1)for a in P for b in P),'constant model')
            classification.append((I,T[2][12]))
    need(classification==[(PORTS[0],1),(PORTS[1],1),(PORTS[2],0),(PORTS[3],0)],
         'exact four symmetry exceptions')
    family_rows=[]
    for I in PORTS:
        chord=int(3 in I)
        eq=source_equations(rawT,I,chord)
        basis=eliminate(eq)
        need(basis is not None and 91-len(basis)==28,'full source dimension')
        T=reduced(rawT,basis);old=old_words(chord);current=T;blue_possible=[]
        for q in old:
            need(sorted(q)==list(P) and all(current[a][b]==1 for a,b in zip(q,q[1:])),
                 'old ledger universally directed')
            extra=eliminate(blue_equations(current,q))
            blue_possible.append(extra is not None)
            current=K(current,q)
        need(current==K(T,P),'old symbolic target')
        need(blue_possible==[False,False,True],'old first two universally black')
        bad=eliminate(eq+blue_equations(K(rawT,P),F))
        need(bad is not None and len(bad)==91,'unique actual failed model')
        literal=reduced(rawT,bad)
        model=symmetry_model(I)
        need(literal==model,'affine and symmetry constructions agree')
        old_flags=literal_check(literal,I,old)
        new_flags=literal_check(literal,I,NEW)
        need([r['blue']for r in old_flags]==[False,False,True],'old hostile flags')
        need(not any(r['blue']for r in new_flags),'new three black moves')
        need([r['endpoint_degrees']for r in new_flags]==[[9,5],[8,7],[6,4]],
             'actual-stage degree certificates')
        B=K(literal,P)
        need(sum(B[4])==sum(B[13])==8,'other-end final degree certificate')
        need(literal[3][5]==0,'symmetry-forced new join')
        family_rows.append({'ports':I,'chord_x_to_c':chord,
                            'full_source_dimension':28,'old_blue_final_dimension':0,
                            'source_rows':rows(literal),'old_flags':old_flags,
                            'repair_words':NEW,'repair_flags':new_flags,
                            'final_word_other_end_degrees':[8,8]})
    # Removing the failed-word-symmetry condition invalidates the new ledger:
    # flip one free P-reflection orbit while preserving the original cell.
    hostile=symmetry_model(PORTS[0])
    for a,b in ((3,5),(8,10)):
        hostile[a][b]^=1;hostile[b][a]^=1
    need(not any(blue_equations(hostile,P)) and hostile[3][5]==1,'perturbation source')
    need(all(hostile[a][b]for a,b in zip(NEW[0],NEW[0][1:])),'hostile first word legal')
    S=K(hostile,NEW[0])
    need(S[3][5]==0 and (3,5)in list(zip(NEW[1],NEW[1][1:])),'missing symmetry join')
    old_report=literal_check(hostile,PORTS[0],old_words(1))
    need(not any(r['blue']for r in old_report),'old ledger covers perturbed source')
    print(json.dumps({'status':'PASS','source':'four full affine cells; exact symmetry exception',
      'classification_port_sets':120,'classification':classification,
      'families':family_rows,'common_symmetry_cycles':[A_CYCLE,C_CYCLE],
      'rectangle_X':sorted(X),'rectangle_Y':sorted(Y),
      'missing_symmetry_hostile':{'source_rows':rows(hostile),'flipped_pairs':[[3,5],[8,10]],
                                  'illegal_repair_join':[3,5],'old_three_black':True},
      'scope':'every completion of four order14 cells; three moves, literal labels; no global connectivity claim'},
      sort_keys=True,indent=2))

if __name__=='__main__':main()

