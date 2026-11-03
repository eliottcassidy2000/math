"""Independent endpoint counts, literal macro replays and actual rank witnesses."""
from itertools import combinations
import json,sympy as s
from effects import HERE,old,junction,need,PREV
from returns import cycle
from singular_return import nr

def endpoints(T,root=None):
    n=len(T);D=[[0]*n for _ in range(1<<n)]
    for v in range(n):
        if root is None or root==v:D[1<<v][v]=1
    for mask in range(1,1<<n):
        for v,c in enumerate(D[mask]):
            if not c:continue
            todo=T[v]&~mask
            while todo:
                bit=todo&-todo;todo-=bit;D[mask|bit][bit.bit_length()-1]+=c
    return D

def values(T,r):
    D=endpoints(T);A=endpoints(T,r);hs=list(map(sum,D));hs[0]=1;ss=list(map(sum,A));full=(1<<len(T))-1;rest=full^(1<<r)
    F=sum(D[m][r]*hs[full^m] for m in range(1,full+1))
    S2=sum(ss[m]*hs[full^m] for m in range(1,full+1))
    P2=sum(hs[m]*hs[full^m] for m in range(full+1))
    D2=sum(hs[m]*hs[rest^m] for m in range(rest+1) if m&rest==m)
    return (hs[-1],ss[-1],D[-1][r],S2,P2,D2),(hs[-1],ss[-1],hs[rest],F)

def main():
    ops={k:s.Matrix(v) for k,v in json.loads((HERE/'effects.json').read_text())['natural_matrices'].items()}
    witnesses=[];basis=[];examined=0
    for T in junction.all_tournaments(5):
        if not junction.strong(T):continue
        for r in range(5):
            U,q,ports,_=junction.literal_seal(T,r,0,(r,next(v for v in range(5) if v!=r)))
            v,z=values(U,q);examined+=1
            need(z==(v[0],v[1],v[1],v[0]),'actual balanced seal')
            if s.Matrix(basis+[v]).rank()>len(basis):
                basis.append(v);witnesses.append({'rows':U,'root':q,'ports':ports,'values':v})
            if len(basis)==6:break
        if len(basis)==6:break
    need(len(basis)==6,'six actual strong order-seven witnesses')
    det=s.Matrix(basis).det();literal=[]
    for source in witnesses[:2]:
        T=tuple(source['rows']);root=source['root'];ports=tuple(source['ports']);v=source['values'];b,t=v[:2]
        for letter in ('3','5','N'):
            out=nr(T,root,ports,b,t) if letter=='N' else cycle(T,root,ports,b,t,int(letter))
            U,q,pp,bp,tp,receipt=out;actual,z=values(U,q)
            expected=tuple(s.Matrix([v])*ops[letter])
            need(actual==expected and z==(bp,tp,tp,bp),'independent literal six-coordinate macro law')
            literal.append({'input_values':v,'letter':letter,'order':len(U),'output_values':actual,'root':q,'ports':pp,'rows':U})
    # Minimum starting-path multiplicity is monotone in a strong ambient tournament.
    monotone=0
    for T in junction.all_tournaments(5):
        if not junction.strong(T):continue
        full=31;rev=tuple(full^row^(1<<i) for i,row in enumerate(T));D=endpoints(rev)
        ambient=min(D[-1])
        for mask in range(1,32):
            value=min(D[mask][r] for r in range(5) if mask>>r&1)
            need(ambient>=value,'all induced-subset minimum-start control');monotone+=1
    out={'status':'FINITE_EXACT_INDEPENDENT_CONTROLS','examined_order7_seals':examined,'rank_witnesses':witnesses,
         'rank_determinant':int(det),'literal_macro_checks':literal,'strong_order5_induced_start_checks':monotone}
    (HERE/'controls.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: six actual order7 witnesses, six independently counted macro outputs, and',monotone,'start-minimum controls')
if __name__=='__main__':main()
