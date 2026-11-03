"""Literal ten-atom bank, profile factors, primitive selector, and one-source repair."""
from pathlib import Path
from itertools import combinations,product
from math import gcd,isqrt,factorial
import json,sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
r=s.Symbol('r');y=s.Symbol('y')
VALUES=(3,49,343,63,5,11,13,35,77,91)
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def rows_from_word(word,ab):
    m=len(word);A=[0]*(m+2)
    for i,x in enumerate(word):
        for j in range(i+1,m):A[i]|=1<<j
        for a in range(2):
            if x>>a&1:A[i]|=1<<(m+a)
            else:A[m+a]|=1<<i
    A[m+(0 if ab else 1)]|=1<<(m+(1 if ab else 0))
    return tuple(A)
def subset_h(A):
    n=len(A);N=1<<n;dp=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,N):
        H[mask]=sum(dp[mask])
        for i,v in enumerate(dp[mask]):
            b=A[i]&~mask
            while v and b:
                t=b&-b;b-=t;dp[mask|t][t.bit_length()-1]+=v
    return H
def strong(A):
    target=(1<<len(A))-1
    for i in range(len(A)):
        seen=1<<i
        while True:
            nxt=seen
            for j in range(len(A)):
                if seen>>j&1:nxt|=A[j]
            if nxt==seen:break
            seen=nxt
        if seen!=target:return False
    return True
def path_poly(A):
    H=subset_h(A);covers=[None]*len(H);covers[0]=[1]
    for mask in range(1,len(H)):
        first=mask&-mask;rest=mask^first;sub=rest;out=[0]*(mask.bit_count()+1)
        while True:
            block=sub|first
            for k,v in enumerate(covers[mask^block]):out[k+1]+=H[block]*v
            if not sub:break
            sub=(sub-1)&rest
        covers[mask]=out
    P=s.expand(sum(v*s.prod(r-j for j in range(k)) for k,v in enumerate(covers[-1])))
    return P,H,covers[-1]
def ordinal(parts):
    n=sum(map(len,parts));out=[];offset=0
    for A in parts:
        later=((1<<n)-1)^((1<<(offset+len(A)))-1)
        out += [(row<<offset)|later for row in A];offset+=len(A)
    return tuple(out)
def converse(A):
    return tuple(((1<<len(A))-1)^row^(1<<i) for i,row in enumerate(A))
def gate(A):
    n=len(A);return tuple(row|(1<<n) for row in A)+(1<<(n+1),(1<<n)-1)
def op(P):
    return s.expand(r*r*P+r*(P.subs(r,r+1)-P.subs(r,r-1)))
def scc(A):
    remaining=set(range(len(A)));components=[]
    # Reachability equivalence, then order components by their uniform arcs.
    reach=[]
    for i in range(len(A)):
        seen={i};stack=[i]
        while stack:
            v=stack.pop()
            for j in range(len(A)):
                if A[v]>>j&1 and j not in seen:seen.add(j);stack.append(j)
        reach.append(seen)
    while remaining:
        i=min(remaining);C=sorted(j for j in remaining if i in reach[j] and j in reach[i])
        components.append(C);remaining.difference_update(C)
    all_components=tuple(components)
    components.sort(key=lambda C:-sum(bool(A[C[0]]>>D[0]&1) for D in all_components if D!=C))
    return components
def extract(A,vertices):
    return tuple(sum(1<<j for j,v in enumerate(vertices) if A[u]>>v&1) for u in vertices)
def main():
    W=json.loads((ROOT/'05-knowledge/results/session156_spectrum_20261103/states/witnesses.json').read_text())
    atoms={};polys={}
    for h in VALUES+(25,):
        A=(126,508,504,240,480,64,384,289,41) if h==343 else rows_from_word(*W[str(h)])
        P,H,covers=path_poly(A)
        need(strong(A) and H[-1]==h and P.subs(r,1)==h,'actual strong atom')
        # Independently compare every polynomial to ordered two-bin subset convolution.
        q=sum(v*H[(len(H)-1)^i] for i,v in enumerate(H))
        need(P.subs(r,2)==q,'profile versus ordered subset convolution')
        atoms[h]={'rows':A,'order':len(A),'H':h,'F2':q,'profile':str(P),'unordered_cover_counts':covers}
        polys[h]=P
    private=[
        y+2,y**3+16*y**2+28*y+4,y**4+28*y**3+148*y**2+152*y+14,
        y**2+14*y+6,y+4,y**2+8*y+2,y**2+8*y+4,
        y**3+12*y**2+20*y+2,y**3+18*y**2+48*y+10,y**3+18*y**2+58*y+14]
    certificates=[]
    for i,j in combinations(range(10),2):
        U,V,g=s.gcdex(private[i],private[j],y);need(g==1,'private factors are coprime')
        denominator=s.ilcm(*[z.q for z in s.Poly(U,y).all_coeffs()+s.Poly(V,y).all_coeffs()])
        U=s.expand(U*denominator);V=s.expand(V*denominator)
        need(s.expand(U*private[i]+V*private[j])==denominator,'integer Bezout certificate')
        certificates.append({'i':i,'j':j,'U':str(U),'V':str(V),'constant':int(denominator)})
    for h,Q in zip(VALUES,private):
        expected=(r**2 if h in (63,5) else r)*Q.subs(y,r*r)
        if h==63:expected*=r*r+2
        need(s.expand(polys[h]-expected)==0,'factorized literal profile')
    loops=[h for h in VALUES if gcd(h,atoms[h]['F2'])>1]
    live=[h for h in VALUES if h not in loops]
    edges=[(a,b) for a,b in combinations(live,2) if gcd(a,atoms[b]['F2'])*gcd(b,atoms[a]['F2'])>1]
    neighbors={a:{b for x,b in edges if x==a}|{x for x,b in edges if b==a} for a in live}
    states={frozenset()};queue=[frozenset()]
    for state in queue:
        for h in live:
            if h in state:continue
            nxt=state|neighbors[h]
            if nxt not in states:states.add(nxt);queue.append(nxt)
    states=sorted(states,key=lambda z:(len(z),sorted(z)));index={F:i for i,F in enumerate(states)}
    dead=len(states);transition=[]
    for F in states:
        transition.append([dead if h not in live or h in F else index[F|neighbors[h]] for h in VALUES])
    transition.append([dead]*10)
    need(len(states)==6,'six live primitive support guards')
    # Reachable product with the parity of the non-5 letter count, plus one dead state.
    phase_states={(0,0)};phase_queue=[(0,0)];phase_table=[]
    for st,b in phase_queue:
        row=[]
        for i,h in enumerate(VALUES):
            dest=transition[st][i];key=None if dest==dead else (dest,b^(h!=5))
            if key is not None and key not in phase_states:phase_states.add(key);phase_queue.append(key)
            row.append(key)
        phase_table.append(row)
    phase_index={key:i for i,key in enumerate(phase_queue)};phase_dead=len(phase_queue)
    phase_table=[[phase_dead if key is None else phase_index[key] for key in row] for row in phase_table]
    phase_table.append([phase_dead]*10)
    accepting={i for i,(st,b) in enumerate(phase_queue) if b==0}
    partition=[int(i in accepting) for i in range(len(phase_table))]
    while True:
        signatures=[(partition[i],tuple(partition[j] for j in row)) for i,row in enumerate(phase_table)]
        labels={v:j for j,v in enumerate(dict.fromkeys(signatures))};new=[labels[v] for v in signatures]
        if len(set(new))==len(set(partition)):break
        partition=new
    need(len(phase_table)==13 and len(set(new))==13,'reachable minimal thirteen-state parity refinement')
    # One-letter distinctions prove minimality of live states; dead is rejection.
    for i,j in combinations(range(6),2):
        need(any((transition[i][a]==dead)!=(transition[j][a]==dead) for a in range(10)),'one-letter state distinction')
    direct_words=0
    for length in range(5):
        for word in product(range(10),repeat=length):
            st=0;h=q=1;parity=0
            for i in word:
                st=transition[st][i];h*=VALUES[i];q*=atoms[VALUES[i]]['F2'];parity^=(VALUES[i]!=5)
            need((st!=dead)==(gcd(h,q)==1),'direct integer primitive predicate')
            if st!=dead:need((parity==0)==((h+q)%3!=0),'direct integer parity predicate')
            direct_words+=1
    dp=[0]*7;dp[0]=1;phase={(0,0):1};count_receipt=[]
    for m in range(13):
        total=sum(dp[:6]);accepted=sum(v for (st,b),v in phase.items() if st!=dead and b==0)
        expected=(2*4**m-2**m+(-3)**m+(-2)**m-(-1)**m)//2
        need(total==2*4**m-2**m and accepted==expected,'all-order word count controls')
        count_receipt.append({'length':m,'primitive_words':total,'primitive_and_not_divisible_by_3_after_G':accepted})
        nxt=[0]*7;np={}
        for st,v in enumerate(dp):
            for dest in transition[st]:nxt[dest]+=v
        for (st,b),v in phase.items():
            for i,h in enumerate(VALUES):
                dest=transition[st][i];key=(dest,b^(h!=5))
                np[key]=np.get(key,0)+v
        dp=nxt;phase=np
    # Explicit value/profile/ordered-source distinction at ten factors.
    pA=polys[3];pC=polys[49];pE=polys[63]
    blocks={h:tuple(atoms[h]['rows']) for h in VALUES}
    baseP=s.expand(pE**2*pC*pA**4);base_h=3969**2;base_q=atoms[63]['F2']**2*atoms[49]['F2']*atoms[3]['F2']**4
    codec=[]
    for m in range(1,4):
        for word in product((0,1),repeat=m):
            component_values=[]
            for bit in word:
                A=[63,63];B=[49,3,3,3,3]
                component_values+=A+B if bit==0 else B+A
            X=ordinal([blocks[h] for h in component_values])
            Z=gate(X);n=len(X)
            need(strong(Z),'strong hidden register')
            need([i for i,row in enumerate(Z) if row.bit_count()==1]==[n],'unique outdegree-one port')
            need([i for i in range(len(Z)) if len(Z)-1-Z[i].bit_count()==1]==[n+1],'unique indegree-one port')
            components=scc(X)
            decoded=[subset_h(extract(X,C))[-1] for C in components]
            need(decoded==component_values,'ordered atomic word recovered')
            # Recovery must use the tournament order, not the original numeric labels.
            relabelled=extract(X,list(reversed(range(n))))
            need([subset_h(extract(relabelled,C))[-1] for C in scc(relabelled)]==component_values,'relabelled SCC recovery')
            Y=ordinal([X,converse(X)]);R=gate(Y)
            N=len(Y);anti=list(range(n,2*n))+list(range(n))+[N+1,N]
            need(all(not (R[i]>>j&1) or (R[anti[j]]>>anti[i]&1) for i in range(N+2) for j in range(N+2) if i!=j),'literal self-converse permutation')
            codec.append({'word':word,'source_order':len(X),'strong_order':len(Z),'self_converse_order':len(R)})
    source25=tuple(atoms[25]['rows']);G25=gate(source25)
    need(subset_h(G25)[-1]==329 and strong(G25),'one new composite source reaches 7*47')
    A=5**6;B=32**6;phase49=[];samples=[]
    for t in range(35):
        H=25*A**t+304*B**t
        need(H%7==0 and (H%49==0)==(t%7==6) and (H%11==0)==(t%5==1) and H%13!=0,'compatible 7/11/13 clocks')
        if H%49==0:phase49.append(t)
        if t<3:samples.append({'t':t,'H':H,'H_over_7':H//7,'order':24*t+8})
    need(all(47%d for d in range(2,isqrt(47)+1)),'prime quotient')
    lifted_phases=[];root=0
    for precision in range(1,9):
        old_mod=7**(precision-1);mod=7**precision
        candidates=[root+j*old_mod for j in range(7)
                    if (25*pow(A,root+j*old_mod,7*mod)+304*pow(B,root+j*old_mod,7*mod))%(7*mod)==0]
        need(len(candidates)==1,'unique compatible 7-adic lift')
        root=candidates[0];lifted_phases.append({'precision':precision,'phase':root,'period':mod})
        for t in (0,1,5,13,101,169):
            big_mod=7**(precision+2)
            ht=(25*pow(A,t,big_mod)+304*pow(B,t,big_mod))%big_mod
            hu=(25*pow(A,t+7**precision,big_mod)+304*pow(B,t+7**precision,big_mod))%big_mod
            need(((hu-ht)//7-5*7**precision)%7**(precision+1)==0,'uniform first-order clock coefficient')
    data={'status':'PROVED_SCOPED','atoms':atoms,'private_factors':list(map(str,private)),
          'private_Bezout_certificates':certificates,'primitive_loops':loops,'primitive_edges':edges,
          'primitive_DFA_states':[sorted(F) for F in states]+['reject'],'transition':transition,
          'primitive_phase_DFA':{'states':phase_queue+['reject'],'transition':phase_table,'accepting':sorted(accepting),'minimal_states':13},
          'direct_integer_word_controls':direct_words,
          'word_counts':count_receipt,'hidden_register_controls':codec,
          'hidden_register_base_H':base_h,'hidden_register_base_F2':base_q,
          'composite25_repair':{'rows':G25,'H':329,'prime_quotient':47,'family':samples,
                                'extra7_phase_mod7':6,'eleven_phase_mod5':1,'thirteen_never_divides':True,'lifted_phases':lifted_phases}}
    (HERE/'bank.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: ten independent profiles; minimal seven-state primitive guard; composite25 repairs 7p')
if __name__=='__main__':main()
