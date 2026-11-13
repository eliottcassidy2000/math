"""Independent two-module transport, factor inverse, and literal controls."""
from pathlib import Path
from functools import lru_cache
from itertools import product,permutations
from math import gcd,isqrt
import json
from response_probe import inventory,pair,need

HERE=Path(__file__).resolve().parent

def tt(n):return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))

def construct(M,N,after):
    blocks=[tt(1),tt(1),M,tt(1),N,tt(1),tt(1)]
    n=sum(map(len,blocks));starts=[];offset=0;rows=[]
    for B in blocks:
        starts.append(offset)
        tail=((1<<n)-1)^((1<<(offset+len(B)))-1)
        rows.extend((row<<offset)|tail for row in B);offset+=len(B)
    wins={0,1,4,5} if after else {0,3,4}
    apex=sum(((1<<len(blocks[i]))-1)<<starts[i] for i in wins)
    return tuple(row|((1<<n) if not(apex>>i&1) else 0) for i,row in enumerate(rows))+(apex,),starts

def state(T,r):
    H,E=inventory(T);mask=len(H)-1
    R=tuple(mask^(1<<v)^row for v,row in enumerate(T));_,S=inventory(R)
    return (H[-1],S[-1][r],H[mask^(1<<r)],sum(E[m][r]*H[mask^m] for m in range(mask+1)))

def birth(T,r,retain):
    n=len(T);U=list(T);U[r]|=1<<n
    U.append(((1<<n)-1)^(1<<r))
    return tuple(U),r if retain else n

def formula(h,p,j,q):
    A=12*p*q+6*h*q+2*p*j-h*j
    Z=12*p*q+2*h*q+6*p*j-h*j
    s0=6*h*q-h*j;q1=(2*p+h)*(2*q+j)
    return (A,s0,s0,A),(Z,q1-h*j,q1,Z+h*j),(3*Z+2*h*j,Z,Z,3*Z+2*h*j)

def run_coefficients(after):
    Q,_=construct(tt(1),tt(1),after)
    out={}
    for obs in ('H','start','delete'):
        C={}
        for r,s in product((1,2),repeat=2):
            rem=[1]*8;rem[2]=r;rem[4]=s
            if obs=='delete':rem[0]=0
            @lru_cache(None)
            def walk(rem,last):
                if not any(rem):return 1
                total=0
                for v,k in enumerate(rem):
                    if k and (last<0 or Q[last]>>v&1):
                        nxt=list(rem);nxt[v]-=1;total+=walk(tuple(nxt),v)
                return total
            last=-1
            if obs=='start':rem[0]-=1;last=0
            C[f'{r},{s}']=walk(tuple(rem),last)
        out[obs]=C
    return out

def arithmetic_candidates(Y):
    """Complete positive nonempty-module arithmetic inverse, not a graph test."""
    out=[]
    if Y<191 or Y%2==0:return out
    for h in range(1,Y//191+1,2):
        for j in range(1,Y//(191*h)+1,2):
            R=Y+4*h*j
            for d in range(1,isqrt(R)+1):
                if R%d:continue
                for x,y in {(d,R//d),(R//d,d)}:
                    if x<13*h or y<15*j or (x-h)%12 or (y-3*j)%12:continue
                    p=(x-h)//6;q=(y-3*j)//6
                    nM=min(p//h-1,p.bit_length()-1)
                    nN=min(q//j-1,q.bit_length()-1)
                    out.append(dict(h=h,p=p,j=j,q=q,module_order_bounds=[nM,nN]))
    return sorted(out,key=lambda x:(x['h'],x['j'],x['p'],x['q']))

def primitive_guard(h,p,j,q):
    return gcd(h,3*p*(2*q+j))==1 and gcd(j,q*(6*p+h))==1

def direct_hp(T):
    return sum(all(T[a]>>b&1 for a,b in zip(word,word[1:])) for word in permutations(range(len(T))))

def main():
    coeff=[run_coefficients(False),run_coefficients(True)]
    expected=[{'H':[63,30,26,12],'start':[11,6,0,0],'delete':[11,6,0,0]},
              {'H':[63,26,30,12],'start':[24,10,10,4],'delete':[25,10,10,4]}]
    for i,C in enumerate(coeff):
        need({obs:list(values.values()) for obs,values in C.items()}==expected[i],'exact run table')
    modules=[tt(1),tt(2),(2,4,1),(2,4,9,3)]
    literal=[]
    for M,N in product(modules,repeat=2):
        h,p=pair(M);j,q=pair(N);expected_states=formula(h,p,j,q)
        T,starts=construct(M,N,False);U,_=construct(M,N,True)
        need(state(T,0)==expected_states[0] and state(U,0)==expected_states[1],'full states')
        V,r=birth(U,0,False);V,r=birth(V,r,True)
        need(state(V,r)==expected_states[2],'literal NR output')
        changed=[(i,k) for i in range(len(T)) for k in range(i+1,len(T)) if ((T[i]^U[i])>>k)&1]
        need(changed==[(starts[i],len(T)-1) for i in (1,3,5)],'three arc support')
        need(T[0]==U[0],'retained root star')
        Y,Z=expected_states[2][:2]
        need(Y==(6*p+h)*(6*q+3*j)-4*h*j,'factored output')
        need(Y==3*expected_states[0][0]+2*(6*(p*j-h*q)+h*j),'compensation law')
        need((gcd(Y,Z)==1)==primitive_guard(h,p,j,q),'primitive iff')
        literal.append(dict(M=M,N=N,h=h,p=p,j=j,q=q,states=expected_states,changed_arcs=changed))
    witness=json.loads((HERE/'response_probe.json').read_text())['nontrivial'][0]
    for name in ('M','N'):
        B=witness[name];need(direct_hp(B['rows'])==B['h'],'independent eight-vertex permutation H')
    M,N=witness['M'],witness['N'];h,p=M['h'],M['p'];j,q=N['h'],N['p']
    states=formula(h,p,j,q)
    need(states[2][0]==3*states[0][0]-2 and primitive_guard(h,p,j,q),'nontrivial negative digit')
    large_before,_=construct(M['rows'],N['rows'],False)
    large_after,_=construct(M['rows'],N['rows'],True)
    large_output,large_root=birth(large_after,0,False)
    large_output,large_root=birth(large_output,large_root,True)
    for T in (large_before,large_after,large_output):
        need(all(not(row>>i&1) and row<(1<<len(T)) for i,row in enumerate(T)) and
             all(((T[i]>>k)&1)+((T[k]>>i)&1)==1 for i in range(len(T)) for k in range(i+1,len(T))),
             'large witness tournament validity without large path DP')
    alignment=(0,1,2,3,4,7,5,6)
    changed=[(i,k) for i in range(8) for k in range(i+1,8)
             if ((M['rows'][i]>>k)&1)!=((N['rows'][alignment[i]]>>alignment[k])&1)]
    need(changed==[(0,1),(1,2),(5,7)],'three-arc relation of the two compensated modules')
    candidate_targets=(7,21,191,197,347,371,881,5299)
    inverse={str(Y):arithmetic_candidates(Y) for Y in candidate_targets}
    need(any(c['h']==c['j']==1 and c['p']==2 and c['q']==4 for c in inverse['347']),'prime inverse positive')
    need(any(c['h']==c['j']==1 and c['p']==4 and c['q']==2 for c in inverse['371']),'seven-prime inverse positive')
    need(not inverse['881'],'complete arithmetic hostile881')
    out=dict(status='PASS',scope='Universal exact quotient-run tables;16 literal ordered module pairs through order16 after NR; module8 permutation check for the compensated negative digit; complete arithmetic inverse for eight named targets',coefficients=coeff,literal=literal,negative_digit=dict(modules=witness,states=states,before_rows=large_before,after_rows=large_after,output_rows=large_output,output_root=large_root,module_alignment=alignment,module_changed_pairs=changed),inverse=inverse)
    (HERE/'polarized_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASS:16 literal module pairs, exact run tables, negative digit and factor inverse')
    print('Arithmetic candidate counts:',{Y:len(c) for Y,c in inverse.items()})

if __name__=='__main__':main()
