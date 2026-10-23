"""Zero-slot marked Fibonacci grammar and exact finite-group block clocks."""
from pathlib import Path
from itertools import product,permutations
from collections import Counter
import argparse,json,math
HERE=Path(__file__).resolve().parent
CHECKS=0
def need(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok:raise ValueError(msg)
H=[1,0,1,1]
while H[-1]<200000:H.append(H[-1]+H[-2])
W=['A','','B','C','CA','CAB']
for k in range(6,len(H)):W.append(W[-1]+W[-2])
U=W[-1]
def zeck(n):
    ans=[]
    for i in range(len(H)-1,2,-1):
        if H[i]<=n:ans.append(i);n-=H[i]
    need(n==0,'greedy incomplete');return ans

def legal(indices,total):
    need(len(indices)==len(set(indices)) and 1 not in indices,'duplicate or selected zero')
    need(all(abs(i-j)>1 for i,j in zip(indices,indices[1:])),'adjacent indices')
    need(sum(H[i] for i in indices)==total,'wrong visible length')

def proper(k):
    need(k>=4,'unit is terminal')
    if k%2:parts=list(range(k-1,1,-2))
    else:parts=list(range(k-1,2,-2))+[0]
    return parts

def expand(k,start,changed):
    if k in(0,2,3):return {0:'A',2:'B',3:'C'}[k]
    parts=proper(k)
    # Each C has a sibling A at the end of an even-index expansion node.
    # Its occurrence address, not just the value1, decides the alternative.
    if k%2==0:
        pos=start+sum(H[i] for i in parts[:-2])
        if pos in changed:parts[-2]=2
    legal(parts,H[k]);text='';offset=start
    for i in parts:text+=expand(i,offset,changed);offset+=H[i]
    return text

def row(n,changed):
    parts=zeck(n)
    if parts and parts[-1]==3 and n in changed:parts[-1]=2
    legal(parts,n);text='';start=1
    for k in parts:text+=expand(k,start,changed);start+=H[k]
    return text

def all_small_reps(N):
    indices=[i for i,h in enumerate(H) if i!=1 and h<=N]
    out={n:set() for n in range(N+1)}
    def visit(t,chosen,total):
        if t==len(indices):out[total].add(tuple(sorted(chosen)));return
        i=indices[t];visit(t+1,chosen,total)
        if total+H[i]<=N and (not chosen or i-chosen[-1]>1):visit(t+1,chosen+[i],total+H[i])
    visit(0,[],0);return out

def grammar_checks():
    supplied='CABCACABCABCACAB'
    need(U.startswith(supplied),'new supplied prefix differs')
    reps=all_small_reps(2000);sizes=Counter()
    for n in range(2001):
        expected={tuple(sorted(zeck(n)))}
        if n>=1:
            z=zeck(n-1);expected.add(tuple(sorted([0]+z)))
            if 3 not in z:expected.add(tuple(sorted([2]+z)))
        if n>=2:
            z=zeck(n-2)
            if 3 not in z:expected.add(tuple(sorted([0,2]+z)))
        need(expected==reps[n],'zero-slot fibre formula');sizes[len(expected)]+=1
    need((0,2) in reps[2],'new zero slot should permit A+B')
    for n in range(2001):
        need(row(n,set())==U[:n],'preferred row mismatch')
    cpos=[i+1 for i,x in enumerate(U[:2000]) if x=='C']
    patterns=[{35},set(cpos),{i for i in cpos if i%5 in(1,2)}]
    for changed in patterns:
        need(all(U[i-1]=='C' for i in changed),'recolouring a non-C')
        target=''.join('B' if i+1 in changed else x for i,x in enumerate(U[:2000]))
        for n in range(2001):need(row(n,changed)==target[:n],'contextual expansion lost shared diagonal')
    short=[i+1 for i,x in enumerate(U[:21]) if x=='C']
    for bits in product((0,1),repeat=len(short)):
        changed={i for i,b in zip(short,bits) if b}
        for n in range(22):
            target=''.join('B' if i+1 in changed else x for i,x in enumerate(U[:n]))
            need(row(n,changed)==target,'short exhaustive recolouring')
    # The zero has a nonzero second register. Every preferred split drops
    # exactly one copy, and each C->B choice drops one more.
    F=[0,1]
    for k in range(2,len(H)):F.append(F[-1]+F[-2])
    nodes=[0]*len(H)
    for k in range(4,len(H)):
        parts=proper(k);nodes[k]=1+sum(nodes[i] for i in parts)
        need(F[k]-sum(F[i] for i in parts)==1,'preferred split zero-charge defect')
        need(nodes[k]==F[k-2],'expansion node count')
    for n in range(2001):
        z=zeck(n);Y=sum(F[i] for i in z);N=sum(nodes[i] for i in z)
        need(Y==2*n-(3*(n+1)-math.isqrt(5*(n+1)**2)-1)//2,'inherited shifted register differs')
        for changed in patterns:
            text=''.join('B' if i+1 in changed else x for i,x in enumerate(U[:n]))
            charge=sum({'A':0,'B':1,'C':2}[x] for x in text)
            need(charge+N+sum(i<=n for i in changed)==Y,'node and recolour debts fail to restore shifted register')
    repaired=row(36,{35})
    need(repaired[:34]==U[:34] and repaired[34:]=='BA','old row35 hostile not repaired')
    return {'zero_never_selected':True,'complete_fibres_through':2000,'fibre_size_histogram_including_zero':sorted(sizes.items()),
            'new_supplied_prefix':supplied,'preferred_prefix_40':U[:40],'old_B35_extension_row36':repaired,
            'zero_charge_accounting_rows':2001,'exhaustive_recolour_subsets_first21':2**len(short),'tested_recolour_patterns_through2000':len(patterns),
            'empty_slot_changes_admissibility':{'new_indices':[0,2],'length':2,'word_descending':'BA'}}

def cycle_periods(G,mul):
    tally=Counter()
    for x,y in product(G,repeat=2):
        a,b=x,y;n=0
        while True:
            a,b=b,mul(a,b);n+=1
            if(a,b)==(x,y):break
            need(n<=len(G)**2,'pair update not a finite permutation')
        tally[n]+=1
    return sorted(tally.items())

def group_checks():
    S=list(permutations(range(3)))
    def compose(p,q):return tuple(p[q[i]] for i in range(3))
    s=cycle_periods(S,compose);need(s==[(1,1),(3,9),(6,18),(8,8)],'S3 periods')
    F=[(a,b) for a in(1,2,4) for b in range(7)]
    A=[[int((j-i)%7 in(1,2,4)) for j in range(7)] for i in range(7)]
    actual={p for p in permutations(range(7)) if all(A[i][j]==A[p[i]][p[j]] for i in range(7) for j in range(i))}
    need(actual=={tuple((a*x+b)%7 for x in range(7)) for a,b in F},'P7 automorphism claim')
    def affine(u,v):return(u[0]*v[0]%7,(u[0]*v[1]+u[1])%7)
    f=cycle_periods(F,affine);need(f==[(1,1),(8,56),(16,48),(56,336)],'P7 group periods')
    # Independent eight-step linear translation matrices, with fixed slope pair.
    matrices=[]
    for a,c in product((1,2,4),repeat=2):
        cols=[]
        for b,d in((1,0),(0,1)):
            x,y=(a,b),(c,d)
            for _ in range(8):x,y=y,affine(x,y)
            need((x[0],y[0])==(a,c),'slopes did not return');cols.append((x[1],y[1]))
        M=[[cols[j][i] for j in range(2)] for i in range(2)]
        if(a,c)==(1,1):need(M==[[6,0],[0,6]],'translation-only matrix')
        else:
            N=[[(M[i][j]-int(i==j))%7 for j in range(2)] for i in range(2)]
            need(any(x for r in N for x in r),'nontrivial unipotent missing')
            need(all(sum(N[i][k]*N[k][j] for k in range(2))%7==0 for i,j in product(range(2),repeat=2)),'square-zero failure')
            need(sum(all(sum(N[i][j]*v[j] for j in range(2))%7==0 for i in range(2)) for v in product(range(7),repeat=2))==7,'fixed line size')
        matrices.append({'slopes':[a,c],'translation_matrix_after8':M})
    # Ordered block products, compared to literal letters for every S3 letter assignment.
    def value(word,letters):
        g=(0,1,2)
        for x in word:g=compose(letters[x],g)
        return g
    for vals in product(S,repeat=3):
        letters=dict(zip('ABC',vals));states={k:value(W[k],letters) for k in(4,5)}
        for k in range(6,15):
            states[k]=compose(states[k-2],states[k-1])
            need(states[k]==value(W[k],letters),'macro loses order or source chart')
    # Counts alone do not preserve the product.
    letters={'A':(1,2,0),'C':(0,2,1)}
    need(value('CA',letters)!=value('AC',letters),'commutative shadow hostile')
    lengths=[1,2]
    for _ in range(1000):lengths.append(sum(lengths[-2:]))
    return {'S3_pair_period_histogram':s,'S3_uniform_period':24,'P7_aut_pair_period_histogram':f,'P7_aut_uniform_period':112,
            'affine_eight_step_matrices':matrices,'S3_letter_assignments_replayed':216,
            'integer_clock_control_bits':lengths[-1].bit_length(),'integer_clock_still_grows':lengths[-1]>lengths[-25]}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    out={'status':'FINITE-EXACT controls for proved scoped rules','grammar':grammar_checks(),'groups':group_checks(),'checks':CHECKS}
    out=json.loads(json.dumps(out));f=HERE/'fibonacci.json'
    if a.write:f.write_text(json.dumps(out,indent=2)+'\n')
    else:need(json.loads(f.read_text())==out,'receipt mismatch')
    print(json.dumps(out,indent=2))
