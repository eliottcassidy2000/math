#!/usr/bin/env python3
"""Verify frozen exact witnesses, integral quotients, exhaustive obstruction, and code.
Python standard library only. All checks remain active under python -O.
Usage: python3 path_charge_code_obstruction_check.py [certificate.json]
"""
import itertools,json,math,sys
from collections import Counter
from pathlib import Path
CHECKS=0

def require(condition,label):
    global CHECKS
    CHECKS+=1
    if not condition:raise ArithmeticError(label)

def determinant(matrix):
    a=[row[:] for row in matrix];n=len(a);previous=1;sign=1
    if not n:return 1
    for k in range(n-1):
        if a[k][k]==0:
            r=next((i for i in range(k+1,n)if a[i][k]),None)
            if r is None:return 0
            a[k],a[r]=a[r],a[k];sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=a[i][j]*pivot-a[i][k]*a[k][j]
                require(numerator%previous==0,'Bareiss exact division')
                a[i][j]=numerator//previous
        for i in range(k+1,n):a[i][k]=0
        previous=pivot
    return sign*a[-1][-1]

def multiplication(p,S):
    n=p-1;M=[[0]*n for _ in range(n)]
    for j in range(n):
        for s in [0]+list(S):
            e=(j+s)%p
            if e==p-1:
                for i in range(n):M[i][j]-=1
            else:M[e][j]+=1
    return M

def prime_divisors(n):
    result=[];d=2
    while d*d<=n:
        if n%d==0:
            result.append(d)
            while n%d==0:n//=d
        d+=1
    if n>1:result.append(n)
    return result

def order(a,p):
    out=1;r=a%p
    while r!=1:out+=1;r=r*a%p
    return out

def poly_remainder(a,f,q):
    a=[x%q for x in a]
    for i in range(len(a)-1,len(f)-2,-1):
        v=a[i]
        for j,x in enumerate(f):a[i-len(f)+1+j]=(a[i-len(f)+1+j]-v*x)%q
    return a[:len(f)-1]

def valid_path(P,p,S):
    return len(P)==p and P[0]==0 and set(P)==set(range(p)) and all((v-u)%p in S for u,v in zip(P,P[1:]))

def verify_ring(row,alphabet):
    p=row['p'];S=set(row['S']);n=p-1;H=row['HNF'];U=row['unimodular_change'];d=row['norm']
    require(len(S)==n//2 and S.isdisjoint({p-s for s in S}),'skew connection set')
    M=multiplication(p,S)
    require(all(sum(M[i][k]*U[k][j]for k in range(n))==H[i][j]for i in range(n)for j in range(n)),'M U=H')
    require(abs(determinant(U))==1,'unimodular change of integral lattice basis')
    require(all(H[i][j]==0 for i in range(n)for j in range(i))and all(H[i][i]>0 for i in range(n)),'upper triangular H positive diagonal')
    require(math.prod(H[i][i]for i in range(n))==d==determinant(M),'exact quotient index')
    K=[[0 if i==j else(1 if (j-i)%p in S else -1)for j in range(p)]for i in range(p)]
    require(determinant([[K[i][j]+int(i==j)for j in range(p)]for i in range(p)])==2**n*d,'independent skew determinant/norm')
    active=[i for i in range(n)if H[i][i]>1];mods=[H[i][i]for i in active]
    require(active==row['active_coordinates']and mods==row['moduli'],'quotient coordinates')
    require(all(H[i][j]==0 for i in active for j in active if i!=j),'active additive diagonal form')
    def reduce(v):
        v=list(v)
        for j in range(n-1,-1,-1):
            q=v[j]//H[j][j]
            for i in range(j+1):v[i]-=q*H[i][j]
        return tuple(v[i]for i in active)
    powers=[reduce([int(i==e)for i in range(n)]if e<n else[-1]*n)for e in range(p)]
    require([list(v)for v in powers]==row['zeta_coordinates'],'all character powers')
    def index(c):
        out=0;stride=1
        for x,q in zip(c,mods):out+=(x%q)*stride;stride*=q
        return out
    def charge(P):return index([sum(i*powers[v][k]for i,v in enumerate(P))for k in range(len(mods))])
    witnessed=0
    for c,w in enumerate(row['witnesses']):
        if w is None:continue
        P=[alphabet.index(v)for v in w]
        require(valid_path(P,p,S),'explicit Hamiltonian witness')
        require(charge(P)==c,'explicit quotient residue witness');witnessed+=1
    require(witnessed==d if row['surjective'] else witnessed<d,'surjectivity certificate scope')
    for ell in prime_divisors(d):
        require(ell!=p and order(ell,p)%2==1,'odd-order norm-prime restriction')
    require(d%p==1,'norm congruence')
    summary={'p':p,'S':sorted(S),'norm':d,'additive_moduli':mods,'witnessed_residues':witnessed,'surjective':row['surjective']}
    if 'exhaustive_counts'in row:
        counts=[0]*d;repaired=[0]*d;nodes=0
        weights=[[tuple(i*v%q for v,q in zip(powers[x],mods))for x in range(p)]for i in range(p)]
        def visit(P,used,C):
            nonlocal nodes
            nodes+=1;i=len(P)
            if i==p:
                counts[index(C)]+=1
                repaired[index([x+y for x,y in zip(C,powers[P[-1]])])]+=1
                return
            for v in range(p):
                if not used>>v&1 and (v-P[-1])%p in S:
                    visit(P+[v],used|(1<<v),tuple((x+y)%q for x,y,q in zip(C,weights[i][v],mods)))
        visit([0],1,(0,)*len(mods))
        require(counts==row['exhaustive_counts'],'full Hamiltonian charge histogram')
        require(repaired==row['exhaustive_repaired_counts'],'full endpoint-repair histogram')
        require(nodes==row['exhaustive_prefix_nodes']and sum(counts)==row['exhaustive_path_count'],'exhaustive traversal sizes')
        require([i for i,c in enumerate(counts)if not c]==row['missing'],'exact missing residues')
        require([i for i,c in enumerate(repaired)if not c]==row['repaired_missing'],'exact repaired missing residues')
        summary.update({'normalized_paths':sum(counts),'prefix_nodes':nodes,'missing':row['missing'],'endpoint_repair_missing':row['repaired_missing'],'fixed_field_fibres':counts[:3],'endpoint_fixed_field_fibres':repaired[:3]})
    return summary

def verify_golay(code):
    f=code['generator_ascending'];h=code['guava_generator_ascending']
    require(f==[2,2,1,2,0,1]and f==[2*x%3 for x in reversed(h)],'reciprocal documented Golay generator')
    require(poly_remainder([1]*11,f,3)==[0]*5,'f divides Phi_11')
    eta=[int(i in {0,2,6,7,8,10})for i in range(11)]
    require(poly_remainder(eta,f,3)==[0]*5,'eta vanishes in the F243 presentation')
    require(order(3,11)==5,'irreducibility via degree of cyclotomic factors')
    rows=[]
    for j in range(6):
        v=[0]*11;v[j:j+6]=f;rows.append(v+[(-sum(v))%3])
    gram=[[sum(a*b for a,b in zip(u,v))%3 for v in rows]for u in rows]
    require(gram==code['extended_generator_gram']==[[0]*6 for _ in range(6)],'extended self-orthogonal Gram matrix')
    words=[]
    for a in itertools.product(range(3),repeat=6):
        v=[sum(a[j]*rows[j][i]for j in range(6))%3 for i in range(11)]
        require((sum(x*x for x in v)+sum(v)**2)%3==0,'code quadratic identity')
        words.append(''.join(map(str,v)))
    require(words==code['words_in_lexicographic_multiplier_order']and len(set(words))==729,'complete 729-word certificate')
    weights=Counter(sum(x!='0'for x in w)for w in words)
    comps=Counter(tuple(w.count(str(i))for i in range(3))for w in words)
    require({str(k):v for k,v in sorted(weights.items())}==code['weight_enumerator'],'weight enumerator')
    require([[list(k),v]for k,v in sorted(comps.items())]==code['composition_enumerator'],'complete composition enumerator')
    require(min(k for k in weights if k)>0 and min(k for k in weights if k)==5,'minimum distance five')
    require([a for a in range(3)if(1+a*a+(1-a)**2)%3==0]==[2],'original fixed-field obstruction')
    require([a for a in range(3)if(1+a*a+(2-a)**2)%3==0]==[1],'endpoint fixed-field obstruction')
    return {'parameters':[11,6,5],'words':len(words),'weight_enumerator':code['weight_enumerator'],'composition_classes':len(comps),'original_allowed_fixed_field_charge':[2],'endpoint_allowed_fixed_field_charge':[1]}

def verify_component(c):
    p=c['p'];S=set(c['S']);ell=c['norm'];z=c['zeta'];powers=[pow(z,i,ell)for i in range(p)]
    require(p==17 and ell==103 and order(z,ell)==p,'degree-one quotient character')
    require(determinant(multiplication(p,S))==ell and(1+sum(powers[s]for s in S))%ell==0,'exact F103 quotient')
    def charge(P):return sum(i*powers[v]for i,v in enumerate(P))%ell
    paths=[];cycle_tests=0
    for index,(parent,i,direction)in enumerate(c['nodes_parent_block_direction']):
        if index==0:
            require([parent,i,direction]==[-1,-1,0],'root marker');P=tuple(c['root'])
        else:
            require(0<=parent<index and 1<=i<p-2 and direction in (1,2),'rotation certificate indices')
            Q=paths[parent];a,b,beta=Q[i:i+3]
            require((a-beta)%p in S,'directed triangle')
            rot=(b,beta,a)if direction==1 else(beta,a,b)
            P=Q[:i]+rot+Q[i+3:]
            increment=(2*powers[a]-powers[b]-powers[beta])if direction==1 else(powers[a]+powers[b]-2*powers[beta])
            require((charge(P)-charge(Q)-increment)%ell==0,'triangle charge increment')
        require(valid_path(P,p,S),'rotation retains boundary arcs and Hamiltonian order');paths.append(P)
        if(-P[-1])%p in S:
            raw=P[1:]+P[:1];normalized=tuple((v-P[1])%p for v in raw)
            require(valid_path(normalized,p,S),'cycle-start rotation path')
            require(charge(raw)==(charge(P)+p)%ell,'cycle-start PLUS p formula')
            require(charge(normalized)==pow(powers[P[1]],-1,ell)*(charge(P)+p)%ell,'normalized cycle-start formula');cycle_tests+=1
    require(len(set(paths))==len(paths),'distinct discovered paths')
    for residue,index in enumerate(c['residue_witness_indices']):
        require(0<=index<len(paths)and charge(paths[index])==residue,'full triangle-component charge coverage')
    P=tuple(c['isolated_path']);require(valid_path(P,p,S)and charge(P)==c['isolated_charge'],'isolated path witness')
    valid_moves=0
    for i in range(p-2):
        a,b,beta=P[i:i+3]
        require((a-beta)%p in S,'every isolated-path consecutive triple is cyclic')
        for rot in ((b,beta,a),(beta,a,b)):
            Q=P[:i]+rot+P[i+3:]
            if all((v-u)%p in S for u,v in zip(Q,Q[1:])):valid_moves+=1
    require(valid_moves==0,'boundary obstruction to every triangle rotation')
    require(tuple((v-P[1])%p for v in P[1:]+P[:1])==P,'cycle-start normalization fixes isolated constant-step path')
    return {'p':p,'norm':ell,'zeta':z,'component_states':len(paths),'witnessed_residues':ell,'verified_cycle_start_moves':cycle_tests,'isolated_charge':charge(P),'isolated_valid_triangle_rotations':valid_moves}

def main():
    if len(sys.argv)>1:path=Path(sys.argv[1])
    else:
        name='path_charge_surjectivity_certificates_20261014.json';here=Path(__file__).resolve()
        path=here.with_name(name)
        if not path.exists():path=here.parents[2]/'05-knowledge'/'results'/name
    data=json.loads(path.read_text());summaries=[]
    for row in data['representatives']:summaries.append(verify_ring(row,data['alphabet']))
    covered=0
    for p in data['scope_primes']:
        expected={frozenset(j+1 if mask>>j&1 else p-j-1 for j in range((p-1)//2))for mask in range(1<<((p-1)//2))}
        observed=set()
        for row in data['representatives']:
            if row['p']!=p:continue
            orbit={frozenset(a*s%p for s in row['S'])for a in range(1,p)}
            require(len(orbit)==row['orbit_size']and observed.isdisjoint(orbit),'disjoint multiplier orbit')
            observed|=orbit
        require(observed==expected,'all labelled prime-cyclic connection sets covered');covered+=len(observed)
    golay=verify_golay(data['golay_code']);component=verify_component(data['p17_triangle_component'])
    require(order(2,1093)==364 and order(2,3511)==1755,'named Wieferich-prime order parity')
    require(3511%4==3 and (3511+1)//4==878 and 878%4==2,'Paley3511 exact two-adic norm exponent')
    result={'status':'PASS','scope_primes':data['scope_primes'],'labelled_connection_sets':covered,'multiplier_representatives':len(summaries),'positive_representatives':sum(r['surjective']for r in summaries),'small_cases':summaries,'golay':golay,'p17_component':component,'wieferich_order_controls':{'ord_1093_2':364,'ord_3511_2':1755,'v2_paley_3511_disc':1755},'active_checks':CHECKS}
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
