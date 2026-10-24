#!/usr/bin/env python3
"""Independent score-layer audit by normalized presentation cubes, not HP DFS."""
from pathlib import Path
from collections import Counter,defaultdict
from itertools import combinations,permutations
from math import factorial
import hashlib,json
HERE=Path(__file__).resolve().parent

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def cyclic(T):
    return sum(T[i][j]*T[j][k]*T[k][i]+T[i][k]*T[k][j]*T[j][i] for i,j,k in combinations(range(len(T)),3))

def boards(n,basepath=False):
    pairs=[(i,j) for i,j in combinations(range(n),2) if not basepath or j-i>1]
    for word in range(1<<len(pairs)):
        T=[[0]*n for _ in range(n)]
        if basepath:
            for i in range(n-1):T[i][i+1]=1
        for k,(i,j) in enumerate(pairs):T[i][j]=(word>>k)&1;T[j][i]=1-T[i][j]
        yield word,tuple(tuple(row) for row in T)

def C(T):
    n=len(T)
    return tuple(tuple(0 if i==j else T[i][j] if abs(i-j)==1 else 1-T[i][j] for j in range(n)) for i in range(n))

def layer(T):
    hist=Counter(map(sum,T));p=len(T)//2
    return tuple(sorted(hist)) if len(hist)==2 and sorted(hist.values())==[p,p] else None

def isomorphic(T,U):
    n=len(T);degreesT=list(map(sum,T));degreesU=list(map(sum,U))
    for perm in permutations(range(n)):
        if any(degreesT[i]!=degreesU[perm[i]] for i in range(n)):continue
        if all(T[i][j]==U[perm[i]][perm[j]] for i,j in combinations(range(n),2)):return True
    return False

def main():
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
    for name,digest in pins.items():need(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,'source SHA256 mismatch: '+name)
    need(pins['sources/candidate-MANIFEST.json']=='e1a91122cbc42f9de3edfc8e650959a1d79f55158e796c007b8dc34fd31432d9','exact frozen score-layer candidate')
    candidate=json.loads((HERE/'sources/candidate-receipt.json').read_text())
    arithmetic=0;adjacent=0;general=0
    # Independently use the closed quadratic in the odd score gap, not choose sums.
    for p in range(1,501):
        layers=[]
        for k in range(1,2*p,2):
            a=(2*p-1-k)//2;b=a+k
            if 2*a<p-1:continue
            arithmetic+=1
            numerator=p*(4*p*p-1-3*k*k)
            need(numerator%12==0,'integral closed c3 quadratic');c=numerator//12
            need(c==(2*p)*(2*p-1)*(2*p-2)//6-p*(a*(a-1)+b*(b-1))//2,'closed quadratic equals score count')
            lo,hi=c-k-1,c+k-1
            need(lo<=c<=hi,'source c3 lies in its own interval')
            need((0 in {-1,-k-1,k-1})==(k==1),'own layer preservation iff balanced necessary gap')
            layers.append((a,k,c,lo,hi))
        layers.sort()
        if p<3:need(len(layers)==1 and layers[0][1]==1,'small feasible p only balanced')
        for left,right in zip(layers,layers[1:]):
            a,k,c,lo,hi=left;aa,kk,cc,ll,hh=right
            need(aa==a+1 and cc-c==p*(k-1),'adjacent cycle gap')
            need(ll-hi==(p-2)*(k-1)>0,'strict complete interval separation');adjacent+=1
        # A second, nonadjacent identity verifies all layer pairs through p=100.
        if p<=100:
            for left,right in combinations(layers,2):
                a,k,c,lo,hi=left;aa,kk,cc,ll,hh=right;d=aa-a
                need(cc-c==p*d*(k-d),'nonadjacent source cycle identity')
                need(ll-hi==(k-d)*(p*d-2)>0,'nonadjacent interval separation');general+=1
    need(arithmetic==candidate['arithmetic_layers_checked']==62750,'exact arithmetic layer total')
    counts=Counter();shifts=Counter();presentation_counts={};loops=Counter();same_score_black=[];hostile=[];partner_records=[]
    for n in (2,4,6):
        p=n//2;counts[f'n{n}_labelled_sources']=sum(layer(T) is not None for _,T in boards(n))
        filtered=0;blue_filtered=0
        for word,T in boards(n,True):
            scores=layer(T)
            if scores is None:continue
            filtered+=1;a,b=scores;k=b-a;degrees=list(map(sum,T));c=cyclic(T);U=C(T);u=cyclic(U)
            need(C(U)==T,'literal C involution; n2 identity permitted')
            need(all(U[i][i+1]==1 for i in range(n-1)),'literal C retains all path arcs')
            need(list(map(sum,U))==[n-1-degrees[v]+(v==0)-(v==n-1) for v in range(n)],'actual full endpoint degree transport')
            jump=u-c
            need(jump==degrees[0]-degrees[-1]-1,'actual literal cyclic jump')
            need(jump in {-1,-k-1,k-1},'actual unrestricted necessary triple list')
            blue=all(T[i][j]!=T[n-1-i][n-1-j] for i,j in combinations(range(n),2))
            if blue:
                blue_filtered+=1
                need(degrees[0]+degrees[-1]==n-1 and degrees[0]!=degrees[-1],'blue complementary endpoint scores')
                need(jump in {-k-1,k-1},'blue necessary two-value list')
            elif degrees[0]==degrees[-1]:same_score_black.append([n,word,jump])
            other=layer(U)
            if other is not None:need(other==scores and k==1 and degrees[0]==b and degrees[-1]==a,'actual layer-union adjacency only balanced high-to-low')
            converse=tuple(tuple(0 if i==j else 1-T[i][j] for j in range(n)) for i in range(n))
            loop=isomorphic(T,U) or isomorphic(converse,U)
            if loop:need(degrees[0]==degrees[-1]+1,'actual full merged loop leg law');loops[n]+=1
            if k>=3:need(not loop,'unbalanced actual merged layer loopfree')
            if n==6 and a==1:
                need(c==2 and c-k-1==-2 and u>=0 and u!=-2,'unattained negative formal value hostile');hostile.append([n,word,c,u])
            shifts[n,a,jump,blue]+=factorial(n)
            partner_records.append({'n':n,'normalized_word':word,'a':a,'b':b,'source_c3':c,'target_c3':u,'blue':blue,'self_or_converse_partner':loop})
        presentation_counts[n]={'normalized_two_score_presentations':filtered,'blue_presentations':blue_filtered,'labelled_path_pairs_per_presentation':factorial(n)}
        counts[f'n{n}_actual_ports']=filtered*factorial(n);counts[f'n{n}_blue_ports']=blue_filtered*factorial(n)
    histogram=[[list(key),value] for key,value in sorted(shifts.items())]
    need(dict(counts)==candidate['actual_universes'],'complete independent labelled-source and weighted actual-port universes')
    need(histogram==candidate['actual_shift_histogram'],'complete actual endpoint shift/blue histogram')
    need(hostile and same_score_black and all(row[2]==-1 for row in same_score_black),'unattained-value and black-third-value hostiles')
    ten=[]
    for a in (3,4):
        k=9-2*a;c=5*(4*25-1-3*k*k)//12;ten.append({'a':a,'b':9-a,'source_c3':c,'possible_blue_target_c3':[c-k-1,c+k-1]})
    need(ten==candidate['order_ten_layers'],'exact order-ten application')
    print(json.dumps({'status':'PASS','scope':'all-size proof reviewed separately; independent p<=500 arithmetic and complete normalized presentation cubes for labelled two-score actual HP ports at n=2,4,6','source_pins':pins,'arithmetic_layers_checked':arithmetic,'adjacent_interval_pairs_checked':adjacent,'all_layer_pairs_through_p100_checked':general,'actual_universes':dict(counts),'actual_shift_histogram':histogram,'presentation_counts':presentation_counts,'normalized_self_or_converse_partner_presentations':sorted(loops.items()),'order_ten_layers':ten,'negative_formal_value_hostiles':hostile,'same_score_black_counterexamples_to_two_value_list':same_score_black,'actual_port_representatives':partner_records},indent=2,sort_keys=True))
if __name__=='__main__':main()
