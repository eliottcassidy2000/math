"""Exact falsification tests and operation-level arithmetic connections."""
from collections import Counter
from itertools import combinations, permutations, product
from pathlib import Path
import json
from experiment import need, rows, ham

HERE=Path(__file__).resolve().parent


def det(a):
    n=len(a);answer=0
    for p in permutations(range(n)):
        value=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        for i in range(n):value*=a[i][p[i]]
        answer+=value
    return answer


def characteristic(out):
    n=len(out);a=[[int(out[i]>>j&1) for j in range(n)] for i in range(n)]
    coeff=[1]
    for k in range(1,n+1):
        coeff.append((-1)**k*sum(det([[a[i][j] for j in idx] for i in idx]) for idx in combinations(range(n),k)))
    # Independent Newton trace recurrence, with exact divisibility.
    power=[[int(i==j) for j in range(n)] for i in range(n)];traces=[0]
    check=[1]
    for k in range(1,n+1):
        power=[[sum(power[i][s]*a[s][j] for s in range(n)) for j in range(n)] for i in range(n)]
        traces.append(sum(power[i][i] for i in range(n)))
        num=-sum(check[k-i]*traces[i] for i in range(1,k+1))
        need(num%k==0,'Newton integer coefficients');check.append(num//k)
    need(coeff==check,'principal-minor and trace polynomials agree')
    return coeff


def projective():
    points=[]
    for v in product(range(3),repeat=3):
        if any(v) and next(x for x in v if x)==1:points.append(v)
    supplied=[(1,0,0),(0,1,0),(0,0,1),(1,1,1),(1,2,1)]
    collinear=[list(c) for c in combinations(range(5),3) if det([supplied[i] for i in c])%3==0]
    anchors=[]
    for anchor in points:
        if anchor in supplied:continue
        ties=[list(c) for c in combinations(range(5),2) if det([supplied[c[0]],supplied[c[1]],anchor])%3==0]
        need(bool(ties),'every proposed anchor has undefined arcs')
        anchors.append({'anchor':anchor,'ties':ties})
    caps=Counter()
    for k in range(3,6):
        caps[k]=sum(all(det([points[i] for i in triple])%3!=0 for triple in combinations(subset,3))
                    for subset in combinations(range(13),k))
    need(caps[4]>0 and caps[5]==0,'no projective pentad without collinear triples')
    return {'points':13,'supplied_collinear_triples_zero_indexed':collinear,'anchors':anchors,'cap_counts':dict(caps)}


def transpose_class_census():
    n=5;edge=list(combinations(range(n),2));top=(1<<10)-1;classes={}
    for mask in range(1024):
        g=rows(5,mask)
        orbit={sum(int(g[p[i]]>>p[j]&1)<<k for k,(i,j) in enumerate(edge)) for p in permutations(range(n))}
        key=min(orbit)
        if key not in classes:classes[key]={'size':len(orbit),'self_converse':(top^mask) in orbit}
    need(len(classes)==12 and sum(v['self_converse'] for v in classes.values())==8,'actual five-vertex self-converse classes')
    return {'unlabelled_classes':12,'self_converse_classes':8,'classes':classes}


def fixed_reversal():
    results=[]
    for n in range(1,7):
        edge=list(combinations(range(n),2));sigma=list(range(n))
        for i in range(0,n-1,2):sigma[i],sigma[i+1]=i+1,i
        found=0
        for mask in range(1<<len(edge)):
            g=rows(n,mask)
            if all((g[sigma[i]]>>sigma[j]&1)!=(g[i]>>j&1) for i,j in edge):found+=1
        need(found==1<<(n*n//4),'fixed reversing involution has exact Boolean dimension')
        results.append({'n':n,'reversing_involution':sigma,'orientations':found})
    return results


def clock_magic():
    cycle=[4,8,3,7,2,6,10,5,9]
    need(all(cycle[(i+1)%9]==2+(cycle[i]-2+4)%9 for i in range(9)),'literal translated modulo-nine clock')
    def period(values):
        return next(k for k in (1,3,9) if all(values[i]==values[(i+k)%9] for i in range(9)))
    need(period([x%3 for x in cycle])==3 and period([x%2 for x in cycle])==9,
         'future parity recovers the whole clock; modulo-three output does not')
    magic=[[5,10,3],[4,6,8],[9,2,7]]
    need(sorted(sum(magic,[]))==sorted(cycle),'same nine labels, explicit arrangement')
    lines=[r for r in magic]+[[magic[i][j] for i in range(3)] for j in range(3)]+[
        [magic[i][i] for i in range(3)],[magic[i][2-i] for i in range(3)]]
    need(all(sum(line)==18 for line in lines),'ordinary magic arrangement')
    totals=[sum(x*x for x in line) for line in lines]
    need(len(set(totals))>1,'squaring does not preserve magic')
    # The lambda notation recovered in the archive is the Lucas cube Lambda_5.
    words=[w for w in range(32) if all(not(w>>i&1 and w>>((i+1)%5)&1) for i in range(5))]
    edges=[(a,b) for a,b in combinations(words,2) if (a^b).bit_count()==1]
    labels={w:sum((w>>i&1)*4**i for i in range(5))%11 for w in words}
    need(len(words)==11 and len(edges)==15 and len(set(labels.values()))==11,'inherited Lucas-five reader')
    return {'cycle':cycle,'modulus':9,'step':4,'magic_array':magic,'magic_line_sum':18,
            'mod3_observer_states':3,'parity_observer_states':9,
            'squared_line_sums':totals,'Lucas5_words':words,'Lucas5_edges':len(edges),'Lucas5_labels':labels}


def main():
    mats=[(6,12,24,17,3),(14,20,24,18,1)]
    result={'matrix_audits':[{'rows':g,'characteristic_descending':characteristic(g),'H':ham(g)} for g in mats],
            'projective':projective(),'self_converse_n5':transpose_class_census(),
            'fixed_reversal_counts':fixed_reversal(),'clock_magic':clock_magic()}
    (HERE/'inspiration_audit.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'polynomials':[x['characteristic_descending'] for x in result['matrix_audits']],
                      'cap_counts':result['projective']['cap_counts'],
                      'square_line_sums':result['clock_magic']['squared_line_sums']}))


if __name__=='__main__':main()
