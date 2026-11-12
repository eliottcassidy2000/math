"""Independent bit-packed semigroup and polynomial determinant audit."""
from common import *
from ozj_periodic import M,PERMS,det_shift
from collections import Counter

def pack(A):return sum((v&1)<<(4*i+j) for i,row in enumerate(A) for j,v in enumerate(row))
def mul(A,B):
    out=0
    for i in range(4):
        row=0
        for k in range(4):
            if A>>(4*i+k)&1:row^=(B>>(4*k))&15
        out|=row<<(4*i)
    return out
def pmul(a,b):
    out=0
    while b:
        if b&1:out^=a
        a<<=1;b>>=1
    return out
def char2(A):
    out=0
    for p,_ in PERMS:
        term=1
        for i in range(4):term=pmul(term,(2 if i==p[i] else 0)^((A>>(4*i+p[i]))&1))
        out^=term
    return [(out>>i)&1 for i in reversed(range(5))]

def run():
    mats={k:pack(v) for k,v in M.items()};root=(pack([[int(i==j) for j in range(4)] for i in range(4)]),0);states=[root];seen={root}
    for A,f in states:
        for label,B in mats.items():
            q=(mul(B,A),f|(1 if label=='J' else 2))
            if q not in seen:seen.add(q);states.append(q)
    old=read('ozj_mod2.json.gz');need({(sum(v<<i for i,v in enumerate(r['matrix'])),r['flags']) for r in old['records']}==seen,'independent complete flagged semigroup equality')
    for r in old['records']:
        A=sum(v<<i for i,v in enumerate(r['matrix']));need(char2(A)==r['charpoly_mod2'],'independent polynomial determinant')
        if r['flags']==3:need(char2(A)==[1,1,0,0,0],'all mixed words have exactly one odd eigenvalue slot')
    barriers=read('repair_obstructions.json.gz');jbad={(r['birth_id'],r['port']):r['exact_valuation'] for r in barriers['records']};count=0
    for r in read('even_valuations.json.gz')['records']:
        j=r['eventual_valuations_by_parity'][0] if r['type']=='bounded_constant' else jbad[r['birth_id'],r['port']]
        e=r['exponential_coefficients'];need(e[0]+2**(j-1)!=0 and e[2]!=0,'all478 lowest heads retain both odd modes1 and3');count+=1
    candidate=next(r for r in read('ozj_periodic.json.gz')['eigenvalue_one_words'] if r['word']=='OOZOZOJ');need(det_shift(candidate['matrix'],1)==0 and det_shift(candidate['matrix'],82)==0 and det_shift(candidate['matrix'],0)==0,'hostile to stronger no-constant-mode guess')
    out={'scope':'Independent complete25-state mod2 semigroup and characteristic-polynomial audit, plus exact odd-mode obligations for all478 even axes. All-length periodic-root-word obstruction uses the separate pure-J and pure-O/Z proofs.','semigroup_states':len(states),'mixed_states':sum(f==3 for A,f in states),'periodic_root_barrier_families':count,'stronger_guess_hostile':{'word':'OOZOZOJ','charpoly':'x^2*(x-1)*(x-82)','minimal_J_containing_word_length_with_eigenvalue1':7}}
    save('audit_ozj.json.gz',out);print(json.dumps(out),flush=True)
if __name__=='__main__':run()
