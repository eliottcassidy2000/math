"""Small exact singular-letter and critical-word grammar checks."""
from itertools import product
from math import gcd
import json
import check as audit


def mul(a,b):
    return tuple(sum(a[3*i+k]*b[3*k+j] for k in range(3)) for i in range(3) for j in range(3))


def rank(m):
    a,b,c,d,e,f,g,h,j=m
    determinant=a*(e*j-f*h)-b*(d*j-f*g)+c*(d*h-e*g)
    minor=any(m[3*r+k]*m[3*s+l]!=m[3*r+l]*m[3*s+k]
              for r in range(3) for s in range(r+1,3)
              for k in range(3) for l in range(k+1,3))
    return 3 if determinant else 2 if minor else 1


def strong(rows):
    for start in range(len(rows)):
        seen=1<<start
        while True:
            more=seen
            for v,row in enumerate(rows):
                if (seen>>v)&1: more|=row
            if more==seen: break
            seen=more
        if seen!=(1<<len(rows))-1: return False
    return True


def main():
    ms=[]
    specs=[]
    for line in audit.SOURCE.read_text().splitlines():
        if line and not line.startswith('#'):
            left,right=line.split('|')
            ms.append(tuple(map(int,left.split())))
            k,bits,w,r=right.split()
            specs.append((int(k),bits,int(w),int(r)))
    C,N=ms[1],ms[5]
    r=(5,2,3)
    ell=N[:3]
    audit.need(ell==(1,4,4),'actual return observer')
    audit.need(mul(N,N)==tuple(a*b for a in r for b in ell),'strong-ray reset')
    X=(2,3,4)
    for a in X:
        ell_a=mul(C,ms[a])[:3]
        for b in X:
            bscale=ms[b][3]
            audit.need(mul(ms[b],ms[a])==tuple(z*v for z in (1,bscale,bscale) for v in ell_a),
                       'transitive-tail ray identity')
        audit.need(mul(N,mul(C,ms[a]))==tuple(z*v for z in r for v in ell_a),
                   'three-letter strong-ray reset')
    ident=(1,0,0,0,1,0,0,0,1)
    minimal=[]
    for length in range(1,4):
        for word in product(range(1,13),repeat=length):
            if any(word[i:i+len(v)]==v for v in minimal for i in range(length-len(v)+1)): continue
            matrix=ident
            for letter in word: matrix=mul(ms[letter],matrix)
            if rank(matrix)==1: minimal.append(word)
    expected=[(a,b) for a in X for b in X]+[(5,5)]+[(a,1,5) for a in X]
    audit.need(set(minimal)==set(expected),'minimal short rank-one patterns')
    checked=critical=terminal_repaired=0
    for length in range(5):
        for word in product(range(1,13),repeat=length):
            state=(1,0,1)
            for letter in word:
                child=audit.apply(ms[letter],state)
                audit.need(gcd(*child)%gcd(*state)==0,'gcd divisibility')
                state=child
            checked+=1
            h=state[0]
            is_critical=audit.prime(h) or (h%7==0 and audit.prime(h//7))
            if not is_critical: continue
            critical+=1
            audit.need(not any(word[i:i+2]==(5,5) for i in range(length-1)),'critical NN')
            audit.need(not any(word[i:i+3] in [(a,1,5) for a in X] for i in range(length-2)), 'critical XCN')
            xx=[i for i in range(length-1) if word[i] in X and word[i+1] in X]
            audit.need(all(i==length-2 for i in xx),'nonterminal critical XX')
            if xx:
                replacement=word[:-1]+(1,)
                other=(1,0,1)
                for letter in replacement: other=audit.apply(ms[letter],other)
                audit.need(other[0]==h,'terminal critical XX repair')
                terminal_repaired+=1
    rows,root,_=audit.compile_word((5,5),specs)
    state=audit.direct_state(rows,root)
    audit.need(strong(rows) and state==(25,10,15),'strong nonprimitive hostile')
    reopening=[]
    for a,target in zip(X,(5,9,17)):
        local_rows,local_root,_=audit.compile_word((a,1),specs)
        local_state=audit.direct_state(local_rows,local_root)
        audit.need(strong(local_rows) and local_state[0]==target,'transient cuts reopen strongly')
        reopening.append({'word':[a+1,2],'rows':local_rows,'root':local_root,'state':local_state})
    iso=(0,4,5,6,1,2,3)
    target_rows=(110,16,2,6,13,30,62)
    audit.need(all(sum(1<<iso[v] for v in range(7) if (rows[u]>>v)&1)==target_rows[iso[u]]
                   for u in range(7)) and iso[root]==1,'absorption endpoint isomorphism')
    seven_rows,seven_root,_=audit.compile_word((6,9,12),specs)
    seven_state=audit.direct_state(seven_rows,seven_root)
    audit.need(strong(seven_rows) and seven_state==(1757,994,1372) and gcd(*seven_state)==7,
               'unattained gcd in a seven-prime state')
    audit.need(audit.prime(251) and seven_state[0]==7*251,'seven-prime target')
    report={'status':'PASS','generator_ranks':[rank(m) for m in ms],
            'minimal_rank_one_patterns_through_length3':[[i+1 for i in w] for w in minimal],
            'r_free_words_through_length4':checked,'critical_words':critical,
            'terminal_XX_replacements':terminal_repaired,
            'transient_cut_reopening':reopening,
            'absorption_endpoint_isomorphism':{'old_to_endpoint':iso,'endpoint_rows':target_rows,'endpoint_root':1},
            'strong_nonprimitive_hostile':{'word':[6,6],'rows':rows,'root':root,'state':state,'gcd':gcd(*state)},
            'unattained_gcd_hostile':{'word':[7,10,13],'rows':seven_rows,'root':seven_root,
                                    'state':seven_state,'gcd':gcd(*seven_state),'H_factorization':[7,251]}}
    output=audit.HERE/'singular_receipt.json'
    output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report))


if __name__=='__main__': main()
