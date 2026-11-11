"""Apply the proved low-defect source classification to the certified prime table."""
from pathlib import Path
from math import isqrt
import json
import audit_defect as audit

HERE=Path(__file__).resolve().parent


def main():
    data=json.loads((HERE/'tail_squares.json').read_text())
    candidates={}
    def add(H,budget,word):
        record={'H':H,'defect':budget,'inports':2,'contact_word':list(word),'order':len(word)+3,
                'source_rows':audit.graph(tuple(word)),'root':len(word)+2}
        if H not in candidates or (record['order'],budget,word)<(candidates[H]['order'],candidates[H]['defect'],tuple(candidates[H]['contact_word'])):
            candidates[H]=record
    for row in data['isolated']:
        if row['square_plus_one']:add(row['H'],row['budget'],tuple(row['word']))
    for row in data['rays']:
        for hit in row['finite_hits']:add(hit['H'],row['budget'],tuple(hit['word']))
    terms=[]
    for line in (HERE/'A002496_roots_through_1000000.txt').read_text().splitlines():
        i,p=map(int,line.split());terms.append((i,p))
    admitted=[];excluded=0
    for i,p in terms:
        if p==2:continue
        n=isqrt(p-1)
        audit.check(n*n+1==p,'table input shape')
        if n and n&(n-1)==0 and n>=8:
            k=2*(n.bit_length()-1)-4
            audit.check(k>=2 and k%2==0 and 16*2**k+1==p,'infinite tail instance')
            add(p,4,(3,)*k+(1,3))
        if p in candidates:
            admitted.append({'index':i,'square_root':n,**candidates[p]})
        else:excluded+=1
    audit.check([x['H'] for x in admitted]==[5,17,37,101,257,4357,15877,65537],'certified table intersection')
    audit.check(len(terms)==54110 and excluded==54101,'complete finite application')
    report={'status':'FINITE-EXACT application of an all-order source exclusion theorem',
        'root_bound':1000000,'certified_prime_terms':len(terms),'exceptional_even_prime':2,
        'admitted_odd_prime_count':len(admitted),'excluded_odd_prime_count':excluded,
        'meaning_of_excluded':'Every rooted tournament with this H and at least two inports has branching defect>=7, including weak sources by the audited SCC-factorization corollary. This excludes all source orders, not only a search box. It does not assert an optimal defect or restrict d1 roots.',
        'explicit_admitted_sources':admitted}
    text=json.dumps(report,indent=2,sort_keys=True)+'\n'
    (HERE/'low_defect_prime_filter.json').write_text(text,encoding='utf-8')
    print(json.dumps({'status':report['status'],'admitted':len(admitted),'excluded':excluded}))


if __name__=='__main__':main()
