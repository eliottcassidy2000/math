"""Independent exact H9 apex-fibre and substitution-response audit.

Universe: all binary cut words, with an ordered transitive base and named apex.
Proof bounds the trimmed H9 core length by4. Finite checks do not supply that
bound. No dependence on inherited tournament-count code.
"""
from itertools import product, permutations
from pathlib import Path
import json

CORES = ('OIII', 'OIOI', 'OOOI')
SUBS = {'A': {'O':'OI','I':'I'}, 'B': {'O':'OI','I':'OIO'},
        'C': {'O':'O','I':'OOI'}}

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def state(word):
    h,p,e = 1,0,0
    for c in word:
        if c == 'O': p,e = 2*p+1,1
        else: h,p,e = h+p+e,2*p,0
    return h,p,e

def rows(word):
    n = len(word)
    ans=[]
    for i,c in enumerate(word):
        ans.append(sum(1 << j for j in range(i+1,n)) + ((1 << n) if c=='I' else 0))
    ans.append(sum(1 << j for j,c in enumerate(word) if c=='O'))
    return tuple(ans)

def literal_h(word):
    adj=rows(word); n=len(adj); size=1<<n
    dp=[[0]*n for _ in range(size)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,size):
        for i in range(n):
            value=dp[mask][i]
            if value:
                choices=adj[i] & (size-1-mask)
                while choices:
                    bit=choices & -choices; choices-=bit
                    dp[mask|bit][bit.bit_length()-1]+=value
    return sum(dp[-1])

def trim(word):
    return word.lstrip('I').rstrip('O')

def substitute(word, branch):
    return ''.join(SUBS[branch][c] for c in word)

def response(word):
    # H(word appended to state h,p,e)=h+alpha*p+beta*e+gamma.
    hp=[0,0,0,0]; pp=[0,1,0,0]; ep=[0,0,1,0]
    for c in word:
        if c=='O':
            pp=[2*x for x in pp]; pp[3]+=1; ep=[0,0,0,1]
        else:
            hp=[a+b+c for a,b,c in zip(hp,pp,ep)]
            pp=[2*x for x in pp]; ep=[0,0,0,0]
    return hp[1:]

def branch_closed(core,a,b,branch):
    w=substitute(core,branch)
    h0,p0,_=state(w); alpha,beta,gamma=response(w)
    need(beta==0 and gamma==h0-1,'core startsO response')
    if branch=='A':
        numerator=2*(3*p0+2)*(4**b-1)+6*b
        need(numerator%9==0,'A divisibility')
        return h0+numerator//9
    if branch=='C':
        numerator=(42*alpha+24)*(8**a-1)+28*a
        need(numerator%49==0,'C divisibility')
        return h0+numerator//49
    numerator=(35*alpha+10)*(8**a-1)+28*a
    need(numerator%49==0,'B left divisibility')
    h=h0+numerator//49
    p=p0+(2**len(w))*5*(8**a-1)//7
    numerator=2*(3*p+2)*(4**b-1)+6*b
    need(numerator%9==0,'B right divisibility')
    return h+numerator//9

def iso(a,b):
    if len(a)!=len(b): return None
    n=len(a)
    for perm in permutations(range(n)):
        if any(a[i].bit_count()!=b[perm[i]].bit_count() for i in range(n)): continue
        if all(bool(a[i]>>j&1)==bool(b[perm[i]]>>perm[j]&1)
               for i in range(n) for j in range(n)):
            return perm
    return None

def main():
    counts=[]; checks=0; dp_checks=0
    for length in range(11):
        words=[''.join(t) for t in product('OI',repeat=length)]
        matches=[w for w in words if state(w)[0]==9]
        need(all(trim(w) in CORES for w in matches),'complete trimmed fibre')
        expected=3*(length-3) if length>=4 else 0
        need(len(matches)==expected,'decorated count')
        counts.append({'vertices':length+1,'decorated_H9':len(matches)})
        for w in words:
            if length<=8:
                need(literal_h(w)==state(w)[0],'literal count recurrence')
                dp_checks+=1
    signatures=[]
    for v in CORES:
        entry={'core':v,'rows':rows(v),'branches':{}}
        for branch in SUBS:
            w=substitute(v,branch)
            entry['branches'][branch]={'word':w,'state':state(w),
                                       'response':response(w)}
            for a in range(9):
                for b in range(9):
                    actual=state(substitute('I'*a+v+'O'*b,branch))[0]
                    need(actual==branch_closed(v,a,b,branch),'closed branch formula')
                    checks+=1
                    if len(substitute('I'*a+v+'O'*b,branch))<=12:
                        need(literal_h(substitute('I'*a+v+'O'*b,branch))==actual,
                             'literal branch graph')
                        dp_checks+=1
        signatures.append(entry)
    need(state('OIOIO')[0]==state('OOOI')[0]==state('OIII')[0]==9,'hostile sources')
    hostile=[state(substitute(w,'A'))[0] for w in ('OIOIO','OOOI','OIII')]
    need(hostile==[107,73,17],'hostile branch outputs')
    chart_switch=iso(rows('OIII'),rows('OOOI'))
    need(chart_switch is not None,'actual apex chart switch')
    need(iso(rows('OIII'),rows('OIOI')) is None,'distinct strong kernels')
    # All-height difference bound in H9 A-branch image: max index M>=7
    # gives absolute difference >=(27*4^M-6*M-96)/9 >30000.
    a_values=[{'core':v,'padding':b,'H':branch_closed(v,0,b,'A')}
              for v in CORES for b in range(7)]
    gap_pairs=[(x,y) for x in a_values for y in a_values if y['H']-x['H']==30000]
    need(not gap_pairs,'no30000 difference in bounded residual')
    need(27*4**7-6*7-96 >9*30000,'all-height difference cutoff')
    output={'status':'FINITE-EXACT controls for proved all-word classification',
      'decorated_counts':counts,'branch_signatures':signatures,
      'branch_formula_controls':checks,'independent_graph_DP_checks':dp_checks,
      'same_graph_chart_switch_OIII_to_OOOI':chart_switch,
      'collision_A_outputs':hostile,
      'A_image_no30000_difference':{'max_index_to_check':6,
          'finite_values':a_values,'pairs':gap_pairs},
      'ordinary_H9_family_generating_function':'2*z^5/(1-z)^2',
      'decorated_H9_family_generating_function':'3*z^5/(1-z)^2'}
    path=Path(__file__).with_name('h9_fibres.json')
    path.write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'result':'PASS','branch_controls':checks,'literal_controls':dp_checks,
                      'output':str(path)}))

if __name__=='__main__': main()
