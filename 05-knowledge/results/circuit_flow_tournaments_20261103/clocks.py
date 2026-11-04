"""Complete induction certificates for the two-gate grammar modulo 8 and 27."""
import json
from math import gcd
from pathlib import Path
from clock_explore import step
from flow_operator import apply,compile_operator,need

HERE=Path(__file__).resolve().parent
OPS={c:compile_operator(s) for c,s in [('L',(16,17,3,7,12)),('U',(16,25,3,5,12))]}


def h(word):
    f=list(range(len(word)+2))
    for letter in word:
        f=[0]+[apply(OPS[letter],r,lambda x:f[x]) for r in range(1,len(f)-1)]
    return f[1]


def full_bank(M):
    states=[tuple(range(M))];words=[''];idx={states[0]:0};edges=[]
    for i,f in enumerate(states):
        row=[]
        for c in 'LU':
            q=step(f,c,M)
            need(q==tuple(apply(OPS[c],r,lambda x:f[x%M])%M for r in range(M)), 'independent exterior compiler modulo clock')
            if q not in idx:
                idx[q]=len(states);states.append(q);words.append(words[i]+c)
            row.append(idx[q])
        edges.append(row)
    return {'M':M,'states':states,'reach_words':words,'edges_LU':edges}


def main():
    eight=full_bank(8);twentyseven=full_bank(27)
    need(len(eight['states'])==5 and len(twentyseven['states'])==162,'complete finite clock bank sizes')
    values=[pow(2,i+12*i*i,27) for i in range(18)]
    need(len(set(values))==18,'clock states separated by current output')
    for f,word,edges in zip(twentyseven['states'],twentyseven['reach_words'],twentyseven['edges_LU']):
        i=values.index(f[1]);need(i==(word.count('L')+14*word.count('U'))%18,'word clock')
        for e,delta in zip(edges,(1,14)):
            need(twentyseven['states'][e][1]==values[(i+delta)%18],'all future transitions mod27')
    for f,word,edges in zip(eight['states'],eight['reach_words'],eight['edges_LU']):
        expected=1 if not word else (3 if word.count('L')%2 else 5)
        need(f[1]==expected,'mod8 word phase')
        a,b=[eight['states'][e][1] for e in edges]
        need((a,b)==({1:(3,5),3:(5,3),5:(3,5)}[f[1]]),'all future transitions mod8')
    combined=[]
    for i in range(18):
        v=next(n for n in range(216) if n%27==values[i] and n%8==(3 if i%2 else 5))
        combined.append(v)
    need(len(set(combined+[1]))==19,'nineteen outputs distinguish nineteen minimal clock states')
    witnesses=[];divisor=0
    for a,b in [('LLU','LUL'),('LLLU','LLUL')]:
        difference=h(a)-h(b);divisor=gcd(divisor,difference)
        witnesses.append({'a':a,'b':b,'counts':[h(a),h(b)],'difference':difference})
    need(divisor==216,'sharp universal order-blind modulus')
    data={'status':'PROVED by complete finite induction certificates','banks':[eight,twentyseven],
          'clock_mod27':values,'clock_mod216_nonempty':combined,'empty_output':1,
          'minimal_observable_states_mod216':19,'sharp_modulus_witnesses':witnesses,
          'scope':'start at singleton; unary L/U only; words inner to outer; no arbitrary sources or named-port attachments'}
    (HERE/'clocks.json').write_text(json.dumps(data,indent=2,sort_keys=True)+'\n',encoding='utf8',newline='\n')
    print('PASS mod8 states5; mod27 states162; sharp clock216, minimal observable states19',flush=True)


if __name__=='__main__':main()
