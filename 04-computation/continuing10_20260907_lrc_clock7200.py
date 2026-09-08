"""Pinned composition of the audited connected7200 word closures."""
from pathlib import Path
from hashlib import sha256
import argparse,json,sys
sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent if HERE.name=='04-computation' else Path('C:/w/s0905')
OUT=ROOT/'05-knowledge/results' if HERE.name=='04-computation' else HERE
parser=argparse.ArgumentParser()
parser.add_argument('--last-b',type=Path,default=ROOT/'05-knowledge/results/continuing10_20260907_lrc_last_b_certificate.json')
args=parser.parse_args()
gates=0
def need(ok,label):
    global gates
    gates+=1
    if not ok:raise ArithmeticError(label)
def read(path,pin):
    raw=path.read_bytes();need(sha256(raw).hexdigest()==pin,'frozen supplier '+path.name)
    return json.loads(raw)
def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':')).encode()
def main():
    old=read(ROOT/'05-knowledge/results/continuing8_20260906_lrc_minimum_tree_certificate.json',
      '580a7c930103aab3bea867ad463a90b0e0208323a90ee95a685ff811a761582d')
    c=next(r for r in old['clocks'] if r['t']==7200)
    need(c['word_count']==76814 and c['survivor_count']==15,'complete inherited7200 universe')
    third=read(ROOT/'05-knowledge/results/continuing10_20260907_lrc_third_wedge_certificate.json',
      '60a6de57a623dbe64a23a77113f140b312f5496247c2253fbe695872740039e9')
    a=read(ROOT/'05-knowledge/results/continuing10_20260907_lrc_last_a_certificate.json',
      'bce97b849036ba89ca8e8d1593b2ab03ff7aa9ecc647bdc72d55bd4962383274')
    b=read(args.last_b,'01381d94aaeaf92ca7dfb37b902c29db16da0aa36c5587b602dd836b5616cabd')
    need(a['E']==116 and b['E']==103,'two literal remaining excesses')
    need(third['remaining_words']==[a['word'],b['word']],'exact remaining-word identity')
    byword={tuple(r['word']):r for r in third['topology']}
    need(len(byword)==15,'all15 supplier rows retained')
    closed=[]
    for d,E,M in c['survivors']:
        t=tuple(d);need(t in byword,'same positional word universe')
        row=byword[t]
        if row['closed']:reason='third-wedge topology'
        elif d==a['word']:reason='last A native zero-arm path'
        elif d==b['word']:reason='last B fresh doubled zero-arm paths'
        else:raise ArithmeticError('unclosed word '+str(d))
        closed.append(dict(word=d,E=E,closure=reason))
    need(sum(r['closure']=='third-wedge topology' for r in closed)==13,'13 plus1 plus1 exact partition')
    need(all(r['credit']>a['E'] for r in a['products']) and len(a['products'])==70,'A complete closing bank')
    need(len(b['wedges'])==2,'B both required path types')
    previous=old['new_scales'];need(len(previous)==7646 and previous==sorted(set(previous)),'old sorted necessary array')
    need(sha256(canonical(previous)).hexdigest()=='8ffc6d14b3883cf7e02c3ab02ddca5339d909a8051411def9096dee83b0aaed7','old semantic array pin')
    need(previous.count(7200)==1,'unique removed clock')
    current=[t for t in previous if t!=7200]
    need(len(current)==7645 and max(current)==11995,'new cardinality and unchanged maximum')
    need(set(previous)-set(current)=={7200} and not set(current)-set(previous),'only authorized mathematical removal')
    cert=dict(status='FINITE-EXACT composition; analytic and audit scope in companion report',
        scope='primitive13, selected6 gcd7200, actual strict-atlas graph on complementary7 connected; no arbitrary entry theorem',
        inherited_word_count=76814,inherited_survivor_count=15,closed_words=closed,
        removed_clocks=[7200],necessary_clocks=current,necessary_clocks_sha256=sha256(canonical(current)).hexdigest(),
        old_count=7646,new_count=7645,maximum=11995)
    (OUT/'continuing10_20260907_lrc_clock7200_certificate.json').write_bytes(canonical(cert)+b'\n')
    print('CONNECTED CLOCK7200 CLOSED: inherited76814 full-profile words;13+1+1 exhaust all15 minimum-tree survivors')
    print('NECESSARY_CLOCKS 7646 ->7645; removed only7200; maximum11995 unchanged')
    print('NEW_ARRAY_SHA256',cert['necessary_clocks_sha256'])
    print('SCOPE selected six-body gcd7200 and actual connected strict complement; general LRC14 and disconnected entry OPEN')
    print('Always-active exact gates:',gates)
if __name__=='__main__':main()
