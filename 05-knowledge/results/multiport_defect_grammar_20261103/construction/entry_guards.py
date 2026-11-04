"""Scoped preparation-locus probe; actual marked rows, no scalar replacement."""
from pathlib import Path
import itertools,json
from algebra import operator,response,strong,need
HERE=Path(__file__).resolve().parent

def main():
    small=0;small_plus4=[]
    for n in range(1,5):
        pairs=list(itertools.combinations(range(n),2))
        for code in range(1<<len(pairs)):
            rows=[0]*n
            for k,(u,v) in enumerate(pairs):
                if code>>k&1:rows[u]|=1<<v
                else:rows[v]|=1<<u
            if not strong(rows):continue
            for w in range(n):
                for a in range(n):
                    item=operator(rows,w,a);small+=1
                    if item['S_a_2']-item['h']==2 and item['L']+item['K']==item['F']:
                        small_plus4.append(item)
    need(not small_plus4,'all labeled strong cores through four exclude this guard')
    bank_path=HERE.parents[1]/'session180_170807_returns_20261104/phases/small_suffix_search.json'
    bank=json.loads(bank_path.read_text())['bank']
    tested=0;plus2=[];plus4=[];both=[]
    for i,core in enumerate(bank):
        guards=[]
        for a in range(len(core['rows'])):
            item=operator(core['rows'],core['port'],a);tested+=1
            sigma=item['S_a_2']-item['h'];tau=item['L']+item['K']-item['F']
            if sigma not in (1,2) or tau:continue
            need(response(item['C2'],item['w'])==(3*item['h']+2*sigma,3*item['F']),'actual digit output')
            item.update(bank_index=i,sigma=sigma,tau=tau)
            (plus2 if sigma==1 else plus4).append(item);guards.append((a,sigma))
        if {sigma for a,sigma in guards}=={1,2}:both.append(dict(bank_index=i,guards=guards))
    need(plus4[0]['rows']==[20,1,2,7,14] and plus4[0]['w']==1 and plus4[0]['a']==2,'canonical smallest found +4 preparation')
    result=dict(status='FINITE_EXACT',small_full_universe={'maximum_order':4,'strong_marked_cases':small,'sigma2_tau0_hits':small_plus4},
        representative_universe={'source':str(bank_path.relative_to(HERE.parents[3])).replace('\\','/'),
            'marked_representatives':len(bank),'auxiliary_choices':tested,
            'scope':'Only the saved actual representative for each h,F pair; not all cores through order six.',
            'sigma1_tau0':plus2,'sigma2_tau0':plus4,'both_guards_same_representative':both})
    (HERE/'entry_guards.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'strong_cases_through4':small,'bank_choices':tested,'plus2':len(plus2),'plus4':len(plus4),'both':len(both)}))

if __name__=='__main__':main()
