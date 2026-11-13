"""Audit every target module among bounded positive-small-defect repairs."""
import json
import sparse_repair as s

def main():
    data=json.loads((s.HERE/'sparse_repair_radius4.json').read_text(encoding='utf-8'));out=[];checks=0
    for hit in data['hits']:
        if not hit['defect']:continue
        child=tuple(hit['rows']);parent=tuple(row&511 for row in child[:9]);port=hit['port']
        for module in range(9):
            if module==port:continue
            checks+=1
            hc=s.exp_coeff(s.module_coeff(child,module));dc=s.exp_coeff(s.module_coeff(parent,module,start=port))
            if hc=={1:23,2:135}:out.append(dict(base=hit,module=module,H=hc,defect=dc))
    (s.HERE/'repair_module_relocation.json').write_text(json.dumps(dict(status='FINITE-EXACT bounded target-module search with exact all-height polynomial certificates',scope='every positive-defect<=2 target found at radius<=4, every one-vertex transitive module other than final root or its sole predecessor',module_checks=checks,hits=out),indent=2)+'\n',encoding='utf-8')
    print('PASS',checks,'target-module choices;',len(out),'exact matching all-height targets')
if __name__=='__main__':main()
