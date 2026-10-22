#!/usr/bin/env python3
import json
from pathlib import Path
from playwright.sync_api import sync_playwright
P=Path(__file__).resolve().parent
with sync_playwright() as p:
 browser=p.chromium.launch(executable_path='/usr/bin/chromium',headless=True,args=['--no-sandbox'])
 page=browser.new_page(viewport={'width':1400,'height':1000},device_scale_factor=1)
 errors=[];requests=[]
 page.on('pageerror',lambda error:errors.append(str(error)))
 page.on('request',lambda req:requests.append(req.url))
 page.goto((P/'atlas.html').as_uri())
 assert page.locator('button.cell').count()==91
 assert page.locator('button.diagonal').count()==13
 assert page.locator('button.leaf').count()==8
 assert '65' in page.locator('#degree-note').inner_text()
 # Every target: exercise actual click and validate inspector/selection.
 data=json.loads((P/'atlas-data.json').read_text())
 for k,t in enumerate(data['targets']):
  page.locator(f'button[data-index="{k}"]').click()
  assert page.locator('#address').inner_text()==f"Address {{{', '.join(map(str,t['root_pair']))}}}"
  assert page.locator('button[aria-pressed="true"]').count()==1
  assert page.locator('#destination').inner_text().endswith(f"{t['total']} captured line"+('s' if t['total']!=1 else ''))
 # Native selector and keyboard.
 page.select_option('#target-select','0');assert '74, 74' in page.locator('#address').inner_text()
 page.locator('#target-select').focus();page.keyboard.press('ArrowDown');assert page.locator('#target-select').input_value()=='1'
 page.locator('button[data-index="89"]').click()
 page.screenshot(path=str(P/'atlas-desktop.png'),full_page=True)
 page.set_viewport_size({'width':390,'height':844});page.screenshot(path=str(P/'atlas-mobile.png'),full_page=True)
 assert page.evaluate('document.documentElement.scrollWidth <= innerWidth')
 assert not errors,errors
 assert all(u.startswith('file:') for u in requests),requests
 report={'passed':True,'browser':'Chromium via Playwright','all_91_cell_clicks':True,'selector_and_keyboard':True,'desktop_viewport':[1400,1000],'mobile_viewport':[390,844],'page_overflow_mobile':False,'javascript_errors':errors,'network_calls':0,'screenshots':['atlas-desktop.png','atlas-mobile.png']}
 (P/'browser-validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));browser.close()
