'use strict';
// A minimal DOM harness validates selection logic, not browser layout/rendering.
const fs=require('fs'),vm=require('vm'),path=require('path'),P=__dirname;
const html=fs.readFileSync(path.join(P,'atlas.html'),'utf8');
const assert=(x,m)=>{if(!x)throw Error(m)};
const nodes=[];
class Element{constructor(tag='div'){this.tag=tag;this.children=[];this.attrs={};this.dataset={};this.value='';this.className='';this.textContent='';this.innerHTML='';nodes.push(this)}append(e){this.children.push(e)}setAttribute(k,v){this.attrs[k]=v}}
const ids={};for(const m of html.matchAll(/\bid="([^"]+)"/g)){assert(!ids[m[1]],'duplicate id');ids[m[1]]=new Element()}
ids['atlas-data'].textContent=html.match(/<script id="atlas-data" type="application\/json">([\s\S]*?)<\/script>/)[1];
const script=html.match(/<script>\s*([\s\S]*?)<\/script>/)[1];
assert(!/fetch\(|XMLHttpRequest|WebSocket|src=["']https?:/.test(html),'network dependency');
const document={getElementById:id=>{assert(ids[id],'missing '+id);return ids[id]},createElement:t=>new Element(t),querySelectorAll:q=>{assert(q==='.cell','unknown selector');return nodes.filter(n=>n.className.split(' ').includes('cell'))}};
const context=vm.createContext({document,console});vm.runInContext(script,context);
const data=JSON.parse(ids['atlas-data'].textContent),buttons=document.querySelectorAll('.cell');assert(buttons.length===91,'button count');
for(let k=0;k<91;k++){buttons[k].onclick();const t=data.targets[k];assert(ids.address.textContent===`Address {${t.root_pair.join(', ')}}`,'address');assert(buttons.filter(b=>b.attrs['aria-pressed']==='true').length===1,'pressed count');assert(ids['degree-note'].textContent.includes(`= ${t.degree}.`),'degree');for(let s of data.sources){const v=t.incidences[s];assert(ids.sources.innerHTML.includes(`${v.blue} blue`),'blue');assert(ids.sources.innerHTML.includes(`${v.black} black`),'black')}}
ids['target-select'].value='0';ids['target-select'].onchange();assert(ids.address.textContent==='Address {74, 74}','selector');
const report={passed:true,method:'Node.js VM with minimal DOM harness; no browser rendering',all_91_click_handlers:true,selector_change_handler:true,unique_ids:true,embedded_data:true,network_dependencies:0,browser_rendering:'NOT VERIFIED: sandbox Chromium socket creation failed; escalated launcher failed before browser startup',keyboard_and_layout:'NOT VERIFIED',screenshots:[]};fs.writeFileSync(path.join(P,'interaction-validation.json'),JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify(report,null,2));
