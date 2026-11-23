// Node.js; also runnable as a self-contained JavaScript script.
"use strict";
const assert = (x, m) => { if (!x) throw Error(m); };
const bitIndex = x => 31 - Math.clz32(x);
const edges = n => Array.from({length:n}, (_,i) =>
  Array.from({length:n-i-1}, (_,k) => [i,i+k+1])).flat();
const flip = (r,[i,j]) => { r[i] ^= 1<<j; r[j] ^= 1<<i; };
function height(r) {
  const n = r.length, F = (1<<n)-1, d = new Uint32Array((1<<n)*n);
  for (let i=0; i<n; i++) d[(1<<i)*n+i] = 1;
  for (let s=1; s<=F; s++) for (let a=s; a; a&=a-1) {
    const i = bitIndex(a&-a), q = d[s*n+i];
    if (q) for (let b=r[i]&(F^s); b; b&=b-1) {
      const j = bitIndex(b&-b);
      d[(s|(1<<j))*n+j] += q;
    }
  }
  let h = 0;
  for (let i=0; i<n; i++) h += d[F*n+i];
  return h;
}
function strong(r) {
  // With an anchored Hamilton path, only reachability back to 0 is needed.
  let t = 1, z = 0;
  while (z !== t) {
    z = t;
    for (let i=1; i<r.length; i++) if (r[i]&t) t |= 1<<i;
  }
  return t === (1<<r.length)-1;
}
function signature(r) {
  const n = r.length, p = [], es = edges(n);
  let best = Infinity;
  function visit(s) {
    if (p.length===n) {
      let x = 0;
      for (const [i,j] of es) x = 2*x+((r[p[i]]>>p[j])&1);
      if (x<best) best = x;
      return;
    }
    for (let i=0; i<n; i++) if (!(s&(1<<i))) {
      p.push(i); visit(s|(1<<i)); p.pop();
    }
  }
  visit(0);
  return best;
}
const cache = new Map();
function scan(n, h, quotient=false) {
  const key = [n,h,quotient].join();
  if (cache.has(key)) return cache.get(key);
  const es = edges(n), free = es.filter(([i,j])=>j>i+1);
  const forward = Array.from({length:n}, (_,i)=>((1<<n)-1)^((1<<(i+1))-1));
  const seen = new Set(), reps = [];
  let max = 0;
  for (let m=0; m<(1<<free.length); m++) {
    const r = forward.slice();
    for (let a=m; a; a&=a-1) flip(r,free[bitIndex(a&-a)]);
    const count = height(r);
    if (count>max) max = count;
    if (count!==h || !strong(r)) continue;
    if (quotient) {
      const sig = signature(r);
      if (seen.has(sig)) continue;
      seen.add(sig);
    }
    reps.push(r);
  }
  const out = {reps,max}; cache.set(key,out); return out;
}
function join(bs) {
  const r = [];
  let l = 0, N = bs.reduce((a,b)=>a+b.length,0);
  for (const b of bs) {
    const suffix = ((1<<N)-1)^((1<<(l+b.length))-1);
    for (const row of b) r.push((row<<l)|suffix);
    l += b.length;
  }
  return r;
}
function padded(blocks,n) {
  const k = n-blocks.reduce((a,b)=>a+b.length,0), out = [];
  if (k<0) return out;
  function distribute(i,k,bs) {
    if (i===blocks.length) {
      out.push(join(bs.concat(Array.from({length:k},()=>[0]))));
      return;
    }
    for (let l=0; l<=k; l++) distribute(i+1,k-l,
      bs.concat(Array.from({length:l},()=>[0]),[blocks[i]]));
  }
  distribute(0,k,[]);
  return out;
}
function cores(n,h,q,expected) {
  const a = scan(n,h,q).reps;
  assert(a.length===expected, [n,h,a.length]);
  return a;
}
assert(scan(4,0).max===5 && scan(5,0).max===15,"order maxima");
const h11 = cores(5,11,false,11).map(a=>[a]);
const h19 = cores(6,19,false,38).map(a=>[a]);
const three = cores(3,3,true,1)[0], five = cores(4,5,true,1)[0];
const thirteen = cores(5,13,true,1)[0];
const h15 = [...cores(5,15,true,2),...cores(6,15,true,2)].map(a=>[a]);
h15.push([three,five],[five,three]);
const h39 = [...cores(6,39,true,0),...cores(7,39,true,2)].map(a=>[a]);
h39.push([three,thirteen],[thirteen,three]);
function audit(bs,h,n,k,shapes,total) {
  let tests = 0, hits = 0;
  const source = bs.flatMap(b=>padded(b,n)), es = edges(n);
  for (const r of source) {
    assert(height(r)===h,"source");
    for (let i=0; i<es.length; i++) {
      flip(r,es[i]);
      if (k===1) { tests++; hits += +(height(r)===h*h); }
      else for (let j=i+1; j<es.length; j++) {
        flip(r,es[j]); tests++; hits += +(height(r)===h*h); flip(r,es[j]);
      }
      flip(r,es[i]);
    }
  }
  assert(source.length===shapes && tests===total && hits===0,[h,k,n,source.length,tests,hits]);
  console.log({h,k,n,shapes,tests,hits});
}
audit(h11,11,7,1,33,693); audit(h11,11,6,2,22,2310);
audit(h15,15,9,1,30,1080); audit(h15,15,7,2,12,2520);
audit(h19,19,8,1,114,3192); audit(h19,19,8,2,114,43092);
audit(h39,39,10,1,20,900); audit(h39,39,9,2,12,7560);
function square(r,e,f,target) {
  const a = height(r);
  flip(r,e); const b = height(r);
  flip(r,f); const d = height(r);
  flip(r,e); const c = height(r);
  flip(r,f);
  assert([a,b,c,d].join()===target.join(),[a,b,c,d]);
  console.log([a,b,c,d]);
}
square([126,100,72,82,102,76,0],[0,2],[0,6],[11,37,111,121]);
square([254,180,248,242,224,192,130,0],[0,1],[0,7],[15,33,231,225]);
square([510,508,504,304,32,64,152,312,112],[0,4],[2,4],[19,311,91,361]);
square([1002,1012,1017,1010,993,960,640,768,576,0],[1,6],[2,9],[39,147,663,1521]);
