
const gcd=(a,b)=>{a=a<0n?-a:a;b=b<0n?-b:b;while(b){[a,b]=[b,a%b]}return a};
class Q{
 constructor(a,b=1n){a=BigInt(a);b=BigInt(b);if(b<0n){a=-a;b=-b}let g=gcd(a,b);this.a=a/g;this.b=b/g}
 add(x){x=q(x);return new Q(this.a*x.b+x.a*this.b,this.b*x.b)}
 neg(){return new Q(-this.a,this.b)}
 sub(x){return this.add(q(x).neg())}
 mul(x){x=q(x);return new Q(this.a*x.a,this.b*x.b)}
 div(x){x=q(x);return new Q(this.a*x.b,this.b*x.a)}
 pow(n){return new Q(this.a**BigInt(n),this.b**BigInt(n))}
 cmp(x){x=q(x);let a=this.a*x.b-x.a*this.b;return a<0n?-1:a>0n?1:0}
 abs(){return this.a<0n?this.neg():this}
 str(){return this.b===1n?""+this.a:this.a+"/"+this.b}
}
const q=x=>x instanceof Q?x:new Q(x), uniq=xs=>[...new Map(xs.map(x=>[x.str(),x])).values()].sort((a,b)=>a.cmp(b));
function triangle(v,t){let x=t.mul(v),r=x.sub(x.a/x.b);return r.cmp(new Q(1,2))<=0?r:q(1).sub(r)}
function build(vs){
 let bs=uniq(vs.flatMap(v=>Array.from({length:2*v+1},(_,k)=>new Q(k,2*v)))),segs=[];
 for(let j=0;j<bs.length-1;j++){
 let a=bs[j],b=bs[j+1],mid=a.add(b).div(2);
 let lines=vs.map(v=>{let xa=triangle(v,a),xb=triangle(v,b),s=xb.sub(xa).div(b.sub(a));return {s,c:xa.sub(s.mul(a))}});
 let cuts=[a,b];
 for(let k=0;k<lines.length;k++)for(let l=k+1;l<lines.length;l++){let u=lines[k],w=lines[l];if(u.s.cmp(w.s)){let t=w.c.sub(u.c).div(u.s.sub(w.s));if(t.cmp(a)>0&&t.cmp(b)<0)cuts.push(t)}}
 cuts=uniq(cuts);
 for(let k=0;k<cuts.length-1;k++){let x=cuts[k],y=cuts[k+1],m=x.add(y).div(2),line=lines.reduce((u,w)=>u.s.mul(m).add(u.c).cmp(w.s.mul(m).add(w.c))<=0?u:w);
 segs.push({x,y,s:line.s,lo:line.s.mul(x).add(line.c),hi:line.s.mul(y).add(line.c)})}
 }
 return segs;
}
function test(vs){
 const segs=build(vs),at=new Map();
 const add=(h,a)=>at.set(h.str(),{h,a:(at.get(h.str())?.a||q(0)).add(a)});
 let moments=Array.from({length:11},()=>q(0));
 for(let e of segs){if(!e.s.a)throw Error("zero-slope");let low=e.lo.cmp(e.hi)<0?e.lo:e.hi,high=e.lo.cmp(e.hi)>0?e.lo:e.hi;add(low,q(-1).div(e.s.abs()));add(high,q(1).div(e.s.abs()));for(let k=0;k<moments.length;k++)moments[k]=moments[k].add(e.hi.pow(k+1).sub(e.lo.pow(k+1)).div(e.s.mul(k+1)))}
 const atoms=[...at.values()].filter(z=>z.a.a).sort((a,b)=>a.h.cmp(b.h));
 for(let k=0;k<moments.length;k++){let rhs=atoms.reduce((s,z)=>s.add(z.a.mul(z.h.pow(k+1))),q(0)).div(k+1);if(rhs.cmp(moments[k]))throw Error("moment mismatch")}
 if(moments[0].cmp(1))throw Error("mass");
 if(atoms.reduce((s,z)=>s.add(z.a),q(0)).a)throw Error("atom mass");
 if(atoms.reduce((s,z)=>s.add(z.a.mul(z.h)),q(0)).cmp(1))throw Error("height normalization");
 let max=segs.reduce((m,e)=>m.cmp(e.lo)>0?m:e.lo,q(0));
 if(atoms.at(-1).h.cmp(max)||atoms.at(-1).a.cmp(0)<=0)throw Error("maximum");
 // independent superlevel lengths at each midpoint between atom heights
 for(let i=0;i<atoms.length-1;i++){let y=atoms[i].h.add(atoms[i+1].h).div(2);
 let direct=segs.reduce((sum,e)=>{let lo=e.lo.cmp(e.hi)<0?e.lo:e.hi,hi=e.lo.cmp(e.hi)>0?e.lo:e.hi;if(y.cmp(lo)<=0)return sum.add(e.y.sub(e.x));if(y.cmp(hi)>=0)return sum;return sum.add(hi.sub(y).div(e.s.abs()))},q(0));
 let formula=atoms.reduce((sum,z)=>z.h.cmp(y)>0?sum.add(z.a.mul(z.h.sub(y))):sum,q(0));if(direct.cmp(formula))throw Error("level mismatch")}
 return {vs,atoms:atoms.map(z=>[z.h.str(),z.a.str()]),moments:moments.map(x=>x.str()),maximum:max.str(),segments:segs.length};
}

const controls=[[1],[2],[1,1],[1,2],[1,3],[1,4],[2,3],[2,4],[1,2,3],[2,3,5],[1,4,7],[3,7,11]];
console.log(JSON.stringify(controls.map(test),null,2));
for(let N=2;N<=40;N+=2){
 const r=test([1,N]);
 if(r.atoms.length-1!==N-1)throw Error("family count "+N);
 const expected=new Map([["0",new Q(-(4*N-2),N).str()]]);
 for(let j=1;j<=N/2;j++)expected.set(new Q(j,N+1).str(),new Q(2*(N+1),N).str());
 for(let j=1;j<N/2;j++)expected.set(new Q(j,N-1).str(),new Q(-2*(N-1),N).str());
 if(r.atoms.some(([h,a])=>expected.get(h)!==a))throw Error("family values "+N);
}
console.log("PASS: controls; even N=2..40; moments k=0..10; superlevel lengths.");

const fs=require('fs');
const census=[];
for(let mask=1;mask<256;mask++){
 const v=Array.from({length:8},(_,i)=>i+1).filter((_,i)=>mask&(1<<i));
 if(v.length<=4)census.push(test(v));
}
for(const v of [Array.from({length:13},(_,i)=>i+1),[...Array.from({length:12},(_,i)=>i+1),182]])census.push(test(v));
fs.writeFileSync(__dirname+'/independent-census.json',JSON.stringify(census,null,2));
let count=0;
for(let a=1;a<=8;a++)for(let b=a;b<=10;b++)for(let c=b;c<=12;c++){
 const x=test([a,b,c]),y=test([2*a,2*b,2*c]);
 if(JSON.stringify(x.atoms)!==JSON.stringify(y.atoms))throw Error('scale mismatch');
 count++;
}
console.log('PASS: 162 original census rows + both large controls + '+count+' triples with doubling.');
