"use strict";
const assert=require("node:assert/strict");
function H(w){let h=1n,x=0n,b=0n;for(const c of w){const d=BigInt(c);if(!d) h+=x+b;x=2n*x+d;b=d;}return h;}
const bits=x=>x.toString(2);
const rc=w=>[...w].reverse().map(d=>d==="0"?"1":"0").join("");
function inv(w){let h=1n; for(let i=0;i<w.length;i++)for(let j=i+1;j<w.length;j++)if(w[i]==="1"&&w[j]==="0")h+=j===i+1?2n:1n<<BigInt(j-i-1);return h;}
function sq(n){if(n<2n)return n;let x=1n<<BigInt(Math.ceil(bits(n).length/2));for(;;){let y=(x+n/x)/2n;if(y>=x)return x;x=y;}}
function pair(N){let x=sq(N/2n);if(x%2n)x++;while(2n*x*x<N)x+=2n;const u=bits(x),v=rc(bits(x+2n));return [u+"0"+v,u+"1"+v];}
function child(h,t){return bits(h-t)+rc(bits(h+t));}
function checkShell(h,k,c){assert(h**c<k&&k+1n<(h+1n)**c);}
for(let i=100n;i<=100000n;i++){
 const w=pair(i),a=H(w[0]),b=H(w[1]);
 assert.equal(b-a,2n);assert(a>i && (b-i)**2n<256n*i);
 if(i<=2000n) assert.deepEqual(w.map(inv),[a,b]);
}
for(let h=5n;h<100000n;h+=2n){
 const a=H(child(h,1n)),b=H(child(h,h===37n?5n:3n));
 checkShell(h,a,2n);checkShell(h,b,2n);assert.notEqual(a,b);
 if(h<1000n){assert.equal(a,inv(child(h,1n)));assert.equal(b,inv(child(h,h===37n?5n:3n)));}
}
let h=63n;
assert.equal(H("1001100"),h);
for(let i=0;i<4;i++){
 const w=pair(h**3n+1n);const v=w.map(H);
 for(const k of v)checkShell(h,k,3n);
 assert.deepEqual(w.map(inv),v);
 if(i<2)console.log("cubic pair",v.map(String),"order",w[0].length+1);
 h=v[0];
}
console.log("square root 5 ->", [1n,3n].map(t=>String(H(child(5n,t)))));
console.log("exception 37 ->", [1n,3n,5n].map(t=>String(H(child(37n,t)))));
console.log("checks passed");
