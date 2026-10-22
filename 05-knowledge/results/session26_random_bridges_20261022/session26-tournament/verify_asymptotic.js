'use strict';
function hstate(w){let h=1,a=0,p=0;for(const b of w){if(b==='1'){a=2*a+1;p=1;}else{h+=a+p;a*=2;p=0;}}return h;}
function F(z){if(z<.25)return 0;if(z<=.5)return 4*z*Math.log(4*z)-4*z+1;if(z<=1)return 4*z-4*z*Math.log(z)-3;return 1;}
let result=[];
for(let m=3;m<=18;m++){let N=2**(m-2),sum=0,bad=0,c=Array(5).fill(0),points=[.25,.3,.5,.75,1],maxError=0;for(let k=0;k<N;k++){let w='1'+k.toString(2).padStart(m-2,'0')+'0';let h=hstate(w);sum+=h;let A=parseInt(w,2),B=parseInt([...w].reverse().map(b=>b==='0'?'1':'0').join(''),2);let err=Math.abs(h-A*B/2**m);if(err>m/2+1)bad++;maxError=Math.max(maxError,err);for(let j=0;j<points.length;j++)if(h/2**m<=points[j])c[j]++;}if(bad||sum/N!==9*2**(m-4)+.5)throw Error('asym '+m);result.push({m,mean:sum/N,maxError,cdfSupAtFivePoints:Math.max(...c.map((n,j)=>Math.abs(n/N-F(points[j]))))});}
console.log(JSON.stringify({status:'PASS',exactClaims:'mean and dyadic product error; all integer operations within IEEE exact range',diagnosticOnly:'five-point CDF differences use floating logarithms',result},null,2));
