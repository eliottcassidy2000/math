const need=(x,m)=>{if(!x)throw Error(m);};
const reps=n=>{let r=[];for(let a=0;a*a<=n/2;a++){let b=Math.floor(Math.sqrt(n-a*a));if(a*a+b*b===n)r.push([a,b]);}return r;};
const D=([a,b])=>[b-a,b+a],phi=([a,b])=>[b*b-a*a,2*a*b,a*a+b*b];
const eq=(a,b)=>JSON.stringify(a)===JSON.stringify(b);
const canon=a=>a.map(x=>x.join(",")).sort();
function sigma(n){let s=0;for(let d=1;d*d<=n;d++)if(n%d===0)s+=d+(d*d===n?0:n/d);return s;}
function eisen(n){const B=Math.floor(Math.sqrt(4*n/3))+1,r=[];for(let a=-B;a<=B;a++)for(let b=-B;b<=B;b++)if(a*a-a*b+b*b===n)r.push([a,b]);return r;}
for(let n=1;n<=10000;n++){let R=reps(n);need(eq(canon(R.map(D)),canon(reps(2*n))),"doubling bijection");for(let v of R){need(eq(D(D(v)),v.map(x=>2*x)),"D squared");let [A,B,C]=phi(v);need(eq(phi(D(v)),[2*B,2*A,2*C]),"Euclid intertwiner");}}
for(let m=1;m<=49;m+=2)for(let k=0;k<=6;k++){let n=2**k*m;need(sigma(n)-n===2**k*(2*sigma(m)-m)-sigma(m),"aliquot");}
for(let n=1;n<=500;n++){need(eq(canon(eisen(4*n)),canon(eisen(n).map(v=>v.map(x=>2*x)))),"Eisenstein quadrupling");if(n%2)need(eisen(2*n).length===0,"Eisenstein inert doubling");}
function splitProduct(m){let product=1;for(let p=2;p*p<=m;p++){let e=0;while(m%p===0){m/=p;e++;}if(p%4===1)product*=2*e+1;}if(m>1&&m%4===1)product*=3;return product;}
for(let m=1;m<=100;m++)for(let k=0;k<=4;k++)need(reps(2**k*m*m).length===(splitProduct(m)+1)/2,"classical count");
const rows=[25,50,100,200,169,338,289,578,9,18].map(n=>({n,representations:reps(n),aliquot:sigma(n)-n,eisensteinShell:eisen(n).length}));
console.log(JSON.stringify({status:"PASS",gaussianRange:10000,eisensteinRange:500,aliquotOddM:49,aliquotMaximumK:6,classicalSquareCountControls:500,rows,scope:"Exact integer finite enumeration, separate from universal elementary proofs."},null,2));