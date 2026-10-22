// Independent exhaustive behavioral-minimization audit for Session 26.
function gcd(a,b){while(b){let c=a%b;a=b;b=c;}return Math.abs(a)}
function check(m,a,b,r,obs=[1,0]){
const key=(x,y,j)=>`${x},${y},${j}`;
let states=[[0,0,0]],ids=new Map([[key(0,0,0),0]]),trans=[];
for(let k=0;k<states.length;k++){let s=states[k];if(s===null){trans[k]=[k,k];continue}let [x,y,j]=s;trans[k]=[];
for(let d=0;d<2;d++){let z=d&&j===r?null:[(y+d*a)%m,(x+y+d*b)%m,d?j+1:0];let ky=z?key(...z):"sink";
if(!ids.has(ky)){ids.set(ky,states.length);states.push(z)}trans[k].push(ids.get(ky))}}
let part=states.map(s=>s?`${(obs[0]*s[0]+obs[1]*s[1])%m}`:"INVALID");
while(true){let map=new Map(),p=states.map((s,i)=>{let k=JSON.stringify([part[i],...trans[i].map(j=>part[j])]);if(!map.has(k))map.set(k,map.size);return map.get(k)});let n=new Set(p).size;if(n===new Set(part).size){part=p;break}part=p}
let h=states.filter(s=>s&&s[2]===0), k=h.filter(s=>{
let [x,y]=s;for(let j=0;j<2;j++){if((obs[0]*x+obs[1]*y)%m)return false;[x,y]=[y,(x+y)%m]}return true}).length;
let pred=(r+1)*h.length/k+1;
return {states:states.length,min:new Set(part).size,pred,h:h.length,k,smith:a||b?(m/gcd(m,gcd(a,b)))*(m/gcd(m,Math.abs(a*a+a*b-b*b)/gcd(a,b))):1}
}
let count=0,fail=[];
for(let m=1;m<=10;m++)for(let a=0;a<m;a++)for(let b=0;b<m;b++)for(let r=1;r<=4;r++)for(let o of [[1,0],[0,0],[1,1],[2,0]]){
let t=check(m,a,b,r,o);count++;if(t.min!==t.pred||t.h!==t.smith)fail.push({m,a,b,r,o,...t})}
console.log(JSON.stringify({count,fail},null,2));
if(fail.length)process.exitCode=1;

