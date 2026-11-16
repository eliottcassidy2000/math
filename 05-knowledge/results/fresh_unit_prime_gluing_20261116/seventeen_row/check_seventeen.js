// Node.js, no dependencies. Run: node check_seventeen.js
function run(log){
  const check=(x,msg)=>{if(!x)throw Error(msg)};
  const pop=new Int8Array(2048);
  for(let s=1;s<pop.length;s++)pop[s]=pop[s>>1]+(s&1);
  const bit=b=>31-Math.clz32(b), pc=s=>pop[s];
  function prime(q){if(q<2)return false; for(let i=2;i*i<=q;i++)if(q%i===0)return false;return true;}
  function dp(a,n,start=-1){
    const F=(1<<n)-1,D=new Int32Array((F+1)*n),H=new Int32Array(F+1);
    H[0]=1;
    for(let i=0;i<n;i++)if(start<0||start===i)D[(1<<i)*n+i]=1;
    for(let S=1;S<=F;S++)for(let b=S;b;b&=b-1){
      let i=bit(b&-b),x=D[S*n+i]; if(!x)continue; H[S]+=x;
      for(let c=a[i]&~S;c;c&=c-1){const v=c&-c;D[(S|v)*n+bit(v)]+=x;}
    }
    return [H,D];
  }
  function strong(a,n){
    const F=(1<<n)-1;
    for(let dir=0;dir<2;dir++){
      let s=1,r=0;
      while(s!==r){r=s;for(let b=r;b;b&=b-1){let i=bit(b&-b);s|=dir?(F^a[i]^(1<<i)):a[i];}}
      if(s!==F)return false;
    }return true;
  }
  let order10=0,strong10=0, padded=[], min=Infinity;
  const rows={};
  function read(a,n){
    if(n===10)order10++;
    if(!strong(a,n))return;
    if(n===10)strong10++;
    const F=(1<<n)-1,[H]=dp(a,n),q=H[F];if(q<17||!prime(q))return;
    const W=[new Int32Array(n*n),new Int32Array(n*n)];
    for(let i=0;i<n;i++){
      const [,D]=dp(a,n,i);
      for(let S=1<<i;S<=F;S++)if(S>>i&1){
        let c=H[F^S], w=W[pc(S)&1];
        for(let b=S;b;b&=b-1){const j=bit(b&-b);w[i*n+j]+=c*D[S*n+j];}
      }
    }
    const tails=[new Int32Array(n),new Int32Array(n)], heads=tails.map(()=>new Int32Array(n));
    for(let p=0;p<2;p++)for(let i=0;i<n;i++)for(let j=0;j<n;j++){
      let w=W[p][i*n+j]; heads[p][i]+=w;tails[p][j]+=w;
    }
    for(let A=1;A<F;A++){
      let B=F^A,d=0;
      for(let b=A;b&&d<=8;b&=b-1){
        const i=bit(b&-b);
        for(let c=B;c;c&=c-1)d+=W[0][i*n+bit(c&-c)];
      }
      if(d>8)continue;check(d>0,"mixed contact must create triangle");
      min=Math.min(min,d);
      const row=rows[q]??={gaps:{},distinct:false,partners:0,cert:null};
      row.gaps[d]=(row.gaps[d]||0)+1;
      if(d<8){
        for(const [side,profiles] of [[B,tails],[A,heads]]){
          let E=0,O=0;
          for(let b=side;b;b&=b-1){let i=bit(b&-b);E+=profiles[0][i];O+=profiles[1][i];}
          check(O>0,"singleton contribution positive");
          const z=(8-d+E)/(O+E);
          if(Number.isInteger(z)&&z>=1&&(z&(z-1))===0)padded.push([q,d,A,E,O]);
        }
      }else{
        const b=a.slice(0,n);
        for(let z=B;z;z&=z-1)b[bit(z&-z)]|=1<<n;
        b.push(A);
        const [Q]=dp(b,n+1);
        check(Q[(F<<1)|1]===q+16,"direct DP mismatch");
        const S=Q.findIndex((x,i)=>x===17&&i!==F);
        if(S>=0){row.partners++;row.distinct=true;if(!row.cert)row.cert=[b,S];}
      }
    }
  }
  function grow(a,tau){
    const n=a.length;
    if(n>=3)read(a,n);
    if(n===10)return;
    const F=(1<<n)-1;
    for(let A=0;A<(1<<(n-1));A++){
      const B=F^A;
      let t=tau;
      for(let b=A;b&&t<=8;b&=b-1)t+=pc(a[bit(b&-b)]&B);
      if(t>8)continue;
      const b=a.slice();
      for(let c=B;c;c&=c-1)b[bit(c&-c)]|=1<<n;
      b.push(A);
      grow(b,t);
    }
  }
  grow([0],0);
  check(order10===100657&&strong10===10300,"order ten counts");
  check(min===5&&padded.length===0,"strict-step or exterior check");
  const expected={17:[{5:51,6:68,7:34,8:119},85],19:[{6:38,7:76,8:38},38],23:[{7:46,8:46},46],29:[{8:58},58]};
  check(Object.keys(rows).length===4,"prime list");
  for(const [q,[gaps,partners]] of Object.entries(expected)){
    const r=rows[q];check(r&&JSON.stringify(r.gaps)===JSON.stringify(gaps)&&partners===r.partners,"table "+q);
  }
  for(const [q,a,S,h] of [
    [17,[94,124,120,112,96,1,32],111,33],
    [19,[30,124,24,48,32,5,61],111,35],
    [23,[94,28,56,48,32,3,62],125,39],
    [29,[190,252,248,16,32,72,153,56],123,45]
  ]){
    const [H]=dp(a,a.length),F=(1<<(a.length-1))-1;
    check(H[S]===17&&S!==F&&H[F]===q&&H[(F<<1)|1]===h,"literal witness");
  }
  const answer={order10,strong10,min,padded,rows};log(answer);return answer;
}
run(x=>console.log(JSON.stringify(x,null,2)));
