import numpy as np
LAM=0.1; d0=LAM/2; SQ=300.0; HD=10.0; L=8; N=24; PLE=3.76; DREF=50.0
def centre(p): p=np.asarray(p,float); return p-(p.max()+p.min())/2
def rescale(p,D): p=centre(p); return p*(D/(p.max()-p.min()))
def nested(N):
    N1=N//2;N2=N-N1
    return np.array(sorted(set(list(range(1,N1+1))+[(i+1)*(N1+1) for i in range(N2)])))[:N]
def endload(N,fm,D):
    nm=int(round(N*fm)); ne=(N-nm)//2; nm=N-2*ne
    left=-D/2+d0*np.arange(ne); right=D/2-d0*np.arange(ne)[::-1]
    mid=d0*(np.arange(nm)-(nm-1)/2)
    return np.sort(np.concatenate([left,mid,right]))
def geo(ap,pts):
    dv=pts-ap; dH=np.abs(dv); d3=np.maximum(np.sqrt(HD**2+dH**2),10.0)
    return d3,np.sin(np.angle(dv))*(dH/d3)
def atoms(pos,d3,s):
    d3=np.atleast_1d(d3)[:,None]; s=np.atleast_1d(s)[:,None]; x=pos[None,:]
    rn=np.sqrt(d3**2+x**2-2*d3*x*s)
    A=(d3/rn)*np.exp(-1j*2*np.pi*rn/LAM)
    return A/np.linalg.norm(A,axis=1,keepdims=True)
CAND={"nested D=7.75":rescale(nested(N)*d0,7.75),
      "dense ULA D=1.15":rescale(np.arange(N)*d0,1.15),
      "end-load 1/4-1/2-1/4 D=2.0":endload(N,0.5,2.0),
      "end-load 1/4-1/2-1/4 D=3.0":endload(N,0.5,3.0),
      "end-load 1/4-1/2-1/4 D=4.0":endload(N,0.5,4.0)}
g=np.arange(1.0,SQ,2.0); GX,GY=np.meshgrid(g,g); GRID=(GX+1j*GY).ravel()
NDEP=6; NTR=150; THR=5.0
res={k:{} for k in CAND}
for snr in (0.0,10.0,20.0):
    for k in CAND: res[k][snr]=[[],[]]
    for dep in range(NDEP):
        rng=np.random.default_rng(500+dep)
        APs=(rng.random(L)+1j*rng.random(L))*SQ
        for name,pos in CAND.items():
            AG=[atoms(pos,*geo(APs[l],GRID)) for l in range(L)]
            ok=0; errs=[]
            for t in range(NTR):
                u=(rng.random()+1j*rng.random())*SQ
                d3u,su=geo(APs,np.array([u])); d3u=d3u.ravel(); su=su.ravel()
                beta=(DREF/d3u)**PLE
                J=np.zeros(len(GRID))
                for l in range(L):
                    a=atoms(pos,d3u[l],su[l])[0]
                    amp=np.sqrt(beta[l]*10**(snr/10)*N)
                    z=amp*a+(rng.normal(size=N)+1j*rng.normal(size=N))/np.sqrt(2)
                    J+=np.abs(AG[l]@z.conj())**2
                e=abs(GRID[np.argmax(J)]-u); errs.append(e)
                if e<THR: ok+=1
            res[name][snr][0].append(100*ok/NTR); res[name][snr][1].append(np.median(errs))
print(f"{NDEP} independent AP deployments x {NTR} users each, 8-AP fusion\n")
for snr in (0.0,10.0,20.0):
    print(f"--- SNR = {snr:.0f} dB ---")
    print(f"{'array':30s} {'correct% (mean+-sd)':>22s} {'median err [m]':>16s}")
    for name in CAND:
        c=np.array(res[name][snr][0]); e=np.array(res[name][snr][1])
        print(f"{name:30s} {c.mean():10.1f} +-{c.std():5.1f}      {e.mean():8.2f} +-{e.std():5.2f}")
    print()
