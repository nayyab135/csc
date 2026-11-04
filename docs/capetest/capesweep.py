import numpy as np
rng=np.random.default_rng(99)
LAM=0.1; d0=LAM/2; SQ=300.0; HD=10.0; L=8; N=24; PLE=3.76; DREF=50.0
def centre(p): p=np.asarray(p,float); return p-(p.max()+p.min())/2
def rescale(p,D): p=centre(p); return p*(D/(p.max()-p.min()))
def nested(N):
    N1=N//2;N2=N-N1
    return np.array(sorted(set(list(range(1,N1+1))+[(i+1)*(N1+1) for i in range(N2)])))[:N]
def endload(N,fm,D):
    nm=int(round(N*fm)); ne=(N-nm)//2; nm=N-2*ne
    if ne==0: return rescale(np.arange(N)*d0,min(D,(N-1)*d0))
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
def psl(pos):
    s=np.linspace(-1,1,20001)
    B=np.abs(np.exp(1j*2*np.pi*np.outer(s,pos)/LAM).sum(1))/len(pos)
    i0=np.argmax(B); i=i0
    while i+1<len(B) and B[i+1]<B[i]: i+=1
    j=i0
    while j-1>=0 and B[j-1]<B[j]: j-=1
    m=B.copy(); m[j:i+1]=0; return m.max()
APs=(rng.random(L)+1j*rng.random(L))*SQ
g=np.arange(1.0,SQ,2.0); GX,GY=np.meshgrid(g,g); GRID=(GX+1j*GY).ravel()
CAND={"nested D=7.75":rescale(nested(N)*d0,7.75),
      "dense ULA D=1.15":rescale(np.arange(N)*d0,1.15)}
for fm in (0.50,0.65,0.80,0.90):
    CAND[f"end-load mid={fm:.2f} D=7.75"]=endload(N,fm,7.75)
for D in (3.0,5.0):
    CAND[f"end-load mid=0.50 D={D:.1f}"]=endload(N,0.50,D)
NTR=400; THR=5.0; SNR=10.0
print(f"8-AP CAPE fusion, SNR={SNR:.0f} dB, {NTR} trials, grid 2 m, threshold {THR} m\n")
print(f"{'array':28s} {'PSL':>6s} {'correct%':>9s} {'median err':>11s} {'err|correct':>12s}")
for name,pos in CAND.items():
    AG=[atoms(pos,*geo(APs[l],GRID)) for l in range(L)]
    ok=0; errs=[]; cond=[]
    for t in range(NTR):
        u=(rng.random()+1j*rng.random())*SQ
        d3u,su=geo(APs,np.array([u])); d3u=d3u.ravel(); su=su.ravel()
        beta=(DREF/d3u)**PLE
        J=np.zeros(len(GRID))
        for l in range(L):
            a=atoms(pos,d3u[l],su[l])[0]
            amp=np.sqrt(beta[l]*10**(SNR/10)*N)
            z=amp*a+(rng.normal(size=N)+1j*rng.normal(size=N))/np.sqrt(2)
            J+=np.abs(AG[l]@z.conj())**2
        e=abs(GRID[np.argmax(J)]-u); errs.append(e)
        if e<THR: ok+=1; cond.append(e)
    print(f"{name:28s} {psl(pos):6.3f} {100*ok/NTR:8.1f}% {np.median(errs):10.2f} m {np.median(cond) if cond else float('nan'):10.2f} m")
