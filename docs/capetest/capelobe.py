import numpy as np
rng=np.random.default_rng(2024)
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
    """3-D distance and effective sine of grid points pts (complex) from AP ap (complex)."""
    dv=pts-ap; dH=np.abs(dv)
    d3=np.sqrt(HD**2+dH**2); d3=np.maximum(d3,10.0)
    s=np.sin(np.angle(dv))*(dH/d3)
    return d3,s

def atoms(pos,d3,s):
    """NUSW unit-norm atoms, rows = points, cols = antennas (matches array_stats)."""
    d3=d3[:,None]; s=s[:,None]; x=pos[None,:]
    rn=np.sqrt(d3**2+x**2-2*d3*x*s)
    A=(d3/rn)*np.exp(-1j*2*np.pi*rn/LAM)
    A/=np.linalg.norm(A,axis=1,keepdims=True)
    return A

# fixed deployment
APs=(rng.random(L)+1j*rng.random(L))*SQ
g=np.arange(1.0,SQ,2.0); GX,GY=np.meshgrid(g,g); GRID=(GX+1j*GY).ravel()
print(f"deployment: L={L} APs, grid {len(GRID)} pts @2 m, N={N}/AP\n")

GEOMS={"dense ULA (D=1.15m)":rescale(np.arange(N)*d0,1.15),
       "nested   (D=7.75m)":rescale(nested(N)*d0,7.75),
       "end-load 1/4-1/2-1/4 (D=7.75m)":endload(N,0.5,7.75)}

def psl(pos):
    s=np.linspace(-1,1,20001)
    B=np.abs(np.exp(1j*2*np.pi*np.outer(s,pos)/LAM).sum(1))/len(pos)
    i0=np.argmax(B); i=i0
    while i+1<len(B) and B[i+1]<B[i]: i+=1
    j=i0
    while j-1>=0 and B[j-1]<B[j]: j-=1
    m=B.copy(); m[j:i+1]=0
    return m.max()

NTR=300; THR=5.0
for gname,pos in GEOMS.items():
    # precompute grid atoms per AP for this geometry
    AG=[];
    for l in range(L):
        d3,s=geo(APs[l],GRID); AG.append(atoms(pos,d3,s))
    print(f"--- {gname}   peak sidelobe={psl(pos):.3f} ---")
    for snr_db in (0.0,10.0):
        for nap in (1,2,4,8):
            ok=0; errs=[]
            for t in range(NTR):
                u=(rng.random()+1j*rng.random())*SQ
                d3u,su=geo(APs,np.array([u])); d3u=d3u.ravel(); su=su.ravel()
                beta=(DREF/d3u)**PLE
                order=np.argsort(-beta); use=order[:nap]
                J=np.zeros(len(GRID))
                for l in use:
                    a=atoms(pos,np.array([d3u[l]]),np.array([su[l]]))[0]
                    amp=np.sqrt(beta[l]*10**(snr_db/10)*N)
                    z=amp*a+(rng.normal(size=N)+1j*rng.normal(size=N))/np.sqrt(2)
                    J+=np.abs(AG[l]@z.conj())**2
                uh=GRID[np.argmax(J)]
                e=abs(uh-u); errs.append(e)
                if e<THR: ok+=1
            print(f"   SNR={snr_db:4.0f}dB  fuse {nap} AP: lobe-correct={100*ok/NTR:5.1f}%  median err={np.median(errs):7.2f} m")
    print()
