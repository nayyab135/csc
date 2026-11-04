import numpy as np
LAM=0.1; d0=LAM/2; SQ=300.0; HD=10.0; L=8; N=24; PLE=3.76; DREF=50.0
def centre(p): p=np.asarray(p,float); return p-(p.max()+p.min())/2
def rescale(p,D): p=centre(p); return p*(D/(p.max()-p.min()))
def nested(N):
    N1=N//2;N2=N-N1
    return np.array(sorted(set(list(range(1,N1+1))+[(i+1)*(N1+1) for i in range(N2)])))[:N]
def endload(N,fm,D):
    nm=int(round(N*fm)); ne=(N-nm)//2; nm=N-2*ne
    return np.sort(np.concatenate([-D/2+d0*np.arange(ne),d0*(np.arange(nm)-(nm-1)/2),D/2-d0*np.arange(ne)[::-1]]))
def hier(N,ncore,D):
    core=d0*(np.arange(ncore)-(ncore-1)/2); no=N-ncore; half=no//2
    return np.sort(np.concatenate([core,np.linspace(-D/2,-D/6,half),np.linspace(D/6,D/2,no-half)]))
def powerlaw(N,D,p):
    M=N//2; h=((np.arange(1,M+1)/M)**p)*D/2
    return np.sort(np.concatenate([-h,h]))[:N]
def geo(ap,pts):
    dv=pts-ap; dH=np.abs(dv); d3=np.maximum(np.sqrt(HD**2+dH**2),10.0)
    return d3,np.sin(np.angle(dv))*(dH/d3)
def atoms(pos,d3,s):
    d3=np.atleast_1d(d3)[:,None]; s=np.atleast_1d(s)[:,None]; x=pos[None,:]
    rn=np.sqrt(d3**2+x**2-2*d3*x*s)
    A=(d3/rn)*np.exp(-1j*2*np.pi*rn/LAM)
    return A/np.linalg.norm(A,axis=1,keepdims=True)
C={"dense ULA D=1.15":rescale(np.arange(N)*d0,1.15),
   "nested D=7.75":rescale(nested(N)*d0,7.75),
   "end-load D=2.0":endload(N,0.5,2.0),
   "R1 hier core16 D=3":hier(N,16,3.0),
   "R1 hier core16 D=6":hier(N,16,6.0),
   "R2 power p=2 D=3":powerlaw(N,3.0,2.0),
   "R2 power p=2 D=6":powerlaw(N,6.0,2.0)}
g=np.arange(1.0,SQ,2.0); GX,GY=np.meshgrid(g,g); GRID=(GX+1j*GY).ravel()
NDEP=5; NTR=150; THR=5.0
for snr in (10.0,20.0):
    print(f"=== SNR={snr:.0f} dB | {NDEP} deployments x {NTR} users | 8-AP CAPE fusion ===")
    print(f"{'geometry':22s} {'correct%':>9s} {'loc err':>9s} {'NMSE':>8s} {'NMSE>0.5':>9s}")
    out={}
    for name,pos in C.items():
        OK=[];E=[];NM=[];BAD=[]
        for dep in range(NDEP):
            rng=np.random.default_rng(1300+dep)
            APs=(rng.random(L)+1j*rng.random(L))*SQ
            AG=[atoms(pos,*geo(APs[l],GRID)) for l in range(L)]
            ok=0
            for t in range(NTR):
                u=(rng.random()+1j*rng.random())*SQ
                d3u,su=geo(APs,np.array([u])); d3u=d3u.ravel(); su=su.ravel()
                beta=(DREF/d3u)**PLE
                J=np.zeros(len(GRID)); At=[]
                for l in range(L):
                    a=atoms(pos,d3u[l],su[l])[0]; At.append(a)
                    amp=np.sqrt(beta[l]*10**(snr/10)*N)
                    z=amp*a+(rng.normal(size=N)+1j*rng.normal(size=N))/np.sqrt(2)
                    J+=np.abs(AG[l]@z.conj())**2
                uh=GRID[np.argmax(J)]; e=abs(uh-u); E.append(e)
                if e<THR: ok+=1
                d3h,sh=geo(APs,np.array([uh])); d3h=d3h.ravel(); sh=sh.ravel()
                num=den=0.0
                for l in range(L):
                    ah=atoms(pos,d3h[l],sh[l])[0]
                    num+=beta[l]*(1-abs(np.vdot(ah,At[l]))**2); den+=beta[l]
                NM.append(num/den); BAD.append(num/den>0.5)
            OK.append(100*ok/NTR)
        print(f"{name:22s} {np.mean(OK):8.1f}% {np.median(E):7.2f} m {np.mean(NM):8.3f} {100*np.mean(BAD):8.1f}%")
    print()
