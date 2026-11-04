import numpy as np
LAM=0.1; d0=LAM/2; SQ=300.0; HD=10.0; L=8; N=24; PLE=3.76; DREF=50.0
def geo(ap,pts):
    dv=pts-ap; dH=np.abs(dv); d3=np.maximum(np.sqrt(HD**2+dH**2),10.0)
    return d3,np.sin(np.angle(dv))*(dH/d3)
def atoms(pos,d3,s):
    d3=np.atleast_1d(d3)[:,None]; s=np.atleast_1d(s)[:,None]; x=pos[None,:]
    rn=np.sqrt(d3**2+x**2-2*d3*x*s)
    A=(d3/rn)*np.exp(-1j*2*np.pi*rn/LAM)
    return A/np.linalg.norm(A,axis=1,keepdims=True)
def hier(N,ncore,D):
    core=d0*(np.arange(ncore)-(ncore-1)/2); no=N-ncore; half=no//2
    out=np.concatenate([np.linspace(-D/2,-D/6,half),np.linspace(D/6,D/2,no-half)])
    allp=np.concatenate([core,out]); o=np.argsort(allp)
    mask=np.concatenate([np.ones(ncore,bool),np.zeros(no,bool)])[o]
    return allp[o],mask
g=np.arange(1.0,SQ,2.0); GX,GY=np.meshgrid(g,g); GRID=(GX+1j*GY).ravel()
NDEP=5; NTR=150; THR=5.0; RAD=25.0
for D in (3.0,6.0,10.0):
    pos,cm=hier(N,16,D)
    print(f"--- R1 hierarchical, core=16 @lambda/2, outriggers to D={D} m ---")
    for snr in (10.0,20.0):
        res={"full-array global":[[],[]], "coarse-to-fine (core then full)":[[],[]]}
        for dep in range(NDEP):
            rng=np.random.default_rng(2200+dep)
            APs=(rng.random(L)+1j*rng.random(L))*SQ
            AGf=[atoms(pos,*geo(APs[l],GRID)) for l in range(L)]
            AGc=[A[:,cm] for A in AGf]
            okF=okC=0; eF=[]; eC=[]
            for t in range(NTR):
                u=(rng.random()+1j*rng.random())*SQ
                d3u,su=geo(APs,np.array([u])); d3u=d3u.ravel(); su=su.ravel()
                beta=(DREF/d3u)**PLE
                Z=[]
                for l in range(L):
                    a=atoms(pos,d3u[l],su[l])[0]
                    amp=np.sqrt(beta[l]*10**(snr/10)*N)
                    Z.append(amp*a+(rng.normal(size=N)+1j*rng.normal(size=N))/np.sqrt(2))
                Jf=np.zeros(len(GRID)); Jc=np.zeros(len(GRID))
                for l in range(L):
                    Jf+=np.abs(AGf[l]@Z[l].conj())**2
                    Jc+=np.abs(AGc[l]@Z[l][cm].conj())**2
                uf=GRID[np.argmax(Jf)]
                uc=GRID[np.argmax(Jc)]                       # coarse: core only (unambiguous)
                sel=np.abs(GRID-uc)<RAD                      # refine with FULL array near coarse
                u2=GRID[sel][np.argmax(Jf[sel])]
                e1=abs(uf-u); e2=abs(u2-u); eF.append(e1); eC.append(e2)
                okF+= e1<THR; okC+= e2<THR
            res["full-array global"][0].append(100*okF/NTR); res["full-array global"][1].append(np.median(eF))
            res["coarse-to-fine (core then full)"][0].append(100*okC/NTR); res["coarse-to-fine (core then full)"][1].append(np.median(eC))
        print(f"  SNR={snr:.0f} dB")
        for k,(o,e) in res.items():
            print(f"    {k:34s} correct={np.mean(o):5.1f}%  median err={np.mean(e):5.2f} m")
    print()
