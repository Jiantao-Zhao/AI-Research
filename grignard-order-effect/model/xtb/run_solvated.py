import json, os, re, subprocess
XTB="/opt/anaconda3/envs/qc/bin/xtb"; HERE=os.path.dirname(os.path.abspath(__file__))
SP=["ch3mgcl_2thf","mgme2_2thf","mgcl2_2thf","mgcl2_4thf","dimer_2thf","thf_free"]
def run(args,tag):
    p=subprocess.run(args,cwd=HERE,capture_output=True,text=True)
    open(os.path.join(HERE,tag+".log"),"w").write(p.stdout+"\n"+p.stderr)
    return p.returncode, p.stdout+p.stderr
def parse(t):
    g=re.search(r"TOTAL FREE ENERGY\s+(-?\d+\.\d+)",t); e=re.search(r"TOTAL ENERGY\s+(-?\d+\.\d+)",t)
    im=re.search(r"# imaginary freq\.\s+(\d+)",t)
    return (float(g.group(1)) if g else None, float(e.group(1)) if e else None, int(im.group(1)) if im else 0)
res={}
for sp in SP:
    opt=os.path.join(HERE,"xtbopt.xyz")
    if os.path.exists(opt): os.remove(opt)
    rc,out=run([XTB,sp+".xyz","--opt","tight","--alpb","thf"],sp+"_opt")
    if rc!=0 or not os.path.exists(opt):
        res[sp]={"error":f"opt rc={rc}","tail":out[-300:]}; print(sp,"OPT FAILED"); continue
    os.replace(opt,os.path.join(HERE,sp+"_opt.xyz"))
    rc2,out2=run([XTB,sp+"_opt.xyz","--ohess","--alpb","thf"],sp+"_freq")
    g,e,n=parse(out2); res[sp]={"G":g,"E":e,"n_imag":n,"rc":rc2}
    print(f"{sp:16s} G={g} E={e} imag={n} rc={rc2}",flush=True)
json.dump(res,open(os.path.join(HERE,"schlenk_solvated.json"),"w"),indent=2)
