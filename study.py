#!/usr/bin/env python3
"""CIR error audit. Synthetic experiments; no market-data calibration.
Run: python study.py --reps 80 --paths 8000 --output-dir results/full
Quick smoke run: python study.py --quick --output-dir results/quick
Implementation checks only: python study.py --verify-only
Existing nonempty output directories are never overwritten.
Seeds isolate scenario/scheme/grid. Standard errors use independent paths.
"""
from __future__ import annotations
import argparse, csv, json, math, platform, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':16, 'axes.labelsize':16,
                     'xtick.labelsize':16, 'ytick.labelsize':16,
                     'axes.titlesize':17, 'legend.fontsize':15.2})

ROOT=Path(__file__).resolve().parent
OUTPUT_NAMES=('replicates.csv','summary.csv','analytic_grid.csv','results.json',
              'cancellation_counterexamples.csv','error_audit.png',
              'price_greek_cancellation.png')
SCENARIOS={
    'regular':dict(kappa=.6,theta=.04,sigma=.1,r0=.03,T=5.),
    'boundary':dict(kappa=.4,theta=.02,sigma=.35,r0=.005,T=5.),
    'severe':dict(kappa=.3,theta=.01,sigma=.6,r0=.001,T=5.)}
SCHEMES=['exact','qe','full_truncation']

def continuous(p, load=1.):
    k,t,s,r,T=(p[x] for x in ['kappa','theta','sigma','r0','T'])
    g=math.sqrt(k*k+2*s*s*load); e=math.expm1(g*T)
    den=(g+k)*e+2*g
    B=2*load*e/den
    logA=(2*k*t/s**2)*(math.log(2*g)+(k+g)*T/2-math.log(den))
    P=math.exp(logA-B*r)
    return P,B,logA

def grid_affine(p,steps,load=1.):
    k,t,s,r,T=(p[x] for x in ['kappa','theta','sigma','r0','T'])
    h=T/steps;q=math.exp(-k*h);c=s*s*(-math.expm1(-k*h))/(4*k)
    d=4*k*t/s**2;a=b=0.
    for _ in range(steps):
        u=b+load*h/2
        a += d/2*math.log1p(2*c*u)
        b = load*h/2+q*u/(1+2*c*u)
    return math.exp(-a-b*r), b, -a

def qe_step(x,p,h,rng):
    k,t,s=p['kappa'],p['theta'],p['sigma'];q=math.exp(-k*h)
    m=t+(x-t)*q
    v=x*s*s*q*(1-q)/k+t*s*s*(1-q)**2/(2*k)
    psi=v/(m*m)
    low=psi<=1.5
    y=np.empty_like(x)
    ps=psi[low];two=2/ps
    b2=two-1+np.sqrt(two*(two-1))
    a=m[low]/(1+b2)
    y[low]=a*(np.sqrt(b2)+rng.standard_normal(a.size))**2
    ps=psi[~low];pointmass=(ps-1)/(ps+1)
    beta=(1-pointmass)/m[~low];u=rng.random(ps.size)
    # Positive continuous exponential component only when U>p.
    yy=np.zeros(ps.size);pos=u>pointmass
    yy[pos]=np.log((1-pointmass[pos])/(1-u[pos]))/beta[pos]
    y[~low]=yy
    return y

def experiment(p,scheme,steps,reps,paths,seed):
    if scheme not in SCHEMES:
        raise ValueError(f'Unknown scheme: {scheme}')
    rng=np.random.default_rng(seed); h=p['T']/steps
    x=np.full((reps,paths),p['r0']); integral=np.zeros_like(x)
    neg=np.zeros(reps);near=np.zeros(reps)
    q=math.exp(-p['kappa']*h)
    c=p['sigma']**2*(1-q)/(4*p['kappa']);d=4*p['kappa']*p['theta']/p['sigma']**2
    start=time.perf_counter()
    for _ in range(steps):
        xp=np.maximum(x,0)
        if scheme=='exact': y=c*rng.noncentral_chisquare(d,q*x/c)
        elif scheme=='qe': y=qe_step(x,p,h,rng)
        else: y=x+p['kappa']*(p['theta']-xp)*h+p['sigma']*np.sqrt(xp*h)*rng.standard_normal(x.shape)
        yp=np.maximum(y,0)
        integral += .5*h*(xp+yp)
        neg += (y<0).mean(axis=1)
        near += (yp<1e-6).mean(axis=1)
        x=y
    values=np.exp(-integral)
    means=values.mean(axis=1); ses=values.std(axis=1,ddof=1)/math.sqrt(paths)
    return means,ses,neg/steps,near/steps,time.perf_counter()-start

def require(condition, message):
    if not condition:
        raise AssertionError(message)

def verify():
    checks={}
    for key,p in SCENARIOS.items():
        P,B,_=continuous(p)
        fine,bfine,_=grid_affine(p,10000)
        require(abs(fine-P)<1e-8, f'{key}: fine-grid price mismatch')
        require(abs(bfine-B)<1e-7, f'{key}: fine-grid sensitivity mismatch')
        for n in [1,5,20]:
            h=p['T']/n
            G,b,_=grid_affine(p,n)
            eps=1e-6
            up=dict(p,r0=p['r0']+eps);dn=dict(p,r0=p['r0']-eps)
            numerical=(grid_affine(up,n)[0]-grid_affine(dn,n)[0])/(2*eps)
            require(abs(numerical+b*G)<1e-8, f'{key}: finite-difference mismatch')
        checks[key]={'continuous_price':P,'fine_grid_abs_error':abs(fine-P)}
    # One-step exact grid closed form is the transform of noncentral chi-square.
    p=SCENARIOS['boundary'];n=10;G,_,_=grid_affine(p,n);G2,_,_=grid_affine(p,n,2)
    m,se,_,_,_=experiment(p,'exact',n,20,10000,91919)
    z=(m.mean()-G)/(m.std(ddof=1)/math.sqrt(len(m)))
    require(abs(z)<5, 'Sample mean differs by more than 5 standard errors')
    ratio=np.mean(se**2)*10000/(G2-G*G)
    require(.90<ratio<1.10, 'Exact-sampler variance check failed')
    checks['exact_sampler']={'pooled_z_to_grid_target':float(z),'estimated_to_exact_variance_ratio':float(ratio)}
    # Check QE conditional moments, for both quadratic and exponential branches.
    checks['qe_moments']=[]
    for x0 in [.0001,.03,.5]:
        rng=np.random.default_rng(19000+round(x0*10000));size=300000;h=.25
        draws=qe_step(np.full(size,x0),p,h,rng);q=math.exp(-p['kappa']*h)
        m=p['theta']+(x0-p['theta'])*q
        v=x0*p['sigma']**2*q*(1-q)/p['kappa']+p['theta']*p['sigma']**2*(1-q)**2/(2*p['kappa'])
        z=(draws.mean()-m)/math.sqrt(v/size)
        require(abs(z)<5, 'Sample mean differs by more than 5 standard errors')
        require(abs(draws.var()/v-1)<.025, 'QE conditional variance check failed')
        checks['qe_moments'].append({'x':x0,'mean_z':float(z),'variance_ratio':float(draws.var()/v)})
    return checks

def save_csv(output_dir,name,rows):
    with (output_dir/name).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def extra_audit(output_dir):
    """Deterministic price/Greek cancellation examples and Wilson intervals."""
    rows=[]
    for name,p in SCENARIOS.items():
        P,B,logA=continuous(p)
        for spy in [1,2,4,12,48]:
            steps=round(p['T']*spy);G,b,loga=grid_affine(p,steps)
            rstar=(loga-logA)/(b-B)
            pp=dict(p,r0=rstar)
            Pc,Bc,_=continuous(pp);Pg,Bg,_=grid_affine(pp,steps)
            require(rstar>0 and abs(Pg-Pc)<1e-13, 'Cancellation check failed')
            rows.append(dict(scenario=name,steps_per_year=spy,rstar=rstar,continuous_price=Pc,grid_price=Pg,
                price_abs_error=abs(Pg-Pc),relative_sensitivity_error=Bg/Bc-1,
                continuous_initial_rate_sensitivity=Bc*Pc,grid_initial_rate_sensitivity=Bg*Pg))
    save_csv(output_dir,'cancellation_counterexamples.csv',rows)
    result=json.loads((output_dir/'results.json').read_text())
    for row in result['summary']:
        n=row['reps'];phat=row['coverage'];z=1.96;den=1+z*z/n
        center=(phat+z*z/(2*n))/den
        half=z*math.sqrt(phat*(1-phat)/n+z*z/(4*n*n))/den
        row['coverage_wilson_lower']=center-half
        row['coverage_wilson_upper']=center+half
    result['cancellation_counterexamples']=rows
    (output_dir/'results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
    save_csv(output_dir,'summary.csv',result['summary'])
    p=SCENARIOS['regular'];steps=5
    x=np.linspace(.02,.06,161);price=[];sensitivity=[]
    for r0 in x:
        pp=dict(p,r0=float(r0));P,B,_=continuous(pp);G,b,_=grid_affine(pp,steps)
        price.append(1e4*(G-P));sensitivity.append(100*(b*G/(B*P)-1))
    rstar=rows[0]['rstar']
    fig,axes=plt.subplots(1,2,figsize=(10,5.2))
    axes[0].plot(100*x,price,color='#086d8f');axes[0].set_ylabel('Price error (bp of unit notional)')
    axes[1].plot(100*x,sensitivity,color='#bc496a');axes[1].set_ylabel('Initial-rate sensitivity error (%)')
    for ax in axes:
        ax.axhline(0,color='#777777',lw=.7);ax.axvline(100*rstar,color='#777777',ls='--',lw=1)
        ax.set_xlabel('Initial short rate (%)');ax.grid(alpha=.15)
    fig.suptitle('A correct price can hide a wrong sensitivity\nExact CIR grid, annual steps',fontsize=18)
    fig.tight_layout(rect=[0,0,1,.91]);fig.savefig(output_dir/'price_greek_cancellation.png',dpi=220);plt.close(fig)

def prepare_output_dir(path):
    """Only an empty destination is allowed; archived study files stay untouched."""
    output_dir=Path(path).expanduser().resolve()
    if output_dir==ROOT or (output_dir.exists() and any(output_dir.iterdir())):
        raise ValueError('Output directory must be empty and cannot be the repository root. '
                         'Choose a new --output-dir; existing results are never overwritten.')
    output_dir.mkdir(parents=True,exist_ok=True)
    return output_dir


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--reps',type=int,default=80,help='independent batches (default: 80)')
    parser.add_argument('--paths',type=int,default=8000,help='paths per batch (default: 8000)')
    parser.add_argument('--output-dir',type=Path,default=Path('results/run'),
                        help='empty output directory (default: results/run)')
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--quick',action='store_true',
                      help='smoke run: 4 batches x 512 paths; overrides --reps/--paths')
    mode.add_argument('--verify-only',action='store_true',
                      help='run implementation checks, print JSON, and write no study files')
    args=parser.parse_args(argv)
    if args.quick:
        args.reps,args.paths=4,512
    if args.reps<2 or args.paths<2:
        parser.error('--reps and --paths must both be at least 2.')
    if args.verify_only:
        print(json.dumps(dict(status='passed',checks=verify(),environment=dict(
            python=platform.python_version(),numpy=np.__version__,
            matplotlib=matplotlib.__version__)),indent=2))
        return
    try:
        output_dir=prepare_output_dir(args.output_dir)
    except (ValueError,OSError) as error:
        parser.error(str(error))
    checks=verify()
    raw=[];summary=[];theory=[]
    for si,(name,p) in enumerate(SCENARIOS.items()):
        P,B,_=continuous(p);dv=1e-4*B*P
        for spy in [2,4,12,48]:
            steps=round(p['T']*spy);h=1/spy;G,b,_=grid_affine(p,steps);G2,_,_=grid_affine(p,steps,2)
            theory.append(dict(scenario=name,steps_per_year=spy,h=h,continuous_price=P,exact_grid_price=G,quadrature_bias_bp=(G-P)*1e4,
                continuous_r0_dv01=dv,exact_grid_r0_dv01=1e-4*b*G,relative_price_error=(G/P-1),relative_r0_dv01_error=b*G/(B*P)-1,
                grid_path_variance=G2-G*G,critical_paths_196=1.96**2*(G2-G*G)/(G-P)**2))
            for mi,scheme in enumerate(SCHEMES):
                means,ses,neg,near,seconds=experiment(p,scheme,steps,args.reps,args.paths,20261011+si*10000+mi*100+spy)
                covered=np.abs(means-P)<=1.96*ses
                grid_covered=np.abs(means-G)<=1.96*ses
                for rep in range(args.reps):
                    raw.append(dict(scenario=name,scheme=scheme,steps_per_year=spy,rep=rep,paths=args.paths,estimate=float(means[rep]),
                        standard_error=float(ses[rep]),covered_continuous=int(covered[rep]),covered_exact_grid=int(grid_covered[rep]),
                        negative_state_fraction=float(neg[rep]),near_zero_fraction=float(near[rep])))
                coverage=float(covered.mean())
                row=dict(scenario=name,scheme=scheme,steps_per_year=spy,paths=args.paths,reps=args.reps,mean_estimate=float(means.mean()),
                    total_bias_bp=float((means.mean()-P)*1e4),bias_mcse_bp=float(means.std(ddof=1)/math.sqrt(args.reps)*1e4),
                    quadrature_bias_bp=(G-P)*1e4,state_bias_bp=float((means.mean()-G)*1e4),
                    coverage=coverage,coverage_mcse=math.sqrt(coverage*(1-coverage)/args.reps),grid_target_coverage=float(grid_covered.mean()),
                    mean_95_half_width_bp=float(1.96*ses.mean()*1e4),negative_state_fraction=float(neg.mean()),near_zero_fraction=float(near.mean()),seconds=seconds)
                summary.append(row);print(name,scheme,spy,'bias bp',round(row['total_bias_bp'],4),'coverage',coverage,flush=True)
    save_csv(output_dir,'replicates.csv',raw);save_csv(output_dir,'summary.csv',summary);save_csv(output_dir,'analytic_grid.csv',theory)
    result=dict(run_mode='quick' if args.quick else 'full',seed=20261011,paths_per_rep=args.paths,reps=args.reps,scenario_parameters=SCENARIOS,checks=checks,summary=summary,analytic_grid=theory,
        environment=dict(python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__),
        caveats=['Synthetic parameters; no market calibration.','Coverage is for normal-approximation MC intervals with fixed parameters.',
        'r0 DV01 is a short-rate initial-state sensitivity, not a parallel yield-curve DV01.','State bias estimate includes finite Monte Carlo uncertainty.'])
    (output_dir/'results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
    fig,axes=plt.subplots(2,3,figsize=(12,8.6),sharex='col')
    colors={'exact':'#086d8f','qe':'#dda137','full_truncation':'#bc496a'}
    for col,name in enumerate(SCENARIOS):
        for scheme in SCHEMES:
            rows=[r for r in summary if r['scenario']==name and r['scheme']==scheme]
            x=[r['steps_per_year'] for r in rows]
            axes[0,col].errorbar(x,[r['total_bias_bp'] for r in rows],yerr=[1.96*r['bias_mcse_bp'] for r in rows],marker='o',label=scheme,color=colors[scheme])
            axes[1,col].plot(x,[100*r['coverage'] for r in rows],'-o',color=colors[scheme])
        axes[0,col].axhline(0,color='#444444',lw=.7);axes[0,col].set_title(name)
        axes[1,col].axhline(95,color='#777777',ls='--',lw=1);axes[1,col].set_ylim(-3,103)
        axes[1,col].set_xlabel('Time steps per year');axes[1,col].set_xscale('log');axes[1,col].set_xticks([2,4,12,48],[2,4,12,48])
        axes[1,col].xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        for ax in axes[:,col]:ax.grid(alpha=.16)
    axes[0,0].set_ylabel('Price bias (bp)');axes[1,0].set_ylabel('Coverage (%)')
    handles,labels=axes[0,0].get_legend_handles_labels()
    label_map={'exact':'Exact','qe':'QE','full_truncation':'Full truncation'}
    axes[0,0].legend(handles,[label_map[s] for s in labels],fontsize=15.2)
    fig.suptitle(f'CIR numerical error audit | {args.reps} independent batches',fontsize=19)
    fig.text(.5,.018,f'Bias bars: 95% intervals for estimated bias across batches.\nEach batch: {args.paths:,} paths; dashed line: nominal 95% MC coverage.',ha='center',fontsize=15.2)
    fig.tight_layout(rect=[0,.088,1,.96]);fig.savefig(output_dir/'error_audit.png',dpi=220);plt.close(fig)
    extra_audit(output_dir)
    print(f'Saved results, CSV and figures to {output_dir}',flush=True)

if __name__=='__main__':main()
