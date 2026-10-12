#!/usr/bin/env python3
"""Reproduce figures/tables: python make_plots.py CASE_DIRECTORY OUTPUT_DIRECTORY"""
import sys,re,csv,hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
src=Path(sys.argv[1]); out=Path(sys.argv[2]); figdir=out/'figures'; tab=out/'tables'
figdir.mkdir(parents=True,exist_ok=True); tab.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.2,'savefig.dpi':220,'svg.fonttype':'none'})
def csvwrite(name,header,rows):
 with (tab/(name+'.csv')).open('w',newline='') as f:
  w=csv.writer(f);w.writerow(header);w.writerows(rows)
def readmon(region,name,kind='normal'):
 base=src/'postProcessing'/region/name; data={}
 for folder in sorted(base.iterdir(),key=lambda p:float(p.name)):
  start=float(folder.name); data={t:v for t,v in data.items() if t<=start}
  for f in sorted(folder.glob('*.dat')):
   for line in f.read_text().splitlines():
    if not line.strip() or line.lstrip().startswith('#'):continue
    p=line.split();t=float(p[0])
    if kind=='extreme':
     if p[1]!='T':continue
     q=re.sub(r'\([^)]*\)','LOCATION',line).split();v=[float(q[2]),float(q[5])]
    elif kind=='patch':
     if p[1]!='fluid_to_solid':continue
     v=list(map(float,p[2:]))
    else:v=list(map(float,p[1:]))
    data[t]=v
 return np.array([[t,*v] for t,v in sorted(data.items())])
def save(fig,name):
 fig.tight_layout();fig.savefig(figdir/(name+'.png'));fig.savefig(figdir/(name+'.svg'));plt.close(fig)
def history(ax):
 ax.set_xlabel('Steady-solver iteration');
 for x in (16000,20000):ax.axvline(x,color='0.5',ls=':',lw=.8)
def plotmon(ax,a,col,label,**kw):ax.plot(a[:,0],a[:,col],label=label,**kw)
mon={}
for key,region,name,kind,cols in [
 ('inlet_mass','fluid','inletMassFlow','normal',['signed_mass_flow_kg_s']),('outlet_mass','fluid','outletMassFlow','normal',['signed_mass_flow_kg_s']),
 ('inlet','fluid','inletAverage','normal',['area_average_p_Pa','area_average_T_K']),('outlet','fluid','outletAverage','normal',['mass_weighted_T_K','mass_weighted_h_J_kg']),
 ('fluid_temperature','fluid','minMaxFluid','extreme',['min_T_K','max_T_K']),('solid_temperature','solid','minMaxSolid','extreme',['min_T_K','max_T_K']),
 ('heated_wall_temperature','solid','topHeatedT','normal',['area_average_T_K']),('interface_heat','fluid','wallHeatFlux','patch',['min_flux_W_m2','max_flux_W_m2','signed_integral_W']),
 ('yplus','fluid','yPlus','patch',['minimum','maximum','average'])]:
 a=readmon(region,name,kind);mon[key]=a;csvwrite(key+'_history',['iteration',*cols],a)
# Restart histories: discard superseded tails and keep checkpoint itself.
records={};provenance=[]
pat=re.compile(r'Solving for (\w+), Initial residual = ([\deE+.-]+), Final residual = ([\deE+.-]+), No Iterations (\d+)')
for log,start in [('log.chtMultiRegionSimpleFoam',0),('log.chtMultiRegionSimpleFoam.restart',16000),('log.chtMultiRegionSimpleFoam.restart2',20000)]:
 records={k:v for k,v in records.items() if k[0]<=start};t=None;region=None;times=[]
 for line in (src/log).read_text().splitlines():
  if line.startswith('Time = '):t=float(line.split('=')[1]);times.append(t)
  if 'Solving for fluid region' in line:region='fluid'
  if 'Solving for solid region' in line:region='solid'
  m=pat.search(line)
  if m and t is not None:records[(t,region,m[1])]=[float(m[2]),float(m[3]),int(m[4]),log]
 provenance.append([log,start,min(times),max(times)])
rows=[[t,r,f,*v] for (t,r,f),v in sorted(records.items())]
csvwrite('residual_history',['iteration','region','field','initial_residual','final_residual','linear_iterations','source_log'],rows)
csvwrite('restart_provenance',['source_log','restart_checkpoint','first_logged_iteration','last_logged_iteration'],provenance)
fig,axs=plt.subplots(2,1,figsize=(10,8),sharex=True)
for ax,fields in zip(axs,[[('fluid',f) for f in ['Ux','Uy','Uz','p_rgh']],[('fluid',f) for f in ['h','k','omega']]+[('solid','h')]]):
 for r,f in fields:
  a=np.array([[t,v[0]] for (t,rr,ff),v in sorted(records.items()) if rr==r and ff==f]);ax.semilogy(a[:,0],a[:,1],lw=.7,label=r+' '+f)
 ax.set_ylabel('Initial equation residual');ax.legend(ncol=4,fontsize=9);history(ax)
axs[0].set_title('Coarse mesh • residual history (restart overlaps removed)');save(fig,'01_residual_convergence')
fig,axs=plt.subplots(2,1,figsize=(10,8),sharex=True)
for key,col,label in [('fluid_temperature',2,'Fluid maximum'),('solid_temperature',2,'Solid maximum'),('heated_wall_temperature',1,'Heated wall area mean'),('outlet',1,'Outlet mass-weighted mean')]:plotmon(axs[0],mon[key],col,label)
axs[0].set_ylabel('Temperature [K]');axs[0].legend(fontsize=9,ncol=2);axs[0].set_title('Coarse mesh • thermal history')
plotmon(axs[1],mon['interface_heat'],3,'Signed interface integral');axs[1].axhline(1350,color='k',ls='--',label='Nominal imposed load: 1350 W');axs[1].set_ylabel('Heat rate [W]');axs[1].legend(fontsize=9)
for ax in axs:history(ax)
save(fig,'02_thermal_history')
mi={r[0]:r[1] for r in mon['inlet_mass']};mo={r[0]:r[1] for r in mon['outlet_mass']}
bal=np.array([[t,mi[t],mo[t],mi[t]+mo[t],100*abs(mi[t]+mo[t])/abs(mi[t])] for t in sorted(mi.keys()&mo.keys()) if mi[t]!=0])
csvwrite('mass_balance_history',['iteration','signed_inlet_kg_s','signed_outlet_kg_s','net_outflow_kg_s','imbalance_percent_of_inflow'],bal)
fig,axs=plt.subplots(2,1,figsize=(10,7),sharex=True);axs[0].plot(bal[:,0],-bal[:,1],label='Inflow magnitude');axs[0].plot(bal[:,0],bal[:,2],ls='--',label='Outflow');axs[0].set_ylabel('Mass flow [kg/s]');axs[0].set_ylim(0,0.014);axs[0].ticklabel_format(axis='y',style='plain',useOffset=False);axs[0].legend();axs[0].set_title('Coarse mesh • boundary mass balance');axs[1].plot(bal[:,0],bal[:,4]);axs[1].set_ylabel('Absolute imbalance [%]')
for ax in axs:history(ax)
save(fig,'03_mass_balance')
axial=[]
for line in (src/'system/planes.txt').read_text().splitlines()[1:]:
 name,x=line.split();b=readmon('fluid','bulk_'+name);p=readmon('fluid','pAvg_'+name);m=readmon('fluid','mdot_'+name)
 assert b[-1,0]==p[-1,0]==m[-1,0]==30000
 axial.append([float(x),*b[-1,1:],p[-1,1],m[-1,1]])
a=np.array(axial);a=np.column_stack([a,(a[0,3]-a[:,3])/1000,a[:,2]-a[0,2]])
csvwrite('axial_profiles',['x_m','mass_weighted_T_K','mass_weighted_h_J_kg','area_average_p_Pa','mass_flow_kg_s','pressure_drop_from_plane0_kPa','enthalpy_rise_from_plane0_J_kg'],a)
fig,axs=plt.subplots(2,1,figsize=(9,7),sharex=True);axs[0].plot(a[:,0],a[:,1],'o-');axs[0].set_ylabel('Mass-weighted temperature [K]');axs[0].set_title('Heated-section axial profiles • iteration 30,000');axs[1].plot(a[:,0],a[:,5],'o-',color='#bc5525');axs[1].set_ylabel('Pressure drop from x = 0.15 m [kPa]');axs[1].set_xlabel('Axial coordinate x [m]');save(fig,'04_axial_profiles')
fig,ax=plt.subplots(figsize=(9,5));y=mon['yplus'];
for col,label in [(1,'Minimum'),(3,'Average (monitor reported)'),(2,'Maximum')]:plotmon(ax,y,col,label,marker='.',ms=3)
ax.set_ylabel('Wall y+ [dimensionless]');ax.set_title('Fluid-to-solid interface • wall y+');history(ax);ax.legend();save(fig,'05_wall_yplus')
q=mon['interface_heat'][-1,3];qadv=a[-1,4]*(a[-1,2]-a[0,2]);gap=1350-q
fig,ax=plt.subplots(figsize=(10,5));labels=['Nominal imposed load','Whole-interface integral','Enthalpy gain: x = 0.15–0.45 m'];vals=[1350,q,qadv];bars=ax.barh(labels,vals,color=['#777777','#267a9d','#b17832']);ax.bar_label(bars,fmt='%.2f W',padding=5);ax.set_xlim(0,1550);ax.set_xlabel('Heat / enthalpy transport rate [W]');ax.set_title('Thermal accounting • different measurement scopes');fig.text(.5,.015,'Not a closed energy balance: inlet enthalpy and all boundary heat flows are not supplied.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.04,1,1]);fig.savefig(figdir/'06_thermal_accounting.png');fig.savefig(figdir/'06_thermal_accounting.svg');plt.close(fig)
summary=[['Final saved iteration',30000,'iteration'],['Fluid cells',192000,'cells'],['Solid cells',120000,'cells'],['Maximum fluid temperature',mon['fluid_temperature'][-1,2],'K'],['Maximum solid temperature',mon['solid_temperature'][-1,2],'K'],['Mass-weighted outlet temperature',mon['outlet'][-1,1],'K'],['Mean heated exterior wall temperature',mon['heated_wall_temperature'][-1,1],'K'],['Interface integral (signed)',q,'W'],['Nominal documented imposed heating',1350,'W'],['Nominal heating minus interface integral',gap,'W'],['Relative difference (nominal load denominator)',100*gap/1350,'%'],['Heated-section pressure drop',a[-1,5],'kPa'],['Heated-section convective enthalpy gain',qadv,'W'],['Wall y+ minimum',y[-1,1],'1'],['Wall y+ average',y[-1,3],'1'],['Wall y+ maximum',y[-1,2],'1']]
csvwrite('final_summary',['quantity','value','unit'],summary)
energy=[['Signed inlet mass flow',bal[-1,1],'kg/s','Outward-normal sign convention'],['Signed outlet mass flow',bal[-1,2],'kg/s','Outward-normal sign convention'],['Net boundary mass flow',bal[-1,3],'kg/s','Zero to printed precision'],['Mass imbalance',bal[-1,4],'%','Relative to inlet magnitude; saved precision only'],['Nominal imposed heat',1350,'W','Documented assumption; current boundary-condition file absent'],['Whole-interface signed integral',q,'W','fluid_to_solid wallHeatFlux monitor'],['Nominal minus interface heat',gap,'W','Unresolved comparison; not a closed system energy residual'],['Relative nominal-interface difference',100*gap/1350,'%','Denominator 1350 W'],['Plane0-to-plane6 enthalpy transport gain',qadv,'W','mdot*(h6-h0); heated-section scope only'],['Full inlet-outlet enthalpy transport gain','not available','W','Inlet enthalpy not monitored; must not assume h_in=0'],['Full energy balance residual','not established','W','Requires all relevant heat/enthalpy fluxes on a consistent control volume']]
csvwrite('mass_energy_balance',['quantity','value','unit','qualification'],energy)
csvwrite('mesh_quality',['region','cells','max_aspect_ratio','max_nonorthogonality_deg','max_skewness','checkMesh'],[['fluid',192000,500.00063,0,4.7521803e-11,'Mesh OK'],['solid',120000,500.00063,0,5.5605468e-11,'Mesh OK']])
manifest={str(f.relative_to(src)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(src.rglob('*')) if f.is_file()};(out/'source_sha256.json').write_text(json.dumps(manifest,indent=2))
md='# Coarse-mesh results\n\nResults extracted from the supplied OpenFOAM logs and postProcessing files. The final restart log ends normally at iteration 30,000. These are steady-solver iterations, not physical seconds.\n\n## Final summary\n\n| Quantity | Value | Unit |\n|---|---:|---|\n'
for k,v,u in summary:md+=f'| {k} | {v:.8g} | {u} |\n'
md+='\n## Mass and thermal accounting\n\n| Quantity | Value | Unit | Qualification |\n|---|---:|---|---|\n'
for k,v,u,n in energy:md+=f'| {k} | {v if isinstance(v,str) else format(v,".8g")} | {u} | {n} |\n'
md+='\n## Axial profiles at iteration 30,000\n\n| x (m) | Bulk T (K) | Bulk h (J/kg) | Mean p (MPa) | Mass flow (kg/s) |\n|---:|---:|---:|---:|---:|\n'
for r in a:md+=f'| {r[0]:.2f} | {r[1]:.5f} | {r[2]:.3f} | {r[3]/1e6:.6f} | {r[4]:.8f} |\n'
md+='''
## Figures

PNG files are high-resolution previews; SVG files remain sharp at any zoom.

'''
for f in sorted(figdir.glob('*.png')):md+=f'![{f.stem}](figures/{f.name})\n\n'
md+='''## Interpretation and limitations

- This is a coarse-mesh demonstration. Medium/fine runs were not completed; mesh independence and validation are not established.
- Residual plots use **initial equation residuals**, not the final residuals of individual linear solves. A stable temperature/heat-rate plateau alone does not establish accuracy or energy closure.
- Restart checkpoints are 16,000 and 20,000 (dotted lines). Superseded history after each checkpoint is discarded; the checkpoint sample is retained when present. No interpolation is used. Monitor sampling intervals differ.
- Mass flow is conserved at the recorded precision at the final iteration; this does not imply an exactly zero mathematical error.
- The nominal 1350 W load is based on the previously documented heating, not reverified from a current boundary-condition file in this archive. Its roughly 10.24% difference from the interface integral remains unresolved.
- The plane0-to-plane6 enthalpy gain uses measured enthalpies and mass flow. Its control volume differs from the whole-interface monitor. It is not a full energy-balance residual. Enthalpy has a reference offset; inlet enthalpy must not be assumed zero.
- y+ is reported as supplied by the OpenFOAM monitor; the average is not reinterpreted as area weighted. Values alone do not validate the turbulence/wall treatment.
- Axial temperature/enthalpy are phi-weighted face-zone averages; pressure is area averaged. These seven stations describe bulk-fluid profiles, not wall-temperature profiles.
- checkMesh reports Mesh OK for both regions, but the maximum aspect ratio is about 500. Mesh checks do not replace discretization-error assessment.

## Files and reproducibility

Copy this directory into `output/coarse/analysis/` in your repository and link this README from the main README. All CSV files include units in column names or an explicit units column. Retain the original solver logs and monitor files separately; source SHA-256 hashes are included to identify the analyzed inputs.

Requires Python 3, NumPy and Matplotlib. Reproduce from the extracted review archive or original case:

```bash
python3 make_plots.py /path/to/methaneCHT /path/to/new-analysis-output
```

The script expects the three named solver logs, `system/planes.txt` and the supplied monitor directories. Nominal heat load (1350 W) and mesh-quality summary are specific to this reviewed case. Change those values only after checking the corresponding inputs/checkMesh log.
'''
(out/'README.md').write_text(md)
print(json.dumps({'summary':summary,'restart_segments':provenance,'residual_rows':len(rows),'figures':len(list(figdir.glob('*.png')))},indent=2))
