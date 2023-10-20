from casm.project import Project, Selection, write_eci
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# construct a 'Selection'
proj = Project()
sel = Selection(proj, 'casm_learn_input', all=False)

# query results into a pandas.DataFrame
comp = 'comp(a)'
Ef = 'formation_energy'
clEf = 'clex(formation_energy)'
hull_dist = 'hull_dist(casm_learn_input,comp)'
cl_hull_dist = 'clex_hull_dist(casm_learn_input,comp)'
config = 'name'
sel.query([comp, Ef, clEf, hull_dist, cl_hull_dist])

# get convex hull configurations, sorted for nice plotting
df = sel.data.sort_values([comp])
print(df)
hull_tol = 1e-6
df_hull = df[df[hull_dist] < hull_tol]
#print(df.to_string())
df_unstable = df[df[hull_dist] >= hull_tol]
print(df_hull[[config,comp,hull_dist,cl_hull_dist]])

cl_df_hull = df[df[cl_hull_dist] < hull_tol]
cl_df_unstable = df[df[cl_hull_dist] >= hull_tol]
print(cl_df_hull[[config,comp,hull_dist,cl_hull_dist]])
df["error"] = df[clEf]-df[Ef]

plt.tick_params(direction='in', which="major", length=7, width=2)
plt.subplots_adjust(hspace=0.1)

plt.rc('font', family='sans-serif')
plt.rc('text', usetex=True)
plt.rcParams['legend.fancybox'] = False
plt.rcParams['legend.loc'] = 'best'
plt.rcParams['legend.numpoints'] = 1
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['legend.framealpha'] = None
plt.rcParams['legend.scatterpoints'] = 2
plt.rcParams['legend.edgecolor'] = 'inherit'
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['lines.markersize'] = 5

fig,(ax1,ax2) = plt.subplots(nrows=1,ncols=2,figsize=(7,3.5))

for i in [ax1,ax2]:
    plt.setp(i.spines.values(), linewidth=1)
    i.tick_params(direction='in', which="major", length=6, width=1)
    i.tick_params(direction='in', which="minor", length=4, width=1)

ax1.set_ylabel(r'Formation Energy (meV/atom)', size=11)
ax1.set_xlabel(r'${\rm x}$ in ${\rm Na_xPb_{1-x}}$', size=11)
ax1.set_xlim([-0.1, 1.1])
ax1.set_ylim([-200, 150])
ax1.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax1.yaxis.set_major_locator(ticker.MultipleLocator(100))
ax1.yaxis.set_minor_locator(ticker.MultipleLocator(50))

ax1.plot(df_hull[comp], df_hull[Ef]*1000, marker = "o", ls = "", lw = 0.0, color = 'Green', label = '')
ax1.plot(df_hull[comp], df_hull[Ef]*1000, marker = "o", ls = "-", lw = 1., color = 'Green', label ="DFT hull")
ax1.scatter(df_unstable['comp(a)'], df_unstable['formation_energy']*1000, s=8, marker = 'o', color='Green',alpha=0.6,label="DFT")

ax1.plot(cl_df_hull[comp], cl_df_hull[clEf]*1000, marker = "o", ls = "", lw = 0.0, color = 'crimson', label = '')
ax1.plot(cl_df_hull[comp], cl_df_hull[clEf]*1000, marker = "o", ls = "-", lw = 1., color = 'crimson', label = "CE hull")
ax1.scatter(cl_df_unstable[comp], cl_df_unstable[clEf]*1000, s=8, marker = '.', color='crimson',alpha=0.6,label="CE")
ax1.hlines(0, -0.1, 1.1, colors = 'black', linestyles = '--', label = '')
ax1.legend()

# plot CE error
ax2.scatter(x=df['error']*1000, y=df[hull_dist]*1000,marker="o",color='lightseagreen',s=10,alpha=0.6)
ax2.set(xlim=(-50,50),ylim=(-5,200),xlabel='Error of CE (meV/atom)',ylabel=r'Energy above Convex Hull (meV/atom)')  
ax2.set_xlabel('Error of CE (meV/atom)',fontsize=12)
ax2.set_ylabel(r'Energy above Convex Hull (meV/atom)')

# draw 2 boxes for showing error boundaries
rect1 = patches.Rectangle((-20,0),40,100,linewidth=0.5,edgecolor="black",facecolor='none',linestyle='--')
rect2 = patches.Rectangle((-10,0),20,50,linewidth=0.5,edgecolor="black",facecolor='none',linestyle='--')   
ax2.add_patch(rect1)
ax2.add_patch(rect2)  
ax2.text(0,110,r'${\pm}$20 meV',color="black",ha='center')
ax2.text(0,60,r'${\pm}$10 meV',color="black",ha='center')
ax2.yaxis.set_minor_locator(ticker.MultipleLocator(50))
ax2.xaxis.set_minor_locator(ticker.MultipleLocator(50))

plt.tight_layout()
plot_title = "final_plot.pdf"
plt.savefig(plot_title, format="pdf", bbox_inches='tight')