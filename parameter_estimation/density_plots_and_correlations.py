import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from scipy import stats
import matplotlib.gridspec as gridspec
import seaborn as sns
from matplotlib.ticker import FixedLocator, FixedFormatter

def get_corr_metrics(res_array, gt_array):
    spman = stats.spearmanr(res_array.flatten(), gt_array[:,:,0].flatten())
    pson = stats.pearsonr(res_array.flatten(), gt_array[:,:,0].flatten())
    return spman, pson


def config_param_dictionary(config_range):
    """
    Create a dictionary with the config number as the keys and the contained
    parameters as values
    """
    my_dict = {}
    for param_config in range(1,config_range+1):
        param_config = str(param_config)
        if param_config == "14":
            # Parameters are ICVF, VCS, D0_intra, D0_extra, kappa
            included_metrics = ['ICVF', 'vCS', 'D0_intra', "D0_extra", "kappa"]
            my_dict.update({param_config: included_metrics})
    return my_dict


def labels_and_title_from_param(param_dict, param, config):
    pname = param_dict[str(config)][param-1]
    if pname == 'ICVF':
        units = ""
        title = "Intracellular\n volume fraction"
    elif pname == 'vCS':
        units = "$(\mu m)$"
        title = "Volume-weighted\n cell size"
    elif pname == 'D0_intra':
        units = "$(\mu m^2 / ms)$"
        title = "Intrinsic intracellular\n diffusivity"
    elif pname == 'D0_extra':
        units = "$(\mu m^2 / ms)$"
        title = "Intrinsic extracellular\n diffusivity"
    elif pname == "kappa":
        units = "$(\mu m/s)$"
        title = "Permeability"
    
    return units, title

substrates = [
    "sub1",
    "sub2",
    "sub3",
    "sub4",
    "sub5",
    "sub6",
    "sub7",
    "sub8",
    "sub9",
    "sub10",
    "sub11",
    "sub12",
    "sub13",
    "sub14",
    "sub15",
    "sub16",
    "sub17",
    "sub18",
]

seq = "PGSEin"
SNR = 50
config = 14
param_dict = config_param_dictionary(config)
save = True
plt.rcParams['font.size'] = 14
cmap = "viridis"

#%%
additional_names = ["$f_{in}$", "vCS$_{cyl}$" , "$D_{0|in}$", "$D_{0|ex}$", "$\kappa$"]
decimals = [2,0,2,2,0]
fig = plt.figure(figsize=(11,8))

gs = gridspec.GridSpec(ncols=6, nrows=2, figure=fig, hspace=0.55, wspace=0.6)

# first row
ax1 = fig.add_subplot(gs[0,:2])
ax2 = fig.add_subplot(gs[0,2:4])
ax3 = fig.add_subplot(gs[0,4:])
# second row
ax4 = fig.add_subplot(gs[1,1:3])
ax5 = fig.add_subplot(gs[1,3:5])

ax = [ax1, ax2, ax3, ax4, ax5]

xlabel = "Ground truth"
ylabel = "Estimation"
ax[0].set_ylabel(ylabel, fontsize=14)
ax[3].set_ylabel(ylabel, fontsize=14)
for par in range(1,len(param_dict[str(config)])+1):
    param = par
    pname = param_dict["14"][param-1]
    # print(param)
    units, title = labels_and_title_from_param(param_dict, param, config)
    one_more_title = additional_names[par-1]
    # to control the ticks better
    decs = decimals[par-1]
    
    res_dir = f"dict_fit/config_{config}"
    gt_dir = f"LeaveOneOut/niftis"
    all_results = []
    all_gts = []
    for sub in substrates:
        all_results.append(nib.load(f"{res_dir}/{seq}_{sub}_out/SNR_{SNR}_par{param}.nii").get_fdata())
        all_gts.append(nib.load(f"{gt_dir}/parameters/{seq}_parameters_SNR{SNR}_{sub}_out.nii").get_fdata()[:,:,:,param-1])

    # Concatenate them
    rsss = np.concatenate(all_results)
    gtss = np.concatenate(all_gts)
    # Correlations
    sp, p = get_corr_metrics(rsss, gtss)

    arr2 = np.array((rsss,gtss[:,:,0]))
    title0 = f"Distribution of {pname}"
    curr_ax = ax[param-1]
    
    a = sns.kdeplot(
                x = arr2[1,:,:].flatten(),
                y = arr2[0,:,:].flatten(), 
                cmap=cmap, 
                fill=True,
                alpha=1,
                levels=100,
                cbar=False,
                ax=curr_ax,
                cut=0,
                bw_method="scott",
                cbar_kws=dict(drawedges=False, boundaries=None),
                thresh=0
                )
    
    ##### Ticks and axes limits
    minx, miny = arr2[1,:,:].flatten().min(), arr2[0,:,:].flatten().min()
    maxx, maxy = arr2[1,:,:].flatten().max(), arr2[0,:,:].flatten().max()
    mins = np.array([minx, miny])
    maxs = np.array([maxx, maxy])
    
    # get the largest min and smallest max to create equal axes
    axrange = np.max(mins), np.min(maxs)
    xlim = (axrange[0], axrange[1])
    ylim = (axrange[0], axrange[1])
    curr_ax.set_xlim(xlim)
    curr_ax.set_ylim(ylim)

    xticks = np.linspace(axrange[0], axrange[1], 5)
    yticks = np.linspace(axrange[0], axrange[1], 5)
    
    curr_ax.xaxis.set_major_locator(FixedLocator(xticks))
    curr_ax.yaxis.set_major_locator(FixedLocator(yticks))
    
    curr_ax.xaxis.set_major_formatter(FixedFormatter([f'{x:.{decs}f}' for x in xticks]))
    curr_ax.yaxis.set_major_formatter(FixedFormatter([f'{y:.{decs}f}' for y in yticks]))
    
    # Identity line
    curr_ax.plot([axrange[0], axrange[1]], [axrange[0], axrange[1]], 'w--', linewidth=1)

    # Title and correlation annotation
    fig.supylabel('Histo-μSim fitting', fontsize=20, fontweight="regular", x=0.04)
    ax[0].set_ylabel(ylabel, fontsize=14)
    if (units and title) != None:
        curr_ax.set_xlabel(xlabel + " " +units, fontsize=13)
        curr_ax.set_title(f"{one_more_title}\n{title}", fontsize=13)
    
    if (sp and p) != None:
        textstr = '\n'.join((
            f"Correlation = {p.correlation: .2f}",
            ))
        props = dict(boxstyle='round', facecolor='white', alpha=1)

        # place a text box in upper left in axes coords
        curr_ax.text(0.05, 0.95, textstr, transform=curr_ax.transAxes,
                 verticalalignment='top', bbox=props)
    
    # Create a new figure and recreate the plot separately
    new_fig, new_ax = plt.subplots(figsize=(4, 4))

    sns.kdeplot(
        x=arr2[1, :, :].flatten(),
        y=arr2[0, :, :].flatten(),
        cmap=cmap,
        fill=True,
        alpha=1,
        levels=30,
        cbar=False,
        ax=new_ax,
        cut=0,
        bw_method="scott",
        thresh=0,
    )

    # Remove white lines in the contours for the new plot
    for acol in new_ax.collections:
        acol.set_edgecolor("face")

    # Set the same limits and labels for the new plot
    new_ax.set_xlim(xlim)
    new_ax.set_ylim(ylim)
    new_ax.set_xlabel("Ground truth", fontsize=18)
    new_ax.set_ylabel("Estimation", fontsize=18)
    if (units and title) != None:
        new_ax.set_xlabel(xlabel + " " +units, fontsize=18)
        new_ax.set_ylabel(ylabel + " " +units, fontsize=18)
        new_ax.set_title(f"{one_more_title}\n{title}", fontsize=18)
    xticks = np.linspace(axrange[0], axrange[1], 4)
    yticks = np.linspace(axrange[0], axrange[1], 5)
    new_ax.xaxis.set_major_locator(FixedLocator(xticks))
    new_ax.yaxis.set_major_locator(FixedLocator(yticks))
    
    new_ax.xaxis.set_major_formatter(FixedFormatter([f'{x:.{decs}f}' for x in xticks]))
    new_ax.yaxis.set_major_formatter(FixedFormatter([f'{y:.{decs}f}' for y in yticks]))
    # Add identity line to the new plot
    new_ax.plot([axrange[0], axrange[1]], [axrange[0], axrange[1]], 'w--', linewidth=1)

    # Add correlation annotation to the new plot
    if (sp and p) is not None:
        new_ax.text(0.05, 0.95, textstr, transform=new_ax.transAxes, verticalalignment='top', bbox=props)

    # Save the individual figure
    new_fig.savefig(f"figs/individual_subplots/{pname}_{seq}.pdf", bbox_inches='tight', format='pdf')
    plt.close(new_fig)
    
    if save:
        plt.savefig(f"figs/{seq}_SNR{SNR}_config{config}_all_params.pdf",dpi=300, bbox_inches='tight', format='pdf')
plt.show()            
