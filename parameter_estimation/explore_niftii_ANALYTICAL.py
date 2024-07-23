import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from sp_funcs import set_rc_params_journal, config_param_dictionary, labels_and_title_from_param, get_corr_metrics
import seaborn as sns
import matplotlib.gridspec as gridspec
import argparse


substrates = ['substrate_1',
              'substrate_2',
              'substrate_3',
              'substrate_4',
              'substrate_5',
              'substrate_6',
              'substrate_7',
              'substrate_8',
              'substrate_9',
              'substrate_10',
              'substrate_11',
              'substrate_12',
              'substrate_13',
              'substrate_14',
              'substrate_15',
              'substrate_16',
              'substrate_17',
              'substrate_18',
              ]


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('true', 'True'):
        return True
    elif v.lower() in ('false', 'False'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


parser = argparse.ArgumentParser(
    prog='explore_niftii_LeaveOneOut_ANALYTICAL',
    description='Generate contour plots of the fitting results of the analytical expression',
    epilog='')
parser.add_argument("SNR", type=int, choices=[
                    20, 50], help="SNR value - 20/50")

args = parser.parse_args()

SNR = args.SNR
print(f"SNR: {SNR}")


set_rc_params_journal(latex=False)
#############################################
seq = "PGSEin"
config = 9
phantom_type = "old_way_cylinders"

fig = plt.figure(figsize=(22, 7))

large = True
if large:
    fig = plt.figure(figsize=(33, 9))
gs = gridspec.GridSpec(ncols=3, nrows=1, figure=fig)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])

ax = [ax1, ax2, ax3]
param_dict = config_param_dictionary(11)
# in the analytical there is not D0ex so we go up to the third parameter of
# config 9
cnt = 0
additional_names = ["$f_{in}$", "vCS$_{cyl}$", "$D_{0|in}$", "$ADC_{ex}$"]
decimals = [2, 0, 2, 2]
for par in range(1, len(param_dict[str(config)])):
    param = par
    pname = param_dict[str(config)][param - 1]
    units, title = labels_and_title_from_param(param_dict, param, config)
    one_more_title = additional_names[cnt]
    decs = decimals[par - 1]

    dir = f"leave_one_out/res_maxlik/{phantom_type}/config_9"
    if param == 1:
        # ICVF
        res_array = nib.load(f"{dir}/SNR_{SNR}_fin.nii").get_fdata()
    elif param == 2:
        # cell size
        res_array = nib.load(f"{dir}/SNR_{SNR}_Lum.nii").get_fdata()
    elif param == 3:
        # D0in
        res_array = nib.load(f"{dir}/SNR_{SNR}_D0um2ms-1.nii").get_fdata()
    gt = nib.load(
        f"leave_one_out/analytical/niftiis/PGSEin_all_params_config_9_SNR_{SNR}.nii")

    gt_arr = gt.get_fdata()[:, :, :, param - 1]
    rsss = res_array
    gtss = gt_arr
    sp, p = get_corr_metrics(rsss, gtss)

    arr2 = np.array((rsss, gtss[:, :, 0]))
    xlabel = "Ground truth"
    ylabel = "Estimation"
    title0 = f"Distribution of {pname}"
    curr_ax = ax[par - 1]
    a = sns.kdeplot(
        x=arr2[1, :, :].flatten(),
        y=arr2[0, :, :].flatten(),
        cmap="viridis",
        fill=True,
        alpha=1,
        levels=100,
        cbar=False,
        ax=curr_ax,
        cut=0,
        bw_method="scott",
        cbar_kws=dict(drawedges=False, boundaries=None),
        label="PASTA",
        thresh=0)

    for acol in a.collections:
        acol.set_edgecolor("face")
    fig.add_axes(a)

    curr_ax.set_aspect('equal', adjustable='box')
    # TICKS and AXES
    xmin, xmax = np.round(arr2[1, :, :].min(), decimals=decs), np.round(
        arr2[1, :, :].max(), decimals=decs)
    ymin, ymax = np.round(arr2[0, :, :].min(), decimals=decs), np.round(
        arr2[0, :, :].max(), decimals=decs)
    xticks = np.round(np.linspace(xmin, xmax, 6), decimals=decs)
    yticks = np.round(np.linspace(ymin, ymax, 6), decimals=decs)
    if par == 1:
        curr_ax.axline((0.027, 0.027), slope=1, color="black", linestyle="--")
    elif par == 2:
        print(xmin, xmax)
        print(ymin, ymax)
        curr_ax.axline((8.43, 8.2), slope=1.02, color="black", linestyle="--")
    else:
        curr_ax.axline((xmin, xmin), slope=1, color="black", linestyle="--")
    curr_ax.set_ylim(xmin, xmax)
    curr_ax.set_yticks(xticks[1:-1])  # setting the ticks the same as x axis
    curr_ax.set_xticks(xticks[1:-1])
    curr_ax.set_xlabel(xlabel)
    curr_ax.set_ylabel(ylabel)
    ax[0].set_ylabel("$\\bf{Classical\ analytical\ fitting}$" + "\n" +
                     f"\n{ylabel}" + " " + units)  # for some reason here \bf needs another \
    if (units and title) != None:
        curr_ax.set_xlabel(xlabel + " " + units)
        curr_ax.set_ylabel(ylabel + " " + units)
        curr_ax.set_title(f"{one_more_title}\n{title}")
    if (sp and p) != None:
        textstr = '\n'.join((
            f"Correlation = {p.correlation: .2f}",
        ))
        props = dict(boxstyle='round', facecolor='white', alpha=1)
        curr_ax.text(0.3, 0.15, textstr, transform=curr_ax.transAxes, fontsize=40,
                     verticalalignment='top', bbox=props)
    cnt += 1
    plt.savefig(f"figs/{seq}_{phantom_type}_SNR_{SNR}_config_{config}_all_params_ANALYTICAL.png",
                dpi=600, transparent=True, bbox_inches='tight', format='png')
