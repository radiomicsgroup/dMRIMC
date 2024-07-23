import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from sp_funcs import set_rc_params_journal, config_param_dictionary, labels_and_title_from_param, get_corr_metrics
import seaborn as sns
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


parser = argparse.ArgumentParser(
    prog='explore_niftii_LeaveOneOut.py',
    description='Generate contour plots of the fitting results',
    epilog='')
parser.add_argument("seq", type=str, choices=[
                    "PGSEin", "PGSEex", "TRSE"], help="protocol name - PGSEin/PGSEex/TRSE")
parser.add_argument("SNR", type=int, choices=[
                    20, 50], help="SNR value - 20/50")
parser.add_argument("config", type=int, choices=[
                    9, 13], help="Configuration argument - 9/13")

args = parser.parse_args()

seq = args.seq
SNR = args.SNR
config = args.config

print(f"Sequence: {seq}")
print(f"SNR: {SNR}")
print(f"Config: {config}")


set_rc_params_journal(latex=False)
#############################################
param_dict = config_param_dictionary(config)

############### Config 9 - Forward model 1###############
if config == 9:
    additional_names = ["$f_{in}$", "vCS$_{cyl}$", "$D_{0|in}$", "$D_{0|ex}$"]
    decimals = [2, 0, 2, 2]
    figsize = (40, 9)
    subplots_num = len(additional_names)
    fig, ax = plt.subplots(1, subplots_num, figsize=figsize,)
    cnt = 0
    fontsz = 40
    for par in range(1, len(param_dict[str(config)]) + 1):
        param = par
        pname = param_dict[str(config)][param - 1]
        units, title = labels_and_title_from_param(
            param_dict, param, config)
        one_more_title = additional_names[cnt]
        decs = decimals[par - 1]

        res_dir = "leave_one_out/res"
        all_results = []
        all_gts = []
        for sub in substrates:
            all_results.append(nib.load(
                f"leave_one_out/res/config_{config}/{seq}_{sub}_out/SNR_{SNR}_par{param}.nii").get_fdata())
            all_gts.append(nib.load(
                f"leave_one_out/niftiis/params/{seq}_params_SNR_{SNR}_config_{config}_{sub}_out_test.nii").get_fdata()[:, :, :, param - 1])

        # Concatenate them
        rsss = np.concatenate(all_results)
        gtss = np.concatenate(all_gts)
        sp, p = get_corr_metrics(rsss, gtss)

        arr2 = np.array((rsss, gtss[:, :, 0]))
        xlabel = "Ground truth"
        ylabel = "Estimation"
        title0 = f"Distribution of {pname}"
        curr_ax = ax.flatten()[par - 1]
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
            thresh=0
        )

        for acol in a.collections:
            acol.set_edgecolor("face")
        fig.add_axes(a)
        # to stop title overlap
        curr_ax.set_aspect("equal")
        fig.subplots_adjust(wspace=0.3)
        curr_ax.set_title(title0)

        # TICKS and AXES
        xmin, xmax = np.round(arr2[1, :, :].min(), decimals=decs), np.round(
            arr2[1, :, :].max(), decimals=decs)
        ymin, ymax = np.round(arr2[0, :, :].min(), decimals=decs), np.round(
            arr2[0, :, :].max(), decimals=decs)
        xticks = np.round(np.linspace(xmin, xmax, 6), decimals=decs)
        yticks = np.round(np.linspace(ymin, ymax, 6), decimals=decs)
        curr_ax.set_yticks(yticks[1:-1])
        curr_ax.set_xticks(xticks[1:-1])

        curr_ax.axline((xticks[1], yticks[1]), slope=1,
                       color="black", linestyle="--")
        curr_ax.set_xlabel(xlabel)
        curr_ax.set_ylabel(ylabel)
        ax[0].set_ylabel("$\\bf{MC-informed \ fitting}$" + "\n" +
                         f"\n{ylabel}" + " ")  # for some reason here \bf needs another \
        if (units and title) != None:
            curr_ax.set_xlabel(xlabel + " " + units)
            curr_ax.set_ylabel(ylabel + " " + units)
            curr_ax.set_title(f"{one_more_title}\n{title}")
        if (sp and p) != None:
            textstr = '\n'.join((
                f"Correlation = {p.correlation: .2f}",
            ))
            props = dict(boxstyle='round', facecolor='white', alpha=1)

            # place a text box in upper left in axes coords
            curr_ax.text(0.05, 0.95, textstr, transform=curr_ax.transAxes, fontsize=fontsz,
                         verticalalignment='top', bbox=props)
        cnt += 1
        plt.savefig(f"figs/{seq}_SNR_{SNR}_config_{config}_all_params.png",
                    dpi=600, transparent=True, bbox_inches='tight', format='png')


############### Config 13 - Forward model 2 ###############
if config == 13:
    additional_names = ["$f_{in}$", "mCS", "varCS",
                        "skewCS", "$D_{0|in}$", "$D_{0|ex}$"]
    figsize = (30, 25)
    subplots_num = len(additional_names)
    fig, ax = plt.subplots(2, 3, figsize=figsize,)
    cnt = 0
    ax_nums = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]
    fontsz = 40
    for par in range(1, len(param_dict[str(config)]) + 1):
        param = par
        pname = param_dict[str(config)][param - 1]
        units, title = labels_and_title_from_param(
            param_dict, param, config)
        one_more_title = additional_names[cnt]

        res_dir = "leave_one_out/res"
        all_results = []
        all_gts = []
        for sub in substrates:
            all_results.append(nib.load(
                f"leave_one_out/res/config_{config}/{seq}_{sub}_out/SNR_{SNR}_par{param}.nii").get_fdata())
            all_gts.append(nib.load(
                f"leave_one_out/niftiis/params/{seq}_params_SNR_{SNR}_config_{config}_{sub}_out_test.nii").get_fdata()[:, :, :, param - 1])

        # Concatenate them
        rsss = np.concatenate(all_results)
        gtss = np.concatenate(all_gts)
        sp, p = get_corr_metrics(rsss, gtss)

        arr2 = np.array((rsss, gtss[:, :, 0]))
        xlabel = "Ground truth"
        ylabel = "Estimation"
        title0 = f"Distribution of {pname}"
        curr_ax = ax.flatten()[par - 1]
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
            thresh=0
        )
        for acol in a.collections:
            acol.set_edgecolor("face")
        fig.add_axes(a)
        # to stop title overlap
        curr_ax.set_aspect("equal")
        fig.subplots_adjust(wspace=0.3)
        curr_ax.set_title(title0)
        curr_ax.axline((arr2[1, :, :].min(), arr2[0, :, :].min()),
                       slope=1, color="black", linestyle="--")
        curr_ax.set_xlabel(xlabel)
        curr_ax.set_ylabel(ylabel)
        if (units and title) != None:
            curr_ax.set_xlabel(xlabel + " " + units)
            curr_ax.set_ylabel(ylabel + " " + units)
            curr_ax.set_title(f"{one_more_title}\n{title}")
        if (sp and p) != None:
            textstr = '\n'.join((
                f"Correlation = {p.correlation: .2f}",
            ))
            props = dict(boxstyle='round', facecolor='white', alpha=1)

            # place a text box in upper left in axes coords
            curr_ax.text(0.05, 0.95, textstr, transform=curr_ax.transAxes, fontsize=fontsz,
                         verticalalignment='top', bbox=props)
        cnt += 1
        plt.savefig(f"figs/{seq}_SNR_{SNR}_config_{config}_all_params.png",
                        dpi=600, transparent=True, bbox_inches='tight', format='png')
