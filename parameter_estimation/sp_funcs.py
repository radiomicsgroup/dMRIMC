import matplotlib.pyplot as plt
from scipy import stats

"""
Helper functions for the signal processing pipeline
"""


def set_rc_params_journal(latex=False):
    """
    Set the matplotlib figure parameters

    Parameters
    ----------
    latex : bool, optional
        Use latex for rendering

    """
    plt.rcParams['font.size'] = 35
    plt.rcParams["font.family"] = "STIXGeneral"
    plt.rcParams["mathtext.fontset"] = 'stix'
    plt.rcParams['figure.figsize'] = [14., 12.0]
    if latex:
        plt.rc('text', usetex=True)
        plt.rcParams['font.size'] = 35
    return


def get_corr_metrics(res_array, gt_array):
    """
    Returns the spearman and pearson correlations

    Parameters
    ----------
    res_array : nd.array
        Array of results from fitting
    gt_array : nd.array
        Array of ground truth of the parameters

    """
    spman = stats.spearmanr(res_array.flatten(), gt_array[:, :, 0].flatten())
    pson = stats.pearsonr(res_array.flatten(), gt_array[:, :, 0].flatten())
    return spman, pson


def config_param_dictionary(config_range):
    """
    Create a dictionary with the config number as the keys and the contained
    parameters as values

    Returns
    -------
    my_dict : dictionary
        Dictionary containing the included metrics


    Parameters
    ----------
    config_range : int
        Up to what configuration to go, best to go 13


    """
    my_dict = {}
    for param_config in range(1, config_range + 1):
        param_config = str(param_config)
        if param_config == "9":
            # Parameters are ICVF, VCS, D0_intra, D0_extra
            included_metrics = ['ICVF', 'vCS', 'D0_intra', 'D0_extra']
            my_dict.update({param_config: included_metrics})
        elif param_config == "13":
            # Parameters are ICVF, mean diameter, var of diameter, skew of diameter, D0_intra, D0_extra
            included_metrics = ['ICVF', "Diam_Mean",
                                "Variance", "Skew", "D0_intra", "D0_extra"]
            my_dict.update({param_config: included_metrics})
    return my_dict


def labels_and_title_from_param(param_dict, param, config):
    """
    Get figure elements based on current configuration of parameters and config

    Parameters
    ----------
    param_dict : dictionary
        Dictionary of the included parameters (see above)
    param : int
        Current parameter in processing
    config : int
        Configuration number

    Returns
    -------
    units : str
        Units of the plot
    title : str
        Title of the plot (i.e name of the parameter)
    """
    pname = param_dict[str(config)][param - 1]
    if pname == 'ICVF':
        units = ""
        title = "Intracellular\n volume fraction"
    elif pname == 'vCS':
        units = "$(\mu m)$"
        title = "Volume-weighted\n cell size"
    elif pname == 'KL_div':
        units = ""
        title = "KL divergence"
    elif pname == 'D0_intra':
        units = "$(\mu m^2 / ms)$"
        title = "Intrinsic intracellular\n diffusivity"
    elif pname == 'D0_extra':
        units = "$(\mu m^2 / ms)$"
        title = "Intrinsic extracellular\n diffusivity"
    elif pname == 'D_normed':
        units = ""
        title = "Normalized distance from center of mass"
    elif pname == 'MC_normed':
        units = ""
        title = "Mean covariance"
    elif pname == 'Diam_IQR':
        units = "$\mu m$"
        title = "IQR of cell size"
    elif pname == 'Gamma_Shape':
        units = ""
        title = "Gamma distribution - Shape a"
    elif pname == 'Gamma_Scale':
        units = "$\mu m$"
        title = "Gamma distribution - Scale b"
    elif pname == 'Lumen_fraction':
        units = ""
        title = "Lumen fraction"
    elif pname == 'Diam_Mean':
        units = "$(\mu m)$"
        title = "Mean cell size"
    elif pname == "(<l6>/<l2>)^1/4":
        units = "$\mu m$"
        title = "Volume-weighted cell size - Cylinders"
    elif pname == "Variance":
        units = "$(\mu m^2)$"
        title = "Variance of cell size"
    elif pname == "Kurtosis":
        units = ""
        title = "Kurtosis of cell size"
    elif pname == "Skew":
        units = ""
        title = "Skew of cell size"

    return units, title
