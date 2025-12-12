import numpy as np

sigs = [
    "reference_signal_arrays/reference_protocol_SIGNALS_with_k0.0um.npy",
    "reference_signal_arrays/reference_protocol_SIGNALS_with_k5.0um.npy",
    "reference_signal_arrays/reference_protocol_SIGNALS_with_k10.0um.npy",
    "reference_signal_arrays/reference_protocol_SIGNALS_with_k15.0um.npy",
    "reference_signal_arrays/reference_protocol_SIGNALS_with_k20.0um.npy",
    "reference_signal_arrays/reference_protocol_SIGNALS_with_k25.0um.npy",
    "reference_signal_arrays/reference_protocol_SIGNALS_with_k30.0um.npy",
    "reference_signal_arrays/reference_protocol_SIGNALS_with_k35.0um.npy",
    "reference_signal_arrays/reference_protocol_SIGNALS_with_k40.0um.npy",
    ]

params = [
    "reference_param_arrays/reference_protocol_PARAMETERS_with_k0.0um.npy",
    "reference_param_arrays/reference_protocol_PARAMETERS_with_k5.0um.npy",
    "reference_param_arrays/reference_protocol_PARAMETERS_with_k10.0um.npy",
    "reference_param_arrays/reference_protocol_PARAMETERS_with_k15.0um.npy",
    "reference_param_arrays/reference_protocol_PARAMETERS_with_k20.0um.npy",
    "reference_param_arrays/reference_protocol_PARAMETERS_with_k25.0um.npy",
    "reference_param_arrays/reference_protocol_PARAMETERS_with_k30.0um.npy",
    "reference_param_arrays/reference_protocol_PARAMETERS_with_k35.0um.npy",
    "reference_param_arrays/reference_protocol_PARAMETERS_with_k40.0um.npy",
    ]

all_signals = []
for sig in sigs:
    all_signals.append(np.load(f"{sig}"))

# Concatenate them
rs = np.concatenate(all_signals)
np.save("all_signals_all_substrates.npy", rs)

all_parameters = []
for param in params:
    all_parameters.append(np.load(f"{param}"))

# Concatenate them
rs = np.concatenate(all_parameters)
np.save("all_parameters_all_substrates.npy", rs)