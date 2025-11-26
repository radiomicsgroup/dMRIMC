
## Histology-informed microstructural diffusion simulations for MRI cancer characterisation — the Histo-μSim framework

<div align="center">
  <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/imgs/diagram22.png" alt="qrcode" width="auto" height="auto">
</div>

This code was developed by Athanasios Grigoriou (<agrigoriou@vhio.net>) and Francesco Grussu (<fgrussu@vhio.net>). **The project that gave rise to these results received the support of a fellowship from ”la Caixa” Foundation (ID 100010434). The fellowship code is "LCF/BQ/PR22/11920010". This study has been funded by Instituto de Salud Carlos III (ISCIII) through the project "PI21/01019" and co-funded by the European Union and by a AEI Severo Ochoa PhD fellowship (PRE2022-102586).**

<div align="center">
  <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/imgs/commbio25.png" alt="commbio" width="auto" height="auto">
</div>

If you find dMRIMC useful, please cite our article:

Athanasios Grigoriou, Carlos Macarro, Marco Palombo, Daniel Navarro-Garcia, Anna Voronova, Kinga Bernatowicz, Ignasi Barba, Alba Escriche, Emanuela Greco, María Abad, Sara Simonetti, Garazi Serna, Richard Mast, Xavier Merino, Núria Roson, Manuel Escobar, Maria Vieito, Paolo Nuciforo, Rodrigo Toledo, Elena Garralda, Roser Sala-Llonch, Els Fieremans, Dmitry S. Novikov, Raquel Perez-Lopez, Francesco Grussu. **"Histology-informed microstructural diffusion simulations for MRI cancer characterisation — the Histo-μSim framework"**. Communications Biology 2025, 8: 1695, doi: [10.1101/2024.07.15.24310280](https://doi.org/10.1038/s42003-025-09096-3).


## General description

This code accompanies our [preprint](https://doi.org/10.1101/2024.07.15.24310280), and provides all the tools needed to implement our proposed _Histo-μSim_ diffusion Magnetic Resonance Imaging (dMRI) technique for cancer imaging. It shows i) how to create realistic substrates from histology that can be used to perform Monte Carlo diffusion simulations; ii) how to synthesise dMRI signals for any protocol of interest given the diffusion simulations; and iii) how to use the synthetic signals to inform the estimation of tissue parameter of interest in any new input dMRI scan.

Here are detailed instructions on:

- generating signals from 2D histology ([tutorial here](https://github.com/radiomicsgroup/dMRIMC/blob/main/manuals/histology_to_signals.md));
- using the synthetic signals to inform microstructural parameter estimation ([tutorial here](https://github.com/radiomicsgroup/dMRIMC/blob/main/manuals/parameter_estimation.md)).

Third-party dependencies are reported in each tutorial.  

Moreover, [here](https://github.com/radiomicsgroup/dMRIMC/tree/main/dictionaries) we share the synthetic dMRI signals and corresponding tissue parameters for the acquisition protocols considered in our [preprint](https://doi.org/10.1101/2024.07.15.24310280).

## License
This repository is distributed under the **Attribution-NonCommercial-ShareAlike 4.0 International license** ([CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0)). Copyright (c) 2024, 2025, Fundació Privada Institut d’Investigació Oncològica de Vall d’Hebron (Vall d'Hebron Institute of Oncology (VHIO), Barcelona, Spain). All rights reserved. Link to license [here](https://github.com/radiomicsgroup/dMRIMC/blob/main/license.txt). 
