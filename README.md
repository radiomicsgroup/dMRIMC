
## A Monte Carlo simulation framework for histology-informed diffusion MRI cancer characterisation and microstructural parameter estimation

<div align="center">
  <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/imgs/diagram22.png" alt="qrcode" width="auto" height="auto">
</div>

This code was developed by Athanasios Grigoriou (<agrigoriou@vhio.net>) and Francesco Grussu (<fgrussu@vhio.net>). **The project that gave rise to these results received the support of a fellowship from ”la Caixa” Foundation (ID 100010434). The fellowship code is "LCF/BQ/PR22/11920010"**.

<div align="center">
  <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/imgs/qr_img_MC_2024_paper.png" alt="qrcode" width="auto" height="auto">
</div>

If you find dMRIMC useful, please cite our preprint:

Athanasios Grigoriou, Carlos Macarro, Marco Palombo, Anna Voronova, Kinga Bernatowicz, Ignasi Barba, Alba Escriche, Emanuela Greco, María Abad, Sara Simonetti, Garazi Serna, Richard Mast, Xavier Merino, Núria Roson, Manuel Escobar, Maria Vieito, Paolo Nuciforo, Rodrigo Toledo, Elena Garralda, Roser Sala-Llonch, Els Fieremans, Dmitry S. Novikov, Raquel Perez-Lopez, and Francesco Grussu. **"A Monte Carlo simulation framework for histology-informed diffusion MRI cancer characterisation and microstructural parameter estimation"**. medRxiv 2024: 2024.07.15.24310280, doi: 10.1101/2024.07.15.24310280. Link to preprint [here](https://www.medrxiv.org/content/10.1101/2024.07.15.24310280v1).

**UNDER CONSTRUCTION -  Note that we are still polishing up the repository so parts are subject to change or could be missing**

## General description

This code accompanies our [preprint](https://www.medrxiv.org/content/10.1101/2024.07.15.24310280v1), sharing the process of creating realistic cancer substrates, running Monte Carlo simulations in them and synthesizing the resulting MRI signal according to any user-specified PGSE protocol. Those signals can the be used to inform microstructural parameter estimation.

Here are detailed instructions on:

- [Generating signals from histology](https://github.com/radiomicsgroup/dMRIMC/blob/main/manuals/histology_to_signals.md)
- [Using synthetic signals to inform microstructural parameter estimation](https://github.com/radiomicsgroup/dMRIMC/blob/main/manuals/parameter_estimation.md)

## License
This repository is distributed under the **Attribution-NonCommercial-ShareAlike 4.0 International license** ([CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0)). Copyright (c) 2024, Fundació Privada Institut d’Investigació Oncològica de Vall d’Hebron (Vall d'Hebron Institute of Oncology (VHIO), Barcelona, Spain). All rights reserved. Link to license [here](https://github.com/radiomicsgroup/dMRIMC/blob/main/license.txt). 


## Dependencies
The code was developed with the following package versions:
- **python** (3.9.17)
- **numpy** (1.25.2)
- **scipy** (1.10.1)
- **pyntcloud** (0.3.1)
- **opencv** (4.7.0)

Blender version `3.6.12` was used for the geometric operations

