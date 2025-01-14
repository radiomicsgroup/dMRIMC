
## Histology-informed microstructural diffusion simulations for MRI cancer characterisation — the Histo-μSim framework

<div align="center">
  <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/imgs/diagram22.png" alt="qrcode" width="auto" height="auto">
</div>

This code was developed by Athanasios Grigoriou (<agrigoriou@vhio.net>) and Francesco Grussu (<fgrussu@vhio.net>). **The project that gave rise to these results received the support of a fellowship from ”la Caixa” Foundation (ID 100010434). The fellowship code is "LCF/BQ/PR22/11920010"**.

<div align="center">
  <img src="https://github.com/radiomicsgroup/dMRIMC/blob/main/imgs/qr_img_MC_2024_paper.png" alt="qrcode" width="auto" height="auto">
</div>

If you find dMRIMC useful, please cite our preprint:

Athanasios Grigoriou, Carlos Macarro, Marco Palombo, Daniel Navarro-Garcia, Anna Voronova, Kinga Bernatowicz, Ignasi Barba, Alba Escriche, Emanuela Greco, María Abad, Sara Simonetti, Garazi Serna, Richard Mast, Xavier Merino, Núria Roson, Manuel Escobar, Maria Vieito, Paolo Nuciforo, Rodrigo Toledo, Elena Garralda, Roser Sala-Llonch, Els Fieremans, Dmitry S. Novikov, Raquel Perez-Lopez, Francesco Grussu. **"Histology-informed microstructural diffusion simulations for MRI cancer characterisation — the Histo-μSim framework"**. medRxiv 2024: 2024.07.15.24310280, doi: 10.1101/2024.07.15.24310280. Link to preprint [here](https://doi.org/10.1101/2024.07.15.24310280).

**UNDER CONSTRUCTION -  Note that we are still polishing up the repository so parts are subject to change or could be missing**

## General description

This code accompanies our [preprint](https://doi.org/10.1101/2024.07.15.24310280), sharing the process of creating realistic cancer substrates, running Monte Carlo simulations in them and synthesizing the resulting MRI signal according to any user-specified PGSE protocol. Those signals can the be used to inform microstructural parameter estimation.

Here are detailed instructions on:

- generating signals from 2D histology ([tutorial here](https://github.com/radiomicsgroup/dMRIMC/blob/main/manuals/histology_to_signals.md));
- using the synthetic signals to inform microstructural parameter estimation ([tutorial here](https://github.com/radiomicsgroup/dMRIMC/blob/main/manuals/parameter_estimation.md)).

Third-party dependencies are reported in each tutorial.  


## License
This repository is distributed under the **Attribution-NonCommercial-ShareAlike 4.0 International license** ([CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0)). Copyright (c) 2024, Fundació Privada Institut d’Investigació Oncològica de Vall d’Hebron (Vall d'Hebron Institute of Oncology (VHIO), Barcelona, Spain). All rights reserved. Link to license [here](https://github.com/radiomicsgroup/dMRIMC/blob/main/license.txt). 
