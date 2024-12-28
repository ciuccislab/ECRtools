# ECRtools: A MATLAB Toolbox for Electrical Conductivity Relaxation (ECR) Analysis and Optimization

**ECRtools** is a powerful and user-friendly MATLAB toolbox designed to streamline the analysis and optimization of Electrical Conductivity Relaxation (ECR) experiments. Developed by Francesco Ciucci at HKUST, this toolbox provides a comprehensive suite of tools for researchers working with ECR measurements, particularly in the field of solid-state ionics.

## Key Features:

*   **GUI Interface:** An intuitive graphical user interface allows for easy input of parameters, visualization of ECR responses, sensitivity analysis, and plotting of confidence regions.
*   **Data Fitting:** Robust fitting algorithms to extract surface exchange coefficients (*k*) and diffusion coefficients (*D*) from ECR data.
*   **Sensitivity Analysis:** Determine the sensitivity of ECR measurements to *k* and *D* parameters, aiding in understanding experimental limitations and parameter identifiability.
*   **Asymptotic Confidence Regions:** Calculate and visualize asymptotic confidence regions for estimated parameters, providing a measure of uncertainty in the fitted values.
*   **Optimal Experimental Design (OED):** Optimize experimental parameters (e.g., sample thickness, measurement time) to maximize information content and improve the accuracy of parameter estimation. Includes both classical and robust OED approaches.
*   **Synthetic Data Generation:** Generate synthetic ECR data with user-defined noise levels for testing and validation of fitting procedures.
*   **Comprehensive Demos:** Seven detailed demos covering all major functionalities, providing a practical guide to using ECRtools.

## Target Audience:

Researchers and engineers working with ECR measurements in various fields, including:

*   Solid-State Ionics
*   Materials Science
*   Electrochemistry
*   Fuel Cells
*   Batteries

## Dependencies:

*   MATLAB (versions > 2014a)
*   [OPTI Toolbox](https://github.com/jonathancurrie/OPTI)

## Installation:

1.  Install the OPTI Toolbox by running `opti_Install.m` from its installation directory.
2.  Download ECRtools and add its folder and subfolders to your MATLAB path.

## How to Cite:

If you use ECRtools in your research, please cite the following publications:

*   Ciucci, F. (2013). Electrical conductivity relaxation measurements: Statistical investigations using sensitivity analysis, optimal experimental design and ECRTOOLS. *Solid State Ionics*, *239*, 28-40. [https://doi.org/10.1016/j.ssi.2013.03.020](https://doi.org/10.1016/j.ssi.2013.03.020)
*   Wan, T. H., Saccoccio, M., Chen, C., & Ciucci, F. (2015). Assessing the identifiability of *k* and *D* in electrical conductivity relaxation via analytical results and nonlinearity estimates. *Solid State Ionics*, *270*, 18-32. [https://doi.org/10.1016/j.ssi.2014.11.026](https://doi.org/10.1016/j.ssi.2014.11.026)

## Support:

For inquiries, please contact [mefrank@ust.hk](mailto:mefrank@ust.hk).

## License:

This project is distributed without any warranty; use it at your own risk.

## Contribute:

Contributions and improvements to ECRtools are welcome! Please feel free to fork the repository and submit pull requests.
