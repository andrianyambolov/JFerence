# JFerence Toolbox

The **JFerence** toolbox is MATLAB toolbox designed for conducting joint inference on impulse responses and forecasts, particularly in Bayesian econometric models. The primary function, `joint_inference`, enables users to compute credible sets, error bands, and central tendencies based on posterior draws, using a range of inference methods and configurations.

### Overview

The `joint_inference` function provides joint inference capabilities for Bayesian time series models by generating credible sets, error bands, and central tendency estimates. The function is highly customizable, allowing users to choose the type of inference, estimation method, loss function, calibration settings, credibility levels, and central tendency measure.

### Example Usage

In the **examples** subfolder, there is a ready-to-run example illustrating how to use the toolbox:

- **`run_me.m`**: This script estimates a standard 3-variable fiscal Vector Autoregression (VAR) model and saves the resulting Impulse Response Functions (IRFs) in `estimation.mat`.

- **`plot_figure.m`**: This file demonstrates how to use `joint_inference` to perform joint inference on the estimated IRFs. It also provides a plotting example to visualize the credible sets and error bands computed by the toolbox.

To use the example:

1. Run `run_me.m` to estimate the fiscal VAR model and generate IRFs.
2. Run `plot_figure.m` to compute joint inference on the IRFs and plot the results.

### Important Notes

- **Calibration**: Calibration methods are available for specific inference methods. For example, `'BDR'` (boundary draw rejection) and `'LQO'` (loss quantile optimization) are compatible with `'min-max'` inference.
- **Loss Functions**: The `LossFunction` field is only used with the `'min-max'` method. If `LossFunction` is not specified, it defaults to `'absolute'`.
- **Central Tendency**: The `CentralTendency` field allows selection of `'median'`, `'mean'`, or `'Bayes'` as the central measure, with `'Bayes'` specifically supported for `'min-max'`.

### Further Information

Further information about the methods, routines and the estimated example is available in the paper *"How to Conduct Joint Bayesian Inference in VAR Models?"*, Andrian Yambolov. 

### License

This toolbox is released under the MIT License.
