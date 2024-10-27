# JFerence Toolbox

The **JFerence** toolbox is MATLAB toolbox designed for conducting joint inference on impulse responses and forecasts, particularly in Bayesian econometric models. The primary function, `joint_inference`, enables users to compute credible sets, error bands, and central tendencies based on posterior draws, using a range of inference methods and configurations.

### Supported Joint Inference Estimators

The JFerence toolbox supports the following joint inference estimators:

- **sup-t**: Proposed by Montiel Olea and Plagborg-Møller (2018) in *Simultaneous confidence bands: Theory, implementation, and an application to SVARs*.

- **min-max**: Introduced by Inoue and Kilian (2022) in *Joint Bayesian Inference About Impulse Responses in VAR Models*. The toolbox further includes calibration routines to adjust error bands, ensuring they achieve the nominal probability content.

- **Bonferroni and Sidak corrections**: These corrections for pointwise intervals are traditionally used as simpler, naive methods for estimating joint error bands.

### Example Usage

In the **examples** subfolder, there is a ready-to-run example illustrating how to use the toolbox:

- **`run_me.m`**: This script estimates a standard 3-variable fiscal Vector Autoregression (VAR) model and saves the resulting Impulse Response Functions (IRFs) in `estimation.mat`.

- **`plot_figure.m`**: This file demonstrates how to use `joint_inference` to perform joint inference on the estimated IRFs. It also provides a plotting example to visualize the credible sets and error bands computed by the toolbox.

To use the example:

1. Run `run_me.m` to estimate the fiscal VAR model and generate IRFs.
2. Run `plot_figure.m` to compute joint inference on the IRFs and plot the results.

### Further Information

Further information about the methods, routines and the estimated example is available in the paper *"How to Conduct Joint Bayesian Inference in VAR Models?"*, Andrian Yambolov. 

### License

This toolbox is released under the MIT License.
