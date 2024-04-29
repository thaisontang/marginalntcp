# marginalntcp
 
We provide source code for the simulation studies and simulated illustration of the proposed methods for the journal article "A marginal model for normal tissue complication probability" by authors Thai-Son Tang, Zhihui Liu, Ali Hosni, John Kim, Olli Saarela.

We recommend running the R scripts using `marginalntcp_code.Rproj` in RStudio to keep the directories consistent. All source code was written using R version 4.1.2 on a macOS computer.

The repository consists of:
* Source code for the simulation studies (`simulation` folder)
* Source code for an illustration of the proposed methods on a simulated dataset (`simulated_illustration` folder)

The `simulation` folder consists of:
- Source code for the simulation study described in the manuscript's main text (`main simulation` folder)
- Source code for the supplementary simulation studies (`supplementary simulation` folder)
- Essential functions required to run all simulation studies (`functions` folder)
- A `data` folder to store estimates from the simulation replicates.

**Instructions** (main simulation): Run the `simulation_macOS.R` file to run the simulation study and produce simulation replicate estimates. Next, run `true_risk.R` to produce the true values. The `simulation_macOS_visualization.R` script produces the perspective and contour plots comparing the estimates with the true values. Results will be stored in the `results` folder. The simulation lasts 19 hours for the bladder DVHs and 12 hours for the skin DVHs.

**Instructions** (supplementary simulations): There are 3 sets of R scripts corresponding to the no confounding, weak confounding and strong confounding simulation scenarios. A main script (`simulation_supplementary_`...) provides code to run the simulation study, a script to produce the true values (`true_risk_`...), and a script to visualize the estimates. Results will be stored in the `results` folder. Each simulation lasts approximately 19 hours for the bladder DVHs and 12 hours for the skin DVHs.

The `simulated_illustration` folder consists of:
- `Simulated Illustration.pdf`: A writeup illustrating the proposed methods on the simulated dataset (stored in `simulated_illustration/data/sim_illustration.RData`).
- The root file `illustration_script.R`.
- A `figures` folder storing the visualized estimates from bootstrap resampling.
- A `data` folder to store estimates from the simulation replicates.

**Instructions**: Run `illustration_script.R`, which contains source code for importing the simulated dataset (`simulated_illustration/data/sim_illustration.RData`), visualizing the dose-volume histograms for the bladder and skin, producing estimates of the proposed methods, producing bootstrap resampled estimates of the proposed methods and their visualizations. Bootstrap resampling lasts 17.1 hours for the bladder DVHs and 26 hours for the skin DVHs.
