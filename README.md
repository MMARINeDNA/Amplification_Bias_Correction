# Amplification_Bias_Correction
Following QAQC of sequence data using the pipeline in the "metabarcoding_QAQC_pipeline" repository, use this step to correct for any primer- or species-specific biases in amplification efficiency.

1. **Step 1**: Use dada2 output and merged sample data sheet (created with the Qmd report) to generate input files for ABC using ABC_input_data.R. _**You will need to edit the top chunk of this code to specify project name, primer, wd (working drive), dada2 output file, and sample data file.**_
2. **Step 2**: Output files from Step 1, as well as two CSVs describing the mock mixture input data (read count and proportion), are formatted and modified in this step and then run in the ABC model from Shelton et al. 2022. The results of the model are used to correct estimates of species included in mock mixtures. 
