# rnaseq_fusion_simulation
For SMC-RNA challenge, code for analyzing datasets of validated fusions and generating synthetic fusion data.


Step 1: Create fusion events, truth file, and RSEM-format fusion reference:
create_fusions_truth.py    

Step 2: Align reads of real sample with RSEM to get alignment model:
run_RSEM_for_model.sh 

Step 3: Adjust estimated isoform values to include expression for fusion genes according to a model.
modify_model_tpm_for_diploid.R

Step 4: Generate reads from diploid and fusion references.
run_RSEM_for_reads.sh
