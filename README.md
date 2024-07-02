# MRI_T1rho_T2star_reliability
A collection of Matlab M-files for reading Philips DICOM MRI T1rho and T2star knee image data, registering the spin lock or echo times using Melastix and elastix, reading OsiriX segmentation CSV files and calculating T1rho or T2star for regions of interest based on tibiofemoral contact.

Please see ImageAnalysisPipeline6.pdf for a flowchart for the sequence of execution of Matlab M-files.

Please see the header comments in the M-files.

Note that the data must be in particular directories.

An additional region of interest on the posterior femoral condyles (not in contact with the tibia) was added to the analysis.  Please see ImageAnalysisPipeline7.pdf for a flowchart for the sequence of execution of Matlab M-files.

M-files for the analysis of the entire cartilage regions (lateral and medial femur, and lateral and medial tibia) were added.  See mri_fitra.m.

Additional M-files for postprocessing the results for publication figures were added.

M-files for analysis of the knee meniscus were added.  See seg_mr_rois.m and mri_mer_fit.m.  The analysis includes eroding the meniscus masks by a pixel to check for edge effects.
