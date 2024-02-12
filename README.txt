If this code is used in a publication, please cite the manuscript:
"CARLA: Adjusted common average referencing for cortico-cortical evoked potential data"
by H Huang, G Ojeda Valencia, NM Gregg, GM Osman, MN Montoya, GA Worrell, KJ Miller, and D Hermes.
A preprint is available currently at https://doi.org/10.48550/arXiv.2310.00185.

Correspondence:
H Huang: huang.harvey@mayo.edu; D Hermes: hermes.dora@mayo.edu

*****

DEPENDENCIES

- mnl_ieegBasics: https://github.com/MultimodalNeuroimagingLab/mnl_ieegBasics
- matmef: https://github.com/MultimodalNeuroimagingLab/matmef
- mnl_seegview: https://github.com/MultimodalNeuroimagingLab/mnl_seegview
- SPM12 (as dependency for mnl_seegview): https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
- vistasoft: https://github.com/vistalab/vistasoft
- Freesurfer v7: https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall
- GIfTI: https://github.com/gllmflndn/gifti

*****

USAGE

a. The data used in this analysis will be available upon manuscript publication on OpenNeuro, in BIDS format. Please download all the data and copy into a folder named "data", in the root level code directory. In the meantime before data release, all results pertaining to simulated CCEPs can still be reproduced by the code alone.

b. Pial, cortical, and subcortical segmentations for each subject were obtained using Freesurfer v7. These were used to localize each stimulation site to a particular tissue type. The relevant Freesurfer outputs for each subject are located in the data subdirectory: data/derivatives/freesurfer

c. All custom analyses were performed in MATLAB R2023a. Step-by-step code blocks and instructions to generate all manuscript figures and results are in "main.m". Please open "main.m" and follow the code sections contained therein.