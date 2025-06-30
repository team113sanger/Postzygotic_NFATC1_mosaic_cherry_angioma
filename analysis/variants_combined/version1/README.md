# Get somatic variant calling results from Figshare

Files with the MAF file contatinig the summary of the results of the variant calling can be found in Figshare Project [here](https://figshare.com/projects/A_post-zygotic_disruptive_germline_NFATC1_variant_in_a_patient_with_segmental_cherry_angiomas/254267). 

The files can be downloaded for plotting in the `analysis/variants_combined/version1` directory. Ensure the `PROJECTDIR` variable is set to the directory where this repository was cloned into. 

The following commands can be used to download the files:

```bash
PROJECTDIR=/lustre/8117_2744_ivo_cherry_angioma_wes
mkdir -p ${PROJECTDIR}/analysis/variants_combined/version1

curl -k -o bundle.zip  https://figshare.com/ndownloader/articles/29437103/versions/1
unzip bundle.zip
```