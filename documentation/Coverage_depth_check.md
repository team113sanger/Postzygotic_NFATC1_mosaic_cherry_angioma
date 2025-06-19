# Coverage depth check 	

To analyse the evenness of coverage depth of the samples, we used a file that contained the coordinates of the baits used in exome captured liftovered to GRChr38 followed by the use of the `samtools depth` command to calculate the coverage depth of each sample. The output was then processed to generate a summary of the coverage depth across all samples.

Samples with a coverage depth of >20X coverage across at least 80% of the baits used were considered for varaint calling. 

## Dependencies

### Required data

### Required software