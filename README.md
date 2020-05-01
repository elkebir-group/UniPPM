# UniPPM

UniPPM is an almost uniform sampler for the Perfect Phylogeny Mixture (PPM) problem, it takes a frequency matrix whoses rows are mixture of clones as input, and sample trees that satisfied the sum condition. UniPPM also support input in the form of reads and solve a maximum likelihood version of the problem (MLPPM). 

![Overview of UniPPM](figure1.png)


## Compliation and Usage instructions

UniPPM is written in Python 3 and requires UniGen:

1. http://www.cs.rice.edu/CS/Verification/Projects/UniGen/

Copy the python files to the path where you build UniGen2 and run following command.


    $python main.py [path/to/input_file.tsv] [output/directory/] [output_file_name] [#number_of_samples] [\lambda (if it is an MLPPM instance)] 

For instance to sample 100 trees from the solutions space of "n7_S81_k1_clustered.tsv" (a PPM instance) you need to run:

    $python main.py input/simulate/n7_S81_k1_clustered.tsv output/ test 100
