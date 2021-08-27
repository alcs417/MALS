# MALS

Cancer Subtype Identification based on Multi-view Subspace Clustering with Adaptive Local Structure Learning

**Method Description**

MALS is a novel similarity-based multi-omics data integration framework to identify cancer subtypes. It integrates multi-view subspace clustering and adaptive local structure learning into a unified framework. Specifically, the model learns a new representation for each omics in the latent subspace and obtains a consensus similarity matrix by adaptively assigning different view weights during the learning process. Therefore, the latent representations and the consensus similarity matrix could be enhanced with each other in an alternative manner. We use Augmented Lagrange Multiplier with Alternating Direction Minimization method to efficiently optimize the proposed objective.

**Software Requirements**

\>= MATLAB 2014b. 

**Usage**

To run our algorithm, please load the script 'processTCGAdata.m' into your MATLAB programming environment and click 'run'. Users can also run the script in standard command-line mode, where you should input the following command to run the algorithm:

matlab -nodisplay -nodesktop -nosplash -r "processTCGAdata;exit;"

All the cancer datasets used in the code can be directly downloaded at http://acgt.cs.tau.ac.il/multi_omic_benchmark/download.html. To run the model, you should first construct a 'struct' data structure in Matlab, where each view represents a type of omics data. We have provided a demo input data "aml.mat" under the './data' directory. 

**Estimation of the Number of Clusters**

In our work, we used rotation cost to estimate the number of clusters given a similarity matrix. We also provided the estimation result by another method, i.e. eigen gap, although we only reported the clustering results obtained by the rotation cost in our manuscript. Users can choose whichever estimation method they prefer. For more details, please refer to 'Similarity network fusion for aggregating data types on a genomic scale'.

**Parameters**

There are three parameters in our algorithm that users can tune according to their own needs, i.e. lambda_R, lambda_A and the number of neighbors k. The default value for k is set to 9. To reset the parameters, simply open the 'processTCGAdata.m' and pass new values for these parameters.

**Input and Output Directories**

To change the input directory, please refer to the 'dataDir' variable in the processTCGAdata.m. For output directory, please refer to the 'outDir' variable in the same script. 

**Interpretation of the Clustering Results**

The output of our model is the cluster labels assignment to all the samples. A sample output could be found at './output/aml.mat'. Users can then perform the survival analysis or other analysis according to their own interests. We also provided the scripts to measure the clustering performance (e.g. ACC, NMI) when the ground truth labels of the samples are given. Please see './funs/ClusteringMeasure.m'.

**Contact**

For any questions regarding our work, please feel free to contact us: ALCS417@sdnu.edu.cn.
