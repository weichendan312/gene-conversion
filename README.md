# gene-conversion
## Citations

- If you use the these methods to identify gene conversion, please cite:

  _Wei et al. (2021) Conversion between 100-million-year-old duplicated genes contributes to rice subspecies divergence. [BMC genomics]()_


## Contents

WCV-I：the Ks values between paralogous and orthologous gene pairs were 
used to infer possible whole gene conversion.

WCV-II: The ratios were used to infer unexpected changes in gene tree 
topology in quartets, depending on whether the paralogues were more 
similar to each other than the orthologues. This is a strict criterion 
used for the detection of whole gene conversion

PCV: A combination of dynamic planning and phylogenetic analysis was used
to document the differences between two aligned bases from paralogous and 
orthologous genes for each genome. In averaged distance arrays, the paralogs 
in each species should be more distant if no PCV was involved.


## Installation

Install ClustalW
```console
wget http://www.clustal.org/download/current/clustalw-2.1.tar.gz
cd clustalw-2.1
./configure
make
```
Install Python2

https://www.python.org/

## Dependencies

Following are a list of third-party python packages that are used by
some routines in the library. These dependencies are _not_ mandatory
since they are only used by a few modules.

- [Bio](http://www.Bio.org)

There are other Python modules here and there in various scripts. The
best way is to install them via `pip install` when you see
`ImportError`.
