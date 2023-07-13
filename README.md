# DERNA

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/derna/README.html#package-derna)

DERNA is a tool that enables the design of RNA sequences based on protein sequences. 
DERNA accepts a protein sequence as input and provides a collection of Pareto optimal solutions consisting of RNA sequences that optimize both minimum free energy and codon adaptation index (CAI). Additionally, DERNA can function as a tool for predicting RNA structures and calculating CAI for given RNA sequences.

## Contents

1. [Installation](#install)
      * [Using conda](#conda) (recommended)
      * [Build from source](#compilation) (alternative)
          * [Dependencies](#dep)
          * [Compilation](#comp)
  2. [Usage instructions](#usage)
      * [Usage Example](#example)
    
        
<a name="install"></a>

## Installation

<a name="conda"></a>

### Using conda

1. Create a new conda environment named "derna" and install dependencies:

   ```bash
   conda create -n derna
   ```

2. Then activate the created environment: `conda activate derna`.
3. Install the package into current environment "derna":

    ```bash
    conda install -c bioconda derna
    ```

<a name="compilation"></a>

### Build from source

<a name="dep"></a>

#### Dependencies

* Recent C++ compiler (C++11)

<a name="comp"></a>
#### Compilation

```
mkdir build
cd build
cmake ..
make
```

<a name="usage"></a>
## Usage instructions

```
-i - <input file path>
-o - <output file path>
-m - model <0,1,-1> , 0 for nussinov, 1 for zuker, -1 for eval
-s - mode <1,2,3>, 1 for mfe, 2 for mfe+cai, 3 for sweep
-l - lambda <[0,1]>
-a - sweep increment <(0,1]>
-r - <input rna file path>
-O - <sweep output csv file name>
-g - <[0,inf)>
-t - threshold tau <(0,1)>
-p - threshold tau2 <(0,1)>
-c - <codon usage table file path>
-d - <energy parameters (model) directory>
```

```
./derna -i <input> -o <output> -m <model> -s <mode> ...
```

```
input: input file path 
output: output file path 
model: integer 0 for Nussinov based model, 1 for Zuker based model, -1 for eval model 
mode: integer 1 for only MFE mode, integer 2 for MFE + CAI mode, integer 3 for lambda swipe mode 
lambda: lambda value for MFE + CAI mode or lambda swipe mode 
incr: increment interval for lambda swipe mode 
swipe: swipe output csv file name 
g: minimal gap allowed in Nussinov based model 
rna: input rna file path for eval model
```

<a name="example"></a>
### Examples

#### Fix $\lambda$ 

`./derna -i ../data/uniprotSeq/P15421.fasta -o P15421_fixed_lambda.txt -m 1 -s 2 -l 0.5`

`cat P15421_fixed_lambda.txt`

```
protein sequence: MYGKIIFVLLLSGIVSISASSTTGVAMHTSTSSSVTKSYISSQTNGITLINWWAMARVIFEVMLVVVGMIILISYCIR

lambda: 0.5
Zuker CAI
Energy: -74.2202
Time taken by DP is : 32sec
lambda: 0.5,O: -7422.02,mfe: -14870,cai: -25.9662,combined: -7422.02
zuker cai bp: (((((((((.....((((((((((.((((((((.(((((((((...))))))))).)))))))).)))))))))))))))))))((((((....((((((((.(((.((((((((((.(((((((((.(((((.((((((((((((((.((((((((((((....)))))).)))))))))))))))))))).))))))))))))))))))))))))))))))))))))))))),size: 234
zuker rna: AUGUAUGGCXXXXXCAUCUUCGUCUUGCUGCUCUCCGGGAUCGUXUCGAUCUCGGCGAGCAGCACGACGGGGGUGGCCAUGCAUACGAGUXXXXGCAGUAGCXUGAXUAAGAGUUAUXUAUCCUCACXGACCAACGGCAUCACCUUGAXAAAUUGGUGGGCGXXGGCCCGCXUAAUUUUCGAGGUGAUGCUGGUGGUCGUGGGGAUGAUAAUUCUUAUCAGCUACUGCAUUCGU.size: 234
zuker cai rna: AUGUAUGGCAAGAUCAUCUUCGUCUUGCUGCUCUCCGGGAUCGUGUCGAUCUCGGCGAGCAGCACGACGGGGGUGGCCAUGCAUACGAGUACCAGCAGUAGCGUGACUAAGAGUUAUAUAUCCUCACAGACCAACGGCAUCACCUUGAUAAAUUGGUGGGCGAUGGCCCGCGUAAUUUUCGAGGUGAUGCUGGUGGUCGUGGGGAUGAUAAUUCUUAUCAGCUACUGCAUUCGU.size: 234
Codon Adaptation Index: 0.716842
Minimum Free Energy: -148.7

```

#### Sweep (default thresholds) 

`./derna -i ../data/uniprotSeq/P15421.fasta -o P15421_sweep.txt -O P15421_sweep -m 1 -s 3`

Estimated time: 10min 

### Evaluate an RNA sequence

`./derna -i ../data/uniprotSeq/P15421.fasta -o P15421_evaluation.txt -r ./data/RNA/P15421_rna.txt -m -1`

`cat P15421_evaluation.txt`

```
protein sequence: MYGKIIFVLLLSGIVSISASSTTGVAMHTSTSSSVTKSYISSQTNGITLINWWAMARVIFEVMLVVVGMIILISYCIR
eval MFE: -148.7
eval CAI: -28.3932
eval standard CAI: 0.694881

```

#### Only consider MFE

`./derna -i ../data/uniprotSeq/P15421.fasta -o P15421_MFE_only.txt -m 1 -s 1`

`cat P15421_MFE_only.txt`

```
protein sequence: MYGKIIFVLLLSGIVSISASSTTGVAMHTSTSSSVTKSYISSQTNGITLINWWAMARVIFEVMLVVVGMIILISYCIR

Zuker
Energy: -148.7
Time taken by DP is : 3sec
Time taken : 3sec
zuker bp:(((((((((.....((((((((((.((((((((.(((((((((...))))))))).)))))))).)))))))))))))))))))((((((....((((((((.(((.((((((((((.(((((((((.(((((.((((((((((((((.((((((((((((....)))))).)))))))))))))))))))).))))))))))))))))))))))))))))))))))))))))), size: 234
zuker rna:AUGUAUGGCXXXXXCAUCUUCGUCCUGCUGCUCUCCGGGAUCGUXUCGAUCUCGGCGAGCAGCACGACGGGGGUGGCCAUGCAUACGAGUXXXXGCAGUAGCXUGAXUAAGAGUUAUXUAUCCUCACXGACCAACGGCAUCACCUUGAXAAAUUGGUGGGCGXXGGCCCGCXUAAUUUUCGAGGUGAUGCUGGUGGUCGUGGGGAUGAUAAUUCUUAUCAGCUACUGCAUUCGU, size: 234
zuker rna:AUGUAUGGCAAAAUCAUCUUCGUCCUGCUGCUCUCCGGGAUCGUUUCGAUCUCGGCGAGCAGCACGACGGGGGUGGCCAUGCAUACGAGUACUAGCAGUAGCGUGACUAAGAGUUAUAUAUCCUCACAGACCAACGGCAUCACCUUGAUAAAUUGGUGGGCGAUGGCCCGCGUAAUUUUCGAGGUGAUGCUGGUGGUCGUGGGGAUGAUAAUUCUUAUCAGCUACUGCAUUCGU, size: 234
zuker cai: 0.694881

```

#### Nussinov based model (Fixed $\lambda$)

`./derna -i ../data/uniprotSeq/P15421.fasta -o P15421_nussinov.txt -m 0 -s 2 -l 0.5 -g 1`

#### Specify Codon Usage Table

`./derna -i ../data/uniprotSeq/P15421.fasta -o P15421_fixed_lambda.txt -m 1 -s 2 -l 0.5 -c ./data/InputFiles/sample_codon_usage.csv`

#### Specify Energy Parameters

`./derna -i ../data/uniprotSeq/P15421.fasta -o P15421_fixed_lambda.txt -m 1 -s 2 -l 0.5 -d ./data/InputFiles/`
