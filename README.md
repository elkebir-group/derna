# DERNA

## Dependencies

* Recent C++ compiler (C++11)

## Compilation instructions

```
mkdir build
cd build
cmake ..
make
```

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


## Examples

### Fix $\lambda$ 

`./derna -i /data/uniprotSeq/P15421.fasta -o P15421_fixed_lambda.txt -m 1 -s 2 -l 0.5`

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

### Sweep (default thresholds) 

`./derna -i ./data/uniprotSeq/P15421.fasta -o P15421_sweep.txt -m 1 -s 3`

Estimated time: 10min 

### Evaluate an RNA sequence

`./derna -i ./data/uniprotSeq/P15421.fasta -o P15421_evaluation.txt -r ./data/RNA/P15421_rna.txt -m -1`

```
protein sequence: MYGKIIFVLLLSGIVSISASSTTGVAMHTSTSSSVTKSYISSQTNGITLINWWAMARVIFEVMLVVVGMIILISYCIR
eval MFE: -148.7
eval CAI: -28.3932
eval standard CAI: 0.694881

```

### Only consider MFE

`./derna -i ./data/uniprotSeq/P15421.fasta -o P15421_MFE_only.txt -m 1 -s 1`

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

### Nussinov based model (Fixed $\lambda$)

`./derna -i ./data/uniprotSeq/P15421.fasta -o P15421_nussinov.txt -m 0 -s 2 -l 0.5 -g 1`

