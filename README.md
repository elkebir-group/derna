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
-s - mode <1,2,3>, 1 for mfe, 2 for mfe+cai, 3 for swipe
-l - lambda <[0,1]>
-a - sweep increment <(0,1]>
-r - <input rna file path>
-O - <swipe output csv file name>
-g - <[0,inf)>
-t - threshold tau <(0,1)>
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

`./derna -i <input> -o <output> -m 1 -s 2 -l 0.5`

### Sweep 

`./derna -i <input> -o <output> -m 1 -s 3 -O <output_csv> (-t 0.0001)`

### Evaluate an RNA sequence

`./derna -i <input> -o <output> -r <rna_sequence> -m -1`

### Only consider MFE

`./derna -i <input> -o <output> -m 1 -s 1`

### Nussinov based model (Fixed $\lambda$)

`./derna -i <input> -o <output> -m 0 -s 2 -l 0.5`

