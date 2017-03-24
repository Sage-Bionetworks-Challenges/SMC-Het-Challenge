# SMC-Het-Challenge

Development of the automated scoring system for the SMC-Het Challenge

## Buiding

Dependancies:
- swig 2.0.12 or greater
- gcc 4.9.2-12 or greater
- OpenMP 4.5 or greater
- python headers
- numpy with headers

Compiling c extensions and wrapping for python:

```bash
cd smc_het_eval
swig -python c_extensions.i
gcc -fopenmp -std=c99 -fpic -c c_extensions.c c_extensions_wrap.c -I /path/to/numpy/core/include -I /path/to/python/include/python2.7
ld -shared c_extensions.o c_extensions_wrap.o -o _c_extensions.so
```

## Evaluator

The evaluator for the SMC-Het Challenges is the python script **SMCScoring.py**.

To run it, you will need to pass in the following arguments..
* **-c** - choose from [1A, 1B, 1C, 2A, 2B, 3A, 3B]
* **--predfiles** - path(s) to the prediction file(s)
* **--truthfiles** - path(s) to the truth file(s)
* **--vcf** - path to the scoring vcf file
* **-o** - desired path to the output file
* **--approx** - (*OPTIONAL*) triggers subsampling mode; requires **2** arguments
  * *sampling fraction* - a float value 0.0 < x < 1.0 that denotes the sampling portion of the full matrix
  * *number of iterations* - integer specifying number of iterations to run
* **--approx_seed** - (*OPTIONAL*) allows you to specify the seed value for the random behaviour of **--approx**

### Examples

#### Running Challenge 1A

```bash
python SMCScoring.py -c 1A --predfiles ./some/path/to/1A_pred.txt --truthfiles ./some/path/to/1A_truth.txt -o ./some/path/to/1A_score.txt
```

#### Running Challenge 3A

```bash
python SMCScoring.py -c 3A --predfiles ./some/path/to/2A_pred.txt ./some/path/to/3A_pred.txt --truthfiles ./some/path/to/2A_truth.txt ./some/path/to/3A_truth.txt --vcf ./some/path/to/scoring.vcf -o ./some/path/to/3A_score.txt
```

#### Running Challenge 2B with Subsampling

```bash
python SMCScoring.py -c 2B --predfiles ./some/path/to/2B_pred.txt.gz --truthfiles ./some/path/to/2B_truth.txt.gz --vcf ./some/path/to/scoring.vcf -o ./some/path/to/2B_score.txt --approx 0.4 10
```
