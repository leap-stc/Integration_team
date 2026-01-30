## Functional Unit Test set up

This example sets up a functional test on Derecho. The functional test is just a place to test your ml model using FTorch in Fortran.
Please file an issue if you run into problems. 

Developed by Adrianna Foster & Linnia Hawkins

### 1) Clone CTSM 
`git clone https://github.com/ESCOMP/CTSM.git CTSM`

I suggest cloning to your work directory $WORK or /glade/work/username/

### 2) Add in some mods
```
cd CTSM
git remote add jedwards https://github.com/jedwards4b/ctsm.git
git fetch jedwards
git checkout ftorch_d1fccec99
./bin/git-fleximod update
cd src/fates
git remote add linnia https://github.com/linniahawkins/fates
git fetch linnia
git checkout ml_example 
```

### 3) Set up your environment
```
export Torch_DIR=/glade/work/jedwards/conda-envs/ml5.6/ 
module load conda
conda activate ctsm_pylib # or some python environment with matplotlib and numpy
```

### 4) Build and test
```
rm -r _build
cd testing
./run_functional_tests.py  -t ml_example
```

Now you can add in your own model and code in ./ml_example/Example.F90


### Updating met forcing 
You may also want to provide your own input data (e.g., 5 years of meteorology) to test your ml model. 
To do this you will create a netcdf and point to it in:
```
CTSM/src/fates/testing/functional_tests.cfg
```
