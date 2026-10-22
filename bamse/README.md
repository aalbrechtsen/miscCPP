# BAMSE

BAMSE (Bayesian Association with Multiple SNPs and Environment) is a program for association studies of multiple SNPs and environmental factors. It is a Bayesian method that can model interactions between a large number of SNPs and environmental risk factors while accounting for multiple testing.

More specifically, BAMSE is a Markov chain Monte Carlo (MCMC) method that identifies sets of SNPs and environmental factors that, when combined, increase disease risk or change the distribution of a quantitative trait:

- Combinations of genetic and environmental factors define **risk sets**. Individuals with genotypes that are members of a risk set have a modified distribution of disease risk or trait value.
- Traits are modelled with normal distributions (quantitative traits) or as binary traits (case-control).
- The Markov chain has a state space of all possible risk sets and the parameters (for example mean and variance) of the trait distribution for each risk set. Its stationary distribution is the posterior density of risk sets and trait distributions, and inference is based on that posterior.
- The chain is simulated with the Metropolis-Hastings algorithm, using reversible jumps to move between risk sets.

This differs from earlier approaches by allowing non-linear models and by handling the multiple testing problem in a computationally and statistically efficient way. The method is developed for unrelated individuals.

The full description of the input format and every option is in [BAMSE_MANUAL.pdf](BAMSE_MANUAL.pdf).

## Install

```bash
git clone https://github.com/aalbrechtsen/miscCPP.git
cd miscCPP/bamse
make
```

This builds the `bamse` program in the same folder. It only needs a C++ compiler.

## Run

BAMSE reads its settings from an options file (`options.txt` by default) and then analyses each input file given on the command line:

```bash
./bamse -o options.txt data1.txt data2.txt
```

Each line of an input file is one individual: the trait value followed by the genotypes for each SNP. The options file sets the number of individuals, SNPs and environmental factors, the priors, the number of MCMC iterations, burn-in and thinning, and whether the trait is quantitative or binary. Each line starts with the value, followed by a description of what it is. See [test_files/options.txt](test_files/options.txt) for an example, and the manual for details.

The results, including the posterior distribution of the number of risk sets, are written to standard output.

## Example

The folder [test_files](test_files) has five simulated data sets with 500 individuals and 20 SNPs, and a matching options file:

```bash
cd test_files
../bamse -o options.txt set1.test > set1.mcmcres
```

The options file runs 100,000 iterations. Lower the number of iterations in `options.txt` for a quick test.

## Cite

Albrechtsen A, Castella S, Andersen G, Hansen T, Pedersen O, Nielsen R (2007). A Bayesian multilocus association method: allowing for higher-order interaction in association studies. *Genetics* 176(2): 1197–1208. https://doi.org/10.1534/genetics.107.071696

## History

Moved here in 2026 from the popgen.dk software wiki. Precompiled Windows and Linux programs from 2007 are still at http://popgen.dk/software/download/bamse/.
