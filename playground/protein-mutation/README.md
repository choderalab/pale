# Protein mutation playground

These directory stores different example systems that will be used for demostrating how to run
free energy simulations using the new `openfe` and `feflow` infrastructure for protein mutations.

This space is intended to be a quick proof of concept repository of notebooks and scripts that
eventually will be transformed into a full user friendly CLI and library.

## Creating the python environment

Since these things are still under heavy development and in "alpha" versions. The environment is
a testing one. This environment can be created using the following command:

```bash
mamba env create -f devtools/conda-envs/test_env.yaml -n pale-dev
```

run from the root directory of your local copy of this repository.
