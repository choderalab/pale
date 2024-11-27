This directory is intended for scripts and files and examples that are not fully developed but
can be useful to have for future reference. Basically an space for quick and dirty scripts and
snippets.

We expect this directory to be removed once we have achieved some stability in the tools and
structure that we want to achieve in this project.


## Boltz prediction script
`boltz_prediction.py` is a script convenient to generate `yaml` files that could be readily
used with `boltz` to make structure predictions. For example you can run it using:

```python
python boltz_tyk2_prediction.py --protein-file protein.pdb --ligand-file ligands.sdf --protein-id MYPROTEIN --output-yaml MYPROTEIN_ligands.yaml
```
