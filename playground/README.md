This directory is intended for scripts and files and examples that are not fully developed but
can be useful to have for future reference. Basically an space for quick and dirty scripts and
snippets.

We expect this directory to be removed once we have achieved some stability in the tools and
structure that we want to achieve in this project.


## Protein mutation
The `protein-mutation` subdirectory stores basic notebooks and scripts for different examples on
systems involving protein mutations. The idea is running free energy predictions using `feflow`
and `openfe` on these systems.

## Boltz prediction script
`boltz_yaml.py` is a script convenient to generate `yaml` files that could be readily
used with `boltz` to make structure predictions. For example you can run it using:

```python
python boltz_yaml.py --protein-file protein.pdb --ligand-file ligands.sdf --protein-id MYPROTEIN --output-yaml MYPROTEIN_ligands.yaml
```

This will generate a `yaml` file in the specified path (`ligand.yaml` by default), that is ready
to be consumed by `boltz predict` command.

## Boltz-1/2 Examples 
Example YAML files for multiple types of use cases is provided in the `boltz-example-yamls` directory

These can all be easily run using the command: 
```
boltz predict [example yaml] --cache [CACHE DIR] --use_msa_server
```
FYI, on our chodera lab installation, our cache dir is `/data1/choderaj/shared/.cache/boltz/`

Example systems include: 
- Monomeric sequences
- Protein + ligand complexes
- Protein + multiple-ligand complexes
- Protein-protein complexes (homodimers, heterodimer, heterooligomers)
