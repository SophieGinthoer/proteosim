# Sophie's Proteosim


This is a python package simulating a bottom-up proteomics workflow.


## Workflow

**Input**:\
protein sequences (FASTA format)

**Processing**:\
protein digestion\
LC RT prediction\
RT window selection\
MS1 m/z calculation\
MS2 peptide fragmentation

**Output**:\
peptide & fragment m/z values\
chromatographic profiles\
simulated MS1/2 spectra


## Installation

```bash
git clone https://github.com/SophieGinthoer/proteosim.git
cd proteosim
pip install -r requirements.txt
pip install .
```

Requires Python 3.10 or newer.


## Modules & functions

### `file_handling`
`read_fasta` reads protein FASTA file into dict of protein IDs & sequences

### `protein_digestion`
`enzyme_cleavage_patterns` defines regular expressions for Trypsin, LysC, LysN & ArgC\
`digest_protein_sequence` digests protein *in silico* & filters peptides by length\
`digest_protein_collection` digests every protein in a collection\
`compute_sequence_coverage` computes percentage of protein covered by peptides

### `liquid_chromatography`
`predict_lc_retention_times` predicts relative RTs\
`plot_retention_time` plots histogram of RTs as chromatographic profile\
`select_retention_time_window` selects peptides eluting within RT window

### `mass_spectra_simulation`
`amino_acid_mass_dalton` masses of aa residue in Da\
`calculate_mol_mass` calculates peptide's molecular mass\
`calculate_mol_mass_collection` calculates molecular masses of peptide list\
`calculate_mz_collection` converts peptide masses to m/z values for given charge\
`plot_spectrum` plots simulated mass spectrum (random intensities)\
`fragment_peptide` generates b- & y-fragments of peptide


## Quick start

```python
import proteosim as ps

proteins = ps.read_fasta("data/sample_proteins.fasta")
peptides = ps.digest_protein_collection(proteins, ps.enzyme_cleavage_patterns["Trypsin"])
```

For complete example, open\
["tutorials/ms_experiment_final.ipynb"](tutorials/ms_experiment_final.ipynb).\
reads FASTA file\
digests proteins with trypsin\
predicts RTs\
selects peptides eluting between 0 & 20\
simulates MS1 & MS2 spectrum of "MATSR".


## Running the tests

```bash
pip install pytest
pytest -v tests
```

Tests also run automatically on GitHub Actions for every push & pull.


## License

MIT