# Enzyme Sequence Reaction Pair (EnzSRP) dataset build script.

This repository provides scripts for creating the enzyme-sequence-reaction-pair (EnzSRP) dataset.
Follow the steps below to replicate our workflow.

If you want to contribute to this project, please refer to [DEVELOPMENT.md](docs/DEVELOPMENT.md).

## 1. Setup environment

### 1.1. Cloning project

```shell
git clone git@github.com:motonuko/enzsrp.git
```

### 1.2 Create environment

```shell
conda env create -f environment.yml
```

### 1.3. Place `.env` file (recommended)

We manage all data paths using the .env file.
These paths will serve as the default paths, simplifying the process of running scripts.
Please create a .env file in the project's root directory:

```shell
touch .env
```

You can also pass custom paths as arguments each time you run the scripts.

## 2. Download files

### 3.1. Download Uniprot database

To show only entries that is reviewed and has catalytic activity(ies), access the following link
https://www.uniprot.org/uniprotkb?query=%28reviewed%3Atrue%29+AND+%28cc_catalytic_activity%3A*%29

or put the following query text to search box in https://www.uniprot.org/

```
(reviewed:true) AND (cc_catalytic_activity:*)
```

In the result page, click 'download' button and download them with 'JSON' format.

In the end, add the downloaded file path to `.env`

```
ENZSRP_ORIGINAL_UNIPROT_REVIEWED_CATALYTIC_ACTIVITY_JSON_FILE="<path-to-downloaded-file>"
ENZSRP_OUTPUT_DIR="<path-to-project-root>/output"
```

If you need exact same original data, please check https://www.uniprot.org/help/synchronization and ftp site to
get original data.

### 3.2. Download Rhea database

Access https://www.rhea-db.org/help/download and download the following files.

- Reactions in `RXN` format
- `rhea-directions.tsv`
- `rhea2metacyc.tsv`

In the end, add the downloaded files paths to `.env`

```
ENZSRP_ORIGINAL_RHEA_DIRECTIONS_FILE="<path-to-downloaded-file-dir>/rhea-directions.tsv"
ENZSRP_ORIGINAL_RHEA2METACYC_FILE="<path-to-downloaded-file-dir>/rhea2metacyc.tsv"
ENZSRP_ORIGINAL_RHEA_RXN_DIR="<path-to-downloaded-dir-parent-dir>/rxn"
```

### 3.3. Download UniParc id mapping

Some catalytic activity information in UniProt contains isoforms. To get isoform sequence we need UniProc data.

Add the downloaded files paths to `.env`. This will be used as default output path in the following steps.

```
ENZSRP_ISOFORM_UNIPARC_ID_MAPPING_FILE="<path-to-downloaded-idmapping_NNNN_NN_NN_isoform_uniparc.json>/"
```

```shell
python bin/create_dataset.py download-isoform-id-uniparc-mapping
```

### 3.4. Download MetaCyc dataset (recommended)

To get reaction direction, we primarily use physiological reaction annotation.
However, many items lacks of physiological annotation.
We complete the lacking by using MetaCyc database data.

https://metacyc.org/

27.0 https://www.metacyc.org/release-notes.shtml

By ignoring all direction undefined activities or
`--use-undirected-rxn` option at enzyme-reaction pair dataset creation script,
that make script to treat all undirected direction as forward direction.

[//]: # (TODO: metacyc lisence)
https://metacyc.org/publications.shtml

In the end, add the downloaded files paths to `.env`

```
ENZSRP_ORIGINAL_METACYC_REACTIONS_DAT_FILE="<path-to-metacyc-reactions.dat-file>"
```

### 3.5. Download M-CSA dataset (optional)

Download from: https://www.ebi.ac.uk/thornton-srv/m-csa/download/

```
ENZSRP_ORIGINAL_MCSA_DATA_DIR="<path-to-m-csa-dataset>"
```

## 3. Create Enzyme Sequence Reaction Pair (EnzSRP) dataset

```shell
python bin/create_dataset.py build-enzyme-reaction-dataset --allow-non-exp-evidence
```

**IMPORTANT:** The values in the `rxn_evidence` and `phy_rxn_evidence` columns of the output file are randomized.
If you need to confirm that the generated file is identical to our reference file,
please use `tests/integration/check_output_file_content.py`.
