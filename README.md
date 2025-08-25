# Enzyme Sequence Reaction Pair (EnzSRP) Dataset Build Scripts.

This repository contains scripts to generate the enzyme sequence reaction pair (EnzSRP) dataset.
Follow the steps below to replicate our workflow.

If you’d like to contribute, please see [DEVELOPMENT.md](docs/DEVELOPMENT.md).

## 1. Environment Setup

### 1.1. Clone the Repository

```shell
git clone git@github.com:motonuko/enzsrp.git
```

### 1.2 Create the Conda Environment

```shell
conda env create -f environment.yml
```

### 1.3. Configure the `.env` File (recommended)

We use a `.env` file to manage data paths.
Defining them here ensures consistent defaults and makes running scripts easier.
Create the file in the project's root directory:

```shell
touch .env
```

You may also override paths by passing them as arguments when running scripts.

## 2. Download Files

### 3.1. Download Entries from UniprotKB

To display only entries that are reviewed and have catalytic activity, use the following link:
https://www.uniprot.org/uniprotkb?query=%28reviewed%3Atrue%29+AND+%28cc_catalytic_activity%3A*%29

Alternatively, you can paste the following query into the search box on https://www.uniprot.org/

```
(reviewed:true) AND (cc_catalytic_activity:*)
```

On the results page, click the Download button and select the JSON format.

Finally, add the path of the downloaded file to your `.env` file.

```
ENZSRP_ORIGINAL_UNIPROT_REVIEWED_CATALYTIC_ACTIVITY_JSON_FILE="<path-to-downloaded-file>"
ENZSRP_OUTPUT_DIR="<path-to-project-root>/output"
```

If you need the exact original data, please refer
to [UniProt synchronization](https://www.uniprot.org/help/synchronization)
and use the FTP site to download the data directly.

### 3.2. Download Rhea Database

Access https://www.rhea-db.org/help/download and download the following files.

- Reactions in `RXN` format
- `rhea-directions.tsv`
- `rhea2metacyc.tsv`

Finally, add the paths of the downloaded files to your `.env` file.

```
ENZSRP_ORIGINAL_RHEA_DIRECTIONS_FILE="<path-to-downloaded-file-dir>/rhea-directions.tsv"
ENZSRP_ORIGINAL_RHEA2METACYC_FILE="<path-to-downloaded-file-dir>/rhea2metacyc.tsv"
ENZSRP_ORIGINAL_RHEA_RXN_DIR="<path-to-downloaded-dir-parent-dir>/rxn"
```

### 3.3. Download UniParc ID Mapping

Some catalytic activity information in UniProt includes isoforms. To obtain the isoform sequences, we need UniProc data.

```shell
python bin/create_dataset.py download-isoform-id-uniparc-mapping
```

Add the paths of the downloaded files to the `.env` file.

```
ENZSRP_ISOFORM_UNIPARC_ID_MAPPING_FILE="<path-to-downloaded-idmapping_NNNN_NN_NN_isoform_uniparc.json>/"
```

### 3.4. Download MetaCyc Dataset (recommended)

To determine the reaction directions, we primarily use physiological reaction annotation.
However, many items lack such information.
We fill these gaps by using data from the MetaCyc database (the `reactions.dat` file).

[MetaCyc](https://metacyc.org/), [27.0](https://www.metacyc.org/release-notes.shtml)


In the end, add the downloaded files paths to `.env`

```
ENZSRP_ORIGINAL_METACYC_REACTIONS_DAT_FILE="<path-to-reactions.dat-file>"
```

## 3. Create Enzyme Sequence Reaction Pair (EnzSRP) Dataset

To generate the EnzSRP dataset, execute the following script.

```shell
python bin/create_dataset.py build-enzyme-reaction-dataset --allow-non-exp-evidence
```

**IMPORTANT:** The values in the `rxn_evidence` and `phy_rxn_evidence` columns of the output file are randomized.
If you need to confirm that the generated file is identical to our reference file,
please use `tests/integration/check_output_file_content.py`.
