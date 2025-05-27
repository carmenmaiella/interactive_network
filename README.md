This network representation is based on information extracted from AlphaBridge, a tool that uses output provided by AlphaFold3, specifically the confusion matrix, to generate various plots and files. These outputs help automate the process of highlighting which interactions in the predicted structure model are confident, making it easier and more straightforward for users to understand.

To learn more, see the GitHub repository: [AlphaBridge GitHub](https://github.com/PDB-REDO/AlphaBridge)

# Running the alphabridge_network.sh

The `alphabridge_network.sh` script runs the `define_interface.py` module from AlphaBridge and uses its output to run the `network_int.py` script. This provides different ways to visualize confident physical interactions between proteins—such as AlphaBridge plots and a network graph, which can be made interactive.
			 
## Installation

Required packages are listed in the `environment.yml` file.

To create and activate the Conda environment:

```bash
conda env create -f environment.yml
conda activate Network
```

## Run the whole pipeline

1. Prepare Your Data

    Download and unzip the AlphaFold3 output ZIP file for a structural prediction of the complex you're interested in.

2. Clone [AlphaBridge repository](https://github.com/PDB-REDO/AlphaBridge)

From the `interactive_network` working directory, run the main script:
```bash
bash alphabridge_network.sh path-to-AlphaFold3-folder path-to-final-output-folder threshold
```
DEMO EXAMPLE
```bash
bash alphabridge_network.sh path-to-AlphaFold3-folder path-to-final-output-folder threshold
```
Threshold: Use one of the following values for threshold:

0.5 → low confidence

0.75 → default confidence

0.9 → high confidence

### Important: CONFIG FILE

The `alphabridge_network.sh` script is designed to use the path to AlphaBridge defined in the `config_alphabridge_path.ini` file.

## Output files

For more information about the output of alphabridge:  https://github.com/PDB-REDO/AlphaBridge

The output includes `.json` files representing two types of networks for each threshold value.

#### Protein Network

- **Node** = protein  
- **Edge** = one or more interfaces predicted by AlphaBridge

#### Interface Network

- **Node** = specific contact region of a protein involved in physical interaction  
  - Nodes from the same protein share the same color
- **Edges**:
  - **Between different proteins** → interaction between contact areas of two different proteins (black edges)
  - **Within the same protein** → highlights nodes from the same protein (edges match node color)


## Interactive Network

A JavaScript script is available that, together with a CSS file, can render an interactive network graph locally using the **D3.js** library.

- **Input**: the `.json` network file
- **Output**: an interactive graph displayed in a local webpage

This allows you to explore physical protein interactions visually and dynamically.
