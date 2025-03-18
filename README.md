This project is based on the information extract from alphabridge, which is a tool that use the output provided by Alphafold3 which  different plots and files useful for the identification of interacting interfaces, if you want to know more about this, here is the link of the github repository: put the link

	running the”mainscriptname”.sh
		“xxx is a script that run the define_interface.py of alphabridge and uses the given output for running the network_int.py script,  providing different way to visualize confident physical interaction between proteins, usign a network, that could be made interactive ”
			
installation

	Packages required to run the find_interface script are listed in the conda environment file environment.yml

Conda environment can be activated with the instruction:
conda env create -f environment.yml

conda activate Network


	RUNNING THE WHOLE  PIPELINE LOCALLY (IF YOU HAVE NOT RUN ALPHABRIDGE YET)
How to run the script
supposedly that you have tun your importat/interesnting/idk  prediction of alphafold3, that could be either:
-structural prediction of all possibile binary interaction of a set of protein of interest
-structural prediction of the complex of your interest 

 Download and unzip the folder from the AlphaFold server.

this pipeline INCLUE running alphabridge, so you need to clone the github repository (put link aain)

moreover, alphabridge allow you to have 3 different level of confidence for the predicted contact interaction, (0.5, 0.75 and 0.9), but for now the network is only provided using only one threshold, so select the one that you preferer between those three options 

To run the script from the interactive_network working directory, use

bash alphabridge_network.sh path-to-AlphaFold3-folder path-to-final-output-folder threshold

DEMO EXAMPLE
bash alphabridge_network.sh path-to-AlphaFold3-folder path-to-final-output-folder threshold

	RUNNING THE NETWORK  PIPELINE LOCALLY (IF YOU HAVE ALREADY RUN ALPHABRIDGE)

How to run the scirpt

To run the script from the interactive_network working directory, use:

for runnng the network of the whole compelx prediction

python3 src/binary_interactionspy -i path-to-AlphaFold3-folder -o path-to-final-output-folder threshold -t threshold

for runnng the network of the whole compelx prediction

python3 src/nework_int.py -i path-to-AlphaFold3-folder -o path-to-final-output-folder threshold -t threshold

The output of alphabridge will be saved in the same folder provided as an input, but the output of the main pipeline will be saved in the desired output folder

IMPORTANT: CONFIG FILE, → the .sh script is written using my personal pathway where I had installed the main alphabridge script, for change this, chang it in the confing file, it will autmatically change the main script 

OUTPUT FILE
for more information about the output of alphabridg, e: put the link

	INTERACTIVE NETWORK: a json file containing 2 different network:
		LOGIC OF THE NETWORK:
			-nodes: is NOT a whole protein but a specific contact area of that protein that is involved in a physical interaction whithin another protein (nodes that belong to the same protein are colored with the same color)
			-edges
-tra diverse proteine → interfaces/physical itneraction between two contact area of different protein (colored in black)
-tra stesse proteine → highlight the nodes that belongs to the same protein (same color of the nodes)
		-NETWORK WITH MORE NODES/NOT LOGIC MERGED → directly take the information form alphabridge and directly display them as a network
		-NETWORK WITH LESS NODES/MERGED → process the information of alphabridge in a way that if different residues form the the same protein are involved in different interfaces/interaction, they will be merged in the same node

HOW CAN YOU VISUALIZE THE NETWORK:
	-d3 web + js + explanation
	-html ???