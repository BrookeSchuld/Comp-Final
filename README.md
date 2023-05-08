# Comp-Final

Introduction: 
The purpose of this code is to tag particles as their most probable type. Takes Monte Carlo generated files from https://portal.nersc.gov/project/dune/data/2x2/simulation/edepsim/ and outputs an array of tagged particles identified by their track id. 

Running: 
To run the program supply the file path. The output file will be titled "Tagged.txt". 
The main code is provided in the "main" folder, with additional functions being called in plotting_functions. 

Reading the Output:
Output is a txt file with two columns- particle ID and the most likely type. 

Analysis: 
An example output file is provided showing both the real particle ID and the tagged ID. The success rate was about 80%. Additional information is provided in the writeup document. 




