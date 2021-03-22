#!/usr/bin/env python3
"""
    EffectorP 3.0: prediction of apoplastic and cytoplasmic effectors in fungi and oomycetes 
    Copyright (C) 2021-2022 Jana Sperschneider  
    This program is free software; you can redistribute it and/or modify  
    it under the terms of the GNU General Public License as published by  
    the Free Software Foundation; either version 3 of the License, or     
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    Contact: jana.sperschneider@anu.edu.au or jana.sperschneider@csiro.au
"""
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
import os
import sys
import functions
import subprocess
import errno
import uuid
import shutil
import tempfile
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# Main Program starts here
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
SCRIPT_PATH = sys.path[0]
# Change the path to WEKA to the appropriate location on your computer if necessary
WEKA_PATH = SCRIPT_PATH  + '/weka-3-8-4/weka.jar'
# Check that the path to the WEKA software exists
path_exists = os.access(WEKA_PATH, os.F_OK)
if not path_exists:
    print()
    print("Path to WEKA software does not exist!")
    print("Check the installation and the given path to the WEKA software %s in EffectorP.py (line 40)." % WEKA_PATH)
    print()
    sys.exit(1)
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# Global variables
EFFECTOR_THRESHOLD = 0.5
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------
commandline = sys.argv[1:]
# -----------------------------------------------------------------------------------------------------------
if commandline:
    FASTA_FILE, FUNGAL_MODE, output_file, effector_output, noneffector_output = functions.scan_arguments(commandline)
    # If no FASTA file was provided with the -i option
    if not FASTA_FILE:
        print()
        print('Please specify a FASTA input file using the -i option!')
        functions.usage()
else:
    functions.usage()
# -----------------------------------------------------------------------------------------------------------
# Temporary folder name identifier that will be used to store results
RESULTS_PATH = tempfile.mkdtemp() + '/'
# -----------------------------------------------------------------------------------------------------------
# Check if FASTA file exists
try:
    open(FASTA_FILE, 'r') 
except OSError as e:
    print("Unable to open FASTA file:", FASTA_FILE)  #Does not exist OR no read permissions
    print("I/O error({0}): {1}".format(e.errno, e.strerror))
    sys.exit(1)
# -----------------------------------------------------------------------------------------------------------
# Try to create folder where results will be stored
try:
    os.mkdir(RESULTS_PATH)
except OSError as exception:        
    if exception.errno != errno.EEXIST:
        raise
# -----------------------------------------------------------------------------------------------------------
# Extract the identifiers and sequences from input FASTA file
ORIGINAL_IDENTIFIERS, SEQUENCES = [], []
for identifier, sequence in functions.SimpleFastaParser(open(FASTA_FILE, 'r')):
  ORIGINAL_IDENTIFIERS.append(identifier)
  SEQUENCES.append(sequence)

SEQUENCES = [seq.upper() for seq in SEQUENCES]
# -----------------------------------------------------------------------------------------------------------
print('-----------------')
print()
print("EffectorP 3.0 is running for", len(ORIGINAL_IDENTIFIERS), "proteins given in FASTA file", FASTA_FILE)
print()
# -----------------------------------------------------------------------------------------------------------
# Write new FASTA file with short identifiers so WEKA has safe input
f_output = RESULTS_PATH + 'short_ids.fasta'
SHORT_IDENTIFIERS = functions.write_FASTA_short_ids(f_output, ORIGINAL_IDENTIFIERS, SEQUENCES)
# -----------------------------------------------------------------------------------------------------------
# Write the WEKA arff file for classification of the input FASTA file
weka_input = RESULTS_PATH + 'weka.arff'    
# Ensembl averaging approach, use seq0,seq1,seq2... as keys in case there are duplicate FASTA identifiers
functions.write_weka_input(weka_input, SHORT_IDENTIFIERS, SEQUENCES)
# -----------------------------------------------------------------------------------------------------------
ensembl_votes_cytoplasmic = {}
ensembl_votes_apoplastic = {}
# -----------------------------------------------------------------------------------------------------------
# Call WEKA models for classification of input FASTA file
# -----------------------------------------------------------------------------------------------------------
print('Ensemble classification')

if FUNGAL_MODE == True:
    ensembl_votes_cytoplasmic = functions.get_model_predictions(WEKA_PATH, RESULTS_PATH, functions.models_bayes_cytoplasmic_fungionly, 'weka.classifiers.bayes.NaiveBayes', ensembl_votes_cytoplasmic, ORIGINAL_IDENTIFIERS, SEQUENCES)
    print('25 percent done...')

    ensembl_votes_cytoplasmic = functions.get_model_predictions(WEKA_PATH, RESULTS_PATH, functions.models_J48_cytoplasmic_fungionly, 'weka.classifiers.trees.J48', ensembl_votes_cytoplasmic, ORIGINAL_IDENTIFIERS, SEQUENCES)
    print('50 percent done...')


if FUNGAL_MODE == False:
    ensembl_votes_cytoplasmic = functions.get_model_predictions(WEKA_PATH, RESULTS_PATH, functions.models_bayes_cytoplasmic, 'weka.classifiers.bayes.NaiveBayes', ensembl_votes_cytoplasmic, ORIGINAL_IDENTIFIERS, SEQUENCES)
    print('25 percent done...')

    ensembl_votes_cytoplasmic = functions.get_model_predictions(WEKA_PATH, RESULTS_PATH, functions.models_J48_cytoplasmic, 'weka.classifiers.trees.J48', ensembl_votes_cytoplasmic, ORIGINAL_IDENTIFIERS, SEQUENCES)
    print('50 percent done...')

ensembl_votes_apoplastic = functions.get_model_predictions(WEKA_PATH, RESULTS_PATH, functions.models_bayes_apoplastic, 'weka.classifiers.bayes.NaiveBayes', ensembl_votes_apoplastic, ORIGINAL_IDENTIFIERS, SEQUENCES)
print('75 percent done...')

ensembl_votes_apoplastic = functions.get_model_predictions(WEKA_PATH, RESULTS_PATH, functions.models_J48_apoplastic, 'weka.classifiers.trees.J48', ensembl_votes_apoplastic, ORIGINAL_IDENTIFIERS, SEQUENCES)
#-------------------------------------------------------------- 
print('All done.')
print()
#--------------------------------------------------------------
# Soft voting on whether a protein is a cytoplasmic/apoplastic effector
#--------------------------------------------------------------
ensemble_predictions, predicted_effectors, predicted_noneffectors, cyto_effectors, apo_effectors, cyto_apo_effectors, apo_cyto_effectors = functions.get_effector_predictions(ORIGINAL_IDENTIFIERS, SEQUENCES, EFFECTOR_THRESHOLD, ensembl_votes_cytoplasmic, ensembl_votes_apoplastic)  
count_cytoplasmic_effectors = len(cyto_effectors) + len(cyto_apo_effectors)
count_apoplastic_effectors = len(apo_effectors) + len(apo_cyto_effectors)
#--------------------------------------------------------------
print(functions.short_output_screen(ensemble_predictions, cyto_effectors, apo_effectors, cyto_apo_effectors, apo_cyto_effectors))

# If user wants the stdout output directed to a specified file
if output_file:
    with open(output_file, 'w') as out:
        # Output predictions for all proteins as tab-delimited table
        output_string = functions.short_output_file(ensemble_predictions, cyto_effectors, apo_effectors, cyto_apo_effectors, apo_cyto_effectors)
        out.writelines(output_string)              
    print('EffectorP results were saved to output file:', output_file) 

print('-----------------')
print(len(ORIGINAL_IDENTIFIERS), 'proteins were provided as input in the following file:', FASTA_FILE)
print('-----------------')
print('Number of predicted effectors:', (count_cytoplasmic_effectors+count_apoplastic_effectors))
print('Number of predicted cytoplasmic effectors:', count_cytoplasmic_effectors)
print('Number of predicted apoplastic effectors:', count_apoplastic_effectors)
print('-----------------')
print(round(100.0*(count_cytoplasmic_effectors+count_apoplastic_effectors)/len(ORIGINAL_IDENTIFIERS),1), 'percent are predicted effectors.')
print(round(100.0*count_cytoplasmic_effectors/len(ORIGINAL_IDENTIFIERS),1), 'percent are predicted cytoplasmic effectors.')
print(round(100.0*count_apoplastic_effectors/len(ORIGINAL_IDENTIFIERS),1), 'percent are predicted apoplastic effectors.')
if FUNGAL_MODE == True:
    print('-----------------')
    print('NOTE: EffectorP was run in fungal mode.')
print('-----------------')
# -----------------------------------------------------------------------------------------------------------
# If the user additionally wants to save the predicted effectors in a provided FASTA file
if effector_output:
    with open(effector_output, 'w') as f_output:
        for effector, prob_cytoplasmic, prob_apoplastic, sequence in predicted_effectors:
            if prob_apoplastic >= EFFECTOR_THRESHOLD and prob_cytoplasmic < EFFECTOR_THRESHOLD:
                f_output.writelines('>' + effector + ' | Apoplastic effector probability: ' + str(prob_apoplastic) + '\n')
                f_output.writelines(sequence + '\n')  
            if prob_apoplastic >= EFFECTOR_THRESHOLD and prob_cytoplasmic >= EFFECTOR_THRESHOLD and prob_apoplastic > prob_cytoplasmic:
                f_output.writelines('>' + effector + ' | Apoplastic effector probability: ' + str(prob_apoplastic) + ' | Cytoplasmic effector probability: ' + str(prob_cytoplasmic) + '\n')
                f_output.writelines(sequence + '\n')                 
            if prob_cytoplasmic >= EFFECTOR_THRESHOLD and prob_apoplastic < EFFECTOR_THRESHOLD:
                f_output.writelines('>' + effector + ' | Cytoplasmic effector probability: ' + str(prob_cytoplasmic) + '\n')
                f_output.writelines(sequence + '\n')  
            if prob_cytoplasmic >= EFFECTOR_THRESHOLD and prob_apoplastic >= EFFECTOR_THRESHOLD and prob_cytoplasmic >= prob_apoplastic:
                f_output.writelines('>' + effector + ' | Cytoplasmic effector probability: ' + str(prob_cytoplasmic) + ' | Apoplastic effector probability: ' + str(prob_apoplastic) + '\n')
                f_output.writelines(sequence + '\n')  

if noneffector_output:
    with open(noneffector_output, 'w') as f_output:
        for effector, prob, sequence in predicted_noneffectors:
            f_output.writelines('>' + effector + ' | Non-effector probability: ' + str(prob) + '\n')
            f_output.writelines(sequence + '\n')  
# -----------------------------------------------------------------------------------------------------------
# Clean up and delete temporary folder that was created
shutil.rmtree(RESULTS_PATH)
# -----------------------------------------------------------------------------------------------------------
try:
    sys.stdout.close()
except:
    pass
try:
    sys.stderr.close()
except:
    pass
# -----------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------