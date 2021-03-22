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
import os
import sys
import subprocess
import io
import getopt
# -----------------------------------------------------------------------------------------------------------
# Global variables
# -----------------------------------------------------------------------------------------------------------
ARFF_HEADER = '''@RELATION effectors
@ATTRIBUTE A NUMERIC
@ATTRIBUTE C NUMERIC
@ATTRIBUTE D NUMERIC
@ATTRIBUTE E NUMERIC
@ATTRIBUTE F NUMERIC
@ATTRIBUTE G NUMERIC
@ATTRIBUTE H NUMERIC
@ATTRIBUTE I NUMERIC
@ATTRIBUTE K NUMERIC
@ATTRIBUTE L NUMERIC
@ATTRIBUTE M NUMERIC
@ATTRIBUTE N NUMERIC
@ATTRIBUTE P NUMERIC
@ATTRIBUTE Q NUMERIC
@ATTRIBUTE R NUMERIC
@ATTRIBUTE S NUMERIC
@ATTRIBUTE T NUMERIC
@ATTRIBUTE V NUMERIC
@ATTRIBUTE W NUMERIC
@ATTRIBUTE Y NUMERIC
@ATTRIBUTE MolecularWeight NUMERIC
@ATTRIBUTE PosCharge NUMERIC
@ATTRIBUTE NegCharge NUMERIC
@ATTRIBUTE Exposed NUMERIC
@ATTRIBUTE Hydrophobicity NUMERIC
@ATTRIBUTE polarity NUMERIC
@ATTRIBUTE flexibility NUMERIC
@ATTRIBUTE aromatic NUMERIC
@ATTRIBUTE polar NUMERIC
@ATTRIBUTE disorder NUMERIC
@ATTRIBUTE Bulky NUMERIC
@ATTRIBUTE Alpha NUMERIC
@ATTRIBUTE Beta NUMERIC
@ATTRIBUTE Coil NUMERIC
@ATTRIBUTE class {effector,non-effector}
@DATA
'''
# -----------------------------------------------------------------------------------------------------------
SCRIPT_PATH = sys.path[0]

models_bayes_cytoplasmic_fungionly = [SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Mycorrhizal_Bayes//trainingdata_iteration85_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Mycorrhizal_Bayes//trainingdata_iteration75_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Mycorrhizal_Bayes//trainingdata_iteration100_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Mycorrhizal_Bayes//trainingdata_iteration41_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Mycorrhizal_Bayes//trainingdata_iteration31_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Pathogens_Bayes//trainingdata_iteration92_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Pathogens_Bayes//trainingdata_iteration16_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Pathogens_Bayes//trainingdata_iteration84_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Pathogens_Bayes//trainingdata_iteration2_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Pathogens_Bayes//trainingdata_iteration96_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Saprophytes_Bayes//trainingdata_iteration34_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Saprophytes_Bayes//trainingdata_iteration64_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Saprophytes_Bayes//trainingdata_iteration88_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Saprophytes_Bayes//trainingdata_iteration23_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Saprophytes_Bayes//trainingdata_iteration32_ratio3.model']

models_J48_cytoplasmic_fungionly = [SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Mycorrhizal_J48//trainingdata_iteration77_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Mycorrhizal_J48//trainingdata_iteration61_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Mycorrhizal_J48//trainingdata_iteration95_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Mycorrhizal_J48//trainingdata_iteration46_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Mycorrhizal_J48//trainingdata_iteration56_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Pathogens_J48//trainingdata_iteration85_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Pathogens_J48//trainingdata_iteration58_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Pathogens_J48//trainingdata_iteration87_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Pathogens_J48//trainingdata_iteration25_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Pathogens_J48//trainingdata_iteration97_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Saprophytes_J48//trainingdata_iteration85_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Saprophytes_J48//trainingdata_iteration13_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Saprophytes_J48//trainingdata_iteration22_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Saprophytes_J48//trainingdata_iteration91_ratio3.model',
SCRIPT_PATH + '/TrainingData_CytoplasmicFungiOnly_Saprophytes_J48//trainingdata_iteration21_ratio3.model']

models_bayes_cytoplasmic = [SCRIPT_PATH + '/TrainingData_Cytoplasmic_Mycorrhizal_Bayes//trainingdata_iteration39_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Mycorrhizal_Bayes//trainingdata_iteration17_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Mycorrhizal_Bayes//trainingdata_iteration76_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Mycorrhizal_Bayes//trainingdata_iteration74_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Mycorrhizal_Bayes//trainingdata_iteration78_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Pathogens_Bayes//trainingdata_iteration79_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Pathogens_Bayes//trainingdata_iteration1_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Pathogens_Bayes//trainingdata_iteration32_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Pathogens_Bayes//trainingdata_iteration70_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Pathogens_Bayes//trainingdata_iteration77_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Saprophytes_Bayes//trainingdata_iteration14_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Saprophytes_Bayes//trainingdata_iteration71_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Saprophytes_Bayes//trainingdata_iteration44_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Saprophytes_Bayes//trainingdata_iteration62_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Saprophytes_Bayes//trainingdata_iteration64_ratio3.model']

models_J48_cytoplasmic = [SCRIPT_PATH + '/TrainingData_Cytoplasmic_Mycorrhizal_J48//trainingdata_iteration36_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Mycorrhizal_J48//trainingdata_iteration14_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Mycorrhizal_J48//trainingdata_iteration63_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Mycorrhizal_J48//trainingdata_iteration60_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Mycorrhizal_J48//trainingdata_iteration29_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Pathogens_J48//trainingdata_iteration67_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Pathogens_J48//trainingdata_iteration92_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Pathogens_J48//trainingdata_iteration78_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Pathogens_J48//trainingdata_iteration1_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Pathogens_J48//trainingdata_iteration34_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Saprophytes_J48//trainingdata_iteration36_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Saprophytes_J48//trainingdata_iteration71_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Saprophytes_J48//trainingdata_iteration67_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Saprophytes_J48//trainingdata_iteration91_ratio3.model',
SCRIPT_PATH + '/TrainingData_Cytoplasmic_Saprophytes_J48//trainingdata_iteration70_ratio3.model']

models_bayes_apoplastic = [SCRIPT_PATH + '/TrainingData_Apoplastic_Animal_Bayes//trainingdata_iteration19_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Animal_Bayes//trainingdata_iteration70_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Animal_Bayes//trainingdata_iteration78_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Animal_Bayes//trainingdata_iteration80_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Animal_Bayes//trainingdata_iteration92_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Pathogens_Bayes//trainingdata_iteration40_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Pathogens_Bayes//trainingdata_iteration31_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Pathogens_Bayes//trainingdata_iteration98_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Pathogens_Bayes//trainingdata_iteration76_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Pathogens_Bayes//trainingdata_iteration53_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Saprophytes_Bayes//trainingdata_iteration61_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Saprophytes_Bayes//trainingdata_iteration79_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Saprophytes_Bayes//trainingdata_iteration29_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Saprophytes_Bayes//trainingdata_iteration91_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Saprophytes_Bayes//trainingdata_iteration32_ratio3.model']

models_J48_apoplastic = [SCRIPT_PATH + '/TrainingData_Apoplastic_Animal_J48//trainingdata_iteration82_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Animal_J48//trainingdata_iteration19_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Animal_J48//trainingdata_iteration46_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Animal_J48//trainingdata_iteration47_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Animal_J48//trainingdata_iteration100_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Pathogens_J48//trainingdata_iteration26_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Pathogens_J48//trainingdata_iteration98_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Pathogens_J48//trainingdata_iteration49_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Pathogens_J48//trainingdata_iteration9_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Pathogens_J48//trainingdata_iteration78_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Saprophytes_J48//trainingdata_iteration25_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Saprophytes_J48//trainingdata_iteration95_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Saprophytes_J48//trainingdata_iteration46_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Saprophytes_J48//trainingdata_iteration92_ratio3.model',
SCRIPT_PATH + '/TrainingData_Apoplastic_Saprophytes_J48//trainingdata_iteration77_ratio3.model']
# -----------------------------------------------------------------------------------------------------------
# Hydrophobicity (Fauchere and Pliska, 1983)
HYDRO_DIC = {
'R': -1.01,
'K': -0.99,  
'D': -0.77,  
'E': -0.64,  
'N': -0.6,
'Q': -0.22,  
'S': -0.04,  
'G': -0.0,  
'H': 0.13,  
'T': 0.26,  
'A': 0.31, 
'P': 0.72,  
'Y': 0.96,  
'V': 1.22,
'C': 1.54,  
'L': 1.7,  
'F': 1.79,  
'I': 1.8,  
'M': 1.23 ,  
'W': 2.25}

# Taken from http://www.cprofiler.org/help.html
# Surface exposure (Janin, 1979), these are free energy values
EXPOSED_DIC = {
'A': 0.3, 
'R': -1.4,
'N': -0.5,
'D': -0.6,  
'C': 0.9,  
'Q': -0.7,  
'E': -0.7,  
'G': 0.3,  
'H': -0.1,  
'I': 0.7,  
'L': 0.5,  
'K': -1.8,  
'M': 0.4,  
'F': 0.5,  
'P': -0.3,  
'S': -0.1,  
'T': -0.2,  
'W': 0.3,  
'Y': -0.4,  
'V': 0.6}

# Flexibility (Vihinen et al., 1994)
FLEX_DIC = {
'A': 0.984, 
'R': 1.008,
'N': 1.048,
'D': 1.068,  
'C': 0.906,  
'Q': 1.037,  
'E': 1.094,  
'G': 1.031,  
'H': 0.950,  
'I': 0.927,  
'L': 0.935,  
'K': 1.102,  
'M': 0.952,  
'F': 0.915,  
'P': 1.049,  
'S': 1.046,  
'T': 0.997,  
'W': 0.904,  
'Y': 0.929,  
'V': 0.931}

# Alpha helix frequency (Nagano, 1973)
ALPHA_DIC = {
'Y': 0.63,  
'P': 0.70,  
'G': 0.72,
'N': 0.77,
'S': 0.78,
'R': 0.83,
'T': 0.87,
'C': 0.94,
'I': 0.94,
'V': 0.97,
'D': 1.00,
'W': 1.06,
'Q': 1.10,
'L': 1.23,
'K': 1.23,
'M': 1.23,
'F': 1.23,
'A': 1.29,
'H': 1.29,
'E': 1.54}

# Beta structure frequency (Nagano, 1973)
BETA_DIC = {
'Y': 1.07,  
'P': 0.75,  
'G': 0.9,
'N': 0.72,
'S': 0.77,
'R': 0.67,
'T': 1.23,
'C': 1.13,
'I': 1.54,
'V': 1.41,
'D': 0.9,
'W': 1.13,
'Q': 1.18,
'L': 1.26,
'K': 0.81,
'M': 1.29,
'F': 1.37,
'A': 0.96,
'H': 0.87,
'E': 0.33}

# Coil propensity (Nagano, 1973)
COIL_DIC = {
    'F' : 0.58,
    'M' : 0.62,
    'L' : 0.63,  
    'A' : 0.72,
    'E' : 0.75,
    'H' : 0.76,
    'I' : 0.8,  
    'Q' : 0.81,  
    'V' : 0.83,  
    'K' : 0.84,  
    'W' : 0.87,  
    'C' : 1.01,  
    'T' : 1.03,  
    'D' : 1.04,
    'R' : 1.33,  
    'S' : 1.34, 
    'G' : 1.35,  
    'Y' : 1.35,  
    'N' : 1.38,  
    'P' : 1.43}

# Polarity (Zimmerman et al., 1968)
POLARITY_DIC = {
'Y': 1.61,  
'P': 1.58,  
'G': 0.0,
'N': 3.38,
'S': 1.67,
'R': 52.0,
'T': 1.66,
'C': 1.48,
'I': 0.13,
'V': 0.13,
'D': 49.7,
'W': 2.1,
'Q': 3.53,
'L': 0.13,
'K': 49.5,
'M': 1.43,
'F': 0.35,
'A': 0.0,
'H': 51.6,
'E': 49.9}

# Disorder propensity (Dunker et al., 2001)
DISORDER_DIC = {
'A': 1.0,
'R': 1.0,
'S': 1.0,
'Q': 1.0,
'E': 1.0,
'G': 1.0,
'K': 1.0,
'P': 1.0,
'D': 0.0,
'H': 0.0,
'M': 0.0,
'T': 0.0,
'N': -1.0,
'C': -1.0,
'I': -1.0,
'L': -1.0,
'F': -1.0,
'W': -1.0,
'Y': -1.0,
'V': -1.0}

# Bulkiness (Zimmerman et al., 1968)
BULKY_DIC = {
'G' : 3.4,
'S' : 9.47, 
'A' : 11.5,     
'D' : 11.68,    
'N' : 12.82,    
'C' : 13.46,    
'E' : 13.57,    
'H' : 13.69,    
'R' : 14.28,
'Q' : 14.45,    
'K' : 15.71,    
'T' : 15.77,
'M' : 16.25,    
'P' : 17.43,
'Y' : 18.03,    
'F' : 19.8, 
'I' : 21.4,
'L' : 21.4,     
'V' : 21.57,    
'W' : 21.67}                                                  

# Charged amino acids, 1 are positively charged residues (K, R); -1 are negatively charged residues (D, E)
CHARGE_DIC = {
    'K' : 1,
    'R' : 1,
    'D' : -1,
    'E' : -1}

# Polarity (Zimmerman et al., 1968)
POLARITY_DIC = {
    'A' : 0.0,
    'G' : 0.0,
    'I' : 0.13,
    'L' : 0.13,
    'V' : 0.13,
    'F' : 0.35,
    'M' : 1.43,
    'C' : 1.48,
    'P' : 1.58,
    'Y' : 1.61,
    'T' : 1.66,
    'S' : 1.67,
    'W' : 2.1,
    'N' : 3.38,
    'Q' : 3.53,
    'K' : 49.5,
    'D' : 49.7,
    'E' : 49.9,
    'H' : 51.6,
    'R' : 52.0}

MOLECULAR_WEIGHT_DIC = {
'A': 71.0788,
'B': 114.5962,
'C': 103.1388,
'D': 115.0886,
'E': 129.1155,
'F': 147.1766,
'G': 57.0519,
'H': 137.1411,
'I': 113.1594,
'J': 113.1594,
'K': 128.1741,
'L': 113.1594,
'M': 131.1926,
'N': 114.1038,
'O': 237.3018,
'P': 97.1167,
'Q': 128.1307,
'R': 156.1875,
'S': 87.0782,
'T': 101.1051,
'U': 150.0388,
'V': 99.1326,
'W': 186.2132,
'X': 118.8860,
'Y': 163.1760,
'Z': 128.6231}    
# -----------------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------------------------------------
def usage():
    """ Function: usage()
        Purpose:  Print helpful information for the user.        
        
        Input:    None.
    
        Return:   Print options for running EffectorP 3.0.       
    """
    print('''
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# EffectorP 3.0: Prediction of apoplastic and cytoplasmic effectors in fungi and oomycetes
# http://effectorp.csiro.au/
# Copyright (C) 2021-2022 Jana Sperschneider.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ''')
    print("Usage for EffectorP 3.0: ")
    print("python EffectorP.py [-options] -i <input_file>")
    print()
    print("where basic options are:")
    print("-f : run in fungal mode")    
    print("-h : show brief help on version and usage")
    print()
    print("options directing output:")
    print("-o <f> : direct tab-delimited output table with predictions to file <f>, not stdout")
    print("-E <f> : save predicted effectors to FASTA file <f>")   
    print("-N <f> : save predicted non-effectors to FASTA file <f>") 
    print()
    print("# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
    print()
    sys.exit()    

    return
# -----------------------------------------------------------------------------------------------------------
def scan_arguments(commandline):
    """ Function: scan_arguments()
        Purpose:  Scan the input options given to the EffectorP program.        
        
        Input:    Input options given by the user.
    
        Return:   Parsed options.
    """
    try:
        opts, args = getopt.getopt(commandline, "hfso:E:N:i:", ["help"])        
    except getopt.GetoptError as err:
        # print help information and exit:
        print(str(err)) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    FASTA_FILE = None
    output_file = None
    effector_output = None
    noneffector_output = None
    FUNGAL_MODE = False

    i_count, o_count, E_count, N_count, P_count = 0, 0, 0, 0, 0
   
    for opt, arg in opts:
        if opt in ("-o"):
            output_file = arg
            o_count += 1
        elif opt in ("-f"):
            FUNGAL_MODE = True            
        elif opt in ("-i"):
            FASTA_FILE = arg
            i_count += 1
        elif opt in ("-E"):
            effector_output = arg
            E_count += 1
        elif opt in ("-N"):
            noneffector_output = arg
            N_count += 1
        elif opt in ("-h", "--help"):
            usage()
        else:
            print()
            print ("Commandline option was supplied that was not recognized!")
            usage()

    if i_count > 1 or o_count > 1 or E_count > 1 or N_count > 1:
       usage()

    return FASTA_FILE, FUNGAL_MODE, output_file, effector_output, noneffector_output
# -----------------------------------------------------------------------------------------------------------
def SimpleFastaParser(handle):
    for line in handle:
        if line[0] == ">":
            title = line[1:].rstrip()
            break

    lines = []
    for line in handle:
        if line[0] == ">":
            yield title, "".join(lines).replace(" ", "").replace("\r", "")
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip())

    yield title, "".join(lines).replace(" ", "").replace("\r", "")
# -----------------------------------------------------------------------------------------------------------
def get_effector_predictions(ORIGINAL_IDENTIFIERS, SEQUENCES, EFFECTOR_THRESHOLD, ensembl_votes_cytoplasmic, ensembl_votes_apoplastic):

    models_cytoplasmic = models_bayes_cytoplasmic + models_J48_cytoplasmic
    models_apoplastic = models_bayes_apoplastic + models_J48_apoplastic

    ensemble_predictions, predicted_effectors, predicted_noneffectors = [], [], []
    cyto_effectors, apo_effectors, cyto_apo_effectors, apo_cyto_effectors = {}, {}, {}, {}

    for index, (ident, seq) in enumerate(zip(ORIGINAL_IDENTIFIERS, SEQUENCES)):

        # Be careful to use this short identifier later, 
        # if all identifiers are equal, predictions will fail if not used
        short_ident = 'protein' + str(index)
        yes_prob_cytoplasmic, no_prob_cytoplasmic = [], []
        yes_prob_apoplastic, no_prob_apoplastic = [], []

        for vote, prob in ensembl_votes_cytoplasmic[short_ident]:

            if vote == 'Non-effector':
                no_prob_cytoplasmic.append(prob)
                yes_prob_cytoplasmic.append(1.0 - prob)                        

            if vote == 'Effector':
                yes_prob_cytoplasmic.append(prob)
                no_prob_cytoplasmic.append(1.0 - prob)        

        for vote, prob in ensembl_votes_apoplastic[short_ident]:

            if vote == 'Non-effector':
                no_prob_apoplastic.append(prob)
                yes_prob_apoplastic.append(1.0 - prob)                        

            if vote == 'Effector':
                yes_prob_apoplastic.append(prob)
                no_prob_apoplastic.append(1.0 - prob) 

        # Soft voting: argmax of the sum of predicted probabilities
        yes_prob_cytoplasmic = round(sum(yes_prob_cytoplasmic)/float(len(models_cytoplasmic)),3)
        no_prob_cytoplasmic = round(sum(no_prob_cytoplasmic)/float(len(models_cytoplasmic)),3)
        yes_prob_apoplastic = round(sum(yes_prob_apoplastic)/float(len(models_apoplastic)),3)
        no_prob_apoplastic = round(sum(no_prob_apoplastic)/float(len(models_apoplastic)),3)     

        cytoplasmic_prediction = False
        apoplastic_prediction = False

        if yes_prob_cytoplasmic >= EFFECTOR_THRESHOLD or yes_prob_apoplastic >= EFFECTOR_THRESHOLD:
            # Is it more likely a cytoplasmic effector
            if yes_prob_cytoplasmic >= yes_prob_apoplastic:

                if yes_prob_apoplastic >= EFFECTOR_THRESHOLD:
                    prediction = 'Cytoplasmic effector (apoplastic effector: ' + str(yes_prob_apoplastic) + ')'
                    prob = yes_prob_cytoplasmic
                    predicted_effectors.append((ident.strip(), yes_prob_cytoplasmic, yes_prob_apoplastic, seq))            
                    cyto_apo_effectors[short_ident] = [yes_prob_cytoplasmic, yes_prob_apoplastic, seq]    
                else:
                    prediction = 'Cytoplasmic effector'
                    prob = yes_prob_cytoplasmic
                    predicted_effectors.append((ident.strip(), yes_prob_cytoplasmic, yes_prob_apoplastic, seq))   
                    cyto_effectors[short_ident] = [yes_prob_cytoplasmic, seq]              

            # Is it more likely an apoplastic effector                
            if yes_prob_apoplastic > yes_prob_cytoplasmic:

                if yes_prob_cytoplasmic >= EFFECTOR_THRESHOLD:
                    prediction = 'Apoplastic effector (cytoplasmic effector: ' + str(yes_prob_cytoplasmic) + ')'
                    prob = yes_prob_apoplastic
                    predicted_effectors.append((ident.strip(), yes_prob_cytoplasmic, yes_prob_apoplastic, seq))     
                    apo_cyto_effectors[short_ident] = [yes_prob_apoplastic, yes_prob_cytoplasmic, seq]     

                else:
                    prediction = 'Apoplastic effector'
                    prob = yes_prob_apoplastic
                    predicted_effectors.append((ident.strip(), yes_prob_cytoplasmic, yes_prob_apoplastic, seq))   
                    apo_effectors[short_ident] = [yes_prob_apoplastic, seq]                                  

        if yes_prob_cytoplasmic < EFFECTOR_THRESHOLD and yes_prob_apoplastic < EFFECTOR_THRESHOLD:
            prediction = 'Non-effector'
            prob = round(min(no_prob_cytoplasmic, no_prob_apoplastic),3)
            predicted_noneffectors.append((ident.strip(), prob, seq))

        ensemble_predictions.append((ident.strip(), prediction, prob, seq))     

    return ensemble_predictions, predicted_effectors, predicted_noneffectors, cyto_effectors, apo_effectors, cyto_apo_effectors, apo_cyto_effectors
# -----------------------------------------------------------------------------------------------------------    
def get_model_predictions(WEKA_PATH, RESULTS_PATH, MODELS, CLASSIFIER, ensembl_votes, ORIGINAL_IDENTIFIERS, SEQUENCES):

    for model in MODELS:
        #--------------------------------------------------------------
        ParamList = ['java', '-cp', WEKA_PATH, CLASSIFIER, '-l', model, '-T', RESULTS_PATH + 'weka.arff', '-p', 'first-last']

        with open(RESULTS_PATH + 'Predictions.txt', 'wb') as out:
            try:
                Process = subprocess.Popen(ParamList, shell=False, stdout=out)
                sts = Process.wait()
                cstdout, cstderr = Process.communicate()

                if Process.returncode:
                    raise Exception("Calling WEKA returned %s"%Process.returncode)
                if cstdout:
                    pass
                elif cstderr:
                    sys.exit()
            except:
                e = sys.exc_info()[1]
                print("Error calling WEKA: %s" % e)
                sys.exit(1)
        #-------------------------------------------------------------- 
        # Parse the WEKA output file
        file_input = RESULTS_PATH + 'Predictions.txt'
        predicted_effectors, predicted_noneffectors, predictions = parse_weka_output(file_input, ORIGINAL_IDENTIFIERS, SEQUENCES)
        
        for index, (ident, prediction, prob, seq) in enumerate(predictions):

            short_ident = 'protein' + str(index)

            if short_ident in ensembl_votes:
                previous_predictions = ensembl_votes[short_ident] 
                ensembl_votes[short_ident] = previous_predictions + [(prediction, prob)]
            else:
                ensembl_votes[short_ident] = [(prediction, prob)]

    return ensembl_votes
# -----------------------------------------------------------------------------------------------------------        
def write_weka_input(weka_input, SHORT_IDENTIFIERS, SEQUENCES):
    """ Function: write_weka_input()
        Purpose:  Given the query identifiers and 
                  protein features, write the input arff file for WEKA. 
              
        Input:    WEKA arff file name, query identifiers.                  
    
        Return:   None. 
    """   
    with open(weka_input, 'w') as f:
        # Create a list of features for each protein
        X = [[] for __ in range(len(SHORT_IDENTIFIERS))]

        for protein_position, TARGET_ID in enumerate(SHORT_IDENTIFIERS):
            TARGET_ID = TARGET_ID.replace('>', '')
            TARGET_ID = TARGET_ID.strip()
            sequence = SEQUENCES[protein_position]

            length = len(sequence)

            # Amino acid frequencies in the sequence
            amino_acid_frequencies = []
            amino_acid_frequencies.append(100.0*sequence.count('A')/length)
            amino_acid_frequencies.append(100.0*sequence.count('C')/length)
            amino_acid_frequencies.append(100.0*sequence.count('D')/length)
            amino_acid_frequencies.append(100.0*sequence.count('E')/length)
            amino_acid_frequencies.append(100.0*sequence.count('F')/length)
            amino_acid_frequencies.append(100.0*sequence.count('G')/length)
            amino_acid_frequencies.append(100.0*sequence.count('H')/length)
            amino_acid_frequencies.append(100.0*sequence.count('I')/length)
            amino_acid_frequencies.append(100.0*sequence.count('K')/length)
            amino_acid_frequencies.append(100.0*sequence.count('L')/length)
            amino_acid_frequencies.append(100.0*sequence.count('M')/length)
            amino_acid_frequencies.append(100.0*sequence.count('N')/length)
            amino_acid_frequencies.append(100.0*sequence.count('P')/length)
            amino_acid_frequencies.append(100.0*sequence.count('Q')/length)
            amino_acid_frequencies.append(100.0*sequence.count('R')/length)
            amino_acid_frequencies.append(100.0*sequence.count('S')/length)
            amino_acid_frequencies.append(100.0*sequence.count('T')/length)
            amino_acid_frequencies.append(100.0*sequence.count('V')/length)
            amino_acid_frequencies.append(100.0*sequence.count('W')/length)
            amino_acid_frequencies.append(100.0*sequence.count('Y')/length)


            molecular_weight = MOLECULAR_WEIGHT(sequence)
            percent_positively_charged, percent_negatively_charged, charge = CHARGE(sequence)
            exposed = EXPOSED(sequence)
            hydrophobicity = HYDROPHOBICITY(sequence)          
            polarity = POLARITY(sequence)
            flexibility = FLEX(sequence)
            aromatic = 100.0*(sequence.count('F') + sequence.count('H') + sequence.count('W') + sequence.count('Y'))/length
            polar = 100.0*(sequence.count('D') + sequence.count('E') + sequence.count('H') + sequence.count('K') + sequence.count('N') + sequence.count('Q') + sequence.count('R') + sequence.count('S') + sequence.count('T') + sequence.count('Z'))/length
            disorder = DISORDER(sequence)
            bulky = BULKY(sequence)
            alpha = ALPHA(sequence)
            beta = BETA(sequence)
            coil = COIL(sequence)

            X[protein_position] = amino_acid_frequencies + [molecular_weight, percent_positively_charged, percent_negatively_charged, exposed] 
            X[protein_position] += [hydrophobicity, polarity, flexibility] + [aromatic, polar, disorder, bulky] + [alpha, beta, coil]

        # Write protein feature data to WEKA arff file
        f.writelines(ARFF_HEADER)
        for index, vector in enumerate(X):
            for feature in vector:
                f.writelines(str(feature) + ',')
            f.writelines('?\n')

    return
# -----------------------------------------------------------------------------------------------------------
def MOLECULAR_WEIGHT(sequence):

    molecular_weight = 0.0

    for aa in sequence:
        if aa.upper() in MOLECULAR_WEIGHT_DIC:
            molecular_weight += MOLECULAR_WEIGHT_DIC[aa.upper()]

    return molecular_weight        
# -----------------------------------------------------------------------------------------------------------
def HYDROPHOBICITY(sequence):

    hydrophobicity = 0

    for aa in sequence:
        if aa.upper() in HYDRO_DIC:
            hydrophobicity += HYDRO_DIC[aa.upper()]

    return hydrophobicity/len(sequence)    
# -----------------------------------------------------------------------------------------------------------
def FLEX(sequence):
    
    flexibility = 0.0
    for aa in sequence:
        if aa.upper() in FLEX_DIC:
            flexibility += FLEX_DIC[aa.upper()]

    return flexibility/len(sequence)    
# -----------------------------------------------------------------------------------------------------------
def CHARGE(sequence):

    positively_charged, negatively_charged, charge = 0, 0, 0

    for aa in sequence:
        if aa.upper() in CHARGE_DIC:
            if CHARGE_DIC[aa.upper()] == 1:
                positively_charged += 1
                charge += 1
            if CHARGE_DIC[aa.upper()] == -1:
                negatively_charged += 1
                charge += -1
        if aa.upper() == 'H':
            charge += 0.5

    return 100.0*(positively_charged)/len(sequence), 100.0*(negatively_charged)/len(sequence), charge/len(sequence)
# -----------------------------------------------------------------------------------------------------------
def POLARITY(sequence):

    polarity = 0

    for aa in sequence:
        if aa.upper() in POLARITY_DIC:
            polarity += POLARITY_DIC[aa.upper()]

    return polarity/len(sequence)
# -----------------------------------------------------------------------------------------------------------
def DISORDER(sequence):
    
    disorder = 0.0
    for aa in sequence:
        if aa.upper() in DISORDER_DIC:
            disorder += DISORDER_DIC[aa.upper()]

    return disorder/len(sequence)
# -----------------------------------------------------------------------------------------------------------
def EXPOSED(sequence):
    
    exposed = 0.0
    for aa in sequence:
        if aa.upper() in EXPOSED_DIC:
            exposed += EXPOSED_DIC[aa.upper()]

    return exposed/len(sequence)
# -----------------------------------------------------------------------------------------------------------
def ALPHA(sequence):

    alpha = 0.0
    for aa in sequence:
        if aa.upper() in ALPHA_DIC:
            alpha += ALPHA_DIC[aa.upper()]

    return alpha/len(sequence)
# -----------------------------------------------------------------------------------------------------------
def BETA(sequence):
    
    beta = 0.0
    for aa in sequence:
        if aa.upper() in BETA_DIC:
            beta += BETA_DIC[aa.upper()]

    return beta/len(sequence)
# -----------------------------------------------------------------------------------------------------------
def COIL(sequence):
    
    coil = 0.0
    for aa in sequence:
        if aa.upper() in COIL_DIC:
            coil += COIL_DIC[aa.upper()]

    return coil/len(sequence)    
# -----------------------------------------------------------------------------------------------------------
def BULKY(sequence):
    
    bulky = 0.0
    for aa in sequence:
        if aa.upper() in BULKY_DIC:
            bulky += BULKY_DIC[aa.upper()]

    return bulky/len(sequence)
# -----------------------------------------------------------------------------------------------------------
def write_FASTA_short_ids(f_output, ORIGINAL_IDENTIFIERS, ORIGINAL_SEQUENCES):
    """ Function: write_FASTA_short_ids()
        Purpose:  Given a list of identifiers and the corresponding list 
                  of sequence, write these to a FASTA file using short
                  identifiers such as protein1, protein2, .... This is 
                  done because some programs like pepstats do not like 
                  long identifier names as input.
              
        Input:    Path to desired FASTA format output file, list of 
                  identifiers and list of corresponding sequences.
    
        Return:   List of short identifiers.
    """

    with open(f_output, 'w') as f:
        SHORT_IDENTIFIERS = []
        # Change identifiers to protein1, protein2, ...
        # and write to temporary file
        SET = zip(ORIGINAL_IDENTIFIERS, ORIGINAL_SEQUENCES)
    
        for index,  (identifier, sequence) in enumerate(SET):
            short_id = '>protein' + str(index)
            SHORT_IDENTIFIERS.append(short_id)
            f.writelines(short_id + '\n')
            f.writelines(sequence + '\n')

    return SHORT_IDENTIFIERS
# -----------------------------------------------------------------------------------------------------------
def parse_weka_output(file_input, ORIGINAL_IDENTIFIERS, SEQUENCES):
    """ Function: parse_weka_output()
        Purpose:  Given the WEKA output file and the query identifiers and sequences, 
                  parse the predicted class for each protein from the WEKA output. 
                  Write the predicted effectors to a FASTA file.
              
        Input:    WEKA output file and the query identifiers and sequences.                  
    
        Return:   The set of predicted effectors only as well as all predictions. 
    """    
    predicted_effectors, predicted_noneffectors, predictions = [], [], []

    with open(file_input) as f:

        content = f.readlines()

        content_start = content.index('    inst#     actual  predicted error prediction (A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,MolecularWeight,PosCharge,NegCharge,Exposed,Hydrophobicity,polarity,flexibility,aromatic,polar,disorder,Bulky,Alpha,Beta,Coil)\n')

        content = content[content_start + 1:]

        for line in content:
            if line.strip():
                position = line.split()[0]
                prediction = line.split()[2]
                prob = float(line.split()[3])        
 
                # WEKA output counts from position 1, our identifiers are counted from zero
                identifier = ORIGINAL_IDENTIFIERS[int(position) - 1]
                sequence = SEQUENCES[int(position) - 1]

                if 'non-eff' in prediction:                               
                    noneffector = identifier.strip()
                    noneffector = noneffector.replace('>', '')  
                    predictions.append((noneffector, 'Non-effector', prob, sequence))
                    predicted_noneffectors.append((noneffector, prob, sequence))
                else:                    
                    effector = identifier.strip()
                    effector = effector.replace('>', '')                                               
                    predictions.append((effector, 'Effector', prob, sequence))
                    # Append predicted effector to list of predicted effectors
                    predicted_effectors.append((effector, prob, sequence))

    return predicted_effectors, predicted_noneffectors, predictions
# -----------------------------------------------------------------------------------------------------------
def short_output_screen(predictions, cyto_effectors, apo_effectors, cyto_apo_effectors, apo_cyto_effectors):
    """ Function: short_output_screen()
        Purpose:  Given the WEKA predictions for each protein, write  
                  string that contains the short output format.
              
        Input:    WEKA predictions for each protein.                  
    
        Return:   String that contains predictions for all proteins as tab-delimited table.
    """
    # Output predictions for all proteins as tab-delimited table

    col_width = max(len(protein) for protein, pred, prob, sequence in predictions) + 1  # padding
    col_width = max(col_width, 10)
    pred_col_width = 20

    short_output_string = "".join("# Identifier".ljust(col_width)) + '\t' 
    short_output_string += "".join("Cytoplasmic effector".ljust(pred_col_width)) + '\t' 
    short_output_string += "".join("Apoplastic effector".ljust(pred_col_width)) + '\t' 
    short_output_string += "".join("Non-effector".ljust(pred_col_width)) + '\t' 
    short_output_string += "".join("Prediction".ljust(pred_col_width)) + '\n' 

    for index, (protein, pred, prob, sequence) in enumerate(predictions):    

        short_ident = 'protein' + str(index)
        if short_ident in cyto_effectors:
            short_output_string += "".join(protein.ljust(col_width)) + '\t' + 'Y' + ' (' + str(prob) + ')           ' + '\t' + "".join('-'.ljust(pred_col_width)) + '\t'
            short_output_string += "".join('-'.ljust(pred_col_width)) + '\t' + 'Cytoplasmic effector' + '\n'            
        elif short_ident in apo_effectors:
            short_output_string += "".join(protein.ljust(col_width)) + '\t' + "".join('-'.ljust(pred_col_width)) + '\t' + 'Y' + ' (' + str(prob) + ')           ' + '\t' 
            short_output_string += "".join('-'.ljust(pred_col_width)) + '\t' + 'Apoplastic effector' + '\n'           
        elif short_ident in cyto_apo_effectors:
            short_output_string += "".join(protein.ljust(col_width)) + '\t' + 'Y' + ' (' + str(cyto_apo_effectors[short_ident][0])+ ')           ' + '\t' 
            short_output_string += 'Y' + ' (' + str(cyto_apo_effectors[short_ident][1]) + ')           ' + '\t' + "".join('-'.ljust(pred_col_width)) + '\t' + 'Cytoplasmic/apoplastic effector' + '\n'  
        elif short_ident in apo_cyto_effectors:
            short_output_string += "".join(protein.ljust(col_width)) + '\t' + 'Y' + ' (' + str(apo_cyto_effectors[short_ident][1])+ ')           ' + '\t' 
            short_output_string += 'Y' + ' (' + str(apo_cyto_effectors[short_ident][0]) + ')           ' + '\t' + "".join('-'.ljust(pred_col_width)) + '\t' + 'Apoplastic/cytoplasmic effector' + '\n'  
        else:
            short_output_string += "".join(protein.ljust(col_width)) + '\t' + "".join('-'.ljust(pred_col_width)) + '\t' + "".join('-'.ljust(pred_col_width)) + '\t' 
            short_output_string += 'Y' + ' (' + str(prob) + ')           ' + '\t' + 'Non-effector' + '\n'  


    return short_output_string
# -----------------------------------------------------------------------------------------------------------
def short_output_file(predictions, cyto_effectors, apo_effectors, cyto_apo_effectors, apo_cyto_effectors):
    """ Function: short_output_file()
        Purpose:  Given the WEKA predictions for each protein, write  
                  string that contains the short output format.
              
        Input:    WEKA predictions for each protein.                  
    
        Return:   String that contains predictions for all proteins as tab-delimited table.
    """
    # Output predictions for all proteins as tab-delimited table


    short_output_string = "# Identifier" + '\t' + "Cytoplasmic effector" + '\t' 
    short_output_string += "Apoplastic effector" + '\t' 
    short_output_string += "Non-effector" + '\t' 
    short_output_string += "Prediction" + '\n' 

    for index, (protein, pred, prob, sequence) in enumerate(predictions):    

        short_ident = 'protein' + str(index)
        
        if short_ident in cyto_effectors:
            short_output_string += protein + '\t' + 'Y' + ' (' + str(prob) + ')' + '\t' + '-' + '\t'
            short_output_string += '-' + '\t' + 'Cytoplasmic effector' + '\n'            
        elif short_ident in apo_effectors:
            short_output_string += protein + '\t' + '-' + '\t' + 'Y' + ' (' + str(prob) + ')' + '\t' 
            short_output_string += '-' + '\t' + 'Apoplastic effector' + '\n'           
        elif short_ident in cyto_apo_effectors:
            short_output_string += protein + '\t' + 'Y' + ' (' + str(cyto_apo_effectors[short_ident][0]) + ')' + '\t' 
            short_output_string += 'Y' + ' (' + str(cyto_apo_effectors[short_ident][1]) + ')' + '\t' + '-' + '\t' + 'Cytoplasmic/apoplastic effector' + '\n'  
        elif short_ident in apo_cyto_effectors:
            short_output_string += protein + '\t' + 'Y' + ' (' + str(apo_cyto_effectors[short_ident][1]) + ')' + '\t' 
            short_output_string += 'Y' + ' (' + str(apo_cyto_effectors[short_ident][0]) + ')' + '\t' + '-' + '\t' + 'Apoplastic/cytoplasmic effector' + '\n'  
        else:
            short_output_string += protein + '\t' + '-' + '\t' + '-' + '\t' 
            short_output_string += 'Y' + ' (' + str(prob) + ')' + '\t' + 'Non-effector' + '\n'  


    return short_output_string    
# -----------------------------------------------------------------------------------------------------------
def long_output(ORIGINAL_IDENTIFIERS, predicted_effectors):
    """ Function: long_output()
        Purpose:  Given the predicted effectors and identifiers for the test set,  
                  write string that contains the long output format.
              
        Input:    Predicted effectors and identifiers of test set.                  
    
        Return:   String that contains list of predicted effectors with posterior probabilites
                  and a short statistic on the percentage of predicted effectors in the test set.
    """
    # Output predicted effectors for long format
    long_output_string = '-----------------\n'
    long_output_string += 'Predicted effectors:\n\n'
    for effector, prob, sequence in predicted_effectors:
        long_output_string += effector + '| Effector probability:' + str(prob) + '\n'

    long_output_string += '-----------------\n\n'
    long_output_string += 'Number of proteins that were tested: ' + str(len(ORIGINAL_IDENTIFIERS)) + '\n' 
    long_output_string += 'Number of predicted effectors: ' + str(len(predicted_effectors)) + '\n' 
    long_output_string += '\n' + '-----------------' + '\n' 
    long_output_string += str(round(100.0*len(predicted_effectors)/len(ORIGINAL_IDENTIFIERS), 1)) + ' percent are predicted to be effectors.'  
    long_output_string += '\n' + '-----------------' + '\n'

    return long_output_string
# -----------------------------------------------------------------------------------------------------------