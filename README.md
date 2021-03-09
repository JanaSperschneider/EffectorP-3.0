### EffectorP-3.0: Prediction of apoplastic and cytoplasmic effectors in fungi and oomycetes

#### What is EffectorP 3.0?

Many fungi and oomycetes species are devasting plant pathogens and secrete effector proteins to facilitate plant infection. Fungi and oomycete pathogens have diverse infection strategies and their effectors do not share sequence homology. However, effectors still have unifying properties: they either localize extracellularly to the plant apoplast or intracellularly to the plant cytoplasm. EffectorP 3.0 exploits this biological signal and uses two machine learning models trained on apoplastic and cytoplasmic effectors, respectively.

For a given set of secreted pathogen proteins, EffectorP 3.0 will predict if a protein is:
* an apoplastic effector
* a cytoplasmic effector
* a dual-localized apoplastic/cytoplasmic effector
* a non-effector

#### What is EffectorP 3.0 not?

EffectorP is **not** a tool for secretome prediction and **not** a tool for bacterial effector prediction. 

EffectorP has been trained to find fungal/oomycete effectors in secretomes, so please run it on a FASTA file of secreted proteins 
to test if they are predicted effectors. It is **essential** 
to use tools such as SignalP, Phobius and TMHMM	to predict first if a protein is likely to be secreted.
Alternatively, high-confidence experimentally determined secretomes instead of computationally predicted secretomes can be submitted to EffectorP. 

#### Installing EffectorP 3.0

EffectorP has been written in Python3 and uses the WEKA 3.8.4 software. To get EffectorP to work on your local machine, you need to unzip the WEKA software, which is already provided in the EffectorP distribution to ensure that a compatible version is used. You also need an installation of Python 3.

0. Download the latest release from this github repo (or alternatively you can clone the github repo and skip step 1).

1. Make sure EffectorP has the permission to execute. Then unpack EffectorP in your desired location
```
tar xvf EffectorP-3.0.tar.gz
chmod -R 755 EffectorP-3.0/
cd EffectorP-3.0
```

2. For WEKA, you need to simply unzip the file weka-3-8-4.zip
```
unzip weka-3-8-4.zip
```
If you are having troube installing WEKA, please see [here](https://www.cs.waikato.ac.nz/~ml/weka/index.html) for help. 

4. Test if EffectorP is working
```
python EffectorP.py -i Effectors.fasta
```

#### EffectorP output format
Run this to get a feel for the output format:
```
python EffectorP.py -i Effectors.fasta
-----------------

EffectorP 3.0 is running for 4 proteins given in FASTA file Effectors.fasta

Ensemble classification
25 percent done...
50 percent done...
75 percent done...
All done.

# Identifier            Cytoplasmic effector    Apoplastic effector     Non-effector            Prediction
AvrSr27 cytoplasmic     Y (0.722)               Y (0.772)               -                       Apoplastic/cytoplasmic effector
AvrLm6                  -                       Y (0.73)                -                       Apoplastic effector
HaCR1 apoplast M4BIN0   Y (0.543)               Y (0.864)               -                       Apoplastic/cytoplasmic effector
PvRXLR53 cytoplasmic    Y (0.788)               -                       -                       Cytoplasmic effector

-----------------
4 proteins were provided as input in the following file: Effectors.fasta
-----------------
Number of predicted effectors: 4
Number of predicted cytoplasmic effectors: 1
Number of predicted apoplastic effectors: 3
-----------------
100.0 percent are predicted effectors.
25.0 percent are predicted cytoplasmic effectors.
75.0 percent are predicted apoplastic effectors.
-----------------

```

EffectorP will return the output as shown in the example above. A summary table will be shown which shows the predictions (effector or non-effector) for each submitted protein.

For a predicted effector, its most likely localization (apoplastic or cytoplasmic) will be returned, with an associated probability.

As a probabilistic classifier (Naive Bayes), EffectorP returns a probability that a tested instance will belong to either the effector or non-effector class and these probabilities are included in the web server output as additional information to researchers. However, in Naive Bayes classification these are known to be only rough estimations and should therefore not be overinterpreted.

We deliberately did not recommend a probability threshold over which a protein would be classified as an effector candidate, as we believe it should remain up to the individual user to interpret their results in the context of additional resources available. For example, a researcher might like to predict the full effector candidate complement using EffectorP and overlay this with in planta expression data to prioritize candidates, whereas in other situations without additional information a list of high-priority candidates as determined by the EffectorP probabilities might be more appropriate. 

#### Citation for EffectorP 

Preprint in preparation...
