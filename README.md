### EffectorP-3.0: Prediction of apoplastic and cytoplasmic effectors in fungi and oomycetes

#### What is EffectorP 3.0?

Many fungi and oomycetes species are devasting plant pathogens and secrete effector proteins to facilitate plant infection. Fungi and oomycete pathogens have diverse infection strategies and their effectors do not share sequence homology. However, effectors still have unifying properties: they either localize extracellularly to the plant apoplast or intracellularly to the plant cytoplasm. EffectorP 3.0 exploits this biological signal and uses two machine learning models trained on apoplastic and cytoplasmic effectors, respectively.

For a given set of secreted pathogen proteins, EffectorP 3.0 will predict if a protein is:
* an apoplastic effector
* a cytoplasmic effector
* a dual-localized apoplastic/cytoplasmic effector
* a non-effector

#### What is EffectorP 3.0 not?

EffectorP is *not* a tool for secretome prediction and *not* a tool for bacterial effector prediction. 

EffectorP has been trained to find fungal/oomycete effectors in secretomes, so please run it on a FASTA file of secreted proteins 
to test if they are predicted effectors. It is recommended 
to use tools such as SignalP, Phobius and TMHMM	to predict first if a protein is likely to be secreted.
Alternatively, high-confidence experimentally determined secretomes instead of computationally predicted secretomes can be submitted to EffectorP. 
