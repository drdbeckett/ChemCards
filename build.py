"""Assemble + validate compounds.json. Each drug carries an expected molecular
formula; we parse the SMILES, compute the formula with RDKit, and refuse to
emit anything that doesn't match (or doesn't parse)."""
import json, sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# (name, smiles, [categories], abbrev, expected_formula or None)
E = None  # no formula check (non-drug)
DATA = [
 # ---- Amino acids (Gly/Glu/Asp also neurotransmitters) ----
 ("Glycine","NCC(=O)O",["Amino acids","Neurotransmitters"],"Gly / G",E),
 ("Alanine","C[C@@H](N)C(=O)O",["Amino acids"],"Ala / A",E),
 ("Valine","CC(C)[C@@H](N)C(=O)O",["Amino acids"],"Val / V",E),
 ("Leucine","CC(C)C[C@@H](N)C(=O)O",["Amino acids"],"Leu / L",E),
 ("Isoleucine","CC[C@H](C)[C@@H](N)C(=O)O",["Amino acids"],"Ile / I",E),
 ("Proline","O=C(O)[C@@H]1CCCN1",["Amino acids"],"Pro / P",E),
 ("Phenylalanine","N[C@@H](Cc1ccccc1)C(=O)O",["Amino acids"],"Phe / F",E),
 ("Tryptophan","N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O",["Amino acids"],"Trp / W",E),
 ("Methionine","N[C@@H](CCSC)C(=O)O",["Amino acids"],"Met / M",E),
 ("Serine","N[C@@H](CO)C(=O)O",["Amino acids"],"Ser / S",E),
 ("Threonine","C[C@@H](O)[C@@H](N)C(=O)O",["Amino acids"],"Thr / T",E),
 ("Cysteine","N[C@@H](CS)C(=O)O",["Amino acids"],"Cys / C",E),
 ("Tyrosine","N[C@@H](Cc1ccc(O)cc1)C(=O)O",["Amino acids"],"Tyr / Y",E),
 ("Asparagine","N[C@@H](CC(N)=O)C(=O)O",["Amino acids"],"Asn / N",E),
 ("Glutamine","N[C@@H](CCC(N)=O)C(=O)O",["Amino acids"],"Gln / Q",E),
 ("Aspartic acid","N[C@@H](CC(=O)O)C(=O)O",["Amino acids","Neurotransmitters"],"Aspartate / Asp / D",E),
 ("Glutamic acid","N[C@@H](CCC(=O)O)C(=O)O",["Amino acids","Neurotransmitters"],"Glutamate / Glu / E",E),
 ("Lysine","NCCCC[C@@H](N)C(=O)O",["Amino acids"],"Lys / K",E),
 ("Arginine","N[C@@H](CCCNC(=N)N)C(=O)O",["Amino acids"],"Arg / R",E),
 ("Histidine","N[C@@H](Cc1c[nH]cn1)C(=O)O",["Amino acids"],"His / H",E),
 # ---- Neurotransmitters (non-amino-acid) ----
 ("Dopamine","NCCc1ccc(O)c(O)c1",["Neurotransmitters"],"DA",E),
 ("Serotonin","NCCc1c[nH]c2ccc(O)cc12",["Neurotransmitters"],"5-HT",E),
 ("Norepinephrine","NC[C@H](O)c1ccc(O)c(O)c1",["Neurotransmitters"],"NE / noradrenaline",E),
 ("Epinephrine","CNC[C@H](O)c1ccc(O)c(O)c1",["Neurotransmitters"],"adrenaline",E),
 ("GABA","NCCCC(=O)O",["Neurotransmitters"],"\u03b3-aminobutyric acid",E),
 ("Acetylcholine","CC(=O)OCC[N+](C)(C)C",["Neurotransmitters"],"ACh",E),
 ("Histamine","NCCc1c[nH]cn1",["Neurotransmitters"],"",E),
 ("Melatonin","COc1ccc2[nH]cc(CCNC(C)=O)c2c1",["Neurotransmitters"],"",E),
 # ---- Nucleic acid bases ----
 ("Adenine","Nc1ncnc2[nH]cnc12",["Nucleic acid bases"],"A",E),
 ("Guanine","Nc1nc2[nH]cnc2c(=O)[nH]1",["Nucleic acid bases"],"G",E),
 ("Cytosine","Nc1cc[nH]c(=O)n1",["Nucleic acid bases"],"C",E),
 ("Thymine","Cc1c[nH]c(=O)[nH]c1=O",["Nucleic acid bases"],"T",E),
 ("Uracil","O=c1cc[nH]c(=O)[nH]1",["Nucleic acid bases"],"U",E),

 # ============ DRUGS ============
 # ---- Antipsychotics ----
 ("Haloperidol","O=C(CCCN1CCC(O)(c2ccc(Cl)cc2)CC1)c1ccc(F)cc1",["Drugs","Antipsychotics"],"Haldol","C21H23ClFNO2"),
 ("Chlorpromazine","CN(C)CCCN1c2ccccc2Sc2ccc(Cl)cc21",["Drugs","Antipsychotics"],"Thorazine / Largactil","C17H19ClN2S"),
 ("Fluphenazine","OCCN1CCN(CCCN2c3ccccc3Sc3ccc(C(F)(F)F)cc32)CC1",["Drugs","Antipsychotics"],"Prolixin / Modecate","C22H26F3N3OS"),
 ("Thioridazine","CN1CCCCC1CCN1c2ccccc2Sc2ccc(SC)cc21",["Drugs","Antipsychotics"],"Mellaril","C21H26N2S2"),
 ("Perphenazine","OCCN1CCN(CCCN2c3ccccc3Sc3ccc(Cl)cc32)CC1",["Drugs","Antipsychotics"],"Trilafon","C21H26ClN3OS"),
 ("Loxapine","CN1CCN(C2=Nc3cc(Cl)ccc3Oc3ccccc32)CC1",["Drugs","Antipsychotics"],"Loxitane / Adasuve","C18H18ClN3O"),
 ("Clozapine","CN1CCN(C2=Nc3cc(Cl)ccc3Nc3ccccc32)CC1",["Drugs","Antipsychotics"],"Clozaril","C18H19ClN4"),
 ("Olanzapine","Cc1cc2c(s1)Nc1ccccc1N=C2N1CCN(C)CC1",["Drugs","Antipsychotics"],"Zyprexa","C17H20N4S"),
 ("Quetiapine","OCCOCCN1CCN(C2=Nc3ccccc3Sc3ccccc32)CC1",["Drugs","Antipsychotics"],"Seroquel","C21H25N3O2S"),
 ("Aripiprazole","O=C1CCc2cc(OCCCCN3CCN(c4cccc(Cl)c4Cl)CC3)ccc2N1",["Drugs","Antipsychotics"],"Abilify","C23H27Cl2N3O2"),
 ("Ziprasidone","O=C1Cc2cc(CCN3CCN(c4nsc5ccccc45)CC3)c(Cl)cc2N1",["Drugs","Antipsychotics"],"Geodon","C21H21ClN4OS"),
 ("Risperidone","Cc1nc2CCCCn2c(=O)c1CCN1CCC(c2noc3cc(F)ccc23)CC1",["Drugs","Antipsychotics"],"Risperdal","C23H27FN4O2"),
 ("Pimozide","O=C1Nc2ccccc2N1C1CCN(CCCC(c2ccc(F)cc2)c2ccc(F)cc2)CC1",["Drugs","Antipsychotics"],"Orap","C28H29F2N3O"),
 ("Sulpiride","CCN1CCCC1CNC(=O)c1cc(S(N)(=O)=O)ccc1OC",["Drugs","Antipsychotics"],"Dogmatil","C15H23N3O4S"),
 ("Amisulpride","CCN1CCCC1CNC(=O)c1cc(S(=O)(=O)CC)c(N)cc1OC",["Drugs","Antipsychotics"],"Solian","C17H27N3O4S"),
 ("Cariprazine","CN(C)C(=O)NC1CCC(CCN2CCN(c3cccc(Cl)c3Cl)CC2)CC1",["Drugs","Antipsychotics"],"Vraylar / Reagila","C21H32Cl2N4O"),
 # ---- SSRIs / SNRIs ----
 ("Fluoxetine","CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1",["Drugs","SSRIs/SNRIs"],"Prozac / Sarafem","C17H18F3NO"),
 ("Sertraline","CNC1CCC(c2ccc(Cl)c(Cl)c2)c2ccccc21",["Drugs","SSRIs/SNRIs"],"Zoloft","C17H17Cl2N"),
 ("Paroxetine","Fc1ccc(C2CCNCC2COc2ccc3c(c2)OCO3)cc1",["Drugs","SSRIs/SNRIs"],"Paxil / Seroxat","C19H20FNO3"),
 ("Citalopram","N#Cc1ccc2c(c1)C(c1ccc(F)cc1)(CCCN(C)C)OC2",["Drugs","SSRIs/SNRIs"],"Celexa / Cipramil","C20H21FN2O"),
 ("Escitalopram","N#Cc1ccc2c(c1)[C@](c1ccc(F)cc1)(CCCN(C)C)OC2",["Drugs","SSRIs/SNRIs"],"Lexapro / Cipralex","C20H21FN2O"),
 ("Fluvoxamine","COCCCC/C(=N\\OCCN)c1ccc(C(F)(F)F)cc1",["Drugs","SSRIs/SNRIs"],"Luvox / Faverin","C15H21F3N2O2"),
 ("Vilazodone","NC(=O)c1cc2cc(N3CCN(CCCCc4c[nH]c5ccc(C#N)cc45)CC3)ccc2o1",["Drugs","SSRIs/SNRIs"],"Viibryd","C26H27N5O2"),
 ("Vortioxetine","Cc1ccc(Sc2ccccc2N2CCNCC2)c(C)c1",["Drugs","SSRIs/SNRIs"],"Trintellix / Brintellix","C18H22N2S"),
 ("Venlafaxine","COc1ccc(C(CN(C)C)C2(O)CCCCC2)cc1",["Drugs","SSRIs/SNRIs"],"Effexor","C17H27NO2"),
 ("Desvenlafaxine","Oc1ccc(C(CN(C)C)C2(O)CCCCC2)cc1",["Drugs","SSRIs/SNRIs"],"Pristiq","C16H25NO2"),
 ("Duloxetine","CNCCC(Oc1cccc2ccccc12)c1cccs1",["Drugs","SSRIs/SNRIs"],"Cymbalta","C18H19NOS"),
 ("Milnacipran","CCN(CC)C(=O)C1(c2ccccc2)CC1CN",["Drugs","SSRIs/SNRIs"],"Savella / Ixel","C15H22N2O"),
 ("Levomilnacipran","CCN(CC)C(=O)[C@]1(c2ccccc2)C[C@@H]1CN",["Drugs","SSRIs/SNRIs"],"Fetzima","C15H22N2O"),
 # ---- Sedatives & anxiolytics ----
 ("Diazepam","CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",["Drugs","Sedatives & anxiolytics"],"Valium / Diastat","C16H13ClN2O"),
 ("Lorazepam","OC1N=C(c2ccccc2Cl)c2cc(Cl)ccc2NC1=O",["Drugs","Sedatives & anxiolytics"],"Ativan","C15H10Cl2N2O2"),
 ("Alprazolam","Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1)=NC2",["Drugs","Sedatives & anxiolytics"],"Xanax","C17H13ClN4"),
 ("Clonazepam","O=C1CN=C(c2ccccc2Cl)c2cc([N+](=O)[O-])ccc2N1",["Drugs","Sedatives & anxiolytics"],"Klonopin / Rivotril","C15H10ClN3O3"),
 ("Temazepam","CN1C(=O)C(O)N=C(c2ccccc2)c2cc(Cl)ccc21",["Drugs","Sedatives & anxiolytics"],"Restoril","C16H13ClN2O2"),
 ("Midazolam","Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",["Drugs","Sedatives & anxiolytics"],"Versed / Dormicum","C18H13ClFN3"),
 ("Oxazepam","OC1N=C(c2ccccc2)c2cc(Cl)ccc2NC1=O",["Drugs","Sedatives & anxiolytics"],"Serax","C15H11ClN2O2"),
 ("Chlordiazepoxide","CNC1=Nc2ccc(Cl)cc2C(c2ccccc2)=[N+]([O-])C1",["Drugs","Sedatives & anxiolytics"],"Librium","C16H14ClN3O"),
 ("Triazolam","Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1Cl)=NC2",["Drugs","Sedatives & anxiolytics"],"Halcion","C17H12Cl2N4"),
 ("Flurazepam","CCN(CC)CCN1C(=O)CN=C(c2ccccc2F)c2cc(Cl)ccc21",["Drugs","Sedatives & anxiolytics"],"Dalmane","C21H23ClFN3O"),
 ("Nitrazepam","O=C1CN=C(c2ccccc2)c2cc([N+](=O)[O-])ccc2N1",["Drugs","Sedatives & anxiolytics"],"Mogadon","C15H11N3O3"),
 ("Clorazepate","O=C1Nc2ccc(Cl)cc2C(c2ccccc2)=NC1C(=O)O",["Drugs","Sedatives & anxiolytics"],"Tranxene","C16H11ClN2O3"),
 ("Zolpidem","Cc1ccc(-c2nc3ccc(C)cn3c2CC(=O)N(C)C)cc1",["Drugs","Sedatives & anxiolytics"],"Ambien / Stilnox","C19H21N3O"),
 ("Zaleplon","CCN(C(C)=O)c1cccc(-c2ccnc3c(C#N)cnn23)c1",["Drugs","Sedatives & anxiolytics"],"Sonata","C17H15N5O"),
 ("Eszopiclone","CN1CCN(C(=O)OC2c3nccnc3C(=O)N2c2ccc(Cl)cn2)CC1",["Drugs","Sedatives & anxiolytics"],"Lunesta","C17H17ClN6O3"),
 ("Phenobarbital","CCC1(c2ccccc2)C(=O)NC(=O)NC1=O",["Drugs","Sedatives & anxiolytics"],"Luminal","C12H12N2O3"),
 ("Pentobarbital","CCCC(C)C1(CC)C(=O)NC(=O)NC1=O",["Drugs","Sedatives & anxiolytics"],"Nembutal","C11H18N2O3"),
 ("Secobarbital","CCCC(C)C1(CC=C)C(=O)NC(=O)NC1=O",["Drugs","Sedatives & anxiolytics"],"Seconal","C12H18N2O3"),
 ("Amobarbital","CCC1(CCC(C)C)C(=O)NC(=O)NC1=O",["Drugs","Sedatives & anxiolytics"],"Amytal","C11H18N2O3"),
 ("Buspirone","O=C1CC2(CCCC2)CC(=O)N1CCCCN1CCN(c2ncccn2)CC1",["Drugs","Sedatives & anxiolytics"],"Buspar","C21H31N5O2"),
 # ---- Opioids (tramadol also SNRI) ----
 ("Morphine","CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5",["Drugs","Opioids"],"MS Contin / Kadian","C17H19NO3"),
 ("Codeine","CN1CC[C@]23c4c5ccc(OC)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5",["Drugs","Opioids"],"","C18H21NO3"),
 ("Oxycodone","CN1CC[C@]23c4c5ccc(OC)c4O[C@H]2C(=O)CC[C@@]3(O)[C@H]1C5",["Drugs","Opioids"],"OxyContin / Roxicodone","C18H21NO4"),
 ("Hydrocodone","CN1CC[C@]23c4c5ccc(OC)c4O[C@H]2C(=O)CC[C@H]3[C@H]1C5",["Drugs","Opioids"],"Hysingla / Norco","C18H21NO3"),
 ("Fentanyl","CCC(=O)N(c1ccccc1)C1CCN(CCc2ccccc2)CC1",["Drugs","Opioids"],"Duragesic / Sublimaze","C22H28N2O"),
 ("Methadone","CCC(=O)C(c1ccccc1)(c1ccccc1)CC(C)N(C)C",["Drugs","Opioids"],"Dolophine / Methadose","C21H27NO"),
 ("Tramadol","CN(C)CC1CCCCC1(O)c1cccc(OC)c1",["Drugs","Opioids","SSRIs/SNRIs"],"Ultram","C16H25NO2"),
 ("Hydromorphone","CN1CC[C@]23c4c5ccc(O)c4O[C@H]2C(=O)CC[C@H]3[C@H]1C5",["Drugs","Opioids"],"Dilaudid","C17H19NO3"),
 ("Naloxone","C=CCN1CC[C@]23c4c5ccc(O)c4O[C@H]2C(=O)CC[C@@]3(O)[C@H]1C5",["Drugs","Opioids"],"Narcan","C19H21NO4"),
 ("Naltrexone","O=C1CC[C@@]2(O)[C@@H]3Cc4ccc(O)c5c4[C@@]2(CCN3CC2CC2)[C@H]1O5",["Drugs","Opioids"],"Vivitrol / Revia","C20H23NO4"),
 ("Diamorphine","CN1CC[C@]23c4c5ccc(OC(C)=O)c4O[C@H]2[C@@H](OC(C)=O)C=C[C@H]3[C@H]1C5",["Drugs","Opioids"],"heroin / diacetylmorphine","C21H23NO5"),
 ("Meperidine","CCOC(=O)C1(c2ccccc2)CCN(C)CC1",["Drugs","Opioids"],"Demerol / pethidine","C15H21NO2"),
 ("Tapentadol","CC[C@@H](c1cccc(O)c1)[C@@H](C)CN(C)C",["Drugs","Opioids"],"Nucynta","C14H23NO"),
 ("Oxymorphone","CN1CC[C@]23c4c5ccc(O)c4O[C@H]2C(=O)CC[C@@]3(O)[C@H]1C5",["Drugs","Opioids"],"Opana","C17H19NO4"),
 ("Loperamide","CN(C)C(=O)C(c1ccccc1)(c1ccccc1)CCN1CCC(O)(c2ccc(Cl)cc2)CC1",["Drugs","Opioids"],"Imodium","C29H33ClN2O2"),
 # ---- Stimulants ----
 ("Caffeine","Cn1cnc2c1c(=O)n(C)c(=O)n2C",["Drugs","Stimulants"],"NoDoz / Vivarin","C8H10N4O2"),
 ("Nicotine","CN1CCC[C@H]1c1cccnc1",["Drugs","Stimulants"],"Nicorette / NicoDerm","C10H14N2"),
 ("Amphetamine","CC(N)Cc1ccccc1",["Drugs","Stimulants"],"Adderall (with dextroamphetamine)","C9H13N"),
 ("Methamphetamine","CNC(C)Cc1ccccc1",["Drugs","Stimulants"],"Desoxyn","C10H15N"),
 ("Methylphenidate","COC(=O)C(c1ccccc1)C1CCCCN1",["Drugs","Stimulants"],"Ritalin / Concerta","C14H19NO2"),
 ("Cocaine","COC(=O)C1C2CCC(CC1OC(=O)c1ccccc1)N2C",["Drugs","Stimulants"],"","C17H21NO4"),
 ("MDMA","CNC(C)Cc1ccc2c(c1)OCO2",["Drugs","Stimulants"],"ecstasy / molly","C11H15NO2"),
 ("Modafinil","NC(=O)CS(=O)C(c1ccccc1)c1ccccc1",["Drugs","Stimulants"],"Provigil","C15H15NO2S"),
 ("Atomoxetine","CNCCC(Oc1ccccc1C)c1ccccc1",["Drugs","Stimulants"],"Strattera","C17H21NO"),
 ("Phentermine","CC(C)(N)Cc1ccccc1",["Drugs","Stimulants"],"Adipex-P","C10H15N"),
 # ---- NSAIDs & analgesics ----
 ("Aspirin","CC(=O)Oc1ccccc1C(=O)O",["Drugs","NSAIDs & analgesics"],"acetylsalicylic acid / ASA / Bayer","C9H8O4"),
 ("Ibuprofen","CC(C)Cc1ccc(C(C)C(=O)O)cc1",["Drugs","NSAIDs & analgesics"],"Advil / Motrin / Nurofen","C13H18O2"),
 ("Acetaminophen","CC(=O)Nc1ccc(O)cc1",["Drugs","NSAIDs & analgesics"],"Paracetamol / Tylenol / Panadol / APAP","C8H9NO2"),
 ("Naproxen","COc1ccc2cc(C(C)C(=O)O)ccc2c1",["Drugs","NSAIDs & analgesics"],"Aleve / Naprosyn","C14H14O3"),
 ("Ketoprofen","CC(C(=O)O)c1cccc(C(=O)c2ccccc2)c1",["Drugs","NSAIDs & analgesics"],"Orudis","C16H14O3"),
 ("Diclofenac","O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl",["Drugs","NSAIDs & analgesics"],"Voltaren","C14H11Cl2NO2"),
 ("Indomethacin","COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1",["Drugs","NSAIDs & analgesics"],"Indocin","C19H16ClNO4"),
 ("Celecoxib","Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1",["Drugs","NSAIDs & analgesics"],"Celebrex","C17H14F3N3O2S"),
 ("Meloxicam","Cc1cnc(NC(=O)C2=C(O)c3ccccc3S(=O)(=O)N2C)s1",["Drugs","NSAIDs & analgesics"],"Mobic","C14H13N3O4S2"),
 ("Piroxicam","CN1S(=O)(=O)c2ccccc2C(O)=C1C(=O)Nc1ccccn1",["Drugs","NSAIDs & analgesics"],"Feldene","C15H13N3O4S"),
 ("Mefenamic acid","Cc1cccc(Nc2ccccc2C(=O)O)c1C",["Drugs","NSAIDs & analgesics"],"Ponstel","C15H15NO2"),
 ("Etoricoxib","Cc1ccc(-c2cc(Cl)cnc2-c2ccc(S(C)(=O)=O)cc2)cn1",["Drugs","NSAIDs & analgesics"],"Arcoxia","C18H15ClN2O2S"),
 # ---- Statins ----
 ("Atorvastatin","CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O",["Drugs","Statins"],"Lipitor / Torvast","C33H35FN2O5"),
 ("Simvastatin","CCC(C)(C)C(=O)O[C@@H]1C[C@@H](C)C=C2C=C[C@H](C)[C@H](CC[C@@H]3C[C@H](O)CC(=O)O3)[C@H]12",["Drugs","Statins"],"Zocor","C25H38O5"),
 ("Rosuvastatin","CC(C)c1nc(N(C)S(C)(=O)=O)nc(-c2ccc(F)cc2)c1/C=C/[C@@H](O)C[C@@H](O)CC(=O)O",["Drugs","Statins"],"Crestor","C22H28FN3O6S"),
 ("Pravastatin","CC[C@H](C)C(=O)O[C@H]1C[C@H](O)C=C2C=C[C@H](C)[C@H](CC[C@@H](O)C[C@@H](O)CC(=O)O)[C@@H]12",["Drugs","Statins"],"Pravachol","C23H36O7"),
 ("Lovastatin","CC[C@H](C)C(=O)O[C@H]1C[C@H](C)C=C2C=C[C@H](C)[C@H](CC[C@@H]3C[C@@H](O)CC(=O)O3)[C@@H]21",["Drugs","Statins"],"Mevacor","C24H36O5"),
 ("Fluvastatin","CC(C)n1c(/C=C/[C@@H](O)C[C@@H](O)CC(=O)O)c(-c2ccc(F)cc2)c2ccccc21",["Drugs","Statins"],"Lescol","C24H26FNO4"),
 ("Pitavastatin","OC(=O)C[C@@H](O)C[C@@H](O)/C=C/c1c(-c2ccc(F)cc2)nc2ccccc2c1C1CC1",["Drugs","Statins"],"Livalo / Zypitamag","C25H24FNO4"),
 # ---- Antihistamines & allergy ----
 ("Diphenhydramine","CN(C)CCOC(c1ccccc1)c1ccccc1",["Drugs","Antihistamines & allergy"],"Benadryl","C17H21NO"),
 ("Loratadine","CCOC(=O)N1CCC(=C2c3ccc(Cl)cc3CCc3cccnc32)CC1",["Drugs","Antihistamines & allergy"],"Claritin","C22H23ClN2O2"),
 ("Desloratadine","Clc1ccc2c(c1)CCc1cccnc1C2=C1CCNCC1",["Drugs","Antihistamines & allergy"],"Clarinex / Neoclarityn","C19H19ClN2"),
 ("Cetirizine","O=C(O)COCCN1CCN(C(c2ccccc2)c2ccc(Cl)cc2)CC1",["Drugs","Antihistamines & allergy"],"Zyrtec / Reactine","C21H25ClN2O3"),
 ("Levocetirizine","O=C(O)COCCN1CCN([C@@H](c2ccccc2)c2ccc(Cl)cc2)CC1",["Drugs","Antihistamines & allergy"],"Xyzal","C21H25ClN2O3"),
 ("Fexofenadine","CC(C)(C(=O)O)c1ccc(C(O)CCCN2CCC(C(O)(c3ccccc3)c3ccccc3)CC2)cc1",["Drugs","Antihistamines & allergy"],"Allegra / Telfast","C32H39NO4"),
 ("Chlorpheniramine","CN(C)CCC(c1ccc(Cl)cc1)c1ccccn1",["Drugs","Antihistamines & allergy"],"Chlor-Trimeton / Piriton","C16H19ClN2"),
 ("Brompheniramine","CN(C)CCC(c1ccc(Br)cc1)c1ccccn1",["Drugs","Antihistamines & allergy"],"Dimetapp","C16H19BrN2"),
 ("Promethazine","CC(CN1c2ccccc2Sc2ccccc21)N(C)C",["Drugs","Antihistamines & allergy"],"Phenergan","C17H20N2S"),
 ("Hydroxyzine","OCCOCCN1CCN(C(c2ccccc2)c2ccc(Cl)cc2)CC1",["Drugs","Antihistamines & allergy"],"Atarax / Vistaril","C21H27ClN2O2"),
 ("Doxylamine","CN(C)CCOC(C)(c1ccccc1)c1ccccn1",["Drugs","Antihistamines & allergy"],"Unisom","C17H22N2O"),
 ("Cyclizine","CN1CCN(C(c2ccccc2)c2ccccc2)CC1",["Drugs","Antihistamines & allergy"],"Marezine","C18H22N2"),
 ("Meclizine","Cc1cccc(CN2CCN(C(c3ccccc3)c3ccc(Cl)cc3)CC2)c1",["Drugs","Antihistamines & allergy"],"Antivert / Bonine","C25H27ClN2"),
 ("Cimetidine","CN/C(=N\\C#N)NCCSCc1nc[nH]c1C",["Drugs","Antihistamines & allergy"],"Tagamet","C10H16N6S"),
 ("Ranitidine","CNC(=C[N+](=O)[O-])NCCSCc1ccc(CN(C)C)o1",["Drugs","Antihistamines & allergy"],"Zantac","C13H22N4O3S"),
 ("Famotidine","N=C(N)Nc1nc(CSCCC(=N)NS(N)(=O)=O)cs1",["Drugs","Antihistamines & allergy"],"Pepcid","C8H15N7O2S3"),
 # ---- Antibiotics ----
 ("Penicillin G","CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O",["Drugs","Antibiotics"],"benzylpenicillin / Pfizerpen","C16H18N2O4S"),
 ("Amoxicillin","CC1(C)S[C@@H]2[C@H](NC(=O)[C@@H](N)c3ccc(O)cc3)C(=O)N2[C@H]1C(=O)O",["Drugs","Antibiotics"],"Amoxil / Moxatag","C16H19N3O5S"),
 ("Ampicillin","CC1(C)S[C@@H]2[C@H](NC(=O)[C@@H](N)c3ccccc3)C(=O)N2[C@H]1C(=O)O",["Drugs","Antibiotics"],"Principen","C16H19N3O4S"),
 ("Cephalexin","CC1=C(C(=O)O)N2C(=O)[C@@H](NC(=O)[C@@H](N)c3ccccc3)[C@H]2SC1",["Drugs","Antibiotics"],"Keflex","C16H17N3O4S"),
 ("Ciprofloxacin","OC(=O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",["Drugs","Antibiotics"],"Cipro","C17H18FN3O3"),
 ("Levofloxacin","C[C@H]1COc2c(N3CCN(C)CC3)c(F)cc3c(=O)c(C(=O)O)cn1c23",["Drugs","Antibiotics"],"Levaquin","C18H20FN3O4"),
 ("Metronidazole","Cc1ncc([N+](=O)[O-])n1CCO",["Drugs","Antibiotics"],"Flagyl","C6H9N3O3"),
 ("Trimethoprim","COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC",["Drugs","Antibiotics"],"Proloprim","C14H18N4O3"),
 ("Sulfamethoxazole","Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1",["Drugs","Antibiotics"],"Gantanol","C10H11N3O3S"),
 ("Sulfanilamide","NS(=O)(=O)c1ccc(N)cc1",["Drugs","Antibiotics"],"","C6H8N2O2S"),
 ("Dapsone","Nc1ccc(S(=O)(=O)c2ccc(N)cc2)cc1",["Drugs","Antibiotics"],"Aczone","C12H12N2O2S"),
 ("Chloramphenicol","O=C(NC(CO)C(O)c1ccc([N+](=O)[O-])cc1)C(Cl)Cl",["Drugs","Antibiotics"],"Chloromycetin","C11H12Cl2N2O5"),
 ("Linezolid","CC(=O)NC[C@H]1CN(c2ccc(N3CCOCC3)c(F)c2)C(=O)O1",["Drugs","Antibiotics"],"Zyvox","C16H20FN3O4"),
 ("Nitrofurantoin","O=C1CN(/N=C/c2ccc([N+](=O)[O-])o2)C(=O)N1",["Drugs","Antibiotics"],"Macrobid / Macrodantin","C8H6N4O5"),
 ("Isoniazid","NNC(=O)c1ccncc1",["Drugs","Antibiotics"],"INH / Nydrazid","C6H7N3O"),
 ("Pyrazinamide","NC(=O)c1cnccn1",["Drugs","Antibiotics"],"PZA","C5H5N3O"),
 # ---- Antifungals ----
 ("Fluconazole","OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F",["Drugs","Antifungals"],"Diflucan","C13H12F2N6O"),
 ("Ketoconazole","CC(=O)N1CCN(c2ccc(OC[C@H]3CO[C@@](Cn4ccnc4)(c4ccc(Cl)cc4Cl)O3)cc2)CC1",["Drugs","Antifungals"],"Nizoral","C26H28Cl2N4O4"),
 ("Clotrimazole","Clc1ccccc1C(c1ccccc1)(c1ccccc1)n1ccnc1",["Drugs","Antifungals"],"Lotrimin / Canesten","C22H17ClN2"),
 ("Miconazole","Clc1ccc(C(Cn2ccnc2)OCc2ccc(Cl)cc2Cl)c(Cl)c1",["Drugs","Antifungals"],"Monistat","C18H14Cl4N2O"),
 ("Econazole","Clc1ccc(C(Cn2ccnc2)OCc2ccc(Cl)cc2)c(Cl)c1",["Drugs","Antifungals"],"Spectazole","C18H15Cl3N2O"),
 ("Voriconazole","C[C@@H](c1ncncc1F)[C@](O)(Cn1cncn1)c1ccc(F)cc1F",["Drugs","Antifungals"],"Vfend","C16H14F3N5O"),
 ("Terbinafine","CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12",["Drugs","Antifungals"],"Lamisil","C21H25N"),
 ("Naftifine","CN(C/C=C/c1ccccc1)Cc1cccc2ccccc12",["Drugs","Antifungals"],"Naftin","C21H21N"),
 ("Butenafine","CN(Cc1ccc(C(C)(C)C)cc1)Cc1cccc2ccccc12",["Drugs","Antifungals"],"Lotrimin Ultra / Mentax","C23H27N"),
 ("Griseofulvin","COC1=CC(=O)C[C@@H](C)[C@]11Oc2c(Cl)c(OC)cc(OC)c2C1=O",["Drugs","Antifungals"],"Grifulvin / Gris-PEG","C17H17ClO6"),
 ("Flucytosine","Nc1nc(=O)[nH]cc1F",["Drugs","Antifungals"],"Ancobon / 5-FC","C4H4FN3O"),
 # ---- Insecticides (not drugs) ----
 ("DDT","ClC(Cl)(Cl)C(c1ccc(Cl)cc1)c1ccc(Cl)cc1",["Insecticides"],"dichlorodiphenyltrichloroethane","C14H9Cl5"),
 ("Malathion","CCOC(=O)CC(SP(=S)(OC)OC)C(=O)OCC",["Insecticides"],"","C10H19O6PS2"),
 ("Parathion","CCOP(=S)(OCC)Oc1ccc([N+](=O)[O-])cc1",["Insecticides"],"","C10H14NO5PS"),
 ("Chlorpyrifos","CCOP(=S)(OCC)Oc1nc(Cl)c(Cl)cc1Cl",["Insecticides"],"Lorsban / Dursban","C9H11Cl3NO3PS"),
 ("Diazinon","CCOP(=S)(OCC)Oc1cc(C)nc(C(C)C)n1",["Insecticides"],"","C12H21N2O3PS"),
 ("Dichlorvos","COP(=O)(OC)OC=C(Cl)Cl",["Insecticides"],"DDVP / Vapona","C4H7Cl2O4P"),
 ("Carbaryl","CNC(=O)Oc1cccc2ccccc12",["Insecticides"],"Sevin","C12H11NO2"),
 ("Carbofuran","CNC(=O)Oc1cccc2c1OC(C)(C)C2",["Insecticides"],"Furadan","C12H15NO3"),
 ("Aldicarb","CSC(C)(C)/C=N/OC(=O)NC",["Insecticides"],"Temik","C7H14N2O2S"),
 ("Permethrin","CC1(C)C(C=C(Cl)Cl)C1C(=O)OCc1cccc(Oc2ccccc2)c1",["Insecticides"],"","C21H20Cl2O3"),
 ("Cypermethrin","CC1(C)C(C=C(Cl)Cl)C1C(=O)OC(C#N)c1cccc(Oc2ccccc2)c1",["Insecticides"],"","C22H19Cl2NO3"),
 ("Deltamethrin","CC1(C)C(C=C(Br)Br)C1C(=O)OC(C#N)c1cccc(Oc2ccccc2)c1",["Insecticides"],"","C22H19Br2NO3"),
 ("Imidacloprid","Clc1ccc(CN2CCN/C2=N\\[N+](=O)[O-])cn1",["Insecticides"],"","C9H10ClN5O2"),
 ("Acetamiprid","Clc1ccc(CN(C)C(C)=NC#N)cn1",["Insecticides"],"","C10H11ClN4"),
 ("Thiamethoxam","[O-][N+](=O)/N=C1\\N(C)COCN1Cc1cnc(Cl)s1",["Insecticides"],"","C8H10ClN5O3S"),
 ("Fipronil","Nc1c(S(=O)C(F)(F)F)c(C#N)nn1-c1c(Cl)cc(C(F)(F)F)cc1Cl",["Insecticides"],"","C12H4Cl2F6N4OS"),
 ("Lindane","ClC1C(Cl)C(Cl)C(Cl)C(Cl)C1Cl",["Insecticides"],"\u03b3-HCH","C6H6Cl6"),
 ("Methoxychlor","COc1ccc(C(c2ccc(OC)cc2)C(Cl)(Cl)Cl)cc1",["Insecticides"],"","C16H15Cl3O2"),
 # ---- Herbicides (not drugs) ----
 ("Glyphosate","OC(=O)CNCP(=O)(O)O",["Herbicides"],"Roundup","C3H8NO5P"),
 ("Atrazine","CCNc1nc(Cl)nc(NC(C)C)n1",["Herbicides"],"","C8H14ClN5"),
 ("Simazine","CCNc1nc(Cl)nc(NCC)n1",["Herbicides"],"","C7H12ClN5"),
 ("2,4-D","OC(=O)COc1ccc(Cl)cc1Cl",["Herbicides"],"2,4-dichlorophenoxyacetic acid","C8H6Cl2O3"),
 ("MCPA","Cc1cc(Cl)ccc1OCC(=O)O",["Herbicides"],"","C9H9ClO3"),
 ("Dicamba","COc1c(Cl)cc(Cl)cc1C(=O)O",["Herbicides"],"","C8H6Cl2O3"),
 ("Paraquat","C[n+]1ccc(-c2cc[n+](C)cc2)cc1",["Herbicides"],"Gramoxone","C12H14N2+2"),
 ("Diquat","c1cc[n+]2c(c1)-c1cccc[n+]1CC2",["Herbicides"],"Reglone","C12H12N2+2"),
 ("Metolachlor","CCc1cccc(C)c1N(C(=O)CCl)C(C)COC",["Herbicides"],"Dual","C15H22ClNO2"),
 ("Alachlor","CCc1cccc(CC)c1N(COC)C(=O)CCl",["Herbicides"],"Lasso","C14H20ClNO2"),
 ("Trifluralin","CCCN(CCC)c1c([N+](=O)[O-])cc(C(F)(F)F)cc1[N+](=O)[O-]",["Herbicides"],"Treflan","C13H16F3N3O4"),
 ("Pendimethalin","CCC(CC)Nc1c([N+](=O)[O-])cc(C)c(C)c1[N+](=O)[O-]",["Herbicides"],"Prowl","C13H19N3O4"),
 ("Glufosinate","CP(=O)(O)CCC(N)C(=O)O",["Herbicides"],"Basta / Liberty","C5H12NO4P"),
 ("Mesotrione","O=C(C1C(=O)CCCC1=O)c1ccc([N+](=O)[O-])c(S(C)(=O)=O)c1",["Herbicides"],"Callisto","C14H13NO7S"),
 ("Bentazon","CC(C)N1C(=O)c2ccccc2NS1(=O)=O",["Herbicides"],"Basagran","C10H12N2O3S"),
 ("Bromoxynil","N#Cc1cc(Br)c(O)c(Br)c1",["Herbicides"],"Buctril","C7H3Br2NO"),
 ("Picloram","Nc1c(Cl)c(Cl)nc(C(=O)O)c1Cl",["Herbicides"],"Tordon","C6H3Cl3N2O2"),
 ("Clopyralid","O=C(O)c1ccc(Cl)nc1Cl",["Herbicides"],"Stinger / Transline","C6H3Cl2NO2"),
 # ---- Others (drugs that don't fit the classes above) ----
 ("Metformin","CN(C)C(=N)NC(N)=N",["Drugs","Others"],"Glucophage / Fortamet / Glumetza","C4H11N5"),
 ("Sildenafil","CCCc1nn(C)c2c1nc([nH]c2=O)-c1cc(S(=O)(=O)N2CCN(C)CC2)ccc1OCC",["Drugs","Others"],"Viagra / Revatio","C22H30N6O4S"),
 ("Omeprazole","COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1",["Drugs","Others"],"Prilosec / Losec","C17H19N3O3S"),
 ("Warfarin","CC(=O)CC(c1ccccc1)C1=C(O)c2ccccc2OC1=O",["Drugs","Others"],"Coumadin / Jantoven","C19H16O4"),
 ("Furosemide","NS(=O)(=O)c1cc(C(=O)O)c(NCc2ccco2)cc1Cl",["Drugs","Others"],"Lasix","C12H11ClN2O5S"),
 ("Lisinopril","NCCCC[C@@H](N[C@@H](CCc1ccccc1)C(=O)O)C(=O)N1CCC[C@H]1C(=O)O",["Drugs","Others"],"Prinivil / Zestril","C21H31N3O5"),
 ("Amlodipine","CCOC(=O)C1=C(COCCN)NC(C)=C(C(=O)OC)C1c1ccccc1Cl",["Drugs","Others"],"Norvasc","C20H25ClN2O5"),
 ("Losartan","CCCCc1nc(Cl)c(CO)n1Cc1ccc(-c2ccccc2-c2nnn[nH]2)cc1",["Drugs","Others"],"Cozaar","C22H23ClN6O"),
 ("Propranolol","CC(C)NCC(O)COc1cccc2ccccc12",["Drugs","Others"],"Inderal","C16H21NO2"),
 ("Metoprolol","COCCc1ccc(OCC(O)CNC(C)C)cc1",["Drugs","Others"],"Lopressor / Toprol","C15H25NO3"),
 ("Atenolol","CC(C)NCC(O)COc1ccc(CC(N)=O)cc1",["Drugs","Others"],"Tenormin","C14H22N2O3"),
 ("Salbutamol","CC(C)(C)NCC(O)c1ccc(O)c(CO)c1",["Drugs","Others"],"albuterol / Ventolin","C13H21NO3"),
 ("Acyclovir","Nc1nc2n(COCCO)cnc2c(=O)[nH]1",["Drugs","Others"],"aciclovir / Zovirax","C8H11N5O3"),
]

out, problems = [], []
for name, smi, cats, abbrev, expect in DATA:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        problems.append((name, "PARSE FAIL", smi)); continue
    formula = rdMolDescriptors.CalcMolFormula(mol)
    base = formula.rstrip("+-").rstrip("0123456789") if formula and formula[-1] in "+-" else formula
    if expect is not None and formula != expect:
        problems.append((name, f"FORMULA {formula} != expected {expect}", smi)); continue
    out.append({"name": name, "smiles": smi, "categories": cats, "abbrev": abbrev})

print(f"OK: {len(out)} / {len(DATA)}")
if problems:
    print("PROBLEMS:")
    for n, msg, smi in problems:
        print(f"  - {n}: {msg}")
    sys.exit(1)

# category counts
from collections import Counter
cc = Counter(cat for c in out for cat in c["categories"])
print("Category counts:", dict(cc))
json.dump(out, open("compounds.json","w"), indent=2, ensure_ascii=False)
print("Wrote compounds.json")
