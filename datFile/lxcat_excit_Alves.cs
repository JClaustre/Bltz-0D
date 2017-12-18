LXCat, www.lxcat.net
Generated on 11 Jul 2017. All rights reserved.

RECOMMENDED REFERENCE FORMAT
- IST-Lisbon database, www.lxcat.net, retrieved on July 11, 2017.
Be aware that some databases and solvers can additionally have instructions how to reference corresponding data.
Please check below in the headers of databases.

CROSS SECTION DATA FORMAT
In downloaded files, each collision process is defined by a block consisting of
1st line
Keyword in capitals indicating the type of the collision. Possible collision types are elastic, effective, excitation,
ionization, or attachment (capital letters required, key words are case sensitive), where "elastic" is used to denote
the elastic momentum transfer cross section and where "effective" denotes the total momentum transfer cross section (sum
of elastic momentum transfer and total inelastic cross sections).  The latter is useful for solving the Boltzmann
equation in the 2-term approximation.
2nd line
Name of the target particle species. This name is a character string, freely chosen by the user, e.g. "Ar". Optionally
for excitation processes, the name of the corresponding excited state can be specified on the same line, separated from
the first name either by arrow "->" (dash + greater than) or by double-head arrow "<->" (less than + dash +
greater than), e.g. "Ar -> Ar*" and "Ar <-> Ar*", respectively. In the later case BOLSIG+ will automatically
define the inverse superelastic process, constructing the superelastic cross-section by detailed balancing, and
considering the indicated excited state as the target. In this case, the ratio of statistical weights must be input in
the 3rd line (see below).  Alternatively, superelastic collisions could be defined explicitly as excitation collisions
with a negative electron energy loss with user input cross sections and species name, "Ar*", for example.
3rd line
For elastic and effective collisions, the ratio of the electron mass to the target particle mass. For excitation or
ionization collisions, the electron energy loss (nominally the threshold energy) in eV. For attachment, the 3rd line is
missing. In case of an excitation process where an excited state has been indicated on the 2nd line using double-head
arrow "<->", the 3rd line must specify also ratio of the statistical weights of the final state to the initial state
as the second parameter in 3rd line this is needed by BOLSIG+ to calculate the de-excitation cross-section.
from 4th line (optionally)
User comments and reference information, maximum 100 lines. The only constraint on format is that these comment lines
must not start with a number.
Finally
Table of the cross section as a function of energy. The table starts and ends by a line of dashes "------" (at least 5),
and has otherwise two numbers per line: the energy in eV and the cross section in m2.

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
DATABASE:         IST-Lisbon database
PERMLINK:         www.lxcat.net/IST-Lisbon
DESCRIPTION:      IST-Lisbon database contains up-to-date electron-neutral scattering cross sections (together with the
                  measured swarm parameters used to validate these data), resulting from the research effort of the
                  Group of Gas Discharges and Gaseous Electronics with IPFN/IST (Instituto de Plasmas e FusÃ£o Nuclear /
                  Instituto Superior TÃ©cnico), Lisbon, Portugal. The data, compiled from the literature, correspond to
                  contributions from different authors (see detailed references in the database). For each gas the
                  database presents a COMPLETE SET of cross sections, validated against measured swarm parameters by
                  solving the two-term homogeneous electron Boltzmann equation. In most cases, predictions are in
                  agreement with measurements within 1â€“20%, for reduced electric fields E/N ~ 1e-4â€“500 Td. To
                  improve predictions at low E/N, some sets need to be completed with rotational cross sections, also
                  available in the database.
CONTACT:          LL Alves and V Guerra
                  e-mail: llalves@@tecnico.ulisboa.pt
HOW TO REFERENCE: L. L. Alves, "The IST-Lisbon database on LXCat" J. Phys. Conf. Series 2014, 565, 1
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

************************************************************************************************************************

COMMENT: L.L. Alves and C.M. Ferreira quotElectron kinetics in weakly ionized helium under DC and HF applied electric
         fieldsquot 1991 J. Phys. D: Appl. Phys. 24 581-592 L.L. Alves, G. Gousset and C.M. Ferreira, quotA
         collisional-radiative model for microwave discharges in helium at low and intermediate pressuresquot 1992 J. Phys. D:
         Appl. Phys. 25 1713-1732. Cross sections were defined for energies up to 100 eV using piecewise linear interpolation,
         and were adjusted (through the multiplication factors defined below) as to yield good agreement between calculated and
         measured swarm parameters. Calculations (i) adopt 500 points energy-grid with constant step-sizes, varying between 0.02
         and 0.2 eV according to the (low/high) reduced-field values considered (ii) neglect the production of secondary
         electrons (born in ionization events) in obtaining the electron energy distribution function. 
         The DIRECT cross sections for the 2-4P1 and 2-3P3 states were obtained from the corresponding TOTAL cross sections
         (reported in the original references, see below), by correcting for the cascade transitions (G. Gousset and C.
         Boulmer-Leborgne 1984 J. Physique 45 689 S. Daviaud 1989 PhD Thesis University Paris-Sud Orsay France). The DIRECT
         cross sections for the nP1 (n>4), nP3 (n>3), nSD1,3 (n=3-7) and nFGHI1,3 (n<7) states were obtained from G.
         Gousset and C. Boulmer-Leborgne 1984 J. Physique 45 689 S. Daviaud 1989 PhD Thesis University Paris-Sud Orsay France.
         The INDIVIDUAL excitation cross sections for metastable 2S3 and 2S1 states were obtained from the corresponding LUMPED
         cross section (reported in the original references), by using the deconvolution parameter suggested by G. Gousset and C.
         Boulmer-Leborgne 1984 J. Physique 45 689 S. Daviaud 1989 PhD Thesis University Paris-Sud Orsay France.

********************************************************** He **********************************************************

EXCITATION
He <-> He(2S3)
 1.982000e+1  3.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(2S3), Excitation
PARAM.:  E = 19.82 eV, g1/g0 = 3, complete set
COMMENT: [e + He(1S1) -> e + He(2S3), Excitation]  Kato and Janev 1992 Atomic and
COMMENT: Plasmaâˆ’Material Interaction Data for Fusion (Supplement to the Journal of Nuclear
COMMENT: Fusion) 3 33 (multiplied by 0.31).
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 1.982000e+1	0.000000e+0
 1.984000e+1	2.387000e-23
 1.990000e+1	5.270000e-23
 2.010000e+1	1.550000e-22
 2.050000e+1	1.395000e-22
 2.100000e+1	1.178000e-22
 2.300000e+1	9.300000e-23
 2.500000e+1	7.750000e-23
 3.000000e+1	5.890000e-23
 4.000000e+1	3.627000e-23
 5.000000e+1	2.468000e-23
 6.000000e+1	1.708000e-23
 7.000000e+1	1.218000e-23
 8.000000e+1	8.928000e-24
 9.000000e+1	6.727000e-24
 1.000000e+2	5.177000e-24
 1.200000e+2	3.224000e-24
 1.500000e+2	1.795000e-24
 2.000000e+2	8.184000e-25
 3.000000e+2	2.644000e-25
 4.000000e+2	1.172000e-25
 5.000000e+2	6.200000e-26
 7.000000e+2	2.381000e-26
 1.000000e+3	8.649000e-27
-----------------------------

EXCITATION
He <-> He(2S1)
 2.062000e+1  1.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(2S1), Excitation
PARAM.:  E = 20.62 eV, g1/g0 = 1, complete set
COMMENT: [e + He(1S1) -> e + He(2S1), Excitation]  de Heer and Jansen 1977 J. Phys. B: At. Mol.
COMMENT: Phys. 10 3741 Belmonte et al 2007  J. Phys. D: Appl. Phys. 40 7343 (multiplied by 0.31).
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.062000e+1	0.000000e+0
 2.100000e+1	7.750000e-23
 2.200000e+1	8.060000e-23
 2.300000e+1	8.370000e-23
 2.400000e+1	8.060000e-23
 2.700000e+1	7.440000e-23
 3.000000e+1	7.130000e-23
 3.500000e+1	6.665000e-23
 4.000000e+1	6.355000e-23
 4.500000e+1	6.045000e-23
 5.000000e+1	5.735000e-23
 6.000000e+1	5.270000e-23
 8.000000e+1	4.650000e-23
 1.000000e+2	4.030000e-23
 2.000000e+2	2.480000e-23
 5.000000e+2	1.240000e-23
 9.000000e+2	7.750000e-24
 1.000000e+3	7.130000e-24
-----------------------------

EXCITATION
He <-> He(2P3)
 2.096000e+1  9.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(2P3), Excitation
PARAM.:  E = 20.96 eV, g1/g0 = 9, complete set
COMMENT: [e + He(1S1) -> e + He(2P3), Excitation]  Ralchenko et al, NIFSâˆ’DATA 59 Research
COMMENT: report NIFSâˆ’DATE serie ISSN 0915âˆ’6364 (Oct. 2000) Belmonte et al 2007  J. Phys. D:
COMMENT: Appl. Phys. 40 7343 (multiplied by 0.6).
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.096000e+1	0.000000e+0
 2.150000e+1	7.200000e-23
 2.200000e+1	9.000000e-23
 2.250000e+1	1.080000e-22
 2.300000e+1	1.170000e-22
 2.400000e+1	1.320000e-22
 2.600000e+1	1.446000e-22
 2.700000e+1	1.482000e-22
 2.800000e+1	1.488000e-22
 2.900000e+1	1.482000e-22
 3.000000e+1	1.458000e-22
 3.500000e+1	1.242000e-22
 4.000000e+1	1.014000e-22
 4.500000e+1	8.220000e-23
 5.000000e+1	6.660000e-23
 6.000000e+1	4.482000e-23
 7.000000e+1	3.132000e-23
 8.000000e+1	2.250000e-23
 9.000000e+1	1.662000e-23
 1.000000e+2	1.260000e-23
 1.200000e+2	7.620000e-24
 1.500000e+2	3.978000e-24
 2.000000e+2	1.662000e-24
 3.000000e+2	4.626000e-25
 4.000000e+2	1.848000e-25
 5.000000e+2	9.060000e-26
 7.000000e+2	3.126000e-26
 1.000000e+3	1.020000e-26
-----------------------------

EXCITATION
He <-> He(2P1)
 2.121800e+1  3.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(2P1), Excitation
PARAM.:  E = 21.218 eV, g1/g0 = 3, complete set
COMMENT: [e + He(1S1) -> e + He(2P1), Excitation] Kato and Janev 1992 Atomic and
COMMENT: Plasmaâˆ’Material Interaction Data for Fusion (Supplement to the Journal of Nuclear
COMMENT: Fusion) 3 33 (multiplied by 1.66).
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.121800e+1	0.000000e+0
 2.130000e+1	4.034000e-23
 2.150000e+1	8.333000e-23
 2.180000e+1	1.260000e-22
 2.200000e+1	1.486000e-22
 2.300000e+1	2.341000e-22
 2.400000e+1	2.988000e-22
 2.500000e+1	3.586000e-22
 2.600000e+1	4.167000e-22
 2.700000e+1	4.748000e-22
 2.800000e+1	5.312000e-22
 2.900000e+1	5.893000e-22
 3.000000e+1	6.457000e-22
 3.500000e+1	9.113000e-22
 4.000000e+1	1.135000e-21
 4.500000e+1	1.313000e-21
 5.000000e+1	1.448000e-21
 6.000000e+1	1.625000e-21
 7.000000e+1	1.716000e-21
 8.000000e+1	1.756000e-21
 9.000000e+1	1.766000e-21
 1.000000e+2	1.756000e-21
 1.200000e+2	1.705000e-21
 1.500000e+2	1.604000e-21
 2.000000e+2	1.433000e-21
 3.000000e+2	1.167000e-21
 4.000000e+2	9.877000e-22
 5.000000e+2	8.582000e-22
 7.000000e+2	6.856000e-22
 1.000000e+3	5.329000e-22
-----------------------------

EXCITATION
He <-> He(3S3)
 2.271900e+1  3.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(3S3), Excitation
PARAM.:  E = 22.719 eV, g1/g0 = 3, complete set
COMMENT: [e + He(1S1) -> e + He(3S3), Excitation]  de Heer 1998 INDC(NDS)âˆ’385, Critically
COMMENT: Assessed Electronâˆ’Impact Excitation Cross Sections for He (11 S), FOMâˆ’Institute for
COMMENT: Atomic and Molecular Physics Amsterdam 1.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.271900e+1	0.000000e+0
 2.272000e+1	9.310000e-23
 2.280000e+1	9.070000e-23
 2.300000e+1	8.590000e-23
 2.330000e+1	8.120000e-23
 2.350000e+1	7.900000e-23
 2.380000e+1	7.640000e-23
 2.400000e+1	7.500000e-23
 2.450000e+1	7.230000e-23
 2.500000e+1	7.010000e-23
 2.600000e+1	6.660000e-23
 2.700000e+1	6.340000e-23
 2.800000e+1	6.050000e-23
 2.900000e+1	5.770000e-23
 3.000000e+1	5.500000e-23
 3.500000e+1	4.300000e-23
 4.000000e+1	3.360000e-23
 4.500000e+1	2.650000e-23
 5.000000e+1	2.120000e-23
 6.000000e+1	1.400000e-23
 7.000000e+1	9.680000e-24
 8.000000e+1	6.960000e-24
 9.000000e+1	5.160000e-24
 1.000000e+2	3.930000e-24
 1.200000e+2	2.420000e-24
 1.500000e+2	1.320000e-24
 2.000000e+2	5.960000e-25
 3.000000e+2	1.880000e-25
 4.000000e+2	8.210000e-26
 5.000000e+2	4.290000e-26
 7.000000e+2	1.600000e-26
 1.000000e+3	5.580000e-27
-----------------------------

EXCITATION
He <-> He(3S1)
 2.291900e+1  1.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(3S1), Excitation
PARAM.:  E = 22.919 eV, g1/g0 = 1, complete set
COMMENT: [e + He(1S1) -> e + He(3S1), Excitation]  de Heer 1998 INDC(NDS)âˆ’385, Critically
COMMENT: Assessed Electronâˆ’Impact Excitation Cross Sections for He (11 S), FOMâˆ’Institute for
COMMENT: Atomic and Molecular Physics Amsterdam 1.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.291900e+1	0.000000e+0
 2.292000e+1	4.160000e-23
 2.300000e+1	4.170000e-23
 2.320000e+1	4.180000e-23
 2.350000e+1	4.200000e-23
 2.370000e+1	4.210000e-23
 2.400000e+1	4.230000e-23
 2.450000e+1	4.250000e-23
 2.550000e+1	4.280000e-23
 2.650000e+1	4.290000e-23
 2.750000e+1	4.290000e-23
 2.850000e+1	4.280000e-23
 2.950000e+1	4.260000e-23
 3.050000e+1	4.240000e-23
 3.500000e+1	4.080000e-23
 4.000000e+1	3.870000e-23
 4.500000e+1	3.660000e-23
 5.000000e+1	3.460000e-23
 6.000000e+1	3.110000e-23
 7.000000e+1	2.830000e-23
 8.000000e+1	2.600000e-23
 9.000000e+1	2.420000e-23
 1.000000e+2	2.260000e-23
 1.200000e+2	2.030000e-23
 1.500000e+2	1.780000e-23
 2.000000e+2	1.520000e-23
 3.000000e+2	1.200000e-23
 4.000000e+2	9.910000e-24
 5.000000e+2	8.390000e-24
 7.000000e+2	6.360000e-24
 1.000000e+3	4.620000e-24
-----------------------------

EXCITATION
He <-> He(3P3)
 2.300900e+1  9.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(3P3), Excitation
PARAM.:  E = 23.009 eV, g1/g0 = 9, complete set
COMMENT: [e + He(1S1) -> e + He(3P3), Excitation]  de Heer 1998 INDC(NDS)âˆ’385, Critically
COMMENT: Assessed Electronâˆ’Impact Excitation Cross Sections for He (11 S), FOMâˆ’Institute for
COMMENT: Atomic and Molecular Physics Amsterdam 1 (multiplied by 0.6).
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.300900e+1	0.000000e+0
 2.301000e+1	1.668000e-23
 2.320000e+1	1.782000e-23
 2.350000e+1	1.950000e-23
 2.400000e+1	2.202000e-23
 2.500000e+1	2.592000e-23
 2.600000e+1	2.862000e-23
 2.700000e+1	3.048000e-23
 2.800000e+1	3.168000e-23
 2.900000e+1	3.234000e-23
 3.000000e+1	3.264000e-23
 3.500000e+1	3.072000e-23
 4.000000e+1	2.676000e-23
 4.500000e+1	2.274000e-23
 5.000000e+1	1.914000e-23
 6.000000e+1	1.356000e-23
 7.000000e+1	9.720000e-24
 8.000000e+1	7.080000e-24
 9.000000e+1	5.268000e-24
 1.000000e+2	3.990000e-24
 1.200000e+2	2.406000e-24
 1.500000e+2	1.248000e-24
 2.000000e+2	5.130000e-25
 3.000000e+2	1.404000e-25
 4.000000e+2	5.568000e-26
 5.000000e+2	2.718000e-26
 7.000000e+2	9.300000e-27
 1.000000e+3	3.012000e-27
-----------------------------

EXCITATION
He <-> He(3D3)
 2.306900e+1  1.500000e+1
SPECIES: e / He
PROCESS: E + He <-> E + He(3D3), Excitation
PARAM.:  E = 23.069 eV, g1/g0 = 15, complete set
COMMENT: [e + He(1S1) -> e + He(3D3), Excitation]  de Heer 1998 INDC(NDS)âˆ’385, Critically
COMMENT: Assessed Electronâˆ’Impact Excitation Cross Sections for He (11 S), FOMâˆ’Institute for
COMMENT: Atomic and Molecular Physics Amsterdam 1.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.306900e+1	0.000000e+0
 2.307000e+1	2.100000e-24
 2.320000e+1	3.890000e-24
 2.330000e+1	6.160000e-24
 2.350000e+1	8.470000e-24
 2.380000e+1	1.010000e-23
 2.400000e+1	1.070000e-23
 2.450000e+1	1.160000e-23
 2.500000e+1	1.210000e-23
 2.600000e+1	1.230000e-23
 2.700000e+1	1.220000e-23
 2.800000e+1	1.190000e-23
 2.900000e+1	1.140000e-23
 3.000000e+1	1.090000e-23
 3.500000e+1	8.150000e-24
 4.000000e+1	5.950000e-24
 4.500000e+1	4.370000e-24
 5.000000e+1	3.260000e-24
 6.000000e+1	1.910000e-24
 7.000000e+1	1.190000e-24
 8.000000e+1	7.800000e-25
 9.000000e+1	5.350000e-25
 1.000000e+2	3.800000e-25
 1.200000e+2	2.100000e-25
 1.500000e+2	1.010000e-25
 2.000000e+2	3.890000e-26
 3.000000e+2	1.020000e-26
 4.000000e+2	4.010000e-27
 5.000000e+2	1.950000e-27
 7.000000e+2	6.670000e-28
 1.000000e+3	2.170000e-28
-----------------------------

EXCITATION
He <-> He(3D1)
 2.306900e+1  5.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(3D1), Excitation
PARAM.:  E = 23.069 eV, g1/g0 = 5, complete set
COMMENT: [e + He(1S1) -> e + He(3D1), Excitation]  de Heer 1998 INDC(NDS)âˆ’385, Critically
COMMENT: Assessed Electronâˆ’Impact Excitation Cross Sections for He (11 S), FOMâˆ’Institute for
COMMENT: Atomic and Molecular Physics Amsterdam 1.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.306900e+1	0.000000e+0
 2.307000e+1	2.400000e-23
 2.310000e+1	2.390000e-23
 2.320000e+1	2.360000e-23
 2.350000e+1	2.280000e-23
 2.370000e+1	2.240000e-23
 2.400000e+1	2.180000e-23
 2.450000e+1	2.090000e-23
 2.500000e+1	2.040000e-23
 2.600000e+1	1.970000e-23
 2.700000e+1	1.960000e-23
 2.800000e+1	1.980000e-23
 2.900000e+1	2.020000e-23
 3.000000e+1	2.070000e-23
 3.500000e+1	2.350000e-23
 4.000000e+1	2.520000e-23
 4.500000e+1	2.580000e-23
 5.000000e+1	2.560000e-23
 6.000000e+1	2.400000e-23
 7.000000e+1	2.190000e-23
 8.000000e+1	1.980000e-23
 9.000000e+1	1.790000e-23
 1.000000e+2	1.630000e-23
 1.200000e+2	1.360000e-23
 1.500000e+2	1.080000e-23
 2.000000e+2	7.860000e-24
 3.000000e+2	5.010000e-24
 4.000000e+2	3.650000e-24
 5.000000e+2	2.860000e-24
 7.000000e+2	1.990000e-24
 1.000000e+3	1.360000e-24
-----------------------------

EXCITATION
He <-> He(3P1)
 2.309000e+1  3.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(3P1), Excitation
PARAM.:  E = 23.09 eV, g1/g0 = 3, complete set
COMMENT: [e + He(1S1) -> e + He(3P1), Excitation]  de Heer 1998 INDC(NDS)âˆ’385, Critically
COMMENT: Assessed Electronâˆ’Impact Excitation Cross Sections for He (11 S), FOMâˆ’Institute for
COMMENT: Atomic and Molecular Physics Amsterdam 1  (multiplied by 1.66).
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.309000e+1	0.000000e+0
 2.310000e+1	1.232000e-23
 2.320000e+1	1.497000e-23
 2.350000e+1	2.241000e-23
 2.370000e+1	2.706000e-23
 2.400000e+1	3.353000e-23
 2.450000e+1	4.333000e-23
 2.500000e+1	5.196000e-23
 2.700000e+1	8.001000e-23
 3.000000e+1	1.144000e-22
 4.000000e+1	2.191000e-22
 5.000000e+1	3.021000e-22
 6.000000e+1	3.569000e-22
 7.000000e+1	3.901000e-22
 8.000000e+1	4.100000e-22
 9.000000e+1	4.183000e-22
 1.000000e+2	4.200000e-22
 1.200000e+2	4.150000e-22
 1.500000e+2	3.951000e-22
 2.000000e+2	3.552000e-22
 3.000000e+2	2.922000e-22
 4.000000e+2	2.473000e-22
 5.000000e+2	2.141000e-22
 7.000000e+2	1.710000e-22
 1.000000e+3	1.328000e-22
-----------------------------

EXCITATION
He <-> He(4S3)
 2.358900e+1  3.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(4S3), Excitation
PARAM.:  E = 23.589 eV, g1/g0 = 3, complete set
COMMENT: [e + He(1S1) -> e + He(4S3), Excitation]  de Heer 1998 INDC(NDS)âˆ’385, Critically
COMMENT: Assessed Electronâˆ’Impact Excitation Cross Sections for He (11 S), FOMâˆ’Institute for
COMMENT: Atomic and Molecular Physics Amsterdam 1.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.358900e+1	0.000000e+0
 2.359000e+1	3.270000e-23
 2.370000e+1	2.920000e-23
 2.380000e+1	2.760000e-23
 2.400000e+1	2.600000e-23
 2.420000e+1	2.530000e-23
 2.450000e+1	2.470000e-23
 2.480000e+1	2.440000e-23
 2.500000e+1	2.430000e-23
 2.550000e+1	2.410000e-23
 2.600000e+1	2.380000e-23
 2.700000e+1	2.340000e-23
 2.800000e+1	2.270000e-23
 2.900000e+1	2.200000e-23
 3.000000e+1	2.120000e-23
 3.500000e+1	1.700000e-23
 4.000000e+1	1.330000e-23
 4.500000e+1	1.050000e-23
 5.000000e+1	8.320000e-24
 6.000000e+1	5.430000e-24
 7.000000e+1	3.710000e-24
 8.000000e+1	2.640000e-24
 9.000000e+1	1.940000e-24
 1.000000e+2	1.470000e-24
 1.200000e+2	8.940000e-25
 1.500000e+2	4.830000e-25
 2.000000e+2	2.150000e-25
 3.000000e+2	6.690000e-26
 4.000000e+2	2.900000e-26
 5.000000e+2	1.510000e-26
 7.000000e+2	5.590000e-27
 1.000000e+3	1.940000e-27
-----------------------------

EXCITATION
He <-> He(4S1)
 2.366900e+1  1.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(4S1), Excitation
PARAM.:  E = 23.669 eV, g1/g0 = 1, complete set
COMMENT: [e + He(1S1) -> e + He(4S1), Excitation]  de Heer 1998 INDC(NDS)âˆ’385, Critically
COMMENT: Assessed Electronâˆ’Impact Excitation Cross Sections for He (11 S), FOMâˆ’Institute for
COMMENT: Atomic and Molecular Physics Amsterdam 1.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.366900e+1	0.000000e+0
 2.367000e+1	1.080000e-23
 2.380000e+1	1.100000e-23
 2.390000e+1	1.110000e-23
 2.400000e+1	1.120000e-23
 2.430000e+1	1.150000e-23
 2.450000e+1	1.170000e-23
 2.500000e+1	1.210000e-23
 2.600000e+1	1.280000e-23
 2.700000e+1	1.330000e-23
 2.800000e+1	1.380000e-23
 2.900000e+1	1.410000e-23
 3.000000e+1	1.430000e-23
 3.100000e+1	1.450000e-23
 3.500000e+1	1.470000e-23
 4.000000e+1	1.440000e-23
 4.500000e+1	1.380000e-23
 5.000000e+1	1.310000e-23
 6.000000e+1	1.180000e-23
 7.000000e+1	1.080000e-23
 8.000000e+1	9.900000e-24
 9.000000e+1	9.220000e-24
 1.000000e+2	8.680000e-24
 1.200000e+2	7.850000e-24
 1.500000e+2	7.000000e-24
 2.000000e+2	6.010000e-24
 3.000000e+2	4.690000e-24
 4.000000e+2	3.820000e-24
 5.000000e+2	3.210000e-24
 7.000000e+2	2.430000e-24
 1.000000e+3	1.770000e-24
-----------------------------

EXCITATION
He <-> He(4P3)
 2.370900e+1  9.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(4P3), Excitation
PARAM.:  E = 23.709 eV, g1/g0 = 9, complete set
COMMENT: [e + He(1S1) -> e + He(4P3), Excitation]  de Heer 1998 INDC(NDS)âˆ’385, Critically
COMMENT: Assessed Electronâˆ’Impact Excitation Cross Sections for He (11 S), FOMâˆ’Institute for
COMMENT: Atomic and Molecular Physics Amsterdam 1 (multiplied by 0.6).
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.370900e+1	0.000000e+0
 2.371000e+1	6.360000e-24
 2.380000e+1	6.600000e-24
 2.390000e+1	6.840000e-24
 2.400000e+1	7.080000e-24
 2.420000e+1	7.560000e-24
 2.450000e+1	8.220000e-24
 2.480000e+1	8.760000e-24
 2.500000e+1	9.120000e-24
 2.550000e+1	9.960000e-24
 2.600000e+1	1.062000e-23
 2.700000e+1	1.164000e-23
 2.800000e+1	1.224000e-23
 2.900000e+1	1.266000e-23
 3.000000e+1	1.290000e-23
 3.500000e+1	1.224000e-23
 4.000000e+1	1.068000e-23
 4.500000e+1	9.060000e-24
 5.000000e+1	7.620000e-24
 6.000000e+1	5.394000e-24
 7.000000e+1	3.882000e-24
 8.000000e+1	2.850000e-24
 9.000000e+1	2.136000e-24
 1.000000e+2	1.626000e-24
 1.200000e+2	9.900000e-25
 1.500000e+2	5.190000e-25
 2.000000e+2	2.154000e-25
 3.000000e+2	5.928000e-26
 4.000000e+2	2.340000e-26
 5.000000e+2	1.140000e-26
 7.000000e+2	3.876000e-27
 1.000000e+3	1.248000e-27
-----------------------------

EXCITATION
He <-> He(4D3)
 2.373900e+1  1.500000e+1
SPECIES: e / He
PROCESS: E + He <-> E + He(4D3), Excitation
PARAM.:  E = 23.739 eV, g1/g0 = 15, complete set
COMMENT: [e + He(1S1) -> e + He(4D3), Excitation]  de Heer 1998 INDC(NDS)âˆ’385, Critically
COMMENT: Assessed Electronâˆ’Impact Excitation Cross Sections for He (11 S), FOMâˆ’Institute for
COMMENT: Atomic and Molecular Physics Amsterdam 1.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.373900e+1	0.000000e+0
 2.374000e+1	1.970000e-24
 2.380000e+1	2.880000e-24
 2.390000e+1	3.640000e-24
 2.400000e+1	4.070000e-24
 2.420000e+1	4.620000e-24
 2.450000e+1	5.150000e-24
 2.480000e+1	5.540000e-24
 2.500000e+1	5.750000e-24
 2.550000e+1	6.170000e-24
 2.600000e+1	6.460000e-24
 2.700000e+1	6.810000e-24
 2.800000e+1	6.920000e-24
 2.900000e+1	6.890000e-24
 3.000000e+1	6.750000e-24
 3.500000e+1	5.440000e-24
 4.000000e+1	4.110000e-24
 4.500000e+1	3.070000e-24
 5.000000e+1	2.310000e-24
 6.000000e+1	1.360000e-24
 7.000000e+1	8.470000e-25
 8.000000e+1	5.540000e-25
 9.000000e+1	3.780000e-25
 1.000000e+2	2.670000e-25
 1.200000e+2	1.450000e-25
 1.500000e+2	6.810000e-26
 2.000000e+2	2.550000e-26
 3.000000e+2	6.390000e-27
 4.000000e+2	2.410000e-27
 5.000000e+2	1.140000e-27
 7.000000e+2	3.750000e-28
 1.000000e+3	1.180000e-28
-----------------------------

EXCITATION
He <-> He(4D1)
 2.373900e+1  5.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(4D1), Excitation
PARAM.:  E = 23.739 eV, g1/g0 = 5, complete set
COMMENT: [e + He(1S1) -> e + He(4D1), Excitation]  de Heer 1998 INDC(NDS)âˆ’385, Critically
COMMENT: Assessed Electronâˆ’Impact Excitation Cross Sections for He (11 S), FOMâˆ’Institute for
COMMENT: Atomic and Molecular Physics Amsterdam 1.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.373900e+1	0.000000e+0
 2.374000e+1	7.890000e-24
 2.380000e+1	7.910000e-24
 2.390000e+1	7.950000e-24
 2.400000e+1	7.990000e-24
 2.430000e+1	8.120000e-24
 2.450000e+1	8.210000e-24
 2.470000e+1	8.310000e-24
 2.500000e+1	8.460000e-24
 2.600000e+1	9.020000e-24
 2.700000e+1	9.590000e-24
 2.800000e+1	1.020000e-23
 2.900000e+1	1.070000e-23
 3.000000e+1	1.120000e-23
 3.500000e+1	1.290000e-23
 4.000000e+1	1.360000e-23
 4.500000e+1	1.360000e-23
 5.000000e+1	1.330000e-23
 6.000000e+1	1.210000e-23
 7.000000e+1	1.090000e-23
 8.000000e+1	9.780000e-24
 9.000000e+1	8.800000e-24
 1.000000e+2	7.960000e-24
 1.200000e+2	6.630000e-24
 1.500000e+2	5.240000e-24
 2.000000e+2	3.840000e-24
 3.000000e+2	2.460000e-24
 4.000000e+2	1.800000e-24
 5.000000e+2	1.420000e-24
 7.000000e+2	9.900000e-25
 1.000000e+3	6.810000e-25
-----------------------------

EXCITATION
He <-> He(4F3)
 2.373900e+1  2.100000e+1
SPECIES: e / He
PROCESS: E + He <-> E + He(4F3), Excitation
PARAM.:  E = 23.739 eV, g1/g0 = 21, complete set
COMMENT: [e + He(1S1) -> e + He(4F3), Excitation]  de Heer 1998 INDC(NDS)âˆ’385, Critically
COMMENT: Assessed Electronâˆ’Impact Excitation Cross Sections for He (11 S), FOMâˆ’Institute for
COMMENT: Atomic and Molecular Physics Amsterdam 1.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.373900e+1	0.000000e+0
 2.374000e+1	3.860000e-24
 2.380000e+1	3.570000e-24
 2.390000e+1	3.180000e-24
 2.400000e+1	2.890000e-24
 2.420000e+1	2.450000e-24
 2.450000e+1	2.030000e-24
 2.480000e+1	1.750000e-24
 2.500000e+1	1.610000e-24
 2.550000e+1	1.350000e-24
 2.600000e+1	1.170000e-24
 2.700000e+1	9.220000e-25
 2.800000e+1	7.600000e-25
 2.900000e+1	6.420000e-25
 3.000000e+1	5.500000e-25
 3.500000e+1	2.910000e-25
 4.000000e+1	1.740000e-25
 4.500000e+1	1.110000e-25
 5.000000e+1	7.440000e-26
 6.000000e+1	3.740000e-26
 7.000000e+1	2.090000e-26
 8.000000e+1	1.260000e-26
 9.000000e+1	8.030000e-27
 1.000000e+2	5.360000e-27
 1.200000e+2	2.660000e-27
 1.500000e+2	1.130000e-27
 2.000000e+2	3.670000e-28
 3.000000e+2	7.500000e-29
 4.000000e+2	2.420000e-29
 5.000000e+2	1.000000e-29
 7.000000e+2	2.640000e-30
 1.000000e+3	6.390000e-31
-----------------------------

EXCITATION
He <-> He(4F1)
 2.373900e+1  7.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(4F1), Excitation
PARAM.:  E = 23.739 eV, g1/g0 = 7, complete set
COMMENT: [e + He(1S1) -> e + He(4F1), Excitation]  de Heer 1998 INDC(NDS)âˆ’385, Critically
COMMENT: Assessed Electronâˆ’Impact Excitation Cross Sections for He (11 S), FOMâˆ’Institute for
COMMENT: Atomic and Molecular Physics Amsterdam 1.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.373900e+1	0.000000e+0
 2.374000e+1	1.760000e-24
 2.380000e+1	1.750000e-24
 2.390000e+1	1.740000e-24
 2.400000e+1	1.720000e-24
 2.430000e+1	1.680000e-24
 2.450000e+1	1.660000e-24
 2.470000e+1	1.630000e-24
 2.500000e+1	1.600000e-24
 2.600000e+1	1.490000e-24
 2.700000e+1	1.390000e-24
 2.800000e+1	1.300000e-24
 2.900000e+1	1.220000e-24
 3.000000e+1	1.140000e-24
 3.500000e+1	8.620000e-25
 4.000000e+1	6.760000e-25
 4.500000e+1	5.470000e-25
 5.000000e+1	4.520000e-25
 6.000000e+1	3.250000e-25
 7.000000e+1	2.460000e-25
 8.000000e+1	1.930000e-25
 9.000000e+1	1.560000e-25
 1.000000e+2	1.290000e-25
 1.200000e+2	9.240000e-26
 1.500000e+2	6.170000e-26
 2.000000e+2	3.730000e-26
 3.000000e+2	1.990000e-26
 4.000000e+2	1.360000e-26
 5.000000e+2	1.040000e-26
 7.000000e+2	7.270000e-27
 1.000000e+3	5.070000e-27
-----------------------------

EXCITATION
He <-> He(4P1)
 2.374000e+1  3.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(4P1), Excitation
PARAM.:  E = 23.74 eV, g1/g0 = 3, complete set
COMMENT: [e + He(1S1) -> e + He(4P1), Excitation]  de Heer 1998 INDC(NDS)âˆ’385, Critically
COMMENT: Assessed Electronâˆ’Impact Excitation Cross Sections for He (11 S), FOMâˆ’Institute for
COMMENT: Atomic and Molecular Physics Amsterdam 1  (multiplied by 1.66).
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.374000e+1	0.000000e+0
 2.380000e+1	2.470000e-24
 2.390000e+1	4.050000e-24
 2.400000e+1	5.560000e-24
 2.450000e+1	1.230000e-23
 2.500000e+1	1.780000e-23
 2.600000e+1	2.620000e-23
 2.700000e+1	3.250000e-23
 2.800000e+1	3.780000e-23
 3.000000e+1	4.660000e-23
 4.000000e+1	9.180000e-23
 5.000000e+1	1.310000e-22
 6.000000e+1	1.570000e-22
 7.000000e+1	1.710000e-22
 8.000000e+1	1.790000e-22
 9.000000e+1	1.830000e-22
 1.000000e+2	1.830000e-22
 1.200000e+2	1.780000e-22
 1.500000e+2	1.680000e-22
 2.000000e+2	1.490000e-22
 3.000000e+2	1.200000e-22
 4.000000e+2	1.000000e-22
 5.000000e+2	8.670000e-23
 7.000000e+2	6.870000e-23
 1.000000e+3	5.310000e-23
-----------------------------

EXCITATION
He <-> He(5S3)
 2.397200e+1  3.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(5S3), Excitation
PARAM.:  E = 23.972 eV, g1/g0 = 3, complete set
COMMENT: [e + He(1S1) -> e + He(5S3), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.397200e+1	0.000000e+0
 2.400000e+1	1.622828e-23
 2.420000e+1	1.405581e-23
 2.450000e+1	1.310830e-23
 2.480000e+1	1.273964e-23
 2.500000e+1	1.259224e-23
 2.550000e+1	1.242155e-23
 2.600000e+1	1.231155e-23
 2.700000e+1	1.206799e-23
 2.800000e+1	1.178064e-23
 2.900000e+1	1.142789e-23
 3.000000e+1	1.104816e-23
 3.500000e+1	8.941359e-24
 4.000000e+1	7.048573e-24
 4.500000e+1	5.579450e-24
 5.000000e+1	4.435841e-24
 6.000000e+1	2.920153e-24
 7.000000e+1	1.996724e-24
 8.000000e+1	1.420788e-24
 9.000000e+1	1.044143e-24
 1.000000e+2	7.905851e-25
 1.200000e+2	4.856297e-25
 1.500000e+2	2.638869e-25
 2.000000e+2	1.187347e-25
 3.000000e+2	3.783982e-26
 4.000000e+2	1.607193e-26
 5.000000e+2	8.292303e-27
 7.000000e+2	3.130804e-27
 1.000000e+3	1.091507e-27
-----------------------------

EXCITATION
He <-> He(5S1)
 2.401100e+1  1.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(5S1), Excitation
PARAM.:  E = 24.011 eV, g1/g0 = 1, complete set
COMMENT: [e + He(1S1) -> e + He(5S1), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.401100e+1	0.000000e+0
 2.430000e+1	5.713379e-24
 2.450000e+1	5.814342e-24
 2.500000e+1	6.051479e-24
 2.600000e+1	6.422814e-24
 2.700000e+1	6.712588e-24
 2.800000e+1	6.964995e-24
 2.900000e+1	7.156681e-24
 3.000000e+1	7.278484e-24
 3.100000e+1	7.379447e-24
 3.500000e+1	7.513824e-24
 4.000000e+1	7.390046e-24
 4.500000e+1	7.104405e-24
 5.000000e+1	6.757502e-24
 6.000000e+1	6.097651e-24
 7.000000e+1	5.579902e-24
 8.000000e+1	5.120539e-24
 9.000000e+1	4.764619e-24
 1.000000e+2	4.482965e-24
 1.200000e+2	4.054986e-24
 1.500000e+2	3.614541e-24
 2.000000e+2	3.105577e-24
 3.000000e+2	2.429737e-24
 4.000000e+2	1.980847e-24
 5.000000e+2	1.665437e-24
 7.000000e+2	1.263778e-24
 1.000000e+3	9.220493e-25
-----------------------------

EXCITATION
He <-> He(5P3)
 2.402800e+1  9.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(5P3), Excitation
PARAM.:  E = 24.028 eV, g1/g0 = 9, complete set
COMMENT: [e + He(1S1) -> e + He(5P3), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.402800e+1	0.000000e+0
 2.420000e+1	3.475927e-24
 2.450000e+1	3.839673e-24
 2.480000e+1	4.175694e-24
 2.500000e+1	4.363557e-24
 2.550000e+1	4.808319e-24
 2.600000e+1	5.204153e-24
 2.700000e+1	5.772480e-24
 2.800000e+1	6.152684e-24
 2.900000e+1	6.399128e-24
 3.000000e+1	6.555859e-24
 3.500000e+1	6.298284e-24
 4.000000e+1	5.552992e-24
 4.500000e+1	4.737826e-24
 5.000000e+1	3.999323e-24
 6.000000e+1	2.852514e-24
 7.000000e+1	2.059528e-24
 8.000000e+1	1.515319e-24
 9.000000e+1	1.137312e-24
 1.000000e+2	8.671788e-25
 1.200000e+2	5.328189e-25
 1.500000e+2	2.817359e-25
 2.000000e+2	1.185396e-25
 3.000000e+2	3.353499e-26
 4.000000e+2	1.295636e-26
 5.000000e+2	6.244644e-27
 7.000000e+2	2.163515e-27
 1.000000e+3	6.985213e-28
-----------------------------

EXCITATION
He <-> He(5D3)
 2.404300e+1  1.500000e+1
SPECIES: e / He
PROCESS: E + He <-> E + He(5D3), Excitation
PARAM.:  E = 24.043 eV, g1/g0 = 15, complete set
COMMENT: [e + He(1S1) -> e + He(5D3), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.404300e+1	0.000000e+0
 2.420000e+1	1.828641e-24
 2.450000e+1	2.347368e-24
 2.480000e+1	2.621725e-24
 2.500000e+1	2.757127e-24
 2.550000e+1	3.019004e-24
 2.600000e+1	3.208933e-24
 2.700000e+1	3.424939e-24
 2.800000e+1	3.522904e-24
 2.900000e+1	3.533368e-24
 3.000000e+1	3.483458e-24
 3.500000e+1	2.845230e-24
 4.000000e+1	2.173880e-24
 4.500000e+1	1.633032e-24
 5.000000e+1	1.232406e-24
 6.000000e+1	7.335844e-25
 7.000000e+1	4.571406e-25
 8.000000e+1	2.989722e-25
 9.000000e+1	2.038916e-25
 1.000000e+2	1.439608e-25
 1.200000e+2	7.902554e-26
 1.500000e+2	3.738091e-26
 2.000000e+2	1.417001e-26
 3.000000e+2	3.646482e-27
 4.000000e+2	1.337999e-27
 5.000000e+2	6.251939e-28
 7.000000e+2	2.095045e-28
 1.000000e+3	6.601656e-29
-----------------------------

EXCITATION
He <-> He(5D1)
 2.404300e+1  5.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(5D1), Excitation
PARAM.:  E = 24.043 eV, g1/g0 = 5, complete set
COMMENT: [e + He(1S1) -> e + He(5D1), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.404300e+1	0.000000e+0
 2.430000e+1	4.088981e-24
 2.450000e+1	4.132632e-24
 2.470000e+1	4.177171e-24
 2.500000e+1	4.250066e-24
 2.600000e+1	4.523362e-24
 2.700000e+1	4.809794e-24
 2.800000e+1	5.111101e-24
 2.900000e+1	5.383913e-24
 3.000000e+1	5.636655e-24
 3.500000e+1	6.527256e-24
 4.000000e+1	6.926709e-24
 4.500000e+1	6.963200e-24
 5.000000e+1	6.829149e-24
 6.000000e+1	6.242118e-24
 7.000000e+1	5.635537e-24
 8.000000e+1	5.065746e-24
 9.000000e+1	4.563074e-24
 1.000000e+2	4.130257e-24
 1.200000e+2	3.446560e-24
 1.500000e+2	2.728168e-24
 2.000000e+2	2.002571e-24
 3.000000e+2	1.286498e-24
 4.000000e+2	9.388031e-25
 5.000000e+2	7.394210e-25
 7.000000e+2	5.166871e-25
 1.000000e+3	3.553838e-25
-----------------------------

EXCITATION
He <-> He(5FG3)
 2.404300e+1  4.800000e+1
SPECIES: e / He
PROCESS: E + He <-> E + He(5FG3), Excitation
PARAM.:  E = 24.043 eV, g1/g0 = 48, complete set
COMMENT: [e + He(1S1) -> e + He(5FG3), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.404300e+1	0.000000e+0
 2.420000e+1	1.644131e-24
 2.450000e+1	1.267710e-24
 2.480000e+1	1.050567e-24
 2.500000e+1	9.524742e-25
 2.550000e+1	7.776065e-25
 2.600000e+1	6.600328e-25
 2.700000e+1	5.156973e-25
 2.800000e+1	4.186780e-25
 2.900000e+1	3.510029e-25
 3.000000e+1	2.995851e-25
 3.500000e+1	1.608061e-25
 4.000000e+1	9.518728e-26
 4.500000e+1	6.052676e-26
 5.000000e+1	4.047777e-26
 6.000000e+1	2.059542e-26
 7.000000e+1	1.145344e-26
 8.000000e+1	6.883884e-27
 9.000000e+1	4.379376e-27
 1.000000e+2	2.918306e-27
 1.200000e+2	1.467484e-27
 1.500000e+2	6.284099e-28
 2.000000e+2	2.077918e-28
 3.000000e+2	4.410830e-29
 4.000000e+2	1.371452e-29
 5.000000e+2	5.582659e-30
 7.000000e+2	1.519541e-30
 1.000000e+3	3.706319e-31
-----------------------------

EXCITATION
He <-> He(5FG1)
 2.404300e+1  1.600000e+1
SPECIES: e / He
PROCESS: E + He <-> E + He(5FG1), Excitation
PARAM.:  E = 24.043 eV, g1/g0 = 16, complete set
COMMENT: [e + He(1S1) -> e + He(5FG1), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.404300e+1	0.000000e+0
 2.430000e+1	8.815893e-25
 2.450000e+1	8.677933e-25
 2.470000e+1	8.557753e-25
 2.500000e+1	8.359562e-25
 2.600000e+1	7.815167e-25
 2.700000e+1	7.292741e-25
 2.800000e+1	6.820211e-25
 2.900000e+1	6.397579e-25
 3.000000e+1	5.993192e-25
 3.500000e+1	4.540248e-25
 4.000000e+1	3.558083e-25
 4.500000e+1	2.876295e-25
 5.000000e+1	2.376145e-25
 6.000000e+1	1.713654e-25
 7.000000e+1	1.295555e-25
 8.000000e+1	1.015789e-25
 9.000000e+1	8.204194e-26
 1.000000e+2	6.780741e-26
 1.200000e+2	4.873978e-26
 1.500000e+2	3.259066e-26
 2.000000e+2	1.973359e-26
 3.000000e+2	1.052895e-26
 4.000000e+2	7.127411e-27
 5.000000e+2	5.429061e-27
 7.000000e+2	3.793626e-27
 1.000000e+3	2.643626e-27
-----------------------------

EXCITATION
He <-> He(5P1)
 2.404600e+1  3.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(5P1), Excitation
PARAM.:  E = 24.046 eV, g1/g0 = 3, complete set
COMMENT: [e + He(1S1) -> e + He(5P1), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.404600e+1	0.000000e+0
 2.450000e+1	4.159851e-24
 2.500000e+1	7.333545e-24
 2.600000e+1	1.200071e-23
 2.700000e+1	1.553896e-23
 2.800000e+1	1.839302e-23
 3.000000e+1	2.300477e-23
 4.000000e+1	4.583129e-23
 5.000000e+1	6.580331e-23
 6.000000e+1	7.937422e-23
 7.000000e+1	8.691765e-23
 8.000000e+1	9.123373e-23
 9.000000e+1	9.346297e-23
 1.000000e+2	9.369600e-23
 1.200000e+2	9.133019e-23
 1.500000e+2	8.633965e-23
 2.000000e+2	7.677994e-23
 3.000000e+2	6.200314e-23
 4.000000e+2	5.171783e-23
 5.000000e+2	4.482085e-23
 7.000000e+2	3.558219e-23
 1.000000e+3	2.752379e-23
-----------------------------

EXCITATION
He <-> He(6S3)
 2.416900e+1  3.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(6S3), Excitation
PARAM.:  E = 24.169 eV, g1/g0 = 3, complete set
COMMENT: [e + He(1S1) -> e + He(6S3), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.416900e+1	0.000000e+0
 2.420000e+1	9.365876e-24
 2.450000e+1	7.900148e-24
 2.480000e+1	7.490377e-24
 2.500000e+1	7.374679e-24
 2.550000e+1	7.215802e-24
 2.600000e+1	7.154792e-24
 2.700000e+1	7.009464e-24
 2.800000e+1	6.864089e-24
 2.900000e+1	6.661616e-24
 3.000000e+1	6.450661e-24
 3.500000e+1	5.244281e-24
 4.000000e+1	4.149395e-24
 4.500000e+1	3.288749e-24
 5.000000e+1	2.618856e-24
 6.000000e+1	1.731121e-24
 7.000000e+1	1.184131e-24
 8.000000e+1	8.425628e-25
 9.000000e+1	6.192243e-25
 1.000000e+2	4.686865e-25
 1.200000e+2	2.892507e-25
 1.500000e+2	1.575971e-25
 2.000000e+2	7.126037e-26
 3.000000e+2	2.295415e-26
 4.000000e+2	9.661240e-27
 5.000000e+2	4.963988e-27
 7.000000e+2	1.890927e-27
 1.000000e+3	6.605792e-28
-----------------------------

EXCITATION
He <-> He(6S1)
 2.419100e+1  1.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(6S1), Excitation
PARAM.:  E = 24.191 eV, g1/g0 = 1, complete set
COMMENT: [e + He(1S1) -> e + He(6S1), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.419100e+1	0.000000e+0
 2.430000e+1	3.250448e-24
 2.450000e+1	3.311525e-24
 2.500000e+1	3.456507e-24
 2.600000e+1	3.677344e-24
 2.700000e+1	3.855254e-24
 2.800000e+1	4.000236e-24
 2.900000e+1	4.122687e-24
 3.000000e+1	4.199043e-24
 3.100000e+1	4.257036e-24
 3.500000e+1	4.344474e-24
 4.000000e+1	4.281864e-24
 4.500000e+1	4.123083e-24
 5.000000e+1	3.925808e-24
 6.000000e+1	3.545688e-24
 7.000000e+1	3.244326e-24
 8.000000e+1	2.978926e-24
 9.000000e+1	2.770606e-24
 1.000000e+2	2.606046e-24
 1.200000e+2	2.357461e-24
 1.500000e+2	2.100986e-24
 2.000000e+2	1.805817e-24
 3.000000e+2	1.414706e-24
 4.000000e+2	1.153888e-24
 5.000000e+2	9.704247e-25
 7.000000e+2	7.372872e-25
 1.000000e+3	5.383756e-25
-----------------------------

EXCITATION
He <-> He(6P3)
 2.420100e+1  9.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(6P3), Excitation
PARAM.:  E = 24.201 eV, g1/g0 = 9, complete set
COMMENT: [e + He(1S1) -> e + He(6P3), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.420100e+1	0.000000e+0
 2.450000e+1	2.099144e-24
 2.480000e+1	2.302462e-24
 2.500000e+1	2.430182e-24
 2.550000e+1	2.692404e-24
 2.600000e+1	2.936888e-24
 2.700000e+1	3.282998e-24
 2.800000e+1	3.525470e-24
 2.900000e+1	3.677743e-24
 3.000000e+1	3.778852e-24
 3.500000e+1	3.654496e-24
 4.000000e+1	3.239619e-24
 4.500000e+1	2.772269e-24
 5.000000e+1	2.344518e-24
 6.000000e+1	1.678674e-24
 7.000000e+1	1.213976e-24
 8.000000e+1	8.941755e-25
 9.000000e+1	6.715967e-25
 1.000000e+2	5.124983e-25
 1.200000e+2	3.163195e-25
 1.500000e+2	1.679634e-25
 2.000000e+2	7.113731e-26
 3.000000e+2	2.038567e-26
 4.000000e+2	7.797844e-27
 5.000000e+2	3.739195e-27
 7.000000e+2	1.307071e-27
 1.000000e+3	4.225447e-28
-----------------------------

EXCITATION
He <-> He(6D3)
 2.420900e+1  1.500000e+1
SPECIES: e / He
PROCESS: E + He <-> E + He(6D3), Excitation
PARAM.:  E = 24.209 eV, g1/g0 = 15, complete set
COMMENT: [e + He(1S1) -> e + He(6D3), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.420900e+1	0.000000e+0
 2.450000e+1	1.223293e-24
 2.480000e+1	1.429324e-24
 2.500000e+1	1.530373e-24
 2.550000e+1	1.704146e-24
 2.600000e+1	1.826159e-24
 2.700000e+1	1.963071e-24
 2.800000e+1	2.032540e-24
 2.900000e+1	2.046518e-24
 3.000000e+1	2.024314e-24
 3.500000e+1	1.664938e-24
 4.000000e+1	1.279374e-24
 4.500000e+1	9.638155e-25
 5.000000e+1	7.284415e-25
 6.000000e+1	4.359608e-25
 7.000000e+1	2.717516e-25
 8.000000e+1	1.777178e-25
 9.000000e+1	1.211699e-25
 1.000000e+2	8.553700e-26
 1.200000e+2	4.720058e-26
 1.500000e+2	2.240368e-26
 2.000000e+2	8.542017e-27
 3.000000e+2	2.225222e-27
 4.000000e+2	8.062364e-28
 5.000000e+2	3.745385e-28
 7.000000e+2	1.266114e-28
 1.000000e+3	3.992229e-29
-----------------------------

EXCITATION
He <-> He(6D1)
 2.420900e+1  5.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(6D1), Excitation
PARAM.:  E = 24.209 eV, g1/g0 = 5, complete set
COMMENT: [e + He(1S1) -> e + He(6D1), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.420900e+1	0.000000e+0
 2.430000e+1	2.346812e-24
 2.450000e+1	2.370274e-24
 2.470000e+1	2.395452e-24
 2.500000e+1	2.434456e-24
 2.600000e+1	2.588482e-24
 2.700000e+1	2.752576e-24
 2.800000e+1	2.923554e-24
 2.900000e+1	3.086606e-24
 3.000000e+1	3.231866e-24
 3.500000e+1	3.753478e-24
 4.000000e+1	3.997279e-24
 4.500000e+1	4.029630e-24
 5.000000e+1	3.958071e-24
 6.000000e+1	3.626779e-24
 7.000000e+1	3.278155e-24
 8.000000e+1	2.949538e-24
 9.000000e+1	2.658359e-24
 1.000000e+2	2.407044e-24
 1.200000e+2	2.010544e-24
 1.500000e+2	1.592742e-24
 2.000000e+2	1.170128e-24
 3.000000e+2	7.528051e-25
 4.000000e+2	5.485842e-25
 5.000000e+2	4.317168e-25
 7.000000e+2	3.020275e-25
 1.000000e+3	2.077279e-25
-----------------------------

EXCITATION
He <-> He(6FGH3)
 2.421000e+1  8.100000e+1
SPECIES: e / He
PROCESS: E + He <-> E + He(6FGH3), Excitation
PARAM.:  E = 24.21 eV, g1/g0 = 81, complete set
COMMENT: [e + He(1S1) -> e + He(6FGH3), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.421000e+1	0.000000e+0
 2.450000e+1	8.423896e-25
 2.480000e+1	6.780260e-25
 2.500000e+1	5.982831e-25
 2.550000e+1	4.767602e-25
 2.600000e+1	4.012283e-25
 2.700000e+1	3.119473e-25
 2.800000e+1	2.514434e-25
 2.900000e+1	2.100317e-25
 3.000000e+1	1.789402e-25
 3.500000e+1	9.671741e-26
 4.000000e+1	5.697392e-26
 4.500000e+1	3.617117e-26
 5.000000e+1	2.416316e-26
 6.000000e+1	1.236661e-26
 7.000000e+1	6.861205e-27
 8.000000e+1	4.117713e-27
 9.000000e+1	2.617355e-27
 1.000000e+2	1.742711e-27
 1.200000e+2	8.819276e-28
 1.500000e+2	3.790995e-28
 2.000000e+2	1.264083e-28
 3.000000e+2	2.729326e-29
 4.000000e+2	8.346665e-30
 5.000000e+2	3.373972e-30
 7.000000e+2	9.313431e-31
 1.000000e+3	2.279450e-31
-----------------------------

EXCITATION
He <-> He(6FGH1)
 2.421000e+1  2.700000e+1
SPECIES: e / He
PROCESS: E + He <-> E + He(6FGH1), Excitation
PARAM.:  E = 24.21 eV, g1/g0 = 27, complete set
COMMENT: [e + He(1S1) -> e + He(6FGH1), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.421000e+1	0.000000e+0
 2.430000e+1	5.177706e-25
 2.450000e+1	5.087868e-25
 2.470000e+1	5.010399e-25
 2.500000e+1	4.913378e-25
 2.600000e+1	4.580376e-25
 2.700000e+1	4.274818e-25
 2.800000e+1	3.997731e-25
 2.900000e+1	3.749116e-25
 3.000000e+1	3.516710e-25
 3.500000e+1	2.666725e-25
 4.000000e+1	2.089101e-25
 4.500000e+1	1.687949e-25
 5.000000e+1	1.394253e-25
 6.000000e+1	1.007074e-25
 7.000000e+1	7.609012e-26
 8.000000e+1	5.963966e-26
 9.000000e+1	4.814991e-26
 1.000000e+2	3.978521e-26
 1.200000e+2	2.864901e-26
 1.500000e+2	1.917007e-26
 2.000000e+2	1.161684e-26
 3.000000e+2	6.198474e-27
 4.000000e+2	4.175509e-27
 5.000000e+2	3.174103e-27
 7.000000e+2	2.217491e-27
 1.000000e+3	1.544674e-27
-----------------------------

EXCITATION
He <-> He(6P1)
 2.421100e+1  3.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(6P1), Excitation
PARAM.:  E = 24.211 eV, g1/g0 = 3, complete set
COMMENT: [e + He(1S1) -> e + He(6P1), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.421100e+1	0.000000e+0
 2.450000e+1	1.748864e-24
 2.500000e+1	3.695667e-24
 2.600000e+1	6.509421e-24
 2.700000e+1	8.653313e-24
 2.800000e+1	1.034823e-23
 3.000000e+1	1.304977e-23
 4.000000e+1	2.616227e-23
 5.000000e+1	3.768984e-23
 6.000000e+1	4.562313e-23
 7.000000e+1	5.010418e-23
 8.000000e+1	5.266970e-23
 9.000000e+1	5.401559e-23
 1.000000e+2	5.422222e-23
 1.200000e+2	5.291293e-23
 1.500000e+2	5.006476e-23
 2.000000e+2	4.458436e-23
 3.000000e+2	3.605491e-23
 4.000000e+2	3.008880e-23
 5.000000e+2	2.607058e-23
 7.000000e+2	2.071715e-23
 1.000000e+3	1.603180e-23
-----------------------------

EXCITATION
He <-> He(7S3)
 2.428500e+1  3.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(7S3), Excitation
PARAM.:  E = 24.285 eV, g1/g0 = 3, complete set
COMMENT: [e + He(1S1) -> e + He(7S3), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.428500e+1	0.000000e+0
 2.450000e+1	5.145552e-24
 2.480000e+1	4.789698e-24
 2.500000e+1	4.687615e-24
 2.550000e+1	4.557541e-24
 2.600000e+1	4.514690e-24
 2.700000e+1	4.423520e-24
 2.800000e+1	4.339628e-24
 2.900000e+1	4.212732e-24
 3.000000e+1	4.083105e-24
 3.500000e+1	3.328101e-24
 4.000000e+1	2.638784e-24
 4.500000e+1	2.092974e-24
 5.000000e+1	1.668160e-24
 6.000000e+1	1.105239e-24
 7.000000e+1	7.561672e-25
 8.000000e+1	5.380412e-25
 9.000000e+1	3.954302e-25
 1.000000e+2	2.992388e-25
 1.200000e+2	1.851590e-25
 1.500000e+2	1.010328e-25
 2.000000e+2	4.580804e-26
 3.000000e+2	1.484166e-26
 4.000000e+2	6.215947e-27
 5.000000e+2	3.186479e-27
 7.000000e+2	1.219747e-27
 1.000000e+3	4.265773e-28
-----------------------------

EXCITATION
He <-> He(7S1)
 2.429800e+1  1.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(7S1), Excitation
PARAM.:  E = 24.298 eV, g1/g0 = 1, complete set
COMMENT: [e + He(1S1) -> e + He(7S1), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.429800e+1	0.000000e+0
 2.430000e+1	2.016868e-24
 2.450000e+1	2.065691e-24
 2.500000e+1	2.156589e-24
 2.600000e+1	2.301127e-24
 2.700000e+1	2.416943e-24
 2.800000e+1	2.507842e-24
 2.900000e+1	2.589215e-24
 3.000000e+1	2.639471e-24
 3.100000e+1	2.675830e-24
 3.500000e+1	2.734471e-24
 4.000000e+1	2.698381e-24
 4.500000e+1	2.600803e-24
 5.000000e+1	2.477858e-24
 6.000000e+1	2.239125e-24
 7.000000e+1	2.048703e-24
 8.000000e+1	1.881731e-24
 9.000000e+1	1.749676e-24
 1.000000e+2	1.645468e-24
 1.200000e+2	1.488586e-24
 1.500000e+2	1.326488e-24
 2.000000e+2	1.140375e-24
 3.000000e+2	8.940776e-25
 4.000000e+2	7.294452e-25
 5.000000e+2	6.135655e-25
 7.000000e+2	4.664927e-25
 1.000000e+3	3.408044e-25
-----------------------------

EXCITATION
He <-> He(7P3)
 2.430400e+1  9.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(7P3), Excitation
PARAM.:  E = 24.304 eV, g1/g0 = 9, complete set
COMMENT: [e + He(1S1) -> e + He(7P3), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.430400e+1	0.000000e+0
 2.450000e+1	1.276359e-24
 2.480000e+1	1.407414e-24
 2.500000e+1	1.487769e-24
 2.550000e+1	1.659950e-24
 2.600000e+1	1.815631e-24
 2.700000e+1	2.046093e-24
 2.800000e+1	2.207106e-24
 2.900000e+1	2.306578e-24
 3.000000e+1	2.374108e-24
 3.500000e+1	2.304953e-24
 4.000000e+1	2.049778e-24
 4.500000e+1	1.757097e-24
 5.000000e+1	1.487587e-24
 6.000000e+1	1.067471e-24
 7.000000e+1	7.726858e-25
 8.000000e+1	5.694917e-25
 9.000000e+1	4.279078e-25
 1.000000e+2	3.266903e-25
 1.200000e+2	2.021545e-25
 1.500000e+2	1.075973e-25
 2.000000e+2	4.573862e-26
 3.000000e+2	1.320045e-26
 4.000000e+2	5.021779e-27
 5.000000e+2	2.401193e-27
 7.000000e+2	8.435121e-28
 1.000000e+3	2.728786e-28
-----------------------------

EXCITATION
He <-> He(7D3)
 2.431000e+1  1.500000e+1
SPECIES: e / He
PROCESS: E + He <-> E + He(7D3), Excitation
PARAM.:  E = 24.31 eV, g1/g0 = 15, complete set
COMMENT: [e + He(1S1) -> e + He(7D3), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.431000e+1	0.000000e+0
 2.450000e+1	6.964449e-25
 2.480000e+1	8.667976e-25
 2.500000e+1	9.311690e-25
 2.550000e+1	1.052883e-24
 2.600000e+1	1.133401e-24
 2.700000e+1	1.229037e-24
 2.800000e+1	1.277626e-24
 2.900000e+1	1.289431e-24
 3.000000e+1	1.277979e-24
 3.500000e+1	1.055444e-24
 4.000000e+1	8.137569e-25
 4.500000e+1	6.140651e-25
 5.000000e+1	4.645036e-25
 6.000000e+1	2.788733e-25
 7.000000e+1	1.738616e-25
 8.000000e+1	1.136972e-25
 9.000000e+1	7.750921e-26
 1.000000e+2	5.470954e-26
 1.200000e+2	3.028035e-26
 1.500000e+2	1.440069e-26
 2.000000e+2	5.508745e-27
 3.000000e+2	1.444881e-27
 4.000000e+2	5.198177e-28
 5.000000e+2	2.406874e-28
 7.000000e+2	8.176705e-29
 1.000000e+3	2.579169e-29
-----------------------------

EXCITATION
He <-> He(7D1)
 2.431000e+1  5.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(7D1), Excitation
PARAM.:  E = 24.31 eV, g1/g0 = 5, complete set
COMMENT: [e + He(1S1) -> e + He(7D1), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.431000e+1	0.000000e+0
 2.450000e+1	1.485063e-24
 2.470000e+1	1.500371e-24
 2.500000e+1	1.524400e-24
 2.600000e+1	1.618997e-24
 2.700000e+1	1.721702e-24
 2.800000e+1	1.828089e-24
 2.900000e+1	1.932730e-24
 3.000000e+1	2.023826e-24
 3.500000e+1	2.354661e-24
 4.000000e+1	2.512981e-24
 4.500000e+1	2.537609e-24
 5.000000e+1	2.494827e-24
 6.000000e+1	2.289392e-24
 7.000000e+1	2.070762e-24
 8.000000e+1	1.864246e-24
 9.000000e+1	1.680773e-24
 1.000000e+2	1.522191e-24
 1.200000e+2	1.272181e-24
 1.500000e+2	1.008292e-24
 2.000000e+2	7.411303e-25
 3.000000e+2	4.772165e-25
 4.000000e+2	3.474708e-25
 5.000000e+2	2.733125e-25
 7.000000e+2	1.913420e-25
 1.000000e+3	1.315970e-25
-----------------------------

EXCITATION
He <-> He(7FGHI3)
 2.431000e+1  1.200000e+2
SPECIES: e / He
PROCESS: E + He <-> E + He(7FGHI3), Excitation
PARAM.:  E = 24.31 eV, g1/g0 = 120, complete set
COMMENT: [e + He(1S1) -> e + He(7FGHI3), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.431000e+1	0.000000e+0
 2.450000e+1	5.811662e-25
 2.480000e+1	4.531068e-25
 2.500000e+1	4.020934e-25
 2.550000e+1	3.136064e-25
 2.600000e+1	2.628429e-25
 2.700000e+1	2.014840e-25
 2.800000e+1	1.617569e-25
 2.900000e+1	1.348401e-25
 3.000000e+1	1.147624e-25
 3.500000e+1	6.227096e-26
 4.000000e+1	3.658300e-26
 4.500000e+1	2.320504e-26
 5.000000e+1	1.549188e-26
 6.000000e+1	7.954782e-27
 7.000000e+1	4.407678e-27
 8.000000e+1	2.643048e-27
 9.000000e+1	1.679199e-27
 1.000000e+2	1.117543e-27
 1.200000e+2	5.675742e-28
 1.500000e+2	2.444902e-28
 2.000000e+2	8.190084e-29
 3.000000e+2	1.784682e-29
 4.000000e+2	5.409126e-30
 5.000000e+2	2.178147e-30
 7.000000e+2	6.058873e-31
 1.000000e+3	1.485650e-31
-----------------------------

EXCITATION
He <-> He(7FGHI1)
 2.431000e+1  4.000000e+1
SPECIES: e / He
PROCESS: E + He <-> E + He(7FGHI1), Excitation
PARAM.:  E = 24.31 eV, g1/g0 = 40, complete set
COMMENT: [e + He(1S1) -> e + He(7FGHI1), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.431000e+1	0.000000e+0
 2.450000e+1	3.238243e-25
 2.470000e+1	3.180021e-25
 2.500000e+1	3.114032e-25
 2.600000e+1	2.905958e-25
 2.700000e+1	2.712332e-25
 2.800000e+1	2.536486e-25
 2.900000e+1	2.378419e-25
 3.000000e+1	2.232666e-25
 3.500000e+1	1.693982e-25
 4.000000e+1	1.326783e-25
 4.500000e+1	1.071702e-25
 5.000000e+1	8.851630e-26
 6.000000e+1	6.399268e-26
 7.000000e+1	4.833297e-26
 8.000000e+1	3.787642e-26
 9.000000e+1	3.057241e-26
 1.000000e+2	2.525743e-26
 1.200000e+2	1.820662e-26
 1.500000e+2	1.218763e-26
 2.000000e+2	7.389012e-27
 3.000000e+2	3.942695e-27
 4.000000e+2	2.648439e-27
 5.000000e+2	2.010893e-27
 7.000000e+2	1.404682e-27
 1.000000e+3	9.782579e-28
-----------------------------

EXCITATION
He <-> He(7P1)
 2.431100e+1  3.000000e+0
SPECIES: e / He
PROCESS: E + He <-> E + He(7P1), Excitation
PARAM.:  E = 24.311 eV, g1/g0 = 3, complete set
COMMENT: [e + He(1S1) -> e + He(7P1), Excitation]  Ralchenko Yu V Janev R K Kato T Fursa D V
COMMENT: Bray I and de Heer F J 2008 At. Data Nucl. Data Tables 94 603.
UPDATED: 2017-05-12 19:55:37
COLUMNS: Energy (eV) | Cross section (m2)
-----------------------------
 2.431100e+1	0.000000e+0
 2.450000e+1	8.305641e-25
 2.500000e+1	2.080933e-24
 2.600000e+1	3.934850e-24
 2.700000e+1	5.321293e-24
 2.800000e+1	6.404980e-24
 3.000000e+1	8.118584e-24
 4.000000e+1	1.633929e-23
 5.000000e+1	2.358719e-23
 6.000000e+1	2.861319e-23
 7.000000e+1	3.147873e-23
 8.000000e+1	3.311992e-23
 9.000000e+1	3.398856e-23
 1.000000e+2	3.414577e-23
 1.200000e+2	3.334384e-23
 1.500000e+2	3.156529e-23
 2.000000e+2	2.813365e-23
 3.000000e+2	2.277061e-23
 4.000000e+2	1.900826e-23
 5.000000e+2	1.646767e-23
 7.000000e+2	1.309379e-23
 1.000000e+3	1.013496e-23
-----------------------------
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx