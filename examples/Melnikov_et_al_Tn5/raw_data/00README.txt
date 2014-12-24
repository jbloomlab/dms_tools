Supplementary Data for Melnikov et al. (2014)

The bundled files contain the post-processed amino acid counts and mutation frequencies on which the analysis in the
main text are based.

Note regarding selection with amikacin: We initially performed two independent rounds of selection (S1 and S2) with all
six aminoglycosides. Analysis of the aminkacin results showed no sign of negative selection at the three lowest
concentrations after using 2.5 ug/mL as the estimated WT MIC. Accordingly, we revised our estimate upwards to 10 ug/mL
and performed two new rounds of selection (S3 and S4). We have included the data from all four rounds.

File names
------------
KKA2_Bkg[12].*     Input/background libraries 1 and 2
KKA2_[S]_[D]_[L].* Post-selection libraries
                     [S] = (S)election round 1, 2, 3 or 4
                     [D] = (D)rug and dilution, eg. Kan12 refers to Kanamycin at 1:2 WT MIC
                     [L] = Input (l)ibrary 1 or 2

File formats
------------
*.aacounts.txt:   Tab-delimited text files containing the number of (high-quality) observations of each
                  amino acid at each position in one sequenced library
*.mutdiff.txt:    Tab-delimited text files containing the observed change in the frequency of non-WT amino acids
                  at each position after one selection:
                    Delta-Mut: (# of non-WT amino acids)/(# of all amino acids)
                    P(Delta-Mut != 1): BHY-corrected p-value (5% FDR) of hypothesis that Delta-Mut in this library
                                       is equal to Delta-Mut in the matched input library (chi^2 with a pseudocount of 1)
*.aadiff.txt:     Tab-delimited text files containing the observed change in the frequency of each non-WT amino acid
                  at each position after one selection:
                    Delta-[AA]: (# of [AA])/(# of all amino acids)
                    P(Delta-[AA] != 1): BHY-corrected p-value (5% FDR) of hypothesis that Delta-[AA] in this library
                                       is equal to Delta-Mut in the matched input library (chi^2 with a pseudocount of 1)

                  Note that the effective lower bounds of Delta-Mut and Delta-[AA] at each position are dependent on
                  the sequence coverage and the input amino acid frequency. Use caution when comparing the absolute values
                  of these statistics, particularly for Delta-[AA] < 0.1. Also note that the P(Delta-[AA] != 1) test is
                  somewhat under-powered due to the relatively low counts for individual substitutions and the large number
                  of tests performed.

Library summary
---------------
  Cov = Median sequence coverage
  Mut = Median observed mutation rate (# non-WT amino acids/# of observed amino acids)
  Neg = Number of positions for which Delta-Mut < 1 at 5% FDR

Library                    Cov   Mut   Neg Notes
KKA2_Bkg1               126881  0.37    NA 
KKA2_Bkg2               142182  0.36    NA 
KKA2_S1_Kan11_L1         65101  0.12   211 
KKA2_S2_Kan11_L2         42046  0.10   207 
KKA2_S1_Kan12_L1         64321  0.22   157 
KKA2_S2_Kan12_L2         25241  0.18   146 
KKA2_S1_Kan14_L1         24380  0.30   116 
KKA2_S2_Kan14_L2         31552  0.29   119 
KKA2_S1_Kan18_L1         86321  0.39   108 
KKA2_S2_Kan18_L2           520  0.34     0 Discarded - Failed sequencing library
KKA2_S1_Neo11_L1         74904  0.09   241 
KKA2_S2_Neo11_L2         42762  0.10   208 
KKA2_S1_Neo12_L1         76258  0.26   142 
KKA2_S2_Neo12_L2         26238  0.27   117 
KKA2_S1_Neo14_L1         78082  0.38   113 
KKA2_S2_Neo14_L2         30496  0.38    91 
KKA2_S1_Neo18_L1         63491  0.42    76 
KKA2_S2_Neo18_L2         29703  0.39    57 
KKA2_S1_G41811_L1        74142  0.11   214 
KKA2_S2_G41811_L2        49791  0.11   195 
KKA2_S1_G41812_L1        63035  0.29   135 
KKA2_S2_G41812_L2        38651  0.30   128 
KKA2_S1_G41814_L1       101415  0.40    99 
KKA2_S2_G41814_L2        26798  0.38    68 
KKA2_S1_G41818_L1        93465  0.37     8 No evidence of selection
KKA2_S2_G41818_L2        36457  0.35     0 No evidence of selection
KKA2_S1_Paro11_L1        47597  0.11   249 
KKA2_S2_Paro11_L2        51104  0.06   249 
KKA2_S1_Paro12_L1        60338  0.18   179 
KKA2_S2_Paro12_L2        55899  0.11   195 
KKA2_S1_Paro14_L1        74210  0.28   134 
KKA2_S2_Paro14_L2        31360  0.26   129 
KKA2_S1_Paro18_L1        83053  0.36   119 
KKA2_S2_Paro18_L2        51721  0.30   125 
KKA2_S1_Ribo11_L1        74303  0.14   196 
KKA2_S2_Ribo11_L2        42213  0.13   184 
KKA2_S1_Ribo12_L1        54583  0.30   132 
KKA2_S2_Ribo12_L2        49357  0.30   120 
KKA2_S1_Ribo14_L1        58106  0.37   117 
KKA2_S2_Ribo14_L2        46939  0.40    97 
KKA2_S1_Ribo18_L1        75475  0.41    94 
KKA2_S2_Ribo18_L2        43658  0.40    88 
KKA2_S1_Ami14_L1         27736  0.31   105 
KKA2_S2_Ami14_L2         54683  0.22   151 
KKA2_S1_Ami18_L1         23711  0.39     0 Discarded - no evidence of selection
KKA2_S2_Ami18_L2         28290  0.36     0 No evidence of selection
KKA2_S1_Ami116_L1        24861  0.38     0 Discarded - no evidence of selection
KKA2_S2_Ami116_L2        34238  0.35     0 No evidence of selection
KKA2_S1_Ami132_L1        30208  0.37     0 Discarded - no evidence of selection
KKA2_S2_Ami132_L2        45092  0.35     0 No evidence of selection
KKA2_S3_Kan11_L1         66108  0.05   210 Shown in Figure 2
KKA2_S3_Kan12_L1        123600  0.22   155 Shown in Figure 2
KKA2_S3_Kan14_L1         63783  0.30   130 Shown in Figure 2
KKA2_S3_Kan18_L1        205099  0.38   113 Shown in Figure 2
KKA2_S3_Ami11_L1         20427  0.04   258 
KKA2_S4_Ami11_L2         59975  0.07   251 
KKA2_S3_Ami12_L1        116441  0.14   183 
KKA2_S4_Ami12_L2         44410  0.30   102 
KKA2_S3_Ami14_L1        168315  0.34    33 
KKA2_S4_Ami14_L2         96649  0.33    29 
KKA2_S3_Ami18_L1        109258  0.34    43 
KKA2_S4_Ami18_L2         47813  0.33    12 
