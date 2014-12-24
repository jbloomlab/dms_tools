import mapmuts.dssp

dssp = mapmuts.dssp.ReadDSSP('1ND4.dssp', 'Tien2013', 'A')

frsa = open('RSAs.txt', 'w')
frsa.write('#Relative solvent accessibilities from a DSSP analysis of the 1ND4 dimer, taking result for chain A. Absolute solvent accessibilities are normalized to relative ones using the maximum solvent accessibilities of Tien et al, PLoS One, 8:e80635.\n#SITE RSA\n')
fss = open('SSs.txt', 'w')
fss.write('#Secondary structures from a DSSP analysis of the 1ND4 dimer A\n#SITE SS\n')
sites = dssp.keys()
sites.sort()
for site in sites:
    fss.write('%d %s\n' % (site, dssp[site]['SS_CLASS']))
    frsa.write('%d %s\n' % (site, dssp[site]['RSA']))
fss.close()
frsa.close()
