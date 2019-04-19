#!/usr/bin/env python

import sys 
import os

output_dir = "erc2-chromatin15state-all-files"

if not os.path.exists(output_dir):
  sys.stderr.write("Creating dir [%s]...\n" % (output_dir))
  os.makedirs(output_dir)

prefix = "/home/cbreeze/for_Alex"

suffix = "_15_coreMarks_mnemonics.bed"

marks = [ '1_TssA',
          '2_TssAFlnk',
          '3_TxFlnk',
          '4_Tx',
          '5_TxWk',
          '6_EnhG',
          '7_Enh',
          '8_ZNF/Rpts',
          '9_Het',
          '10_TssBiv',
          '11_BivFlnk',
          '12_EnhBiv',
          '13_ReprPC',
          '14_ReprPCWk',
          '15_Quies' ]

all = [ 'E001',
        'E002',
        'E003',
        'E004',
        'E005',
        'E006',
        'E007',
        'E008',
        'E009',
        'E010',
        'E011',
        'E012',
        'E013',
        'E014',
        'E015',
        'E016',
        'E017',
        'E018',
        'E019',
        'E020',
        'E021',
        'E022',
        'E023',
        'E024',
        'E025',
        'E026',
        'E027',
        'E028',
        'E029',
        'E030',
        'E031',
        'E032',
        'E033',
        'E034',
        'E035',
        'E036',
        'E037',
        'E038',
        'E039',
        'E040',
        'E041',
        'E042',
        'E043',
        'E044',
        'E045',
        'E046',
        'E047',
        'E048',
        'E049',
        'E050',
        'E051',
        'E052',
        'E053',
        'E054',
        'E055',
        'E056',
        'E057',
        'E058',
        'E059',
        'E061',
        'E062',
        'E063',
        'E065',
        'E066',
        'E067',
        'E068',
        'E069',
        'E070',
        'E071',
        'E072',
        'E073',
        'E074',
        'E075',
        'E076',
        'E077',
        'E078',
        'E079',
        'E080',
        'E081',
        'E082',
        'E083',
        'E084',
        'E085',
        'E086',
        'E087',
        'E088',
        'E089',
        'E090',
        'E091',
        'E092',
        'E093',
        'E094',
        'E095',
        'E096',
        'E097',
        'E098',
        'E099',
        'E100',
        'E101',
        'E102',
        'E103',
        'E104',
        'E105',
        'E106',
        'E107',
        'E108',
        'E109',
        'E110',
        'E111',
        'E112',
        'E113',
        'E114',
        'E115',
        'E116',
        'E117',
        'E118',
        'E119',
        'E120',
        'E121',
        'E122',
        'E123',
        'E124',
        'E125',
        'E126',
        'E127',
        'E128',
        'E129' ]
        
# prefix, suffix, marks, all

for sample in all:
  fns = {}
  fhs = {}
  # set up output file handles for all combinations of per-sample and marks
  for mark in marks:
    fns[mark] = os.path.join(output_dir, "%s_%s.bed" % (sample, mark.replace('/', '-')))
    sys.stderr.write("Setting up output handle to [%s]...\n" % (fns[mark]))
    fhs[mark] = open(fns[mark], "w")
  # split per-sample mnemonics to per-sample, per-mark file
  psm_fn = "%s/%s%s" % (prefix, sample, suffix)
  sys.stderr.write("Reading PSM [%s]...\n" % (psm_fn))
  with open(psm_fn, "r") as psm_fh:
    for line in psm_fh:
      (chr, start, stop, state_call) = line.strip().split('\t')
      fhs[state_call].write('\t'.join([chr, start, stop]) + '\n')
  # close handles
  for mark in marks:
    sys.stderr.write("Closing output handle to [%s]...\n" % (fns[mark]))
    fhs[mark].close()
    fns[mark] = None
    fhs[mark] = None