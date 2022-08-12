# SPDX-FileCopyrightText: 2022 Philipp Basler, Lisa Biermann,
# Margarete Mühlleitner and Jonas Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

import pandas as pd
import numpy as np
import os

'''
The parameters Type, Lambda1 to Lambda5, tanbeta and m12squared should have the
label of the corresponding parameter. If your first column is an index column
you should set HasIndexCol = True otherwise set HasIndexCol = False . With
Seperator you have to tell which seperator your data file is using
(e.g. , \t or space). Your InputFILE will then be saved to OutputFILE
'''

HasIndexCol = True
Seperator = '\t'
InputFileDIR = '/itp/swift/lisab/fast/Data/R2HDM_SAMPLES/Duarte_BSMPT/'
InputFileNAME = 'SampleR2HDM_ScannerS_T2_noSFOEWPT_above9em1.csv'
#InputFileDIR = '../example/'
#InputFileNAME = 'R2HDM_EFT_Input.dat'
OutputFILE = '../example/Input_R2HDMEFT_test.csv'
Type = 'yuktype'
L1 = 'L1'
L2 = 'L2'
L3 = 'L3'
L4 = 'L4'
L5 = 'L5'
TanBeta = 'tbeta'
m12Sq = 'm12sq'
m11Sq = 'm11sq'
m22Sq = 'm22sq'

# EFT-operators
O111111 = 'O111111'
O111122 = 'O111122'
O122111 = 'O122111'
O121211 = 'O121211'
O111112 = 'O111112'
O121221 = 'O121221'
O112212 = 'O112212'
O222222 = 'O222222'
O112222 = 'O112222'
O122122 = 'O122122'
O121222 = 'O121222'
O122222 = 'O122222'
O121212 = 'O121212'

OQU111 = 'OQU111'
OQU112 = 'OQU112'
OQU121 = 'OQU121'
OQU122 = 'OQU122'
OQU211 = 'OQU211'
OQU212 = 'OQU212'
OQU221 = 'OQU221'
OQU222 = 'OQU222'

OQD111 = 'OQD111'
OQD112 = 'OQD112'
OQD121 = 'OQD121'
OQD122 = 'OQD122'
OQD211 = 'OQD211'
OQD212 = 'OQD212'
OQD221 = 'OQD221'
OQD222 = 'OQD222'

OLE111 = 'OLE111'
OLE112 = 'OLE112'
OLE121 = 'OLE121'
OLE122 = 'OLE122'
OLE211 = 'OLE211'
OLE212 = 'OLE212'
OLE221 = 'OLE221'
OLE222 = 'OLE222'

params_dim4 = [Type, L1, L2, L3, L4, L5, m11Sq, m22Sq, m12Sq, TanBeta]
params_dim6_scal = [O111111, O111122, O122111, O121211, O111112, O121221, O112212, O222222, O112222, O122122, O121222, O122222, O121212]
nums = ['111', '112', '121', '122', '211', '212', '221', '222']
params_dim6_up = ['OQU'+num for num in nums]
params_dim6_down = ['OQD'+num for num in nums]
params_dim6_lep = ['OLE'+num for num in nums]
params_dim6_ferm = params_dim6_up + params_dim6_down + params_dim6_lep

InputFILE = os.path.join(InputFileDIR, InputFileNAME)

df = pd.DataFrame()
if HasIndexCol:
    df = pd.read_table(InputFILE, sep=Seperator)
else:
    df = pd.read_table(InputFILE, index_col=False, sep=Seperator)

frontcol = params_dim4 + params_dim6_scal + params_dim6_ferm

Col = [c for c in frontcol if c in df]
ColOtherInF = [c for c in df if c not in frontcol]
ColNotInF = [c for c in frontcol if c not in df]
df = df[Col]
for column in ColNotInF:
    df[column] = 0

df.to_csv(OutputFILE, index=False, sep='\t')
