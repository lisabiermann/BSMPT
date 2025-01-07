#!/usr/bin/env python3

import pandas as pd
import sys

Seperator = '\t'
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

mHl = 'mHl'
mHh = 'mHh'
mA = 'mA'
mHp = 'mHp'
alpha = 'alpha'
v = 'v'

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

def convert(InputFILE, OutputFILE, frontcol):
	df = pd.DataFrame()
	with open(InputFILE, 'r') as file:
		df = pd.read_csv(file, index_col=False, sep=Seperator)

	Col = [c for c in frontcol if c in df]
	ColOtherInF = [c for c in df if c not in frontcol]
	ColNotInF = [c for c in frontcol if c not in df]
	df = df[Col + ColOtherInF]

	for column in ColNotInF:
		df[column] = 0

	with open(OutputFILE, 'w') as file:
		df.to_csv(file, index=True, sep='\t')

if __name__ == "__main__":
	params_dim4 = [mHl, mHh, mA, mHp, TanBeta, m12Sq, alpha, L1, L2, L3, L4, L5, m11Sq, m22Sq, Type, v]
	params_dim6_scal = [O111111, O111122, O122111, O121211, O111112, O121221, O112212, O222222, O112222, O122122, O121222, O122222, O121212]
	nums = ['111', '112', '121', '122', '211', '212', '221', '222']
	params_dim6_up = ['OQU'+num for num in nums]
	params_dim6_down = ['OQD'+num for num in nums]
	params_dim6_lep = ['OLE'+num for num in nums]
	params_dim6_ferm = params_dim6_up + params_dim6_down + params_dim6_lep

	frontcol = params_dim4 + params_dim6_scal + params_dim6_ferm

	convert(sys.argv[1], sys.argv[2], frontcol)
