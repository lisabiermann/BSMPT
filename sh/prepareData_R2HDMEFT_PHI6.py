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

# EFT-operators
O111111 = 'Op6_111111'
O111122 = 'Op6_111122'
O122111 = 'Op6_122111'
O121211 = 'Op6_121211'
O111112 = 'Op6_111112'
O121221 = 'Op6_121221'
O112212 = 'Op6_112212'
O222222 = 'Op6_222222'
O112222 = 'Op6_112222'
O122122 = 'Op6_122122'
O121222 = 'Op6_121222'
O122222 = 'Op6_122222'
O121212 = 'Op6_121212'

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
		df.to_csv(file, index=False, sep='\t')

if __name__ == "__main__":
	params_dim4 = [Type, L1, L2, L3, L4, L5, m11Sq, m22Sq, m12Sq, TanBeta]
	params_dim6_scal = [O111111, O111122, O122111, O121211, O111112, O121221, O112212, O222222, O112222, O122122, O121222, O122222, O121212]
	frontcol = params_dim4 + params_dim6_scal

	convert(sys.argv[1], sys.argv[2], frontcol)
