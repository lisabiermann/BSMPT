# SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas Müller
#
# SPDX-License-Identifier: GPL-3.0-or-later

import pandas as pd
import argparse

####### The parameters Type,Lambda1 to Lambda8, tanbeta, m12squared and v_s should have the label of the
####### corresponding parameter. With Seperator you have to tell which seperator your data file
####### is using (e.g. , \t or space). Your InputFILE will then be saved to OutputFILE.

Seperator='\t'
InputFILE='../example/CN2HDM_Input.dat'
OutputFILE='CN2HDM_Ordered.dat'
Type='yuktype'
Lambda1='L1'
Lambda2='L2'
Lambda3='L3'
Lambda4='L4'
ReLambda5='re_L5'
ImLambda5='im_L5'
Lambda6='L6'
Lambda7='L7'
Lambda8='L8'
Rem12Sq='re_m12sq'
mDM='mDM'
TanBeta='tbeta'
vs='p_vs'


def convert(IndexCol):
	df=pd.DataFrame()
	if IndexCol == 'False':
		df=pd.read_table(InputFILE,index_col=False,sep=Seperator)
	else:
		df=pd.read_table(InputFILE,index_col=int(IndexCol),sep=Seperator)

	frontcol=[Type,Lambda1,Lambda2,Lambda3,Lambda4,ReLambda5,ImLambda5,Lambda6,Lambda7,Lambda8,Rem12Sq,mDM,TanBeta,vs]

	Col = [c for c in frontcol if c in df] + [c for c in df if c not in frontcol]
	df=df[Col]
	df = df.rename(columns={'mDM': 'mDMSq'})
	df['mDMSq'] = df['mDMSq']**2 # loaded value is just mDM

	df.to_csv(OutputFILE,index=False,sep='\t')

parser = argparse.ArgumentParser()
parser.add_argument('-i','--indexcol',help='Column which stores the index of your data', default='False')

if __name__ == "__main__":
	args = parser.parse_args()
	convert(args.indexcol)

