'''
Dixon provides enzyme kinetics with three different modes,
including 1)competitive inhibition, 2) non-competitive inhibition,
and 3) anti-competitive inhibition.
All of them requires arguments x=(Substrates,Inhibitors), Vmax, Km, Ki, 
alpha(for Vmax), and beta(for Km, Ki).
gives return initial rates as array.

This program was written by Kyoichiro Higashi PHD. in Meiji Pharmaceutical University.
2021/01/06
'''
import numpy as np

def make_exog(S,I):
	'''
	arguments:
		S:Substrate in list
		I:Inhibitor(or activator) in list
	return:
		array formatted Substrates-series and inhibitor series
	'''
	return np.array([S*len(I), np.array([[i]*len(S) for i in I]).flatten()])

def competitive(exog, Vmax=1, Km=1, Ki=1):
	'''
	competitive inhibition kinetics
	arguments:
	exog (array): [Substrate,Inhibitor]
	Vmax=1, Km=1, Ki=1
	returns: initial rates'''
	S, I = exog
	v = Vmax * S / (Km * (1 + I/Ki) + S) 
	return v

def noncompetitive(exog, Vmax=1, Km=1, Ki=1, alpha=0):
	''' 
	noncompetitive inhibition kinetics.
	This function treats partial inhibition as well as complete inhibition.
	arguments:
	exog (array): [Substrate,Inhibitor]
	Vmax=1, Km=1, Ki=1, alpha=0
	alpha means ratio of Vmax with maximal inhibitor to Vmax without inhibitor.
	Alpha=0 means Complete inhibition. 
	returns: initial rates'''
	S, I = exog
	v = Vmax*S/(Km+S)*Ki/(Ki+I) + alpha*Vmax*S/(Km+S)*I/(Ki+I) 
	return v

def anticompetitive(exog, Vmax=1, Km=1, Ki=1,alpha=0):
	'''
	anticompetitive (also known as uncompetitive) inhibition kinetics.
	arguments:
	exog (array): [Substrate,Inhibitor]
	Vmax=1, Km=1, Ki=1, alpha=0
	Alpha=0 means Complete inhibition. 
	returns: initial rates'''
	S, I = exog
	v = Vmax * S / (Km + S * (1 + I / Ki)) + alpha * I/Ki * Vmax * S / (Km + S * (1 + I/Ki))
	return v

def mixed(exog, Vmax=1, Km=1, Ki=1, alpha=0.3, beta=3):
	'''
	Mixed type inhibition means random-order binding of inhibitor and substrate. Alpha means the degree of decrease in Vmax, 
	beta means increase in Km.
	exog (array): [Substrate,Inhibitor]
	Vmax=1, Km=1, Ki=1,alpha=0.3,beta=3
	alpha means ratio of inhibited Vmax with maximal inhibitor to Vmax without inhibitor.
	Alpha=0.3 means partial inhibition.
	beta means ratio of Km with inhibitor to Km without inhibitor.
	returns: initial rates	'''
	S, I = exog
	v = Vmax * S / (S + Km * (1 + I/Ki * (1 + S/(Km * beta)))) + \
	Vmax * alpha * S * I / (Ki * beta) / (S + Km * (1 + I/Ki * (1 + S/(Km * beta)))) 
	return v
