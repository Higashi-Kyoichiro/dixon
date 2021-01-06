# dixon
Analysis of enzymatic inhibition and activation with non-linear regression.
# instalation
simply put this dixon.py to your site-packages directory.
# Requirements
* python3
* numpy
# comments
* dixon.py provides classical enzyme kinetic equation with competitive, non-competitive,mixed, and anti-competitive inhibition and activations.
* argments
  * array or list(array([substrate]), array([inhibitor]))
  * parameters for enzyme kinetics:Vmax,Km,Ki,alpha,beta
    * alpha means second Vmax with inhibitor.
    * beta means second Km and second Ki.
* return: array(initial velocity for each (S,I)) 

# technical resources
* competitive inhibition
$v = {Vmax*S} over {(S+ Km*(1+{I} over{Ki}))}$
