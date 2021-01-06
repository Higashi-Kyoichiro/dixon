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
  * v = Vmax*S/(S+ Km*(1 + I/Ki))
* non-competitive inhibition
  * v = Vmax*S/(Km + S)/(1 + I/Ki) + alpha*Vmax*S/(Km + S)*(I/Ki)/(1 + I/Ki))
* mixed inhibition
  * v  = Vmax*S/(S + Km*(1 + I/Ki*(1 + S/(Km*beta)))) + Vmax*alpha*S*I/(Ki*beta)/(S + Km*(1 + I/Ki*(1 + S/(Km*beta)))) 
* anti-competitive inhibition
  * v = Vmax*S/(Km + S*(1 + I/Ki))  + alpha*Vmax*S*I/Ki/(Km + S*(1 + I/Ki))
