--------------------------------------------------------------------------------
# SSVOS_v1_0
--------------------------------------------------------------------------------
Code written by Won-Dong Jang

Contact: Won-Dong Jang, wdjang@mcl.korea.ac.kr

If you want to use this software, please cite:

W.-D. Jang and C.-S. Kim.

Semi-supervised Video Object Segmentation Using Multiple Random Walkers.

In Proceedings of British Machine Vision Conference (BMVC), Sep., (2016)

@inproceedings{jang2016semi-supervised,

 title={Semi-supervised Video Object Segmentation Using Multiple Random Walkers},
 
 author={W.-D. Jang and C.-S. Kim},
 
 booktitle={BMVC},
 
 year={2016}}
 

--------------------------------------------------------------------------------
Quick start
--------------------------------------------------------------------------------
This software was developed under 64-bit Windows with Matlab R2015b. 

There is no guarantee it will run on other operating systems or Matlab versions.

1) CVX optimizer is required. It can be downloaded at: http://cvxr.com/cvx/
⋅⋅* We use the Gurobi optimizer to minimize energy functions. In this case, you need to install Gurobi optimizer additionally.
⋅⋅* The Gurobi optimizer can be replaced by the 'SeDuMi' or 'SDPT3,' which are supported by the default CVX.

2) Please see 'demo.m'

--------------------------------------------------------------------------------
LICENSE
--------------------------------------------------------------------------------
This program is released with a research only license.
