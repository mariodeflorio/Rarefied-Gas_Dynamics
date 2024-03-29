# Rarefied-Gas Dynamics
The Boltzmann equation with the Bhatnagar-Gross-Krook (BGK) collision model has been widely employed to describe the evolution of a gas through a kinetic theory. We propose a new Machine Learning approach to accurately and efficiently learn the solution of a class of problems in the transport theory of rarefied-gas dynamics, employing a particular Physics-Informed Neural Networks method called Extreme Theory of Functional Connections (X-TFC). According to X-TFC, the unknown function is approximated by the Constrained Expressions defined int the Theory of Functional Connections. Constrained Expressions are functionals (e.g., function of functions), which are the sum of a free function and a functional in which the constraints are analytically embedded. The constraints will therefore always be satisfied analytically, no matter what the free-chosen function is. In this work, a single-layer Neural Network trained with Extreme Learning Machine algorithm is employed. The proposed Machine Learning approach learns the solution of the Linear Boundary Value Problem arising from the Thermal Creep, Poiseuille, and Couette flow problems in the BGK approximation between two parallel plates of a channel, for a wide range of Knudsen numbers. The accuracy of X-TFC method is tested and compared with published benchmarks in literature.

For more information, please refer to the following: <br>
(https://github.com/mariodeflorio/Rarefied-Gas_Dynamics/)

<ul>
<li>De Florio, Mario, Enrico Schiassi, Barry D. Ganapol, and Roberto Furfaro. "<a href="https://doi.org/10.1063/5.0046181">Physics-informed neural networks for rarefied-gas dynamics: Thermal creep flow in the Bhatnagar–Gross–Krook approximation</a>." Physics of Fluids 33, no. 4 (2021): 047110.</li>

<li>De Florio, Mario, Enrico Schiassi, Barry D. Ganapol, and Roberto Furfaro. "<a href="https://doi.org/10.1007/s00033-022-01767-z">Physics-Informed Neural Networks for rarefied-gas dynamics: Poiseuille flow in the BGK approximation.</a>." Zeitschrift für angewandte Mathematik und Physik 73, no. 3 (2022): 1-18.</li>
</ul>


## Citation

    @article{de2021physics,
      title={Physics-informed neural networks for rarefied-gas dynamics: Thermal creep flow in the Bhatnagar--Gross--Krook approximation},
      author={De Florio, Mario and Schiassi, Enrico and Ganapol, Barry D and Furfaro, Roberto},
      journal={Physics of Fluids},
      volume={33},
      number={4},
      pages={047110},
      year={2021},
      publisher={AIP Publishing LLC}
    }
    
    @article{de2022physics,
      title={Physics-Informed Neural Networks for rarefied-gas dynamics: Poiseuille flow in the BGK approximation},
      author={De Florio, Mario and Schiassi, Enrico and Ganapol, Barry D and Furfaro, Roberto},
      journal={Zeitschrift f{\"u}r angewandte Mathematik und Physik},
      volume={73},
      number={3},
      pages={1--18},
      year={2022},
      publisher={Springer}
    }
    
    
    
    
    
    
    
    
