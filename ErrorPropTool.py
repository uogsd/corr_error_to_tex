import sympy as sp
from IPython.display import display, Math, Latex
sp.init_printing(use_unicode=True)
import numpy as np

def ptex(output): # print tex
    print(sp.latex(output))
def dtex(output): # display tex
    display(Math(sp.latex(output)))

class ucfunc:
    # Variable explanation:
    
    # "expression" is the function definition in python format:
    # 
    #       Example a*x**2+c+sp.sin(x)
    #
    # "constants" are quantities that have no or negligible uncertainties
    # "dependencies" x_i ... are function arguments f({x_i}) with uncertainties
    # The index "_N" mean numerical and is used for variables that store numerical values for later replacement
    
    # This class uses the formula
    
    #########################################################################
    ### var = g^T V g     with    g_i = (\partial f) / (\partial x_i)     ###
    #########################################################################
    
    # for error propagation
    
    # Use class for example like:
    
    #   > VDiag_N = [[a_uc**2, 0],[0, v_S1_uc**2]]
    #  >  res = ucfunc(DS1, ['a', 'v_S1'], [a, v_S1], ['m'], [m], VDiag_N)       # function from class
    
     
    def __init__(self, expression, dependencies, dependencies_N, constants, constants_N, V_N):
        # define expression and dependencies of function
        self.dependencies_N = dependencies_N
        self.constants_N = constants_N
        self.dim_V = len(dependencies)
        
        # generate symbolic covariance matrix with symbols
        self.V = sp.eye(self.dim_V) 
        self.V_symbols_nonZeroVals = []
        self.V_N_nonZeroVals = []
        for i in range(0, self.dim_V):
            for j in range(0, self.dim_V):
                if V_N[i][j] != 0: 
                    varname = str(f"var{i+1}{j+1}")
                    exec(str(f"tmp = sp.symbols('{varname}', real=True)"))
                    exec(str(f"self.V_symbols_nonZeroVals.append(tmp)")) # for replacement list (to numerical values)
                    exec(str(f"self.V[{i},{j}] = tmp")) # for calcilation
                    self.V_N_nonZeroVals.append(V_N[i][j]) # store numerical values for replacement 
        
        # define symbols for variables
        self.dependencies = []
        for variable in dependencies:
            exec(str(f"{variable} = sp.symbols('{variable}', real=True)"))
            exec(str(f"self.dependencies.append({variable})"))
            
        # define symbols for constants
        self.constants = []
        for variable in constants:
            exec(str(f"{variable} = sp.symbols('{variable}', real=True)"))
            exec(str(f"self.constants.append({variable})"))
        
        # define equation using, the symbol definitions from above are being used now!
        self.expression = eval(expression)
        
        # define g (see formula above)
        g = []
        for variable in self.dependencies:
            g.append(sp.diff(self.expression, variable))
        self.g = sp.Matrix(g)
        
        # calculate variance (see formula above)
        self.var = ((self.g).T*(self.V)*(self.g))[0]
        
        # use result for uncertainty sigma
        self.sigma = sp.sqrt(self.var)
        
        
        
    # Methods for access
    
    def tex(self, x):   # see x in TEX format
        return sp.latex(x)
    
    
    def get_numerical(self, x):            # calculate numerical value of expression part x
        x_N = x
        for index, variable in enumerate(self.dependencies):
            x_N = x_N.subs(variable, self.dependencies_N[index])
        for index, variable in enumerate(self.constants):
            x_N = x_N.subs(variable, self.constants_N[index])
            
        for index, variable in enumerate(self.V_symbols_nonZeroVals):
            x_N = x_N.subs(variable, self.V_N_nonZeroVals[index]) # replace symbol with numerical value
    
        return x_N
    
    
    def var_N(self):
        return self.get_numerical(self.var)
    
    def sigma_N(self):
        return sp.sqrt(self.get_numerical(self.var))
    
    def N(self):
        return str(f"{self.get_numerical(self.expression)} \u00B1 {self.sigma_N()} ")
        
    def formula(self):
        return str(f"{self.tex(self.expression)} \u00B1 {self.tex(self.sigma)}")
    
    def all(self):
        print("Using")
        display(Math("\sigma = g^T V g"))
        print("with")
        display(Math("g_i = \partial_i f(\{x_i\})"))
        print("we obtain")
        display(Math(self.formula()))
        print("\nResult in LaTex - Form:\n")
        print(self.formula())
        print("\n\nNumerical value of result\n")
        print(self.N())