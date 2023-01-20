import sympy as sp
from IPython.display import display, Math, Latex
sp.init_printing(use_unicode=True)
import numpy as np

def ptex(output): # print tex
    print(sp.latex(output))
def dtex(output): # display tex
    display(Math(sp.latex(output)))

class ucfunc:
    # f_py is the function expression in python format
     
    def __init__(self, expression, dependencies, dependencies_N, non_dependencies, non_dependencies_N, V_N):
        # define expression and dependencies of function
        self.dependencies_N = dependencies_N
        self.non_dependencies_N = non_dependencies_N
        self.dim_V = len(dependencies)
        
        self.V = sp.eye(self.dim_V) # correlation matrix
        
        # generate symbolic covariance matrix
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
        
        self.dependencies = []
        for variable in dependencies:
            exec(str(f"{variable} = sp.symbols('{variable}', real=True)"))
            exec(str(f"self.dependencies.append({variable})"))
            
        self.non_dependencies = []
        for variable in non_dependencies:
            exec(str(f"{variable} = sp.symbols('{variable}', real=True)"))
            exec(str(f"self.non_dependencies.append({variable})"))
        
        self.expression = eval(expression)
        
        # define g
        g = []
        for variable in self.dependencies:
            g.append(sp.diff(self.expression, variable))
        self.g = sp.Matrix(g)
        
        # calculate variance formula for corelated uncertainties
        self.var = ((self.g).T*(self.V)*(self.g))[0]
        
    def expression_tex(self):
        return sp.latex(self.expression)
    
    def var_tex(self):
        return sp.latex(self.var)
    
    def sigma_tex(self):
        return sp.latex(sp.sqrt(self.var))
    
    def var_N(self):
        varN = self.var # numerical
        for index, variable in enumerate(self.dependencies):
            varN = varN.subs(variable, self.dependencies_N[index])
        for index, variable in enumerate(self.non_dependencies):
            varN = varN.subs(variable, self.non_dependencies_N[index])
            
        for index, variable in enumerate(self.V_symbols_nonZeroVals):
            varN = varN.subs(variable, self.V_N_nonZeroVals[index]) # replace symbol with numerical value
            
        return varN
    
    def sigma_N(self):
        return np.sqrt(self.var_N())